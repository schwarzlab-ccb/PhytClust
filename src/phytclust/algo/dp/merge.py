"""Combine child clustering tables and record splits for backtracking."""

from dataclasses import dataclass

import numpy as np

_BATCH_MAX_CELLS = 4_000_000

_BATCH_MIN_STATES = 8


def split_offset_dtype(left_state_count: int, right_state_count: int):
    """Return an integer type that can store split offsets and the -1 sentinel."""
    width = min(left_state_count, right_state_count)
    if width <= 127:
        return np.int8
    if width <= 32767:
        return np.int16
    return np.int32


def decode_split_offset(offset: int, k: int, right_state_count: int) -> tuple[int, int]:
    """Recover left and right child states from a stored split offset."""
    left = max(0, k - 1 - (right_state_count - 1)) + offset
    return left, k - 1 - left


def _has_only_one_cluster_state(row) -> bool:
    """Return whether the table has no feasible split beyond one cluster."""
    return len(row) == 2 and not np.isfinite(row[1])


def _merge_with_single_cluster_child(
    left_costs,
    right_costs,
    max_cluster_count,
    costs_out,
    penalized_costs_out,
    split_offsets_out,
    left_cluster_count_sum_of_squares=None,
    cluster_count_sum_of_squares_out=None,
) -> None:
    """Merge a child with only one cluster state by shifting the other table."""
    single_cluster_on_left = _has_only_one_cluster_state(left_costs)
    single_cluster_costs, other_costs = (
        (left_costs, right_costs)
        if single_cluster_on_left
        else (right_costs, left_costs)
    )
    n = min(max_cluster_count, len(other_costs))
    combined_costs = other_costs[:n] + single_cluster_costs[0]
    costs_out[1 : n + 1] = combined_costs
    penalized_costs_out[1 : n + 1] = combined_costs
    if single_cluster_on_left:
        split_offsets_out[1 : n + 1] = 0
    else:
        split_offsets_out[1] = 0
        split_offsets_out[2 : n + 1] = 1
    if cluster_count_sum_of_squares_out is not None:
        cluster_counts = np.arange(1, n + 1)
        a = (
            np.zeros(n, dtype=np.int64)
            if single_cluster_on_left
            else cluster_counts - 1
        )
        cluster_count_sum_of_squares_out[1 : n + 1] = (
            left_cluster_count_sum_of_squares[a] + (cluster_counts - a) ** 2
        )
    infeasible = ~np.isfinite(combined_costs)
    if infeasible.any():
        split_offsets_out[1 : n + 1][infeasible] = -1


def _can_batch_merge(
    left_state_count: int, right_state_count: int, max_cluster_count: int
) -> bool:
    """Return whether the batch calculation fits the grid-size limit."""
    if max_cluster_count < _BATCH_MIN_STATES:
        return False
    column_count = min(left_state_count + right_state_count - 1, max_cluster_count)
    return min(left_state_count, right_state_count) * column_count <= _BATCH_MAX_CELLS


def _choose_balanced_split(
    mask, min_left_state: int, k: int, left_cluster_count_sum_of_squares=None
) -> int:
    """Choose the most evenly distributed split among the candidate states.

    Use stored sums of squared cluster counts when merging a prefix of children.
    For equally balanced splits, choose the lower left state."""
    candidate_indices = np.flatnonzero(mask)
    left_states = min_left_state + candidate_indices
    if left_cluster_count_sum_of_squares is None:
        key = np.abs(2 * left_states + 1 - k)
    else:
        key = left_cluster_count_sum_of_squares[left_states] + (k - left_states) ** 2
    return int(candidate_indices[key.argmin()])


def _store_cluster_count_sum_of_squares(
    cluster_count_sum_of_squares_out, left_cluster_count_sum_of_squares, k: int, a: int
) -> None:
    """Store the combined sum of squared cluster counts for a selected split."""
    if cluster_count_sum_of_squares_out is not None:
        cluster_count_sum_of_squares_out[k] = (
            left_cluster_count_sum_of_squares[a] + (k - a) ** 2
        )


def _merge_all_cluster_counts(
    left_costs,
    right_costs,
    left_outlier_counts,
    right_outlier_counts,
    max_cluster_count,
    relative_cost_tolerance,
    dtype,
    costs_out,
    penalized_costs_out,
    outlier_counts_out,
    split_offsets_out,
    left_cluster_count_sum_of_squares=None,
    cluster_count_sum_of_squares_out=None,
) -> None:
    """Calculate all requested cluster counts using a NumPy grid.

    Each column contains the possible splits for one parent cluster count.
    Use the shorter child table for rows and allocate only requested columns."""
    shorter_state_count, longer_state_count = len(left_costs), len(right_costs)
    right_state_count = longer_state_count
    original_left_state_count = shorter_state_count
    swapped = shorter_state_count > longer_state_count
    if swapped:
        left_costs, right_costs = right_costs, left_costs
        left_outlier_counts, right_outlier_counts = (
            right_outlier_counts,
            left_outlier_counts,
        )
        shorter_state_count, longer_state_count = (
            longer_state_count,
            shorter_state_count,
        )

    column_count = min(shorter_state_count + longer_state_count - 1, max_cluster_count)
    unreachable_outlier_count = np.iinfo(np.int32).max
    track_outlier_counts = (
        left_outlier_counts is not None and right_outlier_counts is not None
    )
    shorter_state_count = min(shorter_state_count, column_count)
    cost_grid = np.full((shorter_state_count, column_count), np.inf, dtype=dtype)
    outlier_count_grid = (
        np.full(
            (shorter_state_count, column_count),
            unreachable_outlier_count,
            dtype=np.int64,
        )
        if track_outlier_counts
        else None
    )
    for i in range(shorter_state_count):
        available_count = min(longer_state_count, column_count - i)
        cost_grid[i, i : i + available_count] = (
            left_costs[i] + right_costs[:available_count]
        )
        if track_outlier_counts:
            outlier_count_grid[i, i : i + available_count] = (
                left_outlier_counts[i] + right_outlier_counts[:available_count]
            )

    minimum_costs = cost_grid.min(axis=0)
    feasible = np.isfinite(minimum_costs)
    if not feasible.any():
        return
    cost_tolerance = np.zeros_like(minimum_costs)
    cost_tolerance[feasible] = relative_cost_tolerance * np.abs(minimum_costs[feasible])
    tied = (
        feasible[None, :]
        & np.isfinite(cost_grid)
        & (cost_grid <= (minimum_costs + cost_tolerance)[None, :])
    )
    if track_outlier_counts:
        min_ns = np.where(tied, outlier_count_grid, unreachable_outlier_count).min(
            axis=0
        )
        candidate_mask = tied & (outlier_count_grid == min_ns[None, :])
    else:
        candidate_mask = tied
    row_indices = np.arange(shorter_state_count)[:, None]
    column_indices = np.arange(column_count)[None, :]
    left_state_indices = (
        (column_indices - row_indices)
        if swapped
        else np.broadcast_to(row_indices, (shorter_state_count, column_count))
    )
    if left_cluster_count_sum_of_squares is None:
        tie_key = np.abs(2 * left_state_indices + 1 - (column_indices + 1)).astype(
            np.int64
        )
    else:
        safe_left = np.clip(left_state_indices, 0, original_left_state_count - 1)
        tie_key = (
            left_cluster_count_sum_of_squares[safe_left]
            + (column_indices + 1 - left_state_indices).astype(np.int64) ** 2
        )
    order = tie_key * (column_count + 1) + left_state_indices
    rows = np.where(candidate_mask, order, np.iinfo(np.int64).max).argmin(axis=0)

    feasible = np.isfinite(minimum_costs)
    feasible_columns = np.flatnonzero(feasible)
    if feasible_columns.size == 0:
        return
    selected_rows = rows[feasible_columns]
    cluster_counts = feasible_columns + 1
    costs_out[cluster_counts] = cost_grid[selected_rows, feasible_columns]
    penalized_costs_out[cluster_counts] = costs_out[cluster_counts]
    if track_outlier_counts:
        outlier_counts_out[cluster_counts] = outlier_count_grid[
            selected_rows, feasible_columns
        ]
    selected_left_states = (
        (feasible_columns - selected_rows) if swapped else selected_rows
    )
    split_offsets_out[cluster_counts] = selected_left_states - np.maximum(
        0, feasible_columns - (right_state_count - 1)
    )
    if cluster_count_sum_of_squares_out is not None:
        cluster_count_sum_of_squares_out[cluster_counts] = (
            left_cluster_count_sum_of_squares[selected_left_states]
            + (cluster_counts - selected_left_states) ** 2
        )


@dataclass(frozen=True)
class MergeSettings:
    """Settings shared by all child-table merges in a clustering run."""

    outlier_detection_enabled: bool
    small_cluster_penalty_enabled: bool
    prioritize_fewer_outliers: bool
    relative_cost_tolerance: float
    dtype: object


def merge_child_tables(
    left_costs,
    right_costs,
    left_penalized_costs,
    right_penalized_costs,
    left_outlier_counts,
    right_outlier_counts,
    max_cluster_count,
    settings,
    costs_out,
    penalized_costs_out,
    outlier_counts_out,
    split_offsets_out,
    left_cluster_count_sum_of_squares=None,
    cluster_count_sum_of_squares_out=None,
) -> None:
    """Combine two child tables and store selected splits for backtracking.

    Child state i represents i + 1 clusters. States i and j produce parent
    state i + j + 1. Leave impossible states at their initialized values.

    With penalties enabled, compare penalized cost, outlier count, then raw
    cost. Otherwise, prioritize outlier count only when requested; compare
    raw cost first in other cases. Resolve remaining ties by distributing
    clusters more evenly among children, then choosing the lower left state."""
    outlier_detection_enabled = settings.outlier_detection_enabled
    small_cluster_penalty_enabled = settings.small_cluster_penalty_enabled
    prioritize_fewer_outliers = settings.prioritize_fewer_outliers
    relative_cost_tolerance = settings.relative_cost_tolerance
    dtype = settings.dtype
    left_state_count, right_state_count = len(left_costs), len(right_costs)

    if (
        not small_cluster_penalty_enabled
        and not outlier_detection_enabled
        and max_cluster_count >= _BATCH_MIN_STATES
        and (
            _has_only_one_cluster_state(left_costs)
            or _has_only_one_cluster_state(right_costs)
        )
    ):
        _merge_with_single_cluster_child(
            left_costs,
            right_costs,
            max_cluster_count,
            costs_out,
            penalized_costs_out,
            split_offsets_out,
            left_cluster_count_sum_of_squares,
            cluster_count_sum_of_squares_out,
        )
        return

    if (
        not small_cluster_penalty_enabled
        and not prioritize_fewer_outliers
        and _can_batch_merge(left_state_count, right_state_count, max_cluster_count)
    ):
        _merge_all_cluster_counts(
            left_costs,
            right_costs,
            left_outlier_counts,
            right_outlier_counts,
            max_cluster_count,
            relative_cost_tolerance,
            dtype,
            costs_out,
            penalized_costs_out,
            outlier_counts_out,
            split_offsets_out,
            left_cluster_count_sum_of_squares,
            cluster_count_sum_of_squares_out,
        )
        return

    for k in range(1, max_cluster_count + 1):
        max_left_state = min(k - 1, left_state_count - 1)
        min_left_state = max(0, k - 1 - (right_state_count - 1))

        if min_left_state > max_left_state:
            continue

        if not outlier_detection_enabled:
            raw_scores = (
                left_costs[min_left_state : max_left_state + 1]
                + right_costs[k - 1 - max_left_state : k - min_left_state][::-1]
            )
            minimum_cost = float(raw_scores.min())
            if not np.isfinite(minimum_cost):
                continue
            cost_tolerance = relative_cost_tolerance * abs(minimum_cost)
            selected_local_index = _choose_balanced_split(
                raw_scores <= minimum_cost + cost_tolerance,
                min_left_state,
                k,
                left_cluster_count_sum_of_squares,
            )
            best_value = raw_scores[selected_local_index]

            costs_out[k] = best_value
            penalized_costs_out[k] = best_value
            split_offsets_out[k] = selected_local_index
            _store_cluster_count_sum_of_squares(
                cluster_count_sum_of_squares_out,
                left_cluster_count_sum_of_squares,
                k,
                min_left_state + selected_local_index,
            )
            continue

        right_start, right_end = k - 1 - max_left_state, k - min_left_state
        raw_scores = (
            left_costs[min_left_state : max_left_state + 1]
            + right_costs[right_start:right_end][::-1]
        )
        outlier_counts = (
            left_outlier_counts[min_left_state : max_left_state + 1]
            + right_outlier_counts[right_start:right_end][::-1]
        )

        if small_cluster_penalty_enabled:
            total_scores = (
                left_penalized_costs[min_left_state : max_left_state + 1]
                + right_penalized_costs[right_start:right_end][::-1]
            )
            minimum_penalized_cost = float(total_scores.min())
            if not np.isfinite(minimum_penalized_cost):
                continue
            cost_tolerance = relative_cost_tolerance * abs(minimum_penalized_cost)
            tied_mask = total_scores <= minimum_penalized_cost + cost_tolerance
            min_n = int(outlier_counts[tied_mask].min())
            candidate_mask = tied_mask & (outlier_counts == min_n)
            minimum_raw_cost = float(np.where(candidate_mask, raw_scores, np.inf).min())
            raw_cost_tolerance = relative_cost_tolerance * abs(minimum_raw_cost)
            best = _choose_balanced_split(
                candidate_mask & (raw_scores <= minimum_raw_cost + raw_cost_tolerance),
                min_left_state,
                k,
                left_cluster_count_sum_of_squares,
            )
            costs_out[k] = raw_scores[best]
            penalized_costs_out[k] = total_scores[best]
            outlier_counts_out[k] = outlier_counts[best]
            split_offsets_out[k] = best
            _store_cluster_count_sum_of_squares(
                cluster_count_sum_of_squares_out,
                left_cluster_count_sum_of_squares,
                k,
                min_left_state + best,
            )
            continue

        if prioritize_fewer_outliers:
            feasible = np.isfinite(raw_scores)
            if not feasible.any():
                continue
            min_n = int(outlier_counts[feasible].min())
            candidate_mask = feasible & (outlier_counts == min_n)
            minimum_raw_cost = float(np.where(candidate_mask, raw_scores, np.inf).min())
            raw_cost_tolerance = relative_cost_tolerance * abs(minimum_raw_cost)
            best = _choose_balanced_split(
                candidate_mask & (raw_scores <= minimum_raw_cost + raw_cost_tolerance),
                min_left_state,
                k,
                left_cluster_count_sum_of_squares,
            )
        else:
            minimum_cost = float(raw_scores.min())
            if not np.isfinite(minimum_cost):
                continue
            cost_tolerance = relative_cost_tolerance * abs(minimum_cost)
            tied_mask = raw_scores <= minimum_cost + cost_tolerance
            if int(np.count_nonzero(tied_mask)) > 1:
                min_n = int(outlier_counts[tied_mask].min())
                best = _choose_balanced_split(
                    tied_mask & (outlier_counts == min_n),
                    min_left_state,
                    k,
                    left_cluster_count_sum_of_squares,
                )
            else:
                best = int(tied_mask.argmax())

        costs_out[k] = raw_scores[best]
        penalized_costs_out[k] = raw_scores[best]
        outlier_counts_out[k] = outlier_counts[best]
        split_offsets_out[k] = best
        _store_cluster_count_sum_of_squares(
            cluster_count_sum_of_squares_out,
            left_cluster_count_sum_of_squares,
            k,
            min_left_state + best,
        )
