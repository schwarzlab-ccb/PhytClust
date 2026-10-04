"""Compare merge implementations against explicit candidate enumeration."""

import numpy as np
import pytest
from phytclust.algo.dp import merge


def reference(
    left,
    right,
    left_total,
    right_total,
    left_counts,
    right_counts,
    cap,
    settings,
    left_balance=None,
):
    outputs = [
        np.full(cap + 1, np.inf),
        np.full(cap + 1, np.inf),
        np.zeros(cap + 1, dtype=int),
        np.full(cap + 1, -1, dtype=int),
    ]
    for k in range(1, cap + 1):
        low = max(0, k - len(right))
        candidates = []
        for i in range(low, min(k, len(left))):
            j = k - 1 - i
            cost = left[i] + right[j]
            total = left_total[i] + right_total[j]
            if not np.isfinite(cost) or not np.isfinite(total):
                continue
            count = (
                left_counts[i] + right_counts[j]
                if settings.outlier_detection_enabled
                else 0
            )
            candidates.append(
                (
                    cost,
                    total,
                    count,
                    ((i + 1) ** 2 if left_balance is None else left_balance[i])
                    + (j + 1) ** 2,
                    i,
                )
            )
        if not candidates:
            continue

        def filter_cost(index):
            nonlocal candidates
            best = min(c[index] for c in candidates)
            candidates = [
                c
                for c in candidates
                if c[index] <= best + settings.relative_cost_tolerance * abs(best)
            ]

        def filter_count():
            nonlocal candidates
            best = min(c[2] for c in candidates)
            candidates = [c for c in candidates if c[2] == best]

        if settings.small_cluster_penalty_enabled:
            filter_cost(1)
            filter_count()
            filter_cost(0)
        elif settings.outlier_detection_enabled and settings.prioritize_fewer_outliers:
            filter_count()
            filter_cost(0)
        else:
            filter_cost(0)
            if settings.outlier_detection_enabled:
                filter_count()
        chosen = min(candidates, key=lambda c: (c[3], c[4]))
        outputs[0][k], outputs[1][k], outputs[2][k], outputs[3][k] = (
            chosen[0],
            chosen[1] if settings.small_cluster_penalty_enabled else chosen[0],
            chosen[2],
            chosen[4] - low,
        )
    return outputs


@pytest.mark.parametrize("lengths", [(2, 12), (12, 2), (9, 13), (13, 9)])
@pytest.mark.parametrize("mode", ["cost", "outliers", "fewer", "penalty"])
@pytest.mark.parametrize("prefix_balance", [False, True])
@pytest.mark.parametrize("batched", [False, True])
def test_merge_matches_reference(monkeypatch, lengths, mode, batched, prefix_balance):
    monkeypatch.setattr(merge, "_BATCH_MIN_STATES", 8 if batched else 10000)
    rng = np.random.default_rng(41)
    left, right = [rng.integers(0, 5, n).astype(float) for n in lengths]
    if lengths[0] == 2:
        left[1] = np.inf
    if lengths[1] == 2:
        right[1] = np.inf
    left_total, right_total = left + 2, right + 1
    left_counts, right_counts = [rng.integers(0, 3, n) for n in lengths]
    settings = merge.MergeSettings(
        mode != "cost", mode == "penalty", mode == "fewer", 0.0, np.float64
    )
    cap = 10
    left_balance = (
        np.arange(len(left), dtype=np.int64) ** 2 + 7 if prefix_balance else None
    )
    balance_output = np.zeros(cap + 1, dtype=np.int64) if prefix_balance else None
    expected = reference(
        left,
        right,
        left_total,
        right_total,
        left_counts,
        right_counts,
        cap,
        settings,
        left_balance=left_balance,
    )
    outputs = [
        np.full(cap + 1, np.inf),
        np.full(cap + 1, np.inf),
        np.zeros(cap + 1, dtype=int),
        np.full(cap + 1, -1, dtype=int),
    ]
    with np.errstate(invalid="raise"):
        merge.merge_child_tables(
            left,
            right,
            left_total,
            right_total,
            left_counts if mode != "cost" else None,
            right_counts if mode != "cost" else None,
            cap,
            settings,
            *outputs,
            left_cluster_count_sum_of_squares=left_balance,
            cluster_count_sum_of_squares_out=balance_output,
        )
    for actual, wanted in zip(outputs, expected):
        np.testing.assert_array_equal(actual, wanted)


@pytest.mark.parametrize("batched", [False, True])
@pytest.mark.parametrize("lengths", [(9, 13), (13, 9)])
@pytest.mark.parametrize("mode", ["cost", "outliers", "fewer", "penalty"])
@pytest.mark.parametrize("cluster_count", [8, 9])
def test_cost_ties_choose_smallest_child_state_difference(
    monkeypatch, batched, lengths, mode, cluster_count
):
    monkeypatch.setattr(merge, "_BATCH_MIN_STATES", 8 if batched else 10000)
    left_costs, right_costs = [np.zeros(length) for length in lengths]
    left_counts, right_counts = [np.zeros(length, dtype=int) for length in lengths]
    settings = merge.MergeSettings(
        mode != "cost", mode == "penalty", mode == "fewer", 0.0, np.float64
    )
    outputs = [
        np.full(cluster_count + 1, np.inf),
        np.full(cluster_count + 1, np.inf),
        np.zeros(cluster_count + 1, dtype=int),
        np.full(cluster_count + 1, -1, dtype=int),
    ]
    merge.merge_child_tables(
        left_costs, right_costs, left_costs, right_costs,
        left_counts if mode != "cost" else None,
        right_counts if mode != "cost" else None,
        cluster_count, settings, *outputs,
    )
    left_state, right_state = merge.decode_split_offset(
        int(outputs[3][cluster_count]), cluster_count, len(right_costs)
    )
    candidates = [
        (i, cluster_count - 1 - i)
        for i in range(len(left_costs))
        if 0 <= cluster_count - 1 - i < len(right_costs)
    ]
    assert (left_state, right_state) == min(
        candidates, key=lambda pair: (abs(pair[0] - pair[1]), pair[0])
    )
