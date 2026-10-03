import numpy as np

from ...exceptions import ConfigurationError
from .merge import split_offset_dtype, MergeSettings, merge_child_tables
from .costs import (
    child_subtree_cluster_cost,
    floating_point_precision,
    branch_support_fraction,
    is_better_partition,
    subtree_has_only_zero_length_branches,
)


def _whole_subtree_cost(node, clustering, child_costs):
    """Return the cost of one cluster containing every child subtree."""
    if not all(np.isfinite(cost) for cost in child_costs):
        return np.inf
    return sum(child_costs) / branch_support_fraction(node, clustering)


def compute_polytomy_dp(
    node,
    clustering,
    global_max_cluster_count,
    min_cluster_size,
    outlier_size_threshold,
    dtype,
):
    """Choose hard or soft mode for a node with more than two children."""
    mode = getattr(clustering, "polytomy_mode", "soft")
    if mode == "hard":
        return compute_polytomy_dp_hard(
            node,
            clustering,
            global_max_cluster_count,
            min_cluster_size,
            outlier_size_threshold,
            dtype,
        )
    if mode != "soft":
        raise ConfigurationError("polytomy_mode must be 'hard' or 'soft'.")
    return compute_polytomy_dp_soft(
        node,
        clustering,
        global_max_cluster_count,
        min_cluster_size,
        outlier_size_threshold,
        dtype,
    )


def compute_polytomy_dp_hard(
    node,
    clustering,
    global_max_cluster_count,
    min_cluster_size,
    outlier_size_threshold,
    dtype,
):
    """Compute costs for a multifurcation without partial sibling groups.

    Allow one cluster containing the whole subtree. Other clusters must
    remain within individual child subtrees. Store each selected child split."""
    children = list(node.clades)
    node_id = clustering.node_to_id[node]
    leaf_count = clustering.num_leaves_per_node[node]
    max_cluster_count = min(leaf_count, global_max_cluster_count)

    outlier_detection_enabled = outlier_size_threshold is not None
    prioritize_fewer_outliers = clustering.outlier.prefer_fewer
    settings = MergeSettings(
        outlier_detection_enabled,
        False,
        prioritize_fewer_outliers,
        8.0 * floating_point_precision(dtype),
        dtype,
    )

    costs = np.full(max_cluster_count + 1, np.inf, dtype=dtype)
    total_costs = costs

    if outlier_detection_enabled:
        outlier_counts = np.zeros(max_cluster_count + 1, dtype=np.int32)

    child_subtree_costs = [
        child_subtree_cluster_cost(child, clustering) for child in children
    ]
    whole_subtree_cost = _whole_subtree_cost(node, clustering, child_subtree_costs)

    clustering.cluster_cost[node] = float(whole_subtree_cost)

    if leaf_count >= min_cluster_size:
        costs[0] = whole_subtree_cost
        total_costs[0] = whole_subtree_cost
    else:
        costs[0] = np.inf
        total_costs[0] = np.inf

    if outlier_detection_enabled:
        outlier_counts[0] = 1 if leaf_count < outlier_size_threshold else 0

    zero_length_tolerance = getattr(clustering, "zero_length_eps", 1e-12)
    no_split_zero = getattr(clustering, "no_split_zero_length", False)
    if no_split_zero and subtree_has_only_zero_length_branches(
        node, clustering, zero_length_tolerance=zero_length_tolerance
    ):
        if outlier_detection_enabled:
            clustering._outlier_counts[node_id] = outlier_counts
        return total_costs, costs, {"mode": "hard", "steps": []}

    first_child = children[0]
    first_id = clustering.node_to_id[first_child]

    prefix_costs = clustering.raw_dp_table[first_id].copy()
    prefix_leaf_count = clustering.num_leaves_per_node[first_child]

    prefix_outlier_counts = (
        clustering._outlier_counts[first_id].copy()
        if outlier_detection_enabled
        else None
    )
    prefix_cluster_balance = (np.arange(len(prefix_costs), dtype=np.int64) + 1) ** 2

    steps = []

    for child in children[1:]:
        child_id = clustering.node_to_id[child]
        child_costs = clustering.raw_dp_table[child_id]
        child_leaf_count = clustering.num_leaves_per_node[child]
        child_outlier_counts = (
            clustering._outlier_counts[child_id] if outlier_detection_enabled else None
        )

        combined_leaf_count = prefix_leaf_count + child_leaf_count
        combined_max_cluster_count = min(combined_leaf_count, global_max_cluster_count)

        combined_costs = np.full(combined_max_cluster_count + 1, np.inf, dtype=dtype)
        combined_outlier_counts = (
            np.zeros(combined_max_cluster_count + 1, dtype=np.int32)
            if outlier_detection_enabled
            else None
        )
        combined_cluster_balance = np.zeros(
            combined_max_cluster_count + 1, dtype=np.int64
        )
        split_offsets = np.full(
            combined_max_cluster_count + 1,
            -1,
            dtype=split_offset_dtype(len(prefix_costs), len(child_costs)),
        )

        merge_child_tables(
            prefix_costs,
            child_costs,
            prefix_costs,
            child_costs,
            prefix_outlier_counts,
            child_outlier_counts,
            combined_max_cluster_count,
            settings,
            combined_costs,
            combined_costs,
            combined_outlier_counts,
            split_offsets,
            left_cluster_count_sum_of_squares=prefix_cluster_balance,
            cluster_count_sum_of_squares_out=combined_cluster_balance,
        )

        steps.append(split_offsets)
        prefix_costs = combined_costs
        prefix_cluster_balance = combined_cluster_balance
        if outlier_detection_enabled:
            prefix_outlier_counts = combined_outlier_counts
        prefix_leaf_count = combined_leaf_count

    last_state = min(max_cluster_count, len(prefix_costs) - 1)
    costs[1 : last_state + 1] = prefix_costs[1 : last_state + 1]
    total_costs[1 : last_state + 1] = prefix_costs[1 : last_state + 1]

    if outlier_detection_enabled:
        outlier_counts[1 : last_state + 1] = prefix_outlier_counts[1 : last_state + 1]
        clustering._outlier_counts[node_id] = outlier_counts

    return total_costs, costs, {"mode": "hard", "steps": steps}


def compute_polytomy_dp_soft(
    node,
    clustering,
    global_max_cluster_count,
    min_cluster_size,
    outlier_size_threshold,
    dtype,
):
    """Compute costs for an unresolved multifurcation.

    Allow groups of two or more child subtrees to form one cluster.
    Individual children may use their own clustering states.

    Store reconstruction choices for selected states. Computation grows
    exponentially with the number of children."""
    children = list(node.clades)
    child_count = len(children)
    max_child_count = getattr(clustering, "soft_polytomy_max_degree", 12)
    if child_count > max_child_count:
        raise ConfigurationError(
            f"Node {node.name!r} has {child_count} children; soft mode allows at most {max_child_count}. "
            f"Increase soft_polytomy_max_degree or use polytomy_mode='hard'."
        )

    node_id = clustering.node_to_id[node]
    leaf_count = clustering.num_leaves_per_node[node]
    max_cluster_count = min(leaf_count, global_max_cluster_count)

    outlier_detection_enabled = outlier_size_threshold is not None
    prioritize_fewer_outliers = clustering.outlier.prefer_fewer
    relative_cost_tolerance = 8.0 * floating_point_precision(dtype)

    costs = np.full(max_cluster_count + 1, np.inf, dtype=dtype)
    total_costs = costs

    if outlier_detection_enabled:
        outlier_counts = np.zeros(max_cluster_count + 1, dtype=np.int32)

    child_subtree_costs = [
        child_subtree_cluster_cost(child, clustering) for child in children
    ]
    child_leaf_counts = [clustering.num_leaves_per_node[child] for child in children]
    whole_subtree_cost = _whole_subtree_cost(node, clustering, child_subtree_costs)

    clustering.cluster_cost[node] = float(whole_subtree_cost)

    if leaf_count >= min_cluster_size:
        costs[0] = whole_subtree_cost
        total_costs[0] = whole_subtree_cost
    else:
        costs[0] = np.inf
        total_costs[0] = np.inf

    if outlier_detection_enabled:
        outlier_counts[0] = 1 if leaf_count < outlier_size_threshold else 0

    zero_length_tolerance = getattr(clustering, "zero_length_eps", 1e-12)
    no_split_zero = getattr(clustering, "no_split_zero_length", False)
    if no_split_zero and subtree_has_only_zero_length_branches(
        node, clustering, zero_length_tolerance=zero_length_tolerance
    ):
        if outlier_detection_enabled:
            clustering._outlier_counts[node_id] = outlier_counts
        return total_costs, costs, {"mode": "soft", "parent": {}}

    all_children_mask = (1 << child_count) - 1

    group_leaf_counts = [0] * (all_children_mask + 1)
    group_costs = [0.0] * (all_children_mask + 1)
    group_feasible = [True] * (all_children_mask + 1)

    for assigned_children_mask in range(1, all_children_mask + 1):
        first_child_bit = assigned_children_mask & -assigned_children_mask
        i = first_child_bit.bit_length() - 1
        previous_mask = assigned_children_mask ^ first_child_bit
        group_leaf_counts[assigned_children_mask] = (
            group_leaf_counts[previous_mask] + child_leaf_counts[i]
        )

        if not group_feasible[previous_mask] or not np.isfinite(child_subtree_costs[i]):
            group_feasible[assigned_children_mask] = False
            group_costs[assigned_children_mask] = np.inf
        else:
            group_costs[assigned_children_mask] = (
                group_costs[previous_mask] + child_subtree_costs[i]
            )

    cost_table = {0: {0: 0.0}}
    outlier_count_table = {0: {0: 0}} if outlier_detection_enabled else None
    cluster_balance_table = {0: {0: 0}}

    reconstruction_choices = {}

    for assigned_children_mask in range(all_children_mask + 1):
        if assigned_children_mask not in cost_table:
            continue

        remaining_children_mask = all_children_mask ^ assigned_children_mask
        if remaining_children_mask == 0:
            continue

        first_remaining_child_bit = remaining_children_mask & -remaining_children_mask
        first_remaining_child_index = first_remaining_child_bit.bit_length() - 1
        child = children[first_remaining_child_index]
        child_id = clustering.node_to_id[child]
        child_costs = clustering.raw_dp_table[child_id]
        child_outlier_counts = (
            clustering._outlier_counts[child_id] if outlier_detection_enabled else None
        )

        feasible_groups = []
        other_remaining_children_mask = (
            remaining_children_mask ^ first_remaining_child_bit
        )
        subset_mask = other_remaining_children_mask
        while subset_mask:
            group_mask = first_remaining_child_bit | subset_mask
            if (
                group_feasible[group_mask]
                and group_leaf_counts[group_mask] >= min_cluster_size
            ):
                group_outlier_count = int(
                    outlier_detection_enabled
                    and group_leaf_counts[group_mask] < outlier_size_threshold
                )
                feasible_groups.append(
                    (group_mask, group_costs[group_mask], group_outlier_count)
                )
            subset_mask = (subset_mask - 1) & other_remaining_children_mask

        for current_cluster_count, current_cost in cost_table[
            assigned_children_mask
        ].items():
            current_outlier_count = (
                outlier_count_table[assigned_children_mask][current_cluster_count]
                if outlier_detection_enabled
                else 0
            )
            current_cluster_balance = cluster_balance_table[assigned_children_mask][
                current_cluster_count
            ]

            for child_state, child_cost in enumerate(child_costs):
                if not np.isfinite(child_cost):
                    continue
                next_cluster_count = current_cluster_count + (child_state + 1)
                if next_cluster_count > max_cluster_count:
                    continue

                candidate_cost = current_cost + child_cost
                candidate_outlier_count = current_outlier_count + (
                    int(child_outlier_counts[child_state])
                    if outlier_detection_enabled
                    else 0
                )
                candidate_cluster_count_sum_of_squares = current_cluster_balance + (
                    child_state + 1
                ) * (child_state + 1)

                next_children_mask = assigned_children_mask | first_remaining_child_bit
                if next_children_mask not in cost_table:
                    cost_table[next_children_mask] = {}
                    cluster_balance_table[next_children_mask] = {}
                    if outlier_detection_enabled:
                        outlier_count_table[next_children_mask] = {}

                best_cost = cost_table[next_children_mask].get(
                    next_cluster_count, np.inf
                )
                best_outlier_count = (
                    outlier_count_table[next_children_mask].get(
                        next_cluster_count, np.iinfo(np.int32).max
                    )
                    if outlier_detection_enabled
                    else 0
                )
                best_cluster_count_sum_of_squares = cluster_balance_table[
                    next_children_mask
                ].get(next_cluster_count, np.iinfo(np.int32).max)

                if is_better_partition(
                    candidate_cost,
                    candidate_outlier_count,
                    best_cost,
                    best_outlier_count,
                    outlier_detection_enabled=outlier_detection_enabled,
                    prioritize_fewer_outliers=prioritize_fewer_outliers,
                    relative_cost_tolerance=relative_cost_tolerance,
                    candidate_cluster_count_sum_of_squares=candidate_cluster_count_sum_of_squares,
                    best_cluster_count_sum_of_squares=best_cluster_count_sum_of_squares,
                ):
                    cost_table[next_children_mask][next_cluster_count] = candidate_cost
                    cluster_balance_table[next_children_mask][
                        next_cluster_count
                    ] = candidate_cluster_count_sum_of_squares
                    if outlier_detection_enabled:
                        outlier_count_table[next_children_mask][
                            next_cluster_count
                        ] = candidate_outlier_count
                    reconstruction_choices[(next_children_mask, next_cluster_count)] = (
                        assigned_children_mask,
                        current_cluster_count,
                        first_remaining_child_bit,
                        child_state,
                    )

            for group_mask, group_cost, group_outlier_count in feasible_groups:
                next_cluster_count = current_cluster_count + 1
                if next_cluster_count <= max_cluster_count:
                    candidate_cost = current_cost + group_cost
                    candidate_outlier_count = (
                        current_outlier_count + group_outlier_count
                    )

                    next_children_mask = assigned_children_mask | group_mask
                    if next_children_mask not in cost_table:
                        cost_table[next_children_mask] = {}
                        cluster_balance_table[next_children_mask] = {}
                        if outlier_detection_enabled:
                            outlier_count_table[next_children_mask] = {}

                    candidate_cluster_count_sum_of_squares = current_cluster_balance + 1

                    best_cost = cost_table[next_children_mask].get(
                        next_cluster_count, np.inf
                    )
                    best_outlier_count = (
                        outlier_count_table[next_children_mask].get(
                            next_cluster_count, np.iinfo(np.int32).max
                        )
                        if outlier_detection_enabled
                        else 0
                    )
                    best_cluster_count_sum_of_squares = cluster_balance_table[
                        next_children_mask
                    ].get(next_cluster_count, np.iinfo(np.int32).max)

                    if is_better_partition(
                        candidate_cost,
                        candidate_outlier_count,
                        best_cost,
                        best_outlier_count,
                        outlier_detection_enabled=outlier_detection_enabled,
                        prioritize_fewer_outliers=prioritize_fewer_outliers,
                        relative_cost_tolerance=relative_cost_tolerance,
                        candidate_cluster_count_sum_of_squares=candidate_cluster_count_sum_of_squares,
                        best_cluster_count_sum_of_squares=best_cluster_count_sum_of_squares,
                    ):
                        cost_table[next_children_mask][
                            next_cluster_count
                        ] = candidate_cost
                        cluster_balance_table[next_children_mask][
                            next_cluster_count
                        ] = candidate_cluster_count_sum_of_squares
                        if outlier_detection_enabled:
                            outlier_count_table[next_children_mask][
                                next_cluster_count
                            ] = candidate_outlier_count
                        reconstruction_choices[
                            (next_children_mask, next_cluster_count)
                        ] = (
                            assigned_children_mask,
                            current_cluster_count,
                            group_mask,
                            None,
                        )

    final_states = cost_table.get(all_children_mask, {})
    for cluster_count, score in final_states.items():
        state_index = cluster_count - 1
        if 0 <= state_index <= max_cluster_count:
            if state_index == 0:
                continue
            costs[state_index] = score
            total_costs[state_index] = score
            if outlier_detection_enabled:
                outlier_counts[state_index] = outlier_count_table[all_children_mask][
                    cluster_count
                ]

    if outlier_detection_enabled:
        clustering._outlier_counts[node_id] = outlier_counts

    return total_costs, costs, {"mode": "soft", "parent": reconstruction_choices}
