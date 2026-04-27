import numpy as np

from ..exceptions import ConfigurationError
from .dp_utils import (
    child_as_atomic_cluster_cost,
    node_support_factor,
    state_better,
    subtree_all_zero,
)


def compute_polytomy_dp(
    node,
    pc,
    max_states_global,
    min_cluster_size,
    outlier_thresh,
    dtype,
    bp_dtype,
):
    mode = getattr(pc, "polytomy_mode", "hard")
    if mode == "hard":
        return compute_polytomy_dp_hard(
            node,
            pc,
            max_states_global,
            min_cluster_size,
            outlier_thresh,
            dtype,
            bp_dtype,
        )
    return compute_polytomy_dp_soft(
        node,
        pc,
        max_states_global,
        min_cluster_size,
        outlier_thresh,
        dtype,
        bp_dtype,
    )


def compute_polytomy_dp_hard(
    node,
    pc,
    max_states_global,
    min_cluster_size,
    outlier_thresh,
    dtype,
    bp_dtype,
):
    """
    HARD polytomy:
      - state 0: all child clades merged into one cluster rooted at `node`
      - states >= 1: no cross-child merge at this node at all
                    (each cluster is entirely inside one child subtree)

    This matches a true hard multifurcation.
    """
    children = list(node.clades)
    node_id = pc.node_to_id[node]
    n_leaves = pc.num_leaves_per_node[node]
    n_states = min(n_leaves, max_states_global)

    use_outlier = outlier_thresh is not None
    prefer_fewer = pc.outlier.prefer_fewer

    raw_array = np.full(n_states + 1, np.inf, dtype=dtype)
    total_array = np.full(n_states + 1, np.inf, dtype=dtype)

    if use_outlier:
        ns_array = np.zeros(n_states + 1, dtype=np.int32)

    # Explicit one-cluster state: all children merged at this node.
    merged_all = 0.0
    feasible_merge_all = True
    for child in children:
        atom_cost = child_as_atomic_cluster_cost(child, pc)
        if not np.isfinite(atom_cost):
            feasible_merge_all = False
            break
        merged_all += atom_cost

    if feasible_merge_all:
        merged_all /= node_support_factor(node, pc)
    else:
        merged_all = np.inf

    pc.cluster_cost[node] = float(merged_all)

    if n_leaves >= min_cluster_size:
        raw_array[0] = merged_all
        total_array[0] = merged_all
    else:
        raw_array[0] = np.inf
        total_array[0] = np.inf

    if use_outlier:
        ns_array[0] = 1 if n_leaves < outlier_thresh else 0

    # Optional: disallow any split if everything below is zero-length.
    zero_eps = getattr(pc, "zero_length_eps", 1e-12)
    no_split_zero = getattr(pc, "no_split_zero_length", False)
    if no_split_zero and subtree_all_zero(node, pc, eps=zero_eps):
        pc._n_small[node_id] = ns_array if use_outlier else None
        return total_array, raw_array, {"mode": "hard", "steps": []}

    # Split-only DP across children.
    # Prefix DP over children, but without allowing merged proper subsets.
    # Therefore only k>=1 states are propagated once we have >1 child.
    first_child = children[0]
    first_id = pc.node_to_id[first_child]

    prefix_raw = pc.raw_dp_table[first_id].copy()
    n_prefix_leaves = pc.num_leaves_per_node[first_child]

    if use_outlier:
        prefix_ns = pc._n_small[first_id].copy()

    steps = []

    for child in children[1:]:
        child_id = pc.node_to_id[child]
        child_raw = pc.raw_dp_table[child_id]
        n_child = pc.num_leaves_per_node[child]

        if use_outlier:
            child_ns = pc._n_small[child_id]

        n_total = n_prefix_leaves + n_child
        n_new_states = min(n_total, max_states_global)

        new_raw = np.full(n_new_states + 1, np.inf, dtype=dtype)
        step_bp = np.full(n_new_states + 1, -1, dtype=bp_dtype)

        if use_outlier:
            new_ns = np.zeros(n_new_states + 1, dtype=np.int32)

        # Note: new_raw[0] stays inf, because with >1 child we do not allow
        # a proper subset of children to be merged in hard mode.
        prefix_len = len(prefix_raw)
        child_len = len(child_raw)

        for k in range(1, n_new_states + 1):
            max_i = min(k - 1, prefix_len - 1)
            min_i = max(0, k - 1 - (child_len - 1))
            if min_i > max_i:
                continue

            # Fast path: no outlier handling. Mirrors the binary-node fast
            # path in dp.py — fused add + argmin, ~4-6x faster per
            # iteration than the previous Python loop.
            if not use_outlier:
                left_slice = prefix_raw[min_i:max_i + 1]
                right_slice = child_raw[k - 1 - max_i:k - min_i][::-1]
                raw_scores = left_slice + right_slice
                best_local = int(raw_scores.argmin())
                best_value = raw_scores[best_local]
                if not np.isfinite(best_value):
                    continue
                new_raw[k] = best_value
                step_bp[k] = min_i + best_local
                continue

            # Outlier-aware path (unchanged semantics).
            best_raw = np.inf
            best_ns = np.iinfo(np.int32).max
            best_i = -1
            for i in range(min_i, max_i + 1):
                j = k - 1 - i
                cand_raw = prefix_raw[i] + child_raw[j]
                cand_ns = int(prefix_ns[i] + child_ns[j])

                if state_better(
                    cand_raw,
                    cand_ns,
                    best_raw,
                    best_ns,
                    use_outlier=use_outlier,
                    prefer_fewer=prefer_fewer,
                ):
                    best_raw = cand_raw
                    best_ns = cand_ns
                    best_i = i

            new_raw[k] = best_raw
            step_bp[k] = best_i
            if best_i >= 0:
                new_ns[k] = best_ns

        steps.append(step_bp)
        prefix_raw = new_raw
        if use_outlier:
            prefix_ns = new_ns
        n_prefix_leaves = n_total

    # Copy split states into the node table.
    upto = min(n_states, len(prefix_raw) - 1)
    raw_array[1 : upto + 1] = prefix_raw[1 : upto + 1]
    total_array[1 : upto + 1] = prefix_raw[1 : upto + 1]

    if use_outlier:
        ns_array[1 : upto + 1] = prefix_ns[1 : upto + 1]
        pc._n_small[node_id] = ns_array

    return total_array, raw_array, {"mode": "hard", "steps": steps}


def compute_polytomy_dp_soft(
    node,
    pc,
    max_states_global,
    min_cluster_size,
    outlier_thresh,
    dtype,
    bp_dtype,
):
    """
    SOFT polytomy:
      - any subset of sibling child clades may form a hidden zero-length clade
      - singleton blocks are handled by the child's own DP state
      - blocks of size >= 2 become one new cluster rooted at a hidden node at
        the polytomy (same distance effect as rooting that group at `node`)

    This is exponential in the node degree.
    """
    children = list(node.clades)
    m = len(children)
    max_deg = getattr(pc, "soft_polytomy_max_degree", 18)
    if m > max_deg:
        raise ConfigurationError(
            f"soft polytomy DP is exponential in node degree; "
            f"node has degree {m}, limit is {max_deg}. "
            f"Increase soft_polytomy_max_degree or use polytomy_mode='hard'."
        )

    node_id = pc.node_to_id[node]
    n_leaves = pc.num_leaves_per_node[node]
    n_states = min(n_leaves, max_states_global)

    use_outlier = outlier_thresh is not None
    prefer_fewer = pc.outlier.prefer_fewer

    raw_array = np.full(n_states + 1, np.inf, dtype=dtype)
    total_array = np.full(n_states + 1, np.inf, dtype=dtype)

    if use_outlier:
        ns_array = np.zeros(n_states + 1, dtype=np.int32)

    # Explicit whole-node one-cluster state
    atom_costs = []
    child_leaf_counts = []
    feasible_merge_all = True
    merged_all = 0.0

    for child in children:
        atom_cost = child_as_atomic_cluster_cost(child, pc)
        atom_costs.append(atom_cost)
        child_leaf_counts.append(pc.num_leaves_per_node[child])

        if not np.isfinite(atom_cost):
            feasible_merge_all = False
        else:
            merged_all += atom_cost

    if feasible_merge_all:
        merged_all /= node_support_factor(node, pc)
    else:
        merged_all = np.inf

    pc.cluster_cost[node] = float(merged_all)

    if n_leaves >= min_cluster_size:
        raw_array[0] = merged_all
        total_array[0] = merged_all
    else:
        raw_array[0] = np.inf
        total_array[0] = np.inf

    if use_outlier:
        ns_array[0] = 1 if n_leaves < outlier_thresh else 0

    zero_eps = getattr(pc, "zero_length_eps", 1e-12)
    no_split_zero = getattr(pc, "no_split_zero_length", False)
    if no_split_zero and subtree_all_zero(node, pc, eps=zero_eps):
        pc._n_small[node_id] = ns_array if use_outlier else None
        return total_array, raw_array, {"mode": "soft", "parent": {}}

    fullmask = (1 << m) - 1

    # Precompute block leaf counts and block costs for subset merges.
    # A block with >=2 children means: take each child as one whole clade,
    # then merge those child-clades into one cluster at a hidden zero-length node.
    block_leaves = [0] * (fullmask + 1)
    block_cost = [0.0] * (fullmask + 1)
    block_feasible = [True] * (fullmask + 1)

    for mask in range(1, fullmask + 1):
        lsb = mask & -mask
        i = lsb.bit_length() - 1
        prev = mask ^ lsb
        block_leaves[mask] = block_leaves[prev] + child_leaf_counts[i]

        if not block_feasible[prev] or not np.isfinite(atom_costs[i]):
            block_feasible[mask] = False
            block_cost[mask] = np.inf
        else:
            block_cost[mask] = block_cost[prev] + atom_costs[i]

    # dp[mask][q] = best raw cost to cover the child subset `mask`
    # with exactly q clusters.
    dp = {0: {0: 0.0}}
    ns_dp = {0: {0: 0}} if use_outlier else None

    # parent[(mask, q)] = (prev_mask, prev_q, block_mask, child_state_idx)
    # If child_state_idx is None => block_mask is a merged multi-child block.
    # Else block_mask must be a singleton and child_state_idx is the child's DP index.
    parent = {}

    for mask in range(fullmask + 1):
        if mask not in dp:
            continue

        rem = fullmask ^ mask
        if rem == 0:
            continue

        # Canonical next child: smallest not-yet-covered child.
        bit_t = rem & -rem
        t = bit_t.bit_length() - 1
        child = children[t]
        child_id = pc.node_to_id[child]
        child_raw = pc.raw_dp_table[child_id]
        child_ns = pc._n_small[child_id] if use_outlier else None

        for q_used, base_raw in dp[mask].items():
            base_ns = ns_dp[mask][q_used] if use_outlier else 0

            # Option 1: handle child t independently via one of its own DP states.
            for s, cr in enumerate(child_raw):
                if not np.isfinite(cr):
                    continue
                q2 = q_used + (s + 1)
                if q2 > n_states + 1:
                    continue

                cand_raw = base_raw + cr
                cand_ns = base_ns + (int(child_ns[s]) if use_outlier else 0)

                mask2 = mask | bit_t
                if mask2 not in dp:
                    dp[mask2] = {}
                    if use_outlier:
                        ns_dp[mask2] = {}

                best_raw = dp[mask2].get(q2, np.inf)
                best_ns = (
                    ns_dp[mask2].get(q2, np.iinfo(np.int32).max) if use_outlier else 0
                )

                if state_better(
                    cand_raw,
                    cand_ns,
                    best_raw,
                    best_ns,
                    use_outlier=use_outlier,
                    prefer_fewer=prefer_fewer,
                ):
                    dp[mask2][q2] = cand_raw
                    if use_outlier:
                        ns_dp[mask2][q2] = cand_ns
                    parent[(mask2, q2)] = (mask, q_used, bit_t, s)

            # Option 2: merge child t with any nonempty subset of the remaining siblings.
            rest = rem ^ bit_t
            sub = rest
            while sub:
                block = bit_t | sub
                if block.bit_count() >= 2:
                    if (
                        block_feasible[block]
                        and block_leaves[block] >= min_cluster_size
                    ):
                        q2 = q_used + 1
                        if q2 <= n_states + 1:
                            cand_raw = base_raw + block_cost[block]
                            cand_ns = base_ns + (
                                1
                                if (
                                    use_outlier and block_leaves[block] < outlier_thresh
                                )
                                else 0
                            )

                            mask2 = mask | block
                            if mask2 not in dp:
                                dp[mask2] = {}
                                if use_outlier:
                                    ns_dp[mask2] = {}

                            best_raw = dp[mask2].get(q2, np.inf)
                            best_ns = (
                                ns_dp[mask2].get(q2, np.iinfo(np.int32).max)
                                if use_outlier
                                else 0
                            )

                            if state_better(
                                cand_raw,
                                cand_ns,
                                best_raw,
                                best_ns,
                                use_outlier=use_outlier,
                                prefer_fewer=prefer_fewer,
                            ):
                                dp[mask2][q2] = cand_raw
                                if use_outlier:
                                    ns_dp[mask2][q2] = cand_ns
                                parent[(mask2, q2)] = (mask, q_used, block, None)

                sub = (sub - 1) & rest

    # Fill node states from subset DP.
    final_states = dp.get(fullmask, {})
    for q, score in final_states.items():
        k = q - 1
        if 0 <= k <= n_states:
            # Keep explicit whole-node state at k=0.
            if k == 0:
                continue
            raw_array[k] = score
            total_array[k] = score
            if use_outlier:
                ns_array[k] = ns_dp[fullmask][q]

    if use_outlier:
        pc._n_small[node_id] = ns_array

    return total_array, raw_array, {"mode": "soft", "parent": parent}
