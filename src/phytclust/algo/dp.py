import numpy as np
from math import ceil
from collections import Counter
from typing import Any
import logging


from ..exceptions import (
    ConfigurationError,
    InvalidKError,
    MissingDPTableError,
    InvalidClusteringError,
)
from .dp_polytomy import compute_polytomy_dp
from .dp_utils import (
    eff_length,
    subtree_all_zero,
    validate_args as _validate_args,
)
from ..validation import (
    validate_and_set_outgroup,
    prune_outgroup,
    resolve_polytomies,
)

logger = logging.getLogger(__name__)


def validate_args(pc) -> None:
    _validate_args(pc)


def prepare_tree(pc) -> None:
    pc.tree, pc.outgroup = validate_and_set_outgroup(
        pc.tree,
        pc.outgroup,
        root_taxon=getattr(pc, "root_taxon", None),
    )
    if not getattr(pc, "optimize_polytomies", True):
        resolve_polytomies(pc.tree)

    # Log polytomy summary
    poly_degrees = [
        len(n.clades)
        for n in pc.tree.find_clades()
        if not n.is_terminal() and len(n.clades) > 2
    ]
    if poly_degrees:
        degree_counts = Counter(poly_degrees)
        summary = ", ".join(
            f"{count} of degree {deg}" for deg, count in sorted(degree_counts.items())
        )
        logger.info(
            "Tree contains %d polytom%s: %s",
            len(poly_degrees),
            "y" if len(poly_degrees) == 1 else "ies",
            summary,
        )
    else:
        logger.debug("Tree is fully bifurcating (no polytomies).")

    # Single postorder pass: build terminal lists and counts together.
    # n.get_terminals() per node would re-walk each subtree (O(N · depth));
    # extending children's lists bottom-up is O(N).
    pc.name_leaves_per_node = {}
    pc.num_leaves_per_node = {}
    for n in pc.tree.find_clades(order="postorder"):
        if n.is_terminal():
            pc.name_leaves_per_node[n] = [n]
            pc.num_leaves_per_node[n] = 1
        else:
            terms: list = []
            cnt = 0
            for child in n.clades:
                terms.extend(pc.name_leaves_per_node[child])
                cnt += pc.num_leaves_per_node[child]
            pc.name_leaves_per_node[n] = terms
            pc.num_leaves_per_node[n] = cnt
    root_cnt = pc.num_leaves_per_node[pc.tree.root]
    pc.num_terminals = root_cnt - 1 if pc.outgroup else root_cnt
    pc.max_k = (
        ceil(pc.num_terminals * pc.max_k_limit)
        if pc.k is None and pc.max_k is None
        else pc.max_k if pc.k is None else None
    )

    pc._tree_wo_outgroup = None
    if pc.outgroup:
        # Newick round-trip is much faster than copy.deepcopy on large trees
        # (deepcopy walks every Bio.Phylo node and cross-reference). Falls
        # back to deepcopy if anything in the tree fails to serialize.
        try:
            from io import StringIO
            from Bio import Phylo as _Phylo

            buf = StringIO()
            _Phylo.write(pc.tree, buf, "newick")
            buf.seek(0)
            pc._tree_wo_outgroup = _Phylo.read(buf, "newick")
        except Exception:
            import copy
            pc._tree_wo_outgroup = copy.deepcopy(pc.tree)
        pc.name_leaves_per_node, pc.num_leaves_per_node = prune_outgroup(
            pc._tree_wo_outgroup, pc.outgroup
        )


def compute_dp_table(pc) -> None:
    tree = pc._tree_wo_outgroup if pc.outgroup else pc.tree
    nodes = list(tree.find_clades(order="postorder"))
    pc.postorder_nodes = nodes
    num_nodes = len(nodes)
    pc.node_to_id = {node: i for i, node in enumerate(nodes)}

    pc.dp_table = [None] * num_nodes  # penalized objective
    pc.raw_dp_table = [None] * num_nodes  # raw WSS
    pc.backptr = [None] * num_nodes
    pc.polytomy_backptr = [None] * num_nodes
    pc.cluster_cost = {}
    pc._root_ties = {}

    dtype = np.float64 if pc.num_terminals > 80_000 else np.float32

    max_states_global = pc.max_k if pc.max_k is not None else pc.num_terminals
    if max_states_global < 1:
        raise InvalidKError("max_k (or implied max_states_global) must be ≥ 1.")

    bp_dtype = np.int16 if max_states_global <= 32767 else np.int32

    min_cluster_size = getattr(pc, "min_cluster_size", 1)
    if min_cluster_size < 1:
        raise ConfigurationError("min_cluster_size must be ≥ 1.")

    outlier_thresh = pc.outlier.size_threshold
    use_outlier = outlier_thresh is not None
    prefer_fewer = pc.outlier.prefer_fewer
    preserve_dp_tables = bool(
        getattr(pc, "preserve_dp_tables", False) or logger.isEnabledFor(logging.DEBUG)
    )

    if use_outlier:
        pc._n_small = [None] * num_nodes

    for node in nodes:
        node_id = pc.node_to_id[node]
        n_leaves = pc.num_leaves_per_node[node]
        n_states = min(n_leaves, max_states_global)

        total_array = np.full(n_states + 1, np.inf, dtype=dtype)
        raw_array = np.full(n_states + 1, np.inf, dtype=dtype)
        backptr_array = np.full((2, n_states + 1), -1, dtype=bp_dtype)

        if node.is_terminal():
            pc.cluster_cost[node] = 0.0

            if n_leaves >= min_cluster_size:
                raw_array[0] = 0.0
            else:
                raw_array[0] = np.inf

            total_array[0] = raw_array[0]
            if use_outlier:
                ns = np.zeros(n_states + 1, dtype=np.int32)
                ns[0] = 1 if n_leaves < outlier_thresh else 0
                pc._n_small[node_id] = ns

            pc.raw_dp_table[node_id] = raw_array
            pc.dp_table[node_id] = total_array
            pc.backptr[node_id] = backptr_array
            continue

        if len(node.clades) > 2 and getattr(pc, "optimize_polytomies", True):
            total_array, raw_array, poly_info = compute_polytomy_dp(
                node,
                pc,
                max_states_global,
                min_cluster_size,
                outlier_thresh,
                dtype,
                bp_dtype,
            )
            pc.dp_table[node_id] = total_array
            pc.raw_dp_table[node_id] = raw_array
            pc.polytomy_backptr[node_id] = poly_info

            if not preserve_dp_tables:
                for child in node.clades:
                    cid = pc.node_to_id[child]
                    pc.dp_table[cid] = None
                    pc.raw_dp_table[cid] = None
                    if use_outlier:
                        pc._n_small[cid] = None
            continue

        left, right = node.clades[0], node.clades[1]
        left_id = pc.node_to_id[left]
        right_id = pc.node_to_id[right]

        left_total = pc.dp_table[left_id]
        right_total = pc.dp_table[right_id]
        left_raw = pc.raw_dp_table[left_id]
        right_raw = pc.raw_dp_table[right_id]

        if left_total is None or right_total is None:
            raise MissingDPTableError(
                "Child DP table missing – compute_dp_table order bug."
            )

        n_left = pc.num_leaves_per_node[left]
        n_right = pc.num_leaves_per_node[right]

        len_left = eff_length(left, pc)
        len_right = eff_length(right, pc)

        raw_one_cluster = (
            left_raw[0] + right_raw[0] + n_left * len_left + n_right * len_right
        )

        if getattr(pc, "use_branch_support", False):
            raw_support = node.confidence if node.confidence is not None else 100.0
            support = max(raw_support / 100.0, pc.min_support)
            raw_one_cluster /= support

        pc.cluster_cost[node] = float(raw_one_cluster)

        if n_leaves >= min_cluster_size:
            raw_array[0] = raw_one_cluster
        else:
            raw_array[0] = np.inf

        backptr_array[0, 0] = 0
        backptr_array[1, 0] = 0

        total_array[0] = raw_array[0]
        if use_outlier:
            left_ns = pc._n_small[left_id]
            right_ns = pc._n_small[right_id]
            ns_array = np.zeros(n_states + 1, dtype=np.int32)
            ns_array[0] = 1 if n_leaves < outlier_thresh else 0

        zero_eps = getattr(pc, "zero_length_eps", 1e-12)
        no_split_zero = getattr(pc, "no_split_zero_length", False)

        if no_split_zero and subtree_all_zero(node, pc, eps=zero_eps):
            # Leave k=0 as the only feasible state; all split states stay inf
            pc.raw_dp_table[node_id] = raw_array
            pc.dp_table[node_id] = total_array
            pc.backptr[node_id] = backptr_array

            if use_outlier:
                pc._n_small[node_id] = ns_array
                if not preserve_dp_tables:
                    pc._n_small[left_id] = None
                    pc._n_small[right_id] = None

            if not preserve_dp_tables:
                pc.dp_table[left_id] = None
                pc.dp_table[right_id] = None
                pc.raw_dp_table[left_id] = None
                pc.raw_dp_table[right_id] = None
            continue

        # Precompute the two slices we need for merging left & right children.
        # Hoisting these out of the k-loop avoids ``np.arange`` per iteration.
        left_len = len(left_raw)
        right_len = len(right_raw)

        for k in range(1, n_states + 1):
            max_i = min(k - 1, left_len - 1)
            min_i = max(0, k - 1 - (right_len - 1))

            if min_i > max_i:
                continue

            # Fast path: no outlier handling. This is the overwhelmingly
            # common case and used to dominate DP time via ``np.isclose``,
            # ``np.flatnonzero``, and ``np.any`` on every k-iteration at
            # every binary node. Here we do a single fused add + argmin,
            # which is ~4-6x faster per iteration in pure NumPy and gives
            # identical results (``np.argmin`` is tie-deterministic and
            # returns the first occurrence of the minimum, matching the
            # previous ``tied[0]`` selection).
            if not use_outlier:
                # Slice left_raw forward and right_raw reversed so that
                # element-wise addition yields the same pairing as
                # ``left_raw[i] + right_raw[k-1-i]`` for i in [min_i, max_i].
                left_slice = left_raw[min_i : max_i + 1]
                right_slice = right_raw[k - 1 - max_i : k - min_i][::-1]
                raw_scores = left_slice + right_slice

                best_local = int(raw_scores.argmin())
                best_value = raw_scores[best_local]
                if not np.isfinite(best_value):
                    continue

                raw_array[k] = best_value
                total_array[k] = best_value
                backptr_array[0, k] = min_i + best_local
                backptr_array[1, k] = k - 1 - (min_i + best_local)
                continue

            # Outlier-aware path (unchanged semantics).
            i_vals = np.arange(min_i, max_i + 1)
            j_vals = k - 1 - i_vals

            raw_scores = left_raw[i_vals] + right_raw[j_vals]

            finite_mask = np.isfinite(raw_scores)
            if not np.any(finite_mask):
                continue
            finite_idx = np.flatnonzero(finite_mask)

            n_sm = left_ns[i_vals] + right_ns[j_vals]

            if prefer_fewer:
                # Lexicographic: minimise outlier count, then raw cost
                n_sm_f = n_sm[finite_idx]
                min_n = int(np.min(n_sm_f))
                candidates = finite_idx[n_sm_f == min_n]
                best = int(candidates[np.argmin(raw_scores[candidates])])
            else:
                # Minimise raw cost, break ties by fewer outliers
                raw_f = raw_scores[finite_idx]
                min_raw = float(np.min(raw_f))
                tied = finite_idx[np.abs(raw_f - min_raw) <= 1e-12]
                if len(tied) > 1:
                    n_sm_t = n_sm[tied]
                    min_n = int(np.min(n_sm_t))
                    still_tied = tied[n_sm_t == min_n]
                    if len(still_tied) > 1 and node is tree.root:
                        pc._root_ties[k + 1] = {
                            "n_solutions": int(len(still_tied)),
                            "total_score": float(min_raw),
                            "outlier_count": int(min_n),
                        }
                    best = int(still_tied[0])
                else:
                    best = int(tied[0])

            raw_array[k] = raw_scores[best]
            total_array[k] = raw_scores[best]
            ns_array[k] = n_sm[best]

            backptr_array[0, k] = i_vals[best]
            backptr_array[1, k] = j_vals[best]

        pc.raw_dp_table[node_id] = raw_array
        pc.dp_table[node_id] = total_array
        pc.backptr[node_id] = backptr_array

        if use_outlier:
            pc._n_small[node_id] = ns_array
            if not preserve_dp_tables:
                pc._n_small[left_id] = None
                pc._n_small[right_id] = None

        if not preserve_dp_tables:
            pc.dp_table[left_id] = None
            pc.dp_table[right_id] = None
            pc.raw_dp_table[left_id] = None
            pc.raw_dp_table[right_id] = None

    pc.postorder_nodes = nodes
    pc._dp_ready = True


def _assign_group(
    group_nodes: list,
    clusters: dict,
    cluster_id: int,
    pc=None,
) -> int:
    """Assign all terminals under each node in group_nodes to cluster_id. Returns the next cluster_id."""
    name_leaves = getattr(pc, "name_leaves_per_node", None) if pc is not None else None
    for group_node in group_nodes:
        # Prefer the pre-cached terminal list (dict lookup) over Biopython's
        # on-the-fly DFS; both are equivalent but the cached path is ~100x
        # faster per call and ``backtrack`` is hot.
        terms = name_leaves.get(group_node) if name_leaves is not None else None
        if terms is None:
            terms = group_node.get_terminals()
        for terminal in terms:
            clusters[terminal] = cluster_id
    return cluster_id + 1


def _expand_prefix_hard(
    pc,
    children: list,
    step_backptrs: list,
    num_children: int,
    state_index: int,
) -> list:
    """
    Iteratively unroll the left-recursive hard-polytomy backtrack chain.

    Returns a list of ("node", node_id, state) tasks in child order.
    """
    pairs: list[tuple[int, int]] = []
    current = state_index
    for i in range(num_children - 1, 0, -1):
        step_bp = step_backptrs[i - 1]
        left = int(step_bp[current])
        if left < 0:
            raise MissingDPTableError(
                "Hard-polytomy back-pointer missing - fatal error."
            )
        right = current - 1 - left
        if right < 0:
            raise MissingDPTableError(
                "Hard-polytomy back-pointer inconsistent - fatal error."
            )
        pairs.append((i, right))
        current = left
    pairs.append((0, current))
    pairs.reverse()
    return [("node", pc.node_to_id[children[i]], s) for i, s in pairs]


def _expand_soft_tasks(
    pc,
    children: list,
    parent_map: dict,
    fullmask: int,
    target_q: int,
) -> list:
    """
    Iteratively walk the soft-polytomy parent chain and reconstruct the task list.

    Returns tasks in forward (root-to-leaves) order.
    """
    chain: list[tuple[int, Any]] = []
    mask, q = fullmask, target_q
    while mask != 0 or q != 0:
        key = (mask, q)
        if key not in parent_map:
            raise MissingDPTableError(
                "Soft-polytomy back-pointer missing - fatal error."
            )
        prev_mask, prev_q, block_mask, child_state = parent_map[key]
        chain.append((block_mask, child_state))
        mask, q = prev_mask, prev_q

    tasks = []
    for block_mask, child_state in reversed(chain):
        bits = [i for i in range(len(children)) if block_mask & (1 << i)]
        if child_state is None:
            tasks.append(("group", [children[i] for i in bits], None))
        else:
            if len(bits) != 1:
                raise MissingDPTableError(
                    "Soft-polytomy singleton back-pointer malformed."
                )
            child = children[bits[0]]
            tasks.append(("node", pc.node_to_id[child], int(child_state)))
    return tasks


def _expand_polytomy_tasks(
    node: Any, pc, c_index: int
) -> list[tuple[str, Any, int | None]]:
    """Return the list of backtrack tasks for a polytomous node at state c_index."""
    info = pc.polytomy_backptr[pc.node_to_id[node]]
    if info is None:
        raise MissingDPTableError("Polytomy back-pointer missing - fatal error.")

    children = list(node.clades)
    mode = info["mode"]

    if mode == "hard":
        return _expand_prefix_hard(pc, children, info["steps"], len(children), c_index)

    if mode == "soft":
        fullmask = (1 << len(children)) - 1
        return _expand_soft_tasks(pc, children, info["parent"], fullmask, c_index + 1)

    raise MissingDPTableError(f"Unknown polytomy back-pointer mode: {mode}")


def backtrack(pc, k: int, *, verbose: bool = False) -> dict[Any, int]:
    if k is None:
        raise InvalidKError("value of k is missing.")
    if k <= 0:
        raise InvalidKError("k must be a positive integer.")
    if not getattr(pc, "_dp_ready", False):
        raise MissingDPTableError(
            "DP table not computed. Call compute_dp_table(pc) first."
        )

    active_tree = pc._tree_wo_outgroup if pc.outgroup else pc.tree
    root = active_tree.root
    root_id = pc.node_to_id[root]

    root_dp = pc.dp_table[root_id]
    if root_dp is None:
        raise MissingDPTableError(
            "DP table at root is missing. Did compute_dp_table(pc) fail?"
        )

    cluster_index = k - 1

    if hasattr(pc, "_root_ties") and k in pc._root_ties:
        info = pc._root_ties[k]
        logger.warning(
            "[phytclust] selected k=%d has %d tied optimal root solutions "
            "(total score=%g, outlier count=%d); using deterministic fallback",
            k,
            info["n_solutions"],
            info["total_score"],
            info["outlier_count"],
        )

    if cluster_index >= len(root_dp):
        max_k_val = getattr(pc, "max_k", None)
        if max_k_val is not None and k > max_k_val:
            raise InvalidKError(
                f"Cannot partition into {k} clusters — max_k is set to {max_k_val}. "
                f"Increase or remove max_k to allow k={k}."
            )
        raise InvalidKError(
            f"Cannot partition into {k} clusters — exceeds tree capacity "
            f"(tree has {pc.num_terminals} terminals)."
        )

    # Warn if solution is infeasible
    min_cluster_size = getattr(pc, "min_cluster_size", 1)
    no_split_zero = getattr(pc, "no_split_zero_length", False)
    if np.isinf(root_dp[cluster_index]):
        reasons = []
        if min_cluster_size > 1:
            reasons.append(f"min_cluster_size={min_cluster_size}")
        if no_split_zero:
            reasons.append("no_split_zero_length=True")
        if not reasons:
            reasons.append("tree structure")
        raise InvalidClusteringError(
            f"Partition into {k} clusters is infeasible given {', '.join(reasons)}."
        )
    clusters: dict[Any, int] = {}
    current_cluster_id = 0

    stack: list[tuple[str, Any, int | None]] = [("node", root_id, cluster_index)]

    while stack:
        item_type, payload, c_index = stack.pop()

        if item_type == "group":
            current_cluster_id = _assign_group(
                payload, clusters, current_cluster_id, pc=pc
            )
            continue

        node_id = payload
        node = pc.postorder_nodes[node_id]
        if c_index is None:
            raise MissingDPTableError("Node task missing c_index during backtrack.")

        if verbose:
            print(f"Visiting node {getattr(node, 'name', '')} with c_index={c_index}")

        if c_index == 0:
            # All leaves in this clade form one cluster. Use the pre-cached
            # terminal list instead of calling ``node.get_terminals()`` — the
            # latter is a full DFS via Biopython and dominates backtrack time
            # on large trees.
            cached_terms = pc.name_leaves_per_node.get(node)
            if cached_terms is None:
                cached_terms = node.get_terminals()
            for t in cached_terms:
                clusters[t] = current_cluster_id
            current_cluster_id += 1
        elif pc.polytomy_backptr[node_id] is not None:
            # Polytomous node: reconstruct child tasks while preserving merged groups.
            tasks = _expand_polytomy_tasks(node, pc, c_index)
            for task in reversed(tasks):
                stack.append(task)
        else:
            # Binary node: unchanged
            left_k = pc.backptr[node_id][0, c_index]
            right_k = pc.backptr[node_id][1, c_index]

            if left_k < 0 or right_k < 0:
                raise MissingDPTableError(
                    "Back-pointer missing - fatal error. Check DP."
                )

            left, right = node.clades[0], node.clades[1]
            stack.append(("node", pc.node_to_id[right], int(right_k)))
            stack.append(("node", pc.node_to_id[left], int(left_k)))

    if current_cluster_id != k:
        raise InvalidClusteringError(
            f"Number of clusters found: {current_cluster_id}, expected: {k}"
        )

    return clusters


def cluster_map(pc, k: int) -> dict[Any, int]:
    """Convenience wrapper — delegates to ``pc.get_clusters(k)``."""
    return pc.get_clusters(k)
