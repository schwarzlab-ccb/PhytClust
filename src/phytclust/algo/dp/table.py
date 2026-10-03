import numpy as np
from math import ceil
from collections import Counter
from typing import Any
import logging
import copy
from numbers import Integral


from ...exceptions import (
    ConfigurationError,
    InvalidKError,
    MissingDPTableError,
    InvalidClusteringError,
)
from .polytomy import compute_polytomy_dp
from .merge import (
    split_offset_dtype,
    MergeSettings,
    merge_child_tables,
    decode_split_offset,
)
from ...utils.traversal import iter_clades, leaf_maps, terminals
from .costs import (
    small_cluster_penalty,
    effective_branch_length,
    small_cluster_penalty_enabled,
    subtree_has_only_zero_length_branches,
    floating_point_precision,
    branch_support_fraction,
    validate_clustering_parameters as _validate_clustering_parameters,
)
from ...validation import (
    validate_and_set_outgroup,
    prune_outgroup,
)

logger = logging.getLogger(__name__)


def validate_clustering_parameters(clustering) -> None:
    """Validate the parameters used to build clustering tables."""
    _validate_clustering_parameters(clustering)


def _copy_tree(tree):
    """Copy a tree and its metadata without recursing through its topology."""
    nodes = list(iter_clades(tree.root))
    memo = {id(node): copy.copy(node) for node in nodes}
    for node in nodes:
        clone = memo[id(node)]
        clone.__dict__ = copy.deepcopy(node.__dict__, memo)
    return copy.deepcopy(tree, memo)


def prepare_tree(clustering) -> None:
    """Normalize the tree, copy and prune its outgroup, and build leaf maps."""
    clustering.tree, clustering.outgroup = validate_and_set_outgroup(
        clustering.tree,
        clustering.outgroup,
        root_taxon=getattr(clustering, "root_taxon", None),
    )

    polytomy_degrees = [
        len(n.clades)
        for n in iter_clades(clustering.tree.root)
        if not n.is_terminal() and len(n.clades) > 2
    ]
    if polytomy_degrees:
        degree_counts = Counter(polytomy_degrees)
        summary = ", ".join(
            f"{count} of degree {deg}" for deg, count in sorted(degree_counts.items())
        )
        logger.info(
            "Tree contains %d polytom%s: %s",
            len(polytomy_degrees),
            "y" if len(polytomy_degrees) == 1 else "ies",
            summary,
        )
    else:
        logger.debug("Tree is fully bifurcating (no polytomies).")

    clustering._tree_wo_outgroup = None
    if clustering.outgroup:
        clustering._tree_wo_outgroup = _copy_tree(clustering.tree)
        clustering.name_leaves_per_node, clustering.num_leaves_per_node = (
            prune_outgroup(clustering._tree_wo_outgroup, clustering.outgroup)
        )
    else:
        clustering.name_leaves_per_node, clustering.num_leaves_per_node = leaf_maps(
            clustering.tree.root
        )

    active_root = (
        clustering._tree_wo_outgroup.root
        if clustering.outgroup
        else clustering.tree.root
    )
    clustering.num_terminals = clustering.num_leaves_per_node[active_root]
    if clustering.k is None and clustering.max_k is None:
        clustering.max_k = ceil(clustering.num_terminals * clustering.max_k_limit)


def compute_dp_table(clustering) -> None:
    """Build clustering costs and backtracking splits, processing children first."""
    tree = clustering._tree_wo_outgroup if clustering.outgroup else clustering.tree
    nodes = list(iter_clades(tree.root, "postorder"))
    clustering.postorder_nodes = nodes
    clustering.dp_children = [tuple(n.clades) for n in nodes]
    num_nodes = len(nodes)
    clustering.node_to_id = {node: i for i, node in enumerate(nodes)}

    clustering.dp_table = [None] * num_nodes
    clustering.raw_dp_table = [None] * num_nodes
    clustering.backptr = [None] * num_nodes
    clustering.polytomy_backptr = [None] * num_nodes
    clustering.cluster_cost = {}
    clustering.tied_optima = {}

    dtype = np.float32 if getattr(clustering, "dp_float32", False) else np.float64
    relative_cost_tolerance = 8.0 * floating_point_precision(dtype)

    dp_cap = getattr(clustering, "_dp_cap", None)
    if dp_cap is not None:
        global_max_cluster_count = dp_cap
    elif clustering.max_k is not None:
        global_max_cluster_count = clustering.max_k
    else:
        global_max_cluster_count = clustering.num_terminals
    if global_max_cluster_count < 1:
        raise InvalidKError("max_k (or implied global_max_cluster_count) must be ≥ 1.")
    clustering._dp_states_cap = global_max_cluster_count

    min_cluster_size = getattr(clustering, "min_cluster_size", 1)
    if min_cluster_size < 1:
        raise ConfigurationError("min_cluster_size must be ≥ 1.")

    outlier_size_threshold = clustering.outlier.size_threshold
    outlier_detection_enabled = outlier_size_threshold is not None
    prioritize_fewer_outliers = clustering.outlier.prefer_fewer
    preserve_dp_tables = bool(
        getattr(clustering, "preserve_dp_tables", False)
        or getattr(clustering, "save_tied_optima", False)
        or logger.isEnabledFor(logging.DEBUG)
    )

    penalty_enabled = small_cluster_penalty_enabled(clustering)
    if penalty_enabled and any(
        (not n.is_terminal()) and len(n.clades) > 2 for n in nodes
    ):
        raise ConfigurationError(
            "outlier penalty (penalty_enabled=True) is supported for "
            "bifurcating trees only; this tree has multifurcations. "
            "Disable the penalty to cluster it."
        )

    merge_settings = MergeSettings(
        outlier_detection_enabled,
        penalty_enabled,
        prioritize_fewer_outliers,
        relative_cost_tolerance,
        dtype,
    )

    if outlier_detection_enabled:
        clustering._outlier_counts = [None] * num_nodes

    for node in nodes:
        node_id = clustering.node_to_id[node]
        leaf_count = clustering.num_leaves_per_node[node]
        max_cluster_count = min(leaf_count, global_max_cluster_count)

        costs = np.full(max_cluster_count + 1, np.inf, dtype=dtype)
        penalized_costs = (
            np.full(max_cluster_count + 1, np.inf, dtype=dtype)
            if penalty_enabled
            else costs
        )

        if node.is_terminal():
            clustering.cluster_cost[node] = 0.0

            if leaf_count >= min_cluster_size:
                costs[0] = 0.0
            else:
                costs[0] = np.inf

            penalized_costs[0] = costs[0]
            if penalty_enabled and np.isfinite(costs[0]):
                penalized_costs[0] = costs[0] + small_cluster_penalty(
                    leaf_count, clustering
                )
            if outlier_detection_enabled:
                leaf_outlier_counts = np.zeros(max_cluster_count + 1, dtype=np.int32)
                leaf_outlier_counts[0] = 1 if leaf_count < outlier_size_threshold else 0
                clustering._outlier_counts[node_id] = leaf_outlier_counts

            clustering.raw_dp_table[node_id] = costs
            clustering.dp_table[node_id] = penalized_costs
            continue

        if len(node.clades) > 2:
            penalized_costs, costs, poly_info = compute_polytomy_dp(
                node,
                clustering,
                global_max_cluster_count,
                min_cluster_size,
                outlier_size_threshold,
                dtype,
            )
            clustering.dp_table[node_id] = penalized_costs if penalty_enabled else costs
            clustering.raw_dp_table[node_id] = costs
            clustering.polytomy_backptr[node_id] = poly_info

            if not preserve_dp_tables:
                for child in node.clades:
                    cid = clustering.node_to_id[child]
                    clustering.dp_table[cid] = None
                    clustering.raw_dp_table[cid] = None
                    if outlier_detection_enabled:
                        clustering._outlier_counts[cid] = None
            continue

        left, right = node.clades[0], node.clades[1]
        left_id = clustering.node_to_id[left]
        right_id = clustering.node_to_id[right]

        left_penalized_costs = clustering.dp_table[left_id]
        right_penalized_costs = clustering.dp_table[right_id]
        left_costs = clustering.raw_dp_table[left_id]
        right_costs = clustering.raw_dp_table[right_id]

        if left_penalized_costs is None or right_penalized_costs is None:
            raise MissingDPTableError(
                f"Child DP table is missing while processing node {node.name!r}."
            )

        split_offsets = np.full(
            max_cluster_count + 1,
            -1,
            dtype=split_offset_dtype(len(left_costs), len(right_costs)),
        )

        left_leaf_count = clustering.num_leaves_per_node[left]
        right_leaf_count = clustering.num_leaves_per_node[right]

        left_branch_length = effective_branch_length(left, clustering)
        right_branch_length = effective_branch_length(right, clustering)

        raw_one_cluster = (
            clustering.cluster_cost[left]
            + clustering.cluster_cost[right]
            + left_leaf_count * left_branch_length
            + right_leaf_count * right_branch_length
        )

        if getattr(clustering, "use_branch_support", False):
            raw_one_cluster /= branch_support_fraction(node, clustering)

        clustering.cluster_cost[node] = float(raw_one_cluster)

        if leaf_count >= min_cluster_size:
            costs[0] = raw_one_cluster
        else:
            costs[0] = np.inf

        split_offsets[0] = 0

        penalized_costs[0] = costs[0]
        if penalty_enabled and np.isfinite(costs[0]):
            penalized_costs[0] = costs[0] + small_cluster_penalty(
                leaf_count, clustering
            )
        if outlier_detection_enabled:
            left_outlier_counts = clustering._outlier_counts[left_id]
            right_outlier_counts = clustering._outlier_counts[right_id]
            outlier_counts = np.zeros(max_cluster_count + 1, dtype=np.int32)
            outlier_counts[0] = 1 if leaf_count < outlier_size_threshold else 0

        zero_eps = getattr(clustering, "zero_length_eps", 1e-12)
        no_split_zero = getattr(clustering, "no_split_zero_length", False)

        if no_split_zero and subtree_has_only_zero_length_branches(
            node, clustering, zero_length_tolerance=zero_eps
        ):
            clustering.raw_dp_table[node_id] = costs
            clustering.dp_table[node_id] = penalized_costs
            clustering.backptr[node_id] = split_offsets

            if outlier_detection_enabled:
                clustering._outlier_counts[node_id] = outlier_counts
                if not preserve_dp_tables:
                    clustering._outlier_counts[left_id] = None
                    clustering._outlier_counts[right_id] = None

            if not preserve_dp_tables:
                clustering.dp_table[left_id] = None
                clustering.dp_table[right_id] = None
                clustering.raw_dp_table[left_id] = None
                clustering.raw_dp_table[right_id] = None
            continue

        merge_child_tables(
            left_costs,
            right_costs,
            left_penalized_costs,
            right_penalized_costs,
            left_outlier_counts if outlier_detection_enabled else None,
            right_outlier_counts if outlier_detection_enabled else None,
            max_cluster_count,
            merge_settings,
            costs,
            penalized_costs,
            outlier_counts if outlier_detection_enabled else None,
            split_offsets,
        )

        clustering.raw_dp_table[node_id] = costs
        clustering.dp_table[node_id] = penalized_costs
        clustering.backptr[node_id] = split_offsets

        if outlier_detection_enabled:
            clustering._outlier_counts[node_id] = outlier_counts
            if not preserve_dp_tables:
                clustering._outlier_counts[left_id] = None
                clustering._outlier_counts[right_id] = None

        if not preserve_dp_tables:
            clustering.dp_table[left_id] = None
            clustering.dp_table[right_id] = None
            clustering.raw_dp_table[left_id] = None
            clustering.raw_dp_table[right_id] = None

    clustering.postorder_nodes = nodes
    clustering._dp_ready = True


def _assign_subtrees_to_cluster(
    group_nodes: list,
    clusters: dict,
    cluster_id: int,
    clustering=None,
) -> int:
    """Assign the leaves of each subtree to one cluster; return the next cluster ID."""
    name_leaves = (
        getattr(clustering, "name_leaves_per_node", None)
        if clustering is not None
        else None
    )
    for group_node in group_nodes:
        terms = name_leaves.get(group_node) if name_leaves is not None else None
        if terms is None:
            terms = terminals(group_node)
        for terminal in terms:
            if terminal in clusters:
                raise InvalidClusteringError(
                    f"Leaf {terminal.name!r} is assigned more than once."
                )
            clusters[terminal] = cluster_id
    return cluster_id + 1


def _state_table_length(clustering, leaf_count: int) -> int:
    """Return the table length, including the one-cluster state."""
    return min(leaf_count, clustering._dp_states_cap) + 1


def _reconstruct_hard_polytomy_tasks(
    clustering,
    children: list,
    step_backptrs: list,
    num_children: int,
    state_index: int,
) -> list:
    """
    Reconstruct child tasks from the stored hard-polytomy splits.

    Returns a list of ("node", node_id, state) tasks in child order.
    """
    pairs: list[tuple[int, int]] = []
    current = state_index
    for i in range(num_children - 1, 0, -1):
        offset = int(step_backptrs[i - 1][current])
        if offset < 0:
            raise MissingDPTableError("Hard-polytomy split is missing.")
        child_len = _state_table_length(
            clustering, clustering.num_leaves_per_node[children[i]]
        )
        left, right = decode_split_offset(offset, current, child_len)
        if right < 0:
            raise MissingDPTableError(
                "Hard-polytomy split contains an invalid child state."
            )
        pairs.append((i, right))
        current = left
    pairs.append((0, current))
    pairs.reverse()
    return [("node", clustering.node_to_id[children[i]], s) for i, s in pairs]


def _reconstruct_soft_polytomy_tasks(
    clustering,
    children: list,
    parent_map: dict,
    all_children_mask: int,
    target_cluster_count: int,
) -> list:
    """
    Reconstruct groups and child tasks from stored soft-polytomy choices.

    Returns tasks in forward (root-to-leaves) order.
    """
    chain: list[tuple[int, Any]] = []
    mask, q = all_children_mask, target_cluster_count
    while mask != 0 or q != 0:
        key = (mask, q)
        if key not in parent_map:
            raise MissingDPTableError("Soft-polytomy reconstruction choice is missing.")
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
            tasks.append(("node", clustering.node_to_id[child], int(child_state)))
    return tasks


def _reconstruct_polytomy_tasks(
    node: Any, clustering, state_index: int
) -> list[tuple[str, Any, int | None]]:
    """Return backtracking tasks for a node with more than two children."""
    node_id = clustering.node_to_id[node]
    info = clustering.polytomy_backptr[node_id]
    if info is None:
        raise MissingDPTableError("Polytomy reconstruction data is missing.")

    children = list(clustering.dp_children[node_id])
    mode = info["mode"]

    if mode == "hard":
        return _reconstruct_hard_polytomy_tasks(
            clustering, children, info["steps"], len(children), state_index
        )

    if mode == "soft":
        all_children_mask = (1 << len(children)) - 1
        return _reconstruct_soft_polytomy_tasks(
            clustering, children, info["parent"], all_children_mask, state_index + 1
        )

    raise MissingDPTableError(f"Unknown polytomy back-pointer mode: {mode}")


def backtrack(clustering, k: int, *, verbose: bool = False) -> dict[Any, int]:
    """Reconstruct a k-cluster partition from stored splits.

    Return a mapping from active leaves to cluster IDs. Reject invalid k,
    impossible partitions, missing splits, and incomplete leaf assignments.
    """
    if k is None:
        raise InvalidKError("value of k is missing.")
    if isinstance(k, (bool, np.bool_)) or not isinstance(k, Integral) or k <= 0:
        raise InvalidKError("k must be a positive integer.")
    if not getattr(clustering, "_dp_ready", False):
        raise MissingDPTableError(
            "DP table not computed. Call compute_dp_table(clustering) first."
        )

    active_tree = (
        clustering._tree_wo_outgroup if clustering.outgroup else clustering.tree
    )
    root = active_tree.root
    root_id = clustering.node_to_id[root]

    root_dp = clustering.dp_table[root_id]
    if root_dp is None:
        raise MissingDPTableError(
            "DP table at root is missing. Did compute_dp_table(clustering) fail?"
        )

    cluster_index = k - 1

    if cluster_index >= len(root_dp):
        if k > clustering.num_terminals:
            raise InvalidKError(
                f"k={k} exceeds the number of leaves ({clustering.num_terminals})."
            )
        max_k_val = getattr(clustering, "max_k", None)
        if max_k_val is not None and k > max_k_val:
            raise InvalidKError(
                f"Cannot partition into {k} clusters — max_k is set to {max_k_val}. "
                f"Increase or remove max_k to allow k={k}."
            )
        raise InvalidKError(
            f"Cannot partition into {k} clusters — exceeds tree capacity "
            f"(tree has {clustering.num_terminals} terminals)."
        )

    min_cluster_size = getattr(clustering, "min_cluster_size", 1)
    no_split_zero = getattr(clustering, "no_split_zero_length", False)
    if not np.isfinite(root_dp[cluster_index]):
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
        item_type, payload, state_index = stack.pop()

        if item_type == "group":
            current_cluster_id = _assign_subtrees_to_cluster(
                payload, clusters, current_cluster_id, clustering=clustering
            )
            continue

        node_id = payload
        node = clustering.postorder_nodes[node_id]
        if state_index is None:
            raise MissingDPTableError("Node task missing state_index during backtrack.")

        if verbose:
            print(
                f"Visiting node {getattr(node, 'name', '')} with state_index={state_index}"
            )

        if state_index == 0:
            cached_terms = clustering.name_leaves_per_node.get(node)
            if cached_terms is None:
                cached_terms = terminals(node)
            for t in cached_terms:
                if t in clusters:
                    raise InvalidClusteringError(
                        f"Leaf {t.name!r} is assigned more than once."
                    )
                clusters[t] = current_cluster_id
            current_cluster_id += 1
        elif clustering.polytomy_backptr[node_id] is not None:
            tasks = _reconstruct_polytomy_tasks(node, clustering, state_index)
            for task in reversed(tasks):
                stack.append(task)
        else:
            offset = int(clustering.backptr[node_id][state_index])
            if offset < 0:
                raise MissingDPTableError(
                    f"Split is missing for node {node.name!r} at state {state_index}."
                )

            left, right = clustering.dp_children[node_id]
            left_k, right_k = decode_split_offset(
                offset,
                state_index,
                _state_table_length(clustering, clustering.num_leaves_per_node[right]),
            )
            stack.append(("node", clustering.node_to_id[right], right_k))
            stack.append(("node", clustering.node_to_id[left], left_k))

    if current_cluster_id != k:
        raise InvalidClusteringError(
            f"Number of clusters found: {current_cluster_id}, expected: {k}"
        )

    if set(clusters) != set(terminals(root)):
        raise InvalidClusteringError(
            "The partition must contain every active leaf exactly once."
        )

    if getattr(clustering, "save_tied_optima", False):
        from .ties import enumerate_tied_optima

        if not hasattr(clustering, "tied_optima") or clustering.tied_optima is None:
            clustering.tied_optima = {}
        clustering.tied_optima[k] = enumerate_tied_optima(
            clustering, k, chosen=clusters
        )

    return clusters


def cluster_map(clustering, k: int) -> dict[Any, int]:
    """Return the partition for k clusters using the clustering instance."""
    return clustering.get_clusters(k)
