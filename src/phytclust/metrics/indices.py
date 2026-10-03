from __future__ import annotations

import statistics
from typing import Any, Optional, List
import numpy as np
from Bio.Phylo.BaseTree import Tree, Clade

from ..utils.traversal import iter_clades, terminals, nonterminals
from ..exceptions import DataError


def colless_index_calc(tree: Tree) -> int:
    """Return the sum of child leaf-count differences at binary nodes.

    Skip nodes with one or more than two children. For a binary tree, this
    is the Colless index; for other trees it is a partial binary-node sum."""
    leaf_count: dict[Any, int] = {}
    colless_sum = 0
    for node in iter_clades(tree.root, "postorder"):
        if node.is_terminal():
            leaf_count[node] = 1
            continue
        leaf_count[node] = sum(leaf_count[c] for c in node.clades)
        if len(node.clades) == 2:
            colless_sum += abs(leaf_count[node.clades[0]] - leaf_count[node.clades[1]])
    return colless_sum


def normalized_colless(tree: Tree) -> float:
    """Divide the binary-node Colless sum by the maximum for a binary tree.

    Return zero for at most two leaves. Skip non-binary nodes; this is not
    a generalized Colless index for multifurcations."""
    colless_sum = colless_index_calc(tree)
    n = sum(1 for _ in terminals(tree.root))
    if n <= 2:
        return 0.0
    return (2.0 * colless_sum) / ((n - 1) * (n - 2))


def calculate_internal_terminal_ratio(tree: Tree) -> float:
    """Return summed internal edge lengths divided by summed terminal lengths.

    Exclude the root length and treat missing lengths as zero. Return
    infinity when terminal lengths sum to zero."""
    internal_len = 0.0
    terminal_len = 0.0
    for node in iter_clades(tree.root):
        if node is tree.root:
            continue
        bl = node.branch_length or 0.0
        if node.is_terminal():
            terminal_len += bl
        else:
            internal_len += bl
    return (internal_len / terminal_len) if terminal_len != 0 else float("inf")


def calculate_int_term_ratio(tree: Tree) -> float:
    """
    Multiply the internal-to-terminal length ratio by n / (2n - 2).

    Here n is the number of leaves. Return infinity for a single leaf.
    """
    ratio = calculate_internal_terminal_ratio(tree)
    n = sum(1 for _ in terminals(tree.root))
    return ratio * (n / ((2 * n) - 2)) if n > 1 else float("inf")


def collect_branch_lengths(node: Clade) -> List[float]:
    """Return present branch lengths below the node, excluding its own length."""
    return [
        float(descendant.branch_length)
        for descendant in iter_clades(node)
        if descendant is not node and descendant.branch_length is not None
    ]


def branch_length_standard_deviation(tree: Tree) -> float:
    """Return the population standard deviation of non-root branch lengths.

    Ignore missing lengths. Return zero when no lengths are available."""
    bl = collect_branch_lengths(tree.root)
    return float(np.std(bl)) if bl else 0.0


def root_to_tip_lengths(tree: Tree) -> List[float]:
    """Return root-to-tip lengths in leaf order, treating missing lengths as zero."""
    lengths = []
    stack = [(tree.root, 0.0)]
    while stack:
        node, distance = stack.pop()
        if node.is_terminal():
            lengths.append(distance)
        else:
            stack.extend(
                (child, distance + float(child.branch_length or 0.0))
                for child in reversed(node.clades)
            )
    return lengths


def root_to_tip_standard_deviation(tree: Tree) -> float:
    """Return the population standard deviation of root-to-tip lengths."""
    vals = root_to_tip_lengths(tree)
    return float(np.std(vals)) if vals else 0.0


def internal_branch_standard_deviation(tree: Tree) -> float:
    """Return the population standard deviation of present internal edge lengths.

    Exclude the root length. Return zero when no lengths are available."""
    vals = [
        float(node.branch_length)
        for node in nonterminals(tree.root)
        if node is not tree.root and node.branch_length is not None
    ]
    return float(np.std(vals)) if vals else 0.0


def terminal_branch_standard_deviation(tree: Tree) -> float:
    """Return the population standard deviation of present terminal edge lengths.

    Exclude the root if it is a leaf. Return zero when no lengths are available."""
    vals = [
        float(node.branch_length)
        for node in terminals(tree.root)
        if node is not tree.root and node.branch_length is not None
    ]
    return float(np.std(vals)) if vals else 0.0


def internal_to_terminal_std_ratio(tree: Tree) -> float:
    """Return internal divided by terminal branch-length standard deviation.

    Return infinity when the terminal standard deviation is zero."""
    v_int = internal_branch_standard_deviation(tree)
    v_term = terminal_branch_standard_deviation(tree)
    return (v_int / v_term) if v_term else float("inf")


def calculate_terminal_contributions(tree: Tree) -> float:
    """
    Percentage of total branch length contributed by terminal edges.
    """
    total = sum(
        (node.branch_length or 0.0)
        for node in iter_clades(tree.root)
        if node is not tree.root
    )
    term = sum(
        (leaf.branch_length or 0.0)
        for leaf in terminals(tree.root)
        if leaf is not tree.root
    )
    return (term / total) * 100.0 if total else 0.0


def sibling_leaf_distances(tree: Tree) -> List[float]:
    """Return distances for pairs of leaves sharing a binary parent."""
    distances = []
    for node in iter_clades(tree.root):
        if len(node.clades) == 2 and all(child.is_terminal() for child in node.clades):
            distances.append(
                sum(float(child.branch_length or 0.0) for child in node.clades)
            )
    return distances


def calculate_variance_of_distances(distances: List[float]) -> Optional[float]:
    """Return sample variance, or None when fewer than two distances are supplied."""
    if distances and len(distances) > 1:
        return float(statistics.variance(distances))
    return None


def calculate_coefficient_of_variation(distances: List[float]) -> Optional[float]:
    """Return sample standard deviation divided by the mean.

    Return None for fewer than two distances or a zero mean."""
    if distances and len(distances) > 1:
        mean = statistics.mean(distances)
        if mean != 0:
            return float(statistics.stdev(distances) / mean)
    return None


def calculate_proportions(tree: Tree, sibling_distances: List[float]) -> List[float]:
    """Divide sibling distances by total non-root branch length.

    Return zeros when the total length is zero."""
    total = sum(
        (node.branch_length or 0.0)
        for node in iter_clades(tree.root)
        if node is not tree.root
    )
    return [(d / total) if total else 0.0 for d in sibling_distances]


def gini_coefficient(proportions_or_tree: List[float] | Tree) -> float:
    """Return the Gini coefficient of finite, non-negative values.

    For a tree, use distances between sibling leaf pairs. Return zero for
    empty input or a zero sum."""
    if isinstance(proportions_or_tree, Tree):
        sibling_distances = sibling_leaf_distances(proportions_or_tree)
        proportions = calculate_proportions(proportions_or_tree, sibling_distances)
    else:
        proportions = list(proportions_or_tree)

    n = len(proportions)
    if n == 0:
        return 0.0
    try:
        values = np.asarray(proportions, dtype=float)
    except (TypeError, ValueError) as exc:
        raise DataError("Gini values must be finite, non-negative numbers.") from exc
    if values.ndim != 1 or not np.all(np.isfinite(values)) or np.any(values < 0):
        raise DataError("Gini values must be finite, non-negative numbers.")
    largest_value = float(values.max())
    if largest_value == 0:
        return 0.0
    sorted_p = sorted(values / largest_value)
    total = sum(sorted_p)
    if total == 0:
        return 0.0
    numerator = sum((2 * i - n - 1) * x for i, x in enumerate(sorted_p, start=1))
    return float(numerator / (n * total))


def variance_ratio(tree: Tree) -> float:
    """Compatibility name for the internal-to-terminal standard-deviation ratio."""
    return internal_to_terminal_std_ratio(tree)


def colless_ratio(tree: Tree) -> float:
    """Backward-compatible normalized Colless ratio in [0, 1] for binary trees."""
    return normalized_colless(tree)


class AlphaIndex:
    """Store leaf spans and edge lengths for repeated alpha calculations.

    The arrays describe the tree when this index was built. Rebuild after
    changes to topology, child order, or branch lengths.
    """

    def __init__(self, tree: Any):
        leaf_pos: dict[int, int] = {}
        node_idx: dict[int, int] = {}
        starts: list[int] = []
        ends: list[int] = []
        lengths: list[float] = []
        for node in iter_clades(tree.root, "postorder"):
            node_idx[id(node)] = len(starts)
            if node.clades:
                starts.append(starts[node_idx[id(node.clades[0])]])
                ends.append(ends[node_idx[id(node.clades[-1])]])
            else:
                leaf_pos[id(node)] = len(leaf_pos)
                starts.append(leaf_pos[id(node)])
                ends.append(leaf_pos[id(node)] + 1)
            lengths.append(float(node.branch_length or 0.0))

        root = node_idx[id(tree.root)]
        edges = np.array([i for i in range(len(starts)) if i != root], dtype=np.int64)
        self.leaf_pos = leaf_pos
        self.n_leaves = len(leaf_pos)
        self.start = np.asarray(starts, dtype=np.int64)[edges]
        self.end = np.asarray(ends, dtype=np.int64)[edges]
        self.length = np.asarray(lengths, dtype=np.float64)[edges]

    def labels(self, cmap: dict[Clade, int]) -> np.ndarray:
        """Cluster labels ``0..m-1`` in leaf order; unassigned leaves get their own."""
        leaf_positions = np.fromiter(
            (self.leaf_pos.get(id(t), -1) for t in cmap),
            dtype=np.int64,
            count=len(cmap),
        )
        cluster_ids = np.fromiter(cmap.values(), dtype=np.int64, count=len(cmap))
        keep = leaf_positions >= 0
        _, dense = np.unique(cluster_ids[keep], return_inverse=True)
        n_clusters = int(dense.max()) + 1 if dense.size else 0
        labels = n_clusters + np.arange(self.n_leaves, dtype=np.int64)
        labels[leaf_positions[keep]] = dense
        return labels

    def within(self, labels: np.ndarray) -> np.ndarray:
        """Per edge: True if it lies inside a cluster rather than on the backbone.

        The edge into ``v`` is inside if all of ``v``'s leaves share one
        cluster and ``v`` does not hold the whole cluster; the node that does
        is the cluster's root and its edge is the stem.
        """
        changes = np.concatenate(([0], np.cumsum(labels[1:] != labels[:-1])))
        single_cluster_subtrees = changes[self.end - 1] == changes[self.start]
        sizes = np.bincount(labels)
        return single_cluster_subtrees & (
            self.end - self.start < sizes[labels[self.start]]
        )


def cluster_alpha(
    tree: Any, cmap: dict[Clade, int], index: Optional[AlphaIndex] = None
) -> dict[str, Any]:
    """
    Alpha = mean extra-cluster branch length / mean intra-cluster branch length.

    Every edge is counted once, as the branch leading into its child node.
    Intra-cluster edges lie below a cluster's root; extra-cluster edges form
    the backbone, including each cluster's stem. A soft-polytomy group, whose
    members are several children of one polytomy, is treated as the clade it
    forms at its implicit zero-length node: edges into its members are intra,
    and the zero-length stem is not an edge of the tree. The root has no
    incoming edge and is not counted, even when it carries a
    ``branch_length`` (as the new root does after outgroup pruning).
    ``None`` is treated as 0. Pass ``index`` to reuse the per-tree arrays
    across clusterings of the same tree.
    """
    if index is None:
        index = AlphaIndex(tree)
    within = index.within(index.labels(cmap))

    intra_count = int(within.sum())
    extra_count = int(within.size - intra_count)
    intra_sum = float(index.length[within].sum())
    extra_sum = float(index.length[~within].sum())

    avg_intra = intra_sum / intra_count if intra_count else float("nan")
    avg_extra = extra_sum / extra_count if extra_count else float("nan")
    if intra_count == 0 or avg_intra == 0:
        alpha = float("inf")
    else:
        alpha = avg_extra / avg_intra

    return {
        "alpha": alpha,
        "avg_extra_branch_length": avg_extra,
        "avg_intra_branch_length": avg_intra,
        "sum_extra_branch_length": extra_sum,
        "sum_intra_branch_length": intra_sum,
        "n_extra_nodes": extra_count,
        "n_intra_nodes": intra_count,
    }


def variance_indices(tree: Tree) -> dict[str, float]:
    """Return standard-deviation metrics under their historical dictionary keys."""
    return {
        "branch_length_variance": branch_length_standard_deviation(tree),
        "total_length_variation": root_to_tip_standard_deviation(tree),
        "internal_variance": internal_branch_standard_deviation(tree),
        "terminal_variance": terminal_branch_standard_deviation(tree),
        "variance_ratio": variance_ratio(tree),
    }


# Compatibility names retain the existing standard-deviation formulas.
calculate_variance_branch_length = branch_length_standard_deviation
total_branch_lengths = root_to_tip_lengths
calculate_total_length_variation = root_to_tip_standard_deviation
calculate_internal_variance = internal_branch_standard_deviation
calculate_terminal_variance = terminal_branch_standard_deviation
variation_ratio = internal_to_terminal_std_ratio
find_siblings = sibling_leaf_distances
