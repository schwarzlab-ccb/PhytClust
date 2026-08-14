from __future__ import annotations

import statistics
from collections import defaultdict
from typing import Any, Optional, List
import numpy as np
from Bio.Phylo.BaseTree import Tree, Clade

from ..utils.traversal import iter_clades


# Colless index (binary only)
def colless_index_calc(tree: Tree) -> int:
    """
    Colless index for *binary* trees: sum over internal nodes of
    |#leaves(left) - #leaves(right)|.
    """
    # O(n) postorder count; get_terminals() per node would re-walk subtrees.
    leaf_count: dict[Any, int] = {}
    colless_sum = 0
    for node in iter_clades(tree.root, "postorder"):
        if node.is_terminal():
            leaf_count[node] = 1
            continue
        leaf_count[node] = sum(leaf_count[c] for c in node.clades)
        if len(node.clades) == 2:
            colless_sum += abs(
                leaf_count[node.clades[0]] - leaf_count[node.clades[1]]
            )
    return colless_sum


def normalized_colless(tree: Tree) -> float:
    """
    Normalized Colless (simple normalization). For non-binary trees,
    we skip non-binary nodes as above.
    """
    colless_sum = colless_index_calc(tree)
    n = tree.count_terminals()
    if n <= 2:
        return 0.0
    return (2.0 * colless_sum) / ((n - 1) * (n - 2))


# Stemmy / Tippy ratio
def calculate_internal_terminal_ratio(tree: Tree) -> float:
    """
    Ratio of internal to terminal branch lengths.
    Treats missing branch lengths as 0.
    """
    internal_len = 0.0
    terminal_len = 0.0
    for node in tree.find_clades():
        bl = node.branch_length or 0.0
        if node.is_terminal():
            terminal_len += bl
        else:
            internal_len += bl
    return (internal_len / terminal_len) if terminal_len != 0 else float("inf")


def calculate_int_term_ratio(tree: Tree) -> float:
    """
    Normalized internal:terminal branch length ratio (simple scaling).
    """
    ratio = calculate_internal_terminal_ratio(tree)
    n = tree.count_terminals()
    return ratio * (n / ((2 * n) - 2)) if n > 1 else float("inf")


# Branch-length variance utilities
def collect_branch_lengths(node: Clade) -> List[float]:
    lengths: List[float] = []
    for child in getattr(node, "clades", []):
        if child.branch_length is not None:
            lengths.append(float(child.branch_length))
        lengths.extend(collect_branch_lengths(child))
    return lengths


def calculate_variance_branch_length(tree: Tree) -> float:
    """
    Standard deviation of all branch lengths in the tree.
    """
    bl = collect_branch_lengths(tree.root)
    return float(np.std(bl)) if bl else 0.0


# Root-to-tip (total) distances
def total_branch_lengths(tree: Tree) -> List[float]:
    """Root-to-tip path length for each terminal."""
    return [float(tree.distance(t)) for t in tree.get_terminals()]


def calculate_total_length_variation(tree: Tree) -> float:
    """Std. dev. of root-to-tip distances."""
    vals = total_branch_lengths(tree)
    return float(np.std(vals)) if vals else 0.0


# Variance on internal vs terminal branches
def calculate_internal_variance(tree: Tree) -> float:
    vals = [
        float(node.branch_length)
        for node in tree.get_nonterminals()
        if node.branch_length is not None
    ]
    return float(np.std(vals)) if vals else 0.0


def calculate_terminal_variance(tree: Tree) -> float:
    vals = [
        float(node.branch_length)
        for node in tree.get_terminals()
        if node.branch_length is not None
    ]
    return float(np.std(vals)) if vals else 0.0


def variation_ratio(tree: Tree) -> float:
    """Internal variance / terminal variance."""
    v_int = calculate_internal_variance(tree)
    v_term = calculate_terminal_variance(tree)
    return (v_int / v_term) if v_term else float("inf")


# Terminal contribution to total length
def calculate_terminal_contributions(tree: Tree) -> float:
    """
    Percentage of total branch length contributed by terminal edges.
    """
    total = sum((node.branch_length or 0.0) for node in tree.find_clades())
    term = sum((leaf.branch_length or 0.0) for leaf in tree.get_terminals())
    return (term / total) * 100.0 if total else 0.0


# Sibling distances (for CV / variance)
def find_siblings(tree: Tree) -> List[float]:
    """
    Distances between siblings in every preterminal clade with exactly 2 terminals.
    """
    dists: List[float] = []
    for clade in tree.find_clades():
        if clade.is_preterminal():
            terminals = clade.get_terminals()
            if len(terminals) == 2:
                dists.append(float(tree.distance(terminals[0], terminals[1])))
    return dists


def calculate_variance_of_distances(distances: List[float]) -> Optional[float]:
    if distances and len(distances) > 1:
        return float(statistics.variance(distances))
    return None


def calculate_coefficient_of_variation(distances: List[float]) -> Optional[float]:
    if distances and len(distances) > 1:
        mean = statistics.mean(distances)
        if mean != 0:
            return float(statistics.stdev(distances) / mean)
    return None


# Gini coefficient on proportions
def calculate_proportions(tree: Tree, sibling_distances: List[float]) -> List[float]:
    total = sum((node.branch_length or 0.0) for node in tree.find_clades())
    return [(d / total) if total else 0.0 for d in sibling_distances]


def gini_coefficient(proportions_or_tree: List[float] | Tree) -> float:
    """
    Gini coefficient of a vector of non-negative proportions.
    Returns 0 if empty or sum is zero.
    """
    if isinstance(proportions_or_tree, Tree):
        sibling_distances = find_siblings(proportions_or_tree)
        proportions = calculate_proportions(proportions_or_tree, sibling_distances)
    else:
        proportions = list(proportions_or_tree)

    n = len(proportions)
    if n == 0:
        return 0.0
    sorted_p = sorted(proportions)
    total = sum(sorted_p)
    if total == 0:
        return 0.0
    numerator = sum((2 * i - n - 1) * x for i, x in enumerate(sorted_p, start=1))
    return float(numerator / (n * total))


def variance_ratio(tree: Tree) -> float:
    """Backward-compatible alias for internal/terminal variance ratio."""
    return variation_ratio(tree)


def colless_ratio(tree: Tree) -> float:
    """Backward-compatible normalized Colless ratio in [0, 1] for binary trees."""
    return normalized_colless(tree)


def _cluster_mrca_ids(tree: Any, cluster_terms: dict[int, list[Clade]]) -> set[int]:
    """ids of every cluster's MRCA, from one postorder pass.

    Bio.Phylo's ``common_ancestor`` re-walks the tree per call (get_path ->
    check_in_path), which costs ~O(N^2) once per cluster and dominated the
    whole run. Counting each cluster's leaves bottom-up finds the same nodes:
    postorder reaches the deepest node holding all of a cluster's leaves first,
    and that node is by definition its MRCA. Children are merged small-into-
    large, so the pass is O(N log N).
    """
    totals = {cid: len(terms) for cid, terms in cluster_terms.items()}
    leaf_cid = {id(t): cid for cid, terms in cluster_terms.items() for t in terms}

    found: dict[int, Clade] = {}
    counts: dict[int, dict[int, int]] = {}

    for node in iter_clades(tree.root, "postorder"):
        if not node.clades:
            cid = leaf_cid.get(id(node))
            sub: dict[int, int] = {} if cid is None else {cid: 1}
        else:
            sub = {}
            for child in node.clades:
                child_sub = counts.pop(id(child), {})
                if len(child_sub) > len(sub):
                    sub, child_sub = child_sub, sub
                for cid, n in child_sub.items():
                    sub[cid] = sub.get(cid, 0) + n
        counts[id(node)] = sub
        for cid, n in sub.items():
            if n == totals[cid] and cid not in found:
                found[cid] = node

    return {id(node) for node in found.values()}


def cluster_alpha(tree: Any, cmap: dict[Clade, int]) -> dict[str, Any]:
    """
    Alpha = mean extra-cluster branch length / mean intra-cluster branch length.

    Intra-cluster nodes are the strict descendants of each cluster's MRCA (the
    MRCA itself is counted as extra-cluster, along with every node outside any
    cluster's MRCA subtree). Each node contributes its ``branch_length`` once;
    ``None`` is treated as 0.
    """
    cluster_terms: dict[int, list[Clade]] = defaultdict(list)
    for term, cid in cmap.items():
        cluster_terms[cid].append(term)

    # A node is "intra" iff it is a strict descendant of some cluster's MRCA.
    # One pre-order walk propagating "an ancestor is an MRCA" downward.
    mrca_ids = _cluster_mrca_ids(tree, cluster_terms)

    intra_sum = 0.0
    intra_count = 0
    extra_sum = 0.0
    extra_count = 0
    stack: list[tuple[Any, bool]] = [(tree.root, False)]
    while stack:
        node, anc_is_mrca = stack.pop()
        bl = float(node.branch_length or 0.0)
        if anc_is_mrca:
            intra_sum += bl
            intra_count += 1
        else:
            extra_sum += bl
            extra_count += 1
        child_flag = anc_is_mrca or (id(node) in mrca_ids)
        for child in node.clades:
            stack.append((child, child_flag))

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
    """Convenience bundle of variance-style tree metrics."""
    return {
        "branch_length_variance": calculate_variance_branch_length(tree),
        "total_length_variation": calculate_total_length_variation(tree),
        "internal_variance": calculate_internal_variance(tree),
        "terminal_variance": calculate_terminal_variance(tree),
        "variance_ratio": variance_ratio(tree),
    }
