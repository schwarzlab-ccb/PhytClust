from __future__ import annotations

from typing import List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from Bio.Phylo.BaseTree import Tree, Clade

from ..exceptions import InvalidTreeError, ConfigurationError
from .traversal import iter_clades, terminals, nonterminals


def get_pairwise_distances(
    tree: Tree,
    mode: str = "terminals",
    as_dataframe: bool = False,
    mrca: Optional[str] = None,
) -> Union[np.ndarray, pd.DataFrame]:
    """
    Return a symmetric matrix of branch-length distances between selected nodes.

    ``mode`` selects leaves, internal nodes, or both (leaves first).
    If ``mrca`` is supplied, select nodes from that subtree, including its
    root when selecting internal nodes. The name must match exactly one node.
    DataFrame labels use node names, with generated labels for unnamed nodes.
    """
    if mode not in {"terminals", "nonterminals", "all"}:
        raise ConfigurationError("mode must be one of {'terminals','nonterminals','all'}")
    clade = tree.root
    if mrca is not None:
        matches = [node for node in iter_clades(tree.root) if node.name == mrca]
        if not matches:
            raise InvalidTreeError(f"No node found with name {mrca}")
        if len(matches) != 1:
            raise InvalidTreeError(f"Node name {mrca!r} is ambiguous ({len(matches)} matches).")
        clade = matches[0]
    if mode == "terminals":
        terms = list(terminals(clade))
    elif mode == "nonterminals":
        terms = list(nonterminals(clade))
    else:
        terms = list(terminals(clade)) + list(nonterminals(clade))

    n = len(terms)
    dist = np.zeros((n, n), dtype=float)

    iu, ju = np.triu_indices(n, k=1)
    for i, j in zip(iu, ju):
        bl = float(tree.distance(terms[i], terms[j]))
        dist[i, j] = bl
        dist[j, i] = bl

    if as_dataframe:
        used_names = {t.name for t in terms if t.name is not None}
        names = []
        for i, node in enumerate(terms):
            name = node.name
            if name is None:
                name = f"node_{i}"
                while name in used_names:
                    name += "_"
                used_names.add(name)
            names.append(name)
        return pd.DataFrame(dist, index=names, columns=names)
    return dist


def get_parent(tree: Tree, child: Clade) -> Optional[Clade]:
    """Return the child's parent, or None for the root or a node outside the tree."""
    path = tree.get_path(child)
    if not path:
        return None
    return path[-2] if len(path) > 1 else tree.root


def count_branches_in_clusters(clusters: dict) -> int:
    """Estimate branch counts assuming each cluster is a rooted binary tree.

    ``clusters`` maps cluster IDs to lists of leaves. Count ``2N - 2``
    branches for each cluster with N > 1 leaves; smaller clusters add zero.
    The estimate assumes two children per internal node.
    """
    branch_count = 0
    for clades in clusters.values():
        N = len(clades)
        if N > 1:
            branch_count += 2 * N - 2
    return branch_count


def find_all_min_indices(arr: List[float]) -> Tuple[List[int], float]:
    """Return all indices equal to the minimum and the minimum value.

    Return ``([], inf)`` for empty input. Reject NaN values.
    """
    if len(arr) == 0:
        return [], float("inf")
    min_value = float("inf")
    min_indices: List[int] = []
    for i, value in enumerate(arr):
        if np.isnan(value):
            raise ValueError("Minimum values must not contain NaN.")
        if value < min_value:
            min_value = value
            min_indices = [i]
        elif value == min_value:
            min_indices.append(i)
    return min_indices, min_value


def rename_internal_nodes(tree: Tree) -> None:
    """
    Rename internal nodes in place to ``internal_1``, ``internal_2``, etc.

    Skip labels already used by leaves.
    """
    internal_node_count = 1
    used_names = {node.name for node in terminals(tree.root)}
    for clade in nonterminals(tree.root):
        while f"internal_{internal_node_count}" in used_names:
            internal_node_count += 1
        clade.name = f"internal_{internal_node_count}"
        used_names.add(clade.name)
        internal_node_count += 1
