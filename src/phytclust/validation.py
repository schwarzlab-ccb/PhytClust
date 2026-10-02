from __future__ import annotations

import logging
import math
import string
from collections import deque
from typing import Any

from Bio.Phylo.BaseTree import Clade, Tree

from .exceptions import InvalidTreeError
from .utils.traversal import LeafSpans, iter_clades, leaf_maps, nonterminals, terminals

logger = logging.getLogger(__name__)


def _named_clade(tree: Tree, name: str) -> Clade:
    """Resolve a literal, unique node label; reject missing or ambiguous labels."""
    matches = [node for node in iter_clades(tree.root) if node.name == name]
    if not matches:
        raise InvalidTreeError(f"Node '{name}' not found in the tree.")
    if len(matches) != 1:
        raise InvalidTreeError(f"Node name '{name}' is ambiguous ({len(matches)} matches).")
    return matches[0]


def root_tree_at_taxon(tree: Tree, root_taxon: str | None) -> Tree:
    """
    Root tree at a specified taxon or use midpoint rooting if root_taxon is "midpoint".

    Parameters
    ----------
    tree : Tree
        The phylogenetic tree (modified in-place).
    root_taxon : str, optional
        Name of taxon to root on. If "midpoint", use midpoint rooting.
        If None, tree is used as-is (assumed already rooted).

    Returns
    -------
    tree : Tree
        The rooted tree.
    """
    if root_taxon is None:
        return tree

    if root_taxon == "midpoint":
        if len(list(terminals(tree.root))) < 2:
            raise InvalidTreeError("Midpoint rooting requires at least two leaves.")
        ensure_branch_lengths(tree)
        for node in iter_clades(tree.root):
            if node.branch_length is None:
                node.branch_length = 0.0
        tree.root_at_midpoint()
        return tree

    root_clade = _named_clade(tree, root_taxon)
    tree.root_with_outgroup(root_clade)
    return tree


def validate_and_set_outgroup(
    tree: Tree,
    outgroup: str | None,
    root_taxon: str | None = None,
) -> tuple[Tree, str | None]:
    """
    Optionally root the tree, validate presence of outgroup (if provided),
    collapse unprotected single-child nodes, and normalize names/branch lengths.
    The tree is modified in place. Requested node labels must be unique.
    Returns (tree, outgroup).
    """
    if outgroup is not None:
        _named_clade(tree, outgroup)
    if root_taxon != "midpoint":
        ensure_branch_lengths(tree)
    tree = root_tree_at_taxon(tree, root_taxon)
    merge_single_child_clades(tree, outgroup=outgroup)
    rename_nodes(tree, outgroup)
    return tree, outgroup


def prune_outgroup(
    tree: Tree, outgroup: str | None
) -> tuple[LeafSpans, dict[Any, int]]:
    """
    Remove the outgroup from ``tree`` and return the leaf maps that remain.
    Tree is modified in-place.
    If the outgroup is at root with two children, we keep the sibling as the new root.
    Otherwise remove the outgroup subtree and collapse its parent if unary,
    preserving distances between remaining leaves.
    """
    if outgroup is None:
        return leaf_maps(tree.root)

    outgroup_clade = _named_clade(tree, outgroup)

    if outgroup_clade is tree.root:
        raise InvalidTreeError("The entire tree cannot be used as an outgroup.")

    if tree.root and len(tree.root.clades) == 2 and outgroup_clade in tree.root.clades:
        sibling = (
            tree.root.clades[0]
            if tree.root.clades[1] is outgroup_clade
            else tree.root.clades[1]
        )
        tree.root = sibling
    else:
        path = tree.get_path(outgroup_clade)
        parent = tree.root if len(path) == 1 else path[-2]
        parent.clades.remove(outgroup_clade)
        if len(parent.clades) == 1:
            child = parent.clades[0]
            if parent is tree.root:
                tree.root = child
            else:
                child.branch_length = (child.branch_length or 0.0) + (
                    parent.branch_length or 0.0
                )
                grandparent = tree.root if len(path) == 2 else path[-3]
                grandparent.clades[grandparent.clades.index(parent)] = child

    return leaf_maps(tree.root)


def is_outgroup_valid(tree: Tree, outgroup: str) -> bool:
    """True if any clade in the tree has name == outgroup."""
    return any(node.name == outgroup for node in iter_clades(tree.root))


def rename_nodes(tree: Tree, outgroup: str | None = None) -> None:
    """
    Give unique names to all nodes
    Keep the outgroup name as-is if unique

    Unnamed nodes receive generated ``internal_node_*`` labels, Nodes carrying a name that another node also uses are suffixed ``_1``, ``_2``, ... with a warning.
    """
    ordered = list(nonterminals(tree.root)) + list(terminals(tree.root))

    existing_names = {
        name for name in (getattr(n, "name", None) for n in ordered) if name
    }
    outgroup_count = (
        sum(1 for n in ordered if getattr(n, "name", None) == outgroup)
        if outgroup
        else 0
    )

    node_names = set([outgroup]) if outgroup else set()
    internal_node_counter = 0

    for node in ordered:
        name = getattr(node, "name", None)

        if not name or (outgroup and name == outgroup and outgroup_count > 1):
            while True:
                new_name = "internal_node_" + (
                    string.ascii_uppercase[internal_node_counter % 26]
                    + str(internal_node_counter // 26)
                )
                internal_node_counter += 1
                if new_name not in node_names and new_name not in existing_names:
                    node.name = new_name
                    break

        elif (not outgroup or name != outgroup) and name in node_names:
            suffix = 1
            base = name
            new_name = f"{base}_{suffix}"
            while new_name in node_names or new_name in existing_names:
                suffix += 1
                new_name = f"{base}_{suffix}"
            logger.warning(
                "Node name '%s' is duplicated; renaming to '%s'. "
                "Output rows will reference the renamed label.",
                base,
                new_name,
            )
            node.name = new_name

        node_names.add(node.name)


def merge_single_child_clades(tree: Tree, outgroup: str | None = None) -> None:
    """
    Collapse chains of single-child clades by summing branch lengths.
    Preserve nodes carrying the requested outgroup label and their parents.
    The surviving node retains its original confidence and other metadata.
    Warn when a different, non-missing child confidence is discarded.
    """
    queue: deque[Clade] = deque([tree.root])
    while queue:
        clade = queue.popleft()
        while len(clade.clades) == 1:
            child = clade.clades[0]
            if outgroup is not None and outgroup in (clade.name, child.name):
                break
            if child.confidence is not None and child.confidence != clade.confidence:
                logger.warning(
                    "Collapsing single-child clade %r into parent %r: "
                    "retaining parent confidence %r; discarding child confidence %r.",
                    child.name,
                    clade.name,
                    clade.confidence,
                    child.confidence,
                )
            clade.name = getattr(child, "name", clade.name)
            clade.branch_length = (clade.branch_length or 0.0) + (
                child.branch_length or 0.0
            )
            clade.clades = child.clades
        queue.extend(clade.clades)


def ensure_branch_lengths(tree: Tree) -> None:
    """
    Normalize branch lengths before clustering.

    Non-numeric and non-finite lengths (including the root length) are rejected.

    - Negative lengths are clamped to 0.0 with a warning; a negative branch
      has no meaning for the clustering objective, and 0.0 is the neutral
      value the DP already handles.
    - If no branch has a positive length (all None/0/negative), every
      non-root branch is set to 1.0 (unweighted tree).
    - If lengths are mixed (some present, some missing/zero/clamped), they
      are left as-is, but the non-positive branches are treated as length
      0.0 by the DP, which may yield zero within-cluster dispersion and
      undefined (inf) alpha.
    """
    nodes = list(iter_clades(tree.root))
    for node in nodes:
        if node.branch_length is not None:
            try:
                finite = math.isfinite(node.branch_length)
            except (TypeError, ValueError, OverflowError):
                finite = False
            if not finite:
                raise InvalidTreeError(
                    f"Branch length for node {node.name!r} must be a finite number."
                )
    clades = [cl for cl in nodes if cl is not tree.root]

    n_missing = 0
    n_negative = 0
    for cl in clades:
        length = cl.branch_length or 0.0
        if length < 0.0:
            n_negative += 1
            cl.branch_length = 0.0
        if not length > 0.0:
            n_missing += 1

    if n_negative:
        logger.warning(
            "%d of %d branches have a negative length; they are clamped to "
            "0.0, which can produce zero within-cluster dispersion and "
            "undefined (inf) alpha for the affected clusters.",
            n_negative,
            len(clades),
        )

    if n_missing == 0:
        return

    if n_missing == len(clades):
        logger.warning(
            "Tree has no branch lengths. "
            "PhytClust will assume all branches have length 1.0."
        )
        for cl in clades:
            cl.branch_length = 1.0
        return

    logger.warning(
        "%d of %d branches have no positive length; they are treated as "
        "length 0.0, which can produce zero within-cluster dispersion and "
        "undefined (inf) alpha for the affected clusters. Provide complete "
        "branch lengths to avoid this.",
        n_missing,
        len(clades),
    )
