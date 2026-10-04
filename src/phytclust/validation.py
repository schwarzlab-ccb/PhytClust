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
    Prepare the tree for clustering and return ``(tree, outgroup)``.

    Check the outgroup name, normalize branch lengths, and root the tree
    if requested. Collapse single-child nodes except the outgroup and its
    parent, then give each node a unique name.

    Modify the tree in place. Root and outgroup names must each match one node.
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
    When replacing a parent with its remaining child, keep the parent's
    confidence and warn if a different child confidence is discarded.
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
        _retain_parent_confidence(tree.root, sibling)
        tree.root = sibling
    else:
        path = tree.get_path(outgroup_clade)
        parent = tree.root if len(path) == 1 else path[-2]
        parent.clades.remove(outgroup_clade)
        if len(parent.clades) == 1:
            child = parent.clades[0]
            _retain_parent_confidence(parent, child)
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
    Give each node a unique name, preserving a unique outgroup name.

    Assign ``internal_node_*`` labels to unnamed nodes. Add numeric
    suffixes to duplicate names and log a warning for each rename.
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


def _retain_parent_confidence(parent: Clade, child: Clade) -> None:
    """Copy the parent's confidence to the child, warning if support is lost."""
    if child.confidence is not None and child.confidence != parent.confidence:
        logger.warning(
            "Collapsing single-child clade %r into parent %r: "
            "retaining parent confidence %r; discarding child confidence %r.",
            child.name,
            parent.name,
            parent.confidence,
            child.confidence,
        )
    child.confidence = parent.confidence


def merge_single_child_clades(tree: Tree, outgroup: str | None = None) -> None:
    """
    Collapse chains of single-child clades by summing branch lengths.
    Preserve nodes carrying the requested outgroup label and their parents.
    Keep the parent's confidence and other metadata. Log a warning if the
    child has a different confidence value that is not None.
    """
    queue: deque[Clade] = deque([tree.root])
    while queue:
        clade = queue.popleft()
        while len(clade.clades) == 1:
            child = clade.clades[0]
            if outgroup is not None and outgroup in (clade.name, child.name):
                break
            _retain_parent_confidence(clade, child)
            clade.name = getattr(child, "name", clade.name)
            clade.branch_length = (clade.branch_length or 0.0) + (
                child.branch_length or 0.0
            )
            clade.clades = child.clades
        queue.extend(clade.clades)


def ensure_branch_lengths(tree: Tree) -> None:
    """
    Normalize branch lengths before clustering. Zero lengths are valid.

    Reject non-numeric lengths, NaN, and infinity, including at the root.

    - Set negative non-root branch lengths to 0.0 and log a warning.
    - If no non-root branch has a positive length, set all non-root branch
      lengths to 1.0, treating the tree as unweighted.
    - Otherwise, keep positive lengths and leave missing lengths as None.
      Clustering treats missing and zero lengths as zero, which can produce
      zero within-cluster dispersion and infinite alpha values.
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

    n_nonpositive = 0
    n_missing = 0
    n_negative = 0
    for cl in clades:
        if cl.branch_length is None:
            n_missing += 1
        length = cl.branch_length or 0.0
        if length < 0.0:
            n_negative += 1
            cl.branch_length = 0.0
        if not length > 0.0:
            n_nonpositive += 1

    if n_negative:
        logger.warning(
            "%d of %d branches have a negative length; they are set to "
            "0.0, which can produce zero within-cluster dispersion and "
            "undefined (inf) alpha for the affected clusters.",
            n_negative,
            len(clades),
        )

    if n_nonpositive == 0:
        return

    if n_nonpositive == len(clades):
        logger.warning(
            "All non-root branch lengths are zero or missing. "
            "PhytClust will treat this tree as unweighted, using length 1.0."
        )
        for cl in clades:
            cl.branch_length = 1.0
        return

    if n_missing:
        logger.warning(
            "%d of %d branches have missing lengths; clustering treats them as 0.0.",
            n_missing,
            len(clades),
        )
    logger.info(
        "%d of %d branches have zero or missing lengths. Zero-length branches "
        "are valid; clusters with zero dispersion can have infinite alpha values.",
        n_nonpositive,
        len(clades),
    )
