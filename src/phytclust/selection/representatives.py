from __future__ import annotations

from typing import Any, Optional, Tuple
import heapq
from numbers import Integral
import math

import numpy as np

from Bio.Phylo.BaseTree import Tree

from ..exceptions import ConfigurationError, InvalidTreeError
from ..utils.traversal import iter_clades


def _leaf_name_lookup(tree):
    """Require unique named leaves and finite non-negative non-root lengths."""
    lookup = {}
    for node in iter_clades(tree.root):
        if node is not tree.root and node.branch_length is not None:
            try:
                valid = math.isfinite(node.branch_length) and node.branch_length >= 0
            except (TypeError, ValueError, OverflowError):
                valid = False
            if not valid:
                raise InvalidTreeError(
                    f"Branch length for node {node.name!r} must be finite and non-negative."
                )
        if node.is_terminal():
            if not isinstance(node.name, str) or not node.name:
                raise InvalidTreeError(
                    "Representative selection requires named leaves."
                )
            if node.name in lookup:
                raise InvalidTreeError(f"Leaf name {node.name!r} is duplicated.")
            lookup[node.name] = node
    return lookup


def _tree_paths(tree):
    """Build parent and level lookups in one traversal."""
    parents, levels = {tree.root: None}, {tree.root: 0}
    for node in iter_clades(tree.root):
        for child in node.clades:
            parents[child] = node
            levels[child] = levels[node] + 1
    return parents, levels


def _descendant_distances(root):
    """Accumulate distances below a node without subtracting large root depths."""
    distances = {root: 0.0}
    for node in iter_clades(root):
        for child in node.clades:
            distances[child] = distances[node] + float(child.branch_length or 0.0)
    return distances


def _common_ancestor(leaves, parents, levels):
    """Find the common ancestor using parent links and node levels."""
    ancestor = leaves[0]
    for leaf in leaves[1:]:
        candidate = leaf
        while levels[ancestor] > levels[candidate]:
            ancestor = parents[ancestor]
        while levels[candidate] > levels[ancestor]:
            candidate = parents[candidate]
        while ancestor is not candidate:
            ancestor, candidate = parents[ancestor], parents[candidate]
    return ancestor


def map_terminal_to_internal(tree: Tree) -> dict[str, list]:
    """Map leaf names to their internal ancestors, from root to parent."""
    _leaf_name_lookup(tree)
    terminal_to_internal: dict[str, list] = {}
    stack = [(tree.root, [])]
    while stack:
        clade, path = stack.pop()
        if not clade.clades:
            name = getattr(clade, "name", None)
            if name is not None:
                terminal_to_internal[name] = path
        else:
            new_path = path + [clade]
            for child in clade.clades:
                stack.append((child, new_path))
    return terminal_to_internal


def _get_output_maximizing_pd(ranked_nodes: list[Tuple[str, float]]) -> Tuple[str, str]:
    sum_distances = sum(distance for _, distance in ranked_nodes)
    node_names = [name for name, _ in ranked_nodes]
    maximizing_pd_output = f"Sum of selected distance scores: {sum_distances}"
    chosen_leaves_output = f"Chosen leaves: {', '.join(node_names)}"
    return maximizing_pd_output, chosen_leaves_output


def maximize_pd(tree: Tree, num_species: Optional[int] = None) -> list[Tuple[str, str]]:
    """Rank leaves by summed distance to all leaves and return one text summary.

    This compatibility name does not calculate phylogenetic diversity."""
    ranked_nodes = rank_terminal_nodes(
        tree, num_species=num_species, mode="maximize", distance_ref="all"
    )
    outputs: list[Tuple[str, str]] = []
    maximizing_pd_output, chosen_leaves_output = _get_output_maximizing_pd(ranked_nodes)
    outputs.append((maximizing_pd_output, chosen_leaves_output))
    return outputs


def _sum_distances_to_members(
    tree: Tree, members: set[str] | None = None
) -> dict[Any, float]:
    """Sum distances between selected leaves in two iterative tree passes.

    When members is None, include every leaf. Return scores only for members.
    The first pass counts members and sums descendant distances; the second
    propagates distance sums across parent-child edges."""
    root = tree.root
    member_counts: dict[Any, int] = {}
    descendant_distance_sums: dict[Any, float] = {}
    for node in iter_clades(tree.root, "postorder"):
        if node.is_terminal():
            member_counts[node] = int(members is None or node.name in members)
            descendant_distance_sums[node] = 0.0
        else:
            child_member_count = 0
            child_distance_sum = 0.0
            for child in node.clades:
                branch_length = float(child.branch_length or 0.0)
                child_member_count += member_counts[child]
                child_distance_sum += (
                    descendant_distance_sums[child]
                    + member_counts[child] * branch_length
                )
            member_counts[node] = child_member_count
            descendant_distance_sums[node] = child_distance_sum

    member_count = member_counts[root]
    distance_sums: dict[Any, float] = {root: descendant_distance_sums[root]}
    for node in iter_clades(tree.root):
        for child in node.clades:
            branch_length = float(child.branch_length or 0.0)
            distance_sums[child] = (
                distance_sums[node]
                + (member_count - 2 * member_counts[child]) * branch_length
            )

    return {
        n: distance_sums[n]
        for n in member_counts
        if n.is_terminal() and (members is None or n.name in members)
    }


def _sum_distances_to_all_leaves(tree: Tree) -> dict[Any, float]:
    """Sum each terminal's distances to every terminal in the tree."""
    return _sum_distances_to_members(tree)


def compute_species_distance(
    tree: Tree,
    terminals: list[str],
    *,
    distance_ref: str = "all",
) -> dict[str, float]:
    """
    Compute a distance measure for each terminal.

    distance_ref:
        - 'all': sum of distances from this terminal to all other terminals
          (O(n) via the reroot technique)
        - 'mrca': distance from this terminal to the MRCA of all 'terminals'

    Returns: dict {terminal_name: distance_value}
    """
    lookup = _leaf_name_lookup(tree)
    missing = set(terminals).difference(lookup)
    if missing:
        raise ConfigurationError(f"Species absent from the tree: {sorted(missing)}")
    distances: dict[str, float] = {}

    if distance_ref == "all":
        by_node = _sum_distances_to_all_leaves(tree)
        wanted = set(terminals)
        for node, dist in by_node.items():
            name = getattr(node, "name", None)
            if name in wanted:
                distances[name] = float(dist)

    elif distance_ref == "mrca":
        if not terminals:
            return {}
        parents, levels = _tree_paths(tree)
        mrca = _common_ancestor([lookup[name] for name in terminals], parents, levels)
        local_distances = _descendant_distances(mrca)
        for name in terminals:
            distances[name] = local_distances[lookup[name]]

    else:
        raise ConfigurationError("distance_ref must be 'all' or 'mrca'")

    return distances


def rank_terminal_nodes(
    tree: Tree,
    num_species: Optional[int] = None,
    *,
    mode: str = "maximize",
    distance_ref: str = "all",
) -> list[Tuple[str, float]]:
    """Return leaf names and fixed distance scores in ranked order.

    Maximize selects highest scores; minimize selects lowest scores.
    Equal scores choose the alphabetically first name. Return at most
    num_species leaves, or all leaves when it is None."""
    if num_species is not None and (
        isinstance(num_species, (bool, np.bool_))
        or not isinstance(num_species, Integral)
        or num_species < 1
    ):
        raise ConfigurationError("num_species must be a positive integer or None.")
    terminals = list(_leaf_name_lookup(tree))
    distance_scores = compute_species_distance(
        tree, terminals, distance_ref=distance_ref
    )

    if mode not in {"maximize", "minimize"}:
        raise ConfigurationError("mode must be 'maximize' or 'minimize'")

    if mode == "maximize":
        queue = [(-distance_scores[t], t) for t in terminals]
    else:
        queue = [(distance_scores[t], t) for t in terminals]
    heapq.heapify(queue)

    selected: list[Tuple[str, float]] = []
    seen = set()
    while queue and (num_species is None or len(selected) < num_species):
        priority, terminal = heapq.heappop(queue)
        if terminal in seen:
            continue
        seen.add(terminal)
        actual = -priority if mode == "maximize" else priority
        selected.append((terminal, float(actual)))

    return selected


_STRATEGY_MAP = {
    "central": ("minimize", "mrca"),
    "divergent": ("maximize", "mrca"),
    "medoid": ("minimize", "all"),
}


def select_representative_species(
    tree: Tree,
    clusters: dict[str, int],
    *,
    strategy: str = "central",
    mode: Optional[str] = None,
    distance_ref: Optional[str] = None,
) -> list[str]:
    """Select one member of each cluster.

    Central selects the leaf closest to the cluster common ancestor; divergent
    selects the farthest. Medoid minimizes summed distances to cluster members.
    Equal scores choose the alphabetically first name.

    Override strategy defaults with mode (minimize/maximize) and distance_ref
    (mrca/all). Names must identify unique leaves in the tree."""
    if mode is None or distance_ref is None:
        if strategy not in _STRATEGY_MAP:
            raise ConfigurationError(
                "strategy must be 'central', 'divergent', or 'medoid'"
            )
        strategy_mode, strategy_distance_reference = _STRATEGY_MAP[strategy]
        mode = strategy_mode if mode is None else mode
        distance_ref = (
            strategy_distance_reference if distance_ref is None else distance_ref
        )
    if mode not in {"maximize", "minimize"}:
        raise ConfigurationError("mode must be 'maximize' or 'minimize'")
    if distance_ref not in {"all", "mrca"}:
        raise ConfigurationError("distance_ref must be 'all' or 'mrca'")
    lookup = _leaf_name_lookup(tree)
    missing = set(clusters).difference(lookup)
    if missing:
        raise ConfigurationError(
            f"Cluster contains species absent from the tree: {sorted(missing)}"
        )
    parents, levels = _tree_paths(tree)
    by_cluster: dict[int, list[str]] = {}
    for species_name, cluster_id in clusters.items():
        by_cluster.setdefault(cluster_id, []).append(species_name)

    representatives: list[str] = []

    for cluster_id, species_list in by_cluster.items():
        if len(species_list) == 1:
            representatives.append(species_list[0])
            continue

        mrca = _common_ancestor(
            [lookup[name] for name in species_list], parents, levels
        )

        sub_tree = Tree(root=mrca, rooted=True)

        members = set(species_list)
        if distance_ref == "all":
            distances = {
                node.name: distance
                for node, distance in _sum_distances_to_members(
                    sub_tree, members
                ).items()
            }
        elif distance_ref == "mrca":
            local_distances = _descendant_distances(mrca)
            distances = {name: local_distances[lookup[name]] for name in species_list}
        else:
            raise ConfigurationError("distance_ref must be 'all' or 'mrca'")
        if mode not in {"maximize", "minimize"}:
            raise ConfigurationError("mode must be 'maximize' or 'minimize'")
        missing = members.difference(distances)
        if missing:
            raise ConfigurationError(
                f"Cluster contains species absent from the tree: {sorted(missing)}"
            )
        best = min(
            members,
            key=lambda name: (
                -distances[name] if mode == "maximize" else distances[name],
                name,
            ),
        )
        representatives.append(best)

    return representatives
