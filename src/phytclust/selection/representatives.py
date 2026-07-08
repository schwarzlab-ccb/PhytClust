from __future__ import annotations

from typing import Any, Optional, Tuple
import heapq

from Bio.Phylo.BaseTree import Tree

from ..exceptions import ConfigurationError


#  terminal -> path of internal nodes
def map_terminal_to_internal(tree: Tree) -> dict[str, list]:
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


# PD pick
def _get_output_maximizing_pd(ranked_nodes: list[Tuple[str, float]]) -> Tuple[str, str]:
    sum_distances = sum(distance for _, distance in ranked_nodes)
    node_names = [name for name, _ in ranked_nodes]
    maximizing_pd_output = f"Maximizing PD to get {sum_distances}"
    chosen_leaves_output = f"Chosen leaves: {', '.join(node_names)}"
    return maximizing_pd_output, chosen_leaves_output


def maximize_pd(tree: Tree, num_species: Optional[int] = None) -> list[Tuple[str, str]]:
    """
    Convenience wrapper: greedy selection by 'all' distances, maximizing sum.
    Returns human-readable summaries for each incremental selection set.
    """
    ranked_nodes = rank_terminal_nodes(
        tree, num_species=num_species, mode="maximize", distance_ref="all"
    )
    outputs: list[Tuple[str, str]] = []
    maximizing_pd_output, chosen_leaves_output = _get_output_maximizing_pd(ranked_nodes)
    outputs.append((maximizing_pd_output, chosen_leaves_output))
    return outputs


# Distances used for ranking
def _sum_distances_to_all_leaves(tree: Tree) -> dict[Any, float]:
    """For every terminal, the sum of path distances to all other terminals.

    Computed in O(n) with the standard reroot technique (two passes) instead of
    O(n^2) pairwise ``tree.distance`` calls (each itself a path trace):

    - postorder: ``cnt[node]`` = #leaves below, ``down[node]`` = sum of distances
      from node to leaves below it;
    - preorder: ``f[node]`` = sum of distances from node to *all* leaves, via
      ``f[child] = f[node] + (L - 2*cnt[child]) * edge(child)``.

    A leaf's distance to itself is 0, so ``f[leaf]`` is exactly its sum of
    distances to every other leaf.
    """
    root = tree.root
    cnt: dict[Any, int] = {}
    down: dict[Any, float] = {}
    for node in tree.find_clades(order="postorder"):
        if node.is_terminal():
            cnt[node] = 1
            down[node] = 0.0
        else:
            c = 0
            d = 0.0
            for child in node.clades:
                w = float(child.branch_length or 0.0)
                c += cnt[child]
                d += down[child] + cnt[child] * w
            cnt[node] = c
            down[node] = d

    L = cnt[root]
    f: dict[Any, float] = {root: down[root]}
    for node in tree.find_clades(order="preorder"):
        for child in node.clades:
            w = float(child.branch_length or 0.0)
            f[child] = f[node] + (L - 2 * cnt[child]) * w

    return {n: f[n] for n in cnt if n.is_terminal()}


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
        mrca = tree.common_ancestor(terminals)
        for t in terminals:
            distances[t] = float(tree.distance(t, mrca))

    else:
        raise ConfigurationError("distance_ref must be 'all' or 'mrca'")

    return distances


# Greedy ranking
def rank_terminal_nodes(
    tree: Tree,
    num_species: Optional[int] = None,
    *,
    mode: str = "maximize",
    distance_ref: str = "all",
) -> list[Tuple[str, float]]:
    """
    Rank/Select terminals by a simple scalar "distance" measure.

    mode:
        - 'maximize': pick highest distance first (greedy)
        - 'minimize': pick lowest distance first

    Returns a list of (terminal_name, distance_value) in chosen order.
    """
    terminals = [t.name for t in tree.get_terminals() if t.name]
    base = compute_species_distance(tree, terminals, distance_ref=distance_ref)

    if mode not in {"maximize", "minimize"}:
        raise ConfigurationError("mode must be 'maximize' or 'minimize'")

    if mode == "maximize":
        queue = [(-base[t], t) for t in terminals]
    else:
        queue = [(base[t], t) for t in terminals]
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


# One representative per cluster
_STRATEGY_MAP = {
    # strategy -> (mode, distance_ref)
    "central": ("minimize", "mrca"),   # leaf closest to the cluster MRCA
    "divergent": ("maximize", "mrca"),  # leaf farthest from the cluster MRCA
    "medoid": ("minimize", "all"),      # min total distance to cluster members
}


def select_representative_species(
    tree: Tree,
    clusters: dict[str, int],
    *,
    strategy: str = "central",
    mode: Optional[str] = None,
    distance_ref: Optional[str] = None,
) -> list[str]:
    """
    Pick a single representative species per cluster.

    clusters: mapping {species_name -> cluster_id}

    strategy (the simple knob — choose what "representative" means):
        - "central"   : leaf *closest* to the cluster MRCA (default)
        - "divergent" : leaf *farthest* from the cluster MRCA (most divergent)
        - "medoid"    : leaf minimising total distance to all cluster members
                        (O(n) per cluster via the reroot technique)

    Advanced: pass ``mode`` ("minimize"/"maximize") and/or ``distance_ref``
    ("mrca"/"all") to override the strategy mapping directly.
    """
    if mode is None or distance_ref is None:
        if strategy not in _STRATEGY_MAP:
            raise ConfigurationError(
                "strategy must be 'central', 'divergent', or 'medoid'"
            )
        s_mode, s_ref = _STRATEGY_MAP[strategy]
        mode = mode or s_mode
        distance_ref = distance_ref or s_ref
    # group by cluster
    by_cluster: dict[int, list[str]] = {}
    for sp, cid in clusters.items():
        by_cluster.setdefault(cid, []).append(sp)

    reps: list[str] = []

    for cid, species_list in by_cluster.items():
        if len(species_list) == 1:
            reps.append(species_list[0])
            continue

        # rank within each cluster by the chosen criterion
        # Use the full tree; MRCA is computed from the species_list subset
        mrca = tree.common_ancestor(species_list)

        # Build a temporary subtree rooted at MRCA for distance calc context
        sub_tree = Tree(root=mrca, rooted=True)

        ranked = rank_terminal_nodes(
            sub_tree,
            num_species=1,
            mode=mode,
            distance_ref=distance_ref,
        )
        reps.append(ranked[0][0] if ranked else species_list[0])

    return reps
