"""Enumerate partitions tied for the optimal cost, up to a configured limit.

The DP selects one partition using cost, outlier, and cluster-allocation
preferences. List that partition first, followed by other cost ties.
Enumeration requires the child DP tables to be preserved.
"""

from typing import Any, Optional
from numbers import Integral

import numpy as np

from ...exceptions import ConfigurationError, InvalidKError, MissingDPTableError
from .costs import child_subtree_cluster_cost, floating_point_precision
from ...utils.traversal import terminals


class _OptimalPartitionEnumerator:
    def __init__(self, clustering, max_solutions: int):
        self.clustering = clustering
        self.max_solutions = max_solutions
        self.truncated = False
        self.partition_cache: dict[tuple, list] = {}
        self._soft_cache: dict[int, tuple] = {}
        self._hard_cache: dict[int, list] = {}

        root = (
            clustering._tree_wo_outgroup if clustering.outgroup else clustering.tree
        ).root
        root_table = clustering.dp_table[clustering.node_to_id[root]]
        if root_table is None:
            raise MissingDPTableError("DP table at root is missing.")
        self.relative_cost_tolerance = 8.0 * floating_point_precision(root_table.dtype)

    def table(self, node_id: int):
        arr = self.clustering.dp_table[node_id]
        if arr is None:
            raise MissingDPTableError(
                "Child DP tables were freed. save_tied_optima requires the "
                "tables to be preserved; recompute the DP with it enabled."
            )
        return arr

    def find_tied_cost_indices(self, costs) -> list[int]:
        """Positions within tie tolerance of the finite minimum of ``costs``."""
        costs = np.asarray(costs, dtype=float)
        finite = np.isfinite(costs)
        if not finite.any():
            return []
        best = float(costs[finite].min())
        cost_tolerance = self.relative_cost_tolerance * abs(best)
        return [int(i) for i in np.flatnonzero(costs <= best + cost_tolerance)]

    def append_partition_combinations(
        self, out: list, prefix_partitions: list, suffix_partitions: list
    ) -> bool:
        """Append the cross product to ``out``. Returns True if the max_solutions was hit."""
        for a in prefix_partitions:
            for b in suffix_partitions:
                if len(out) >= self.max_solutions:
                    self.truncated = True
                    return True
                out.append(a + b)
        return False

    def enumerate_state_partitions(self, node_id: int, state_index: int) -> list:
        key = (node_id, state_index)
        cached = self.partition_cache.get(key)
        if cached is not None:
            return cached

        own = self.table(node_id)
        if state_index >= len(own) or not np.isfinite(own[state_index]):
            self.partition_cache[key] = []
            return []

        node = self.clustering.postorder_nodes[node_id]
        if state_index == 0:
            result = [((node,),)]
        elif self.clustering.polytomy_backptr[node_id] is not None:
            result = self.polytomy(node, node_id, state_index)
        else:
            result = self.binary(node_id, state_index)

        self.partition_cache[key] = result
        return result

    def binary(self, node_id: int, k: int) -> list:
        left, right = self.clustering.dp_children[node_id]
        left_id = self.clustering.node_to_id[left]
        right_id = self.clustering.node_to_id[right]
        left_dp, right_dp = self.table(left_id), self.table(right_id)

        out: list = []
        for i, j in self.find_tied_child_states(left_dp, right_dp, k):
            if self.append_partition_combinations(
                out,
                self.enumerate_state_partitions(left_id, i),
                self.enumerate_state_partitions(right_id, j),
            ):
                break
        return out

    def find_tied_child_states(
        self, left_dp, right_dp, k: int
    ) -> list[tuple[int, int]]:
        """Tied (i, j) child-state pairs with i + j + 1 == k, as in the DP."""
        max_i = min(k - 1, len(left_dp) - 1)
        min_i = max(0, k - 1 - (len(right_dp) - 1))
        if min_i > max_i:
            return []
        costs = left_dp[min_i : max_i + 1] + right_dp[k - 1 - max_i : k - min_i][::-1]
        return [
            (min_i + p, k - 1 - (min_i + p)) for p in self.find_tied_cost_indices(costs)
        ]

    def polytomy(self, node, node_id: int, state_index: int) -> list:
        mode = self.clustering.polytomy_backptr[node_id]["mode"]
        degree = len(self.clustering.dp_children[node_id])
        if mode == "hard":
            return self.enumerate_hard_prefix(node, node_id, degree - 1, state_index)
        if mode == "soft":
            return self.enumerate_soft_groups(
                node, node_id, (1 << degree) - 1, state_index + 1
            )
        raise MissingDPTableError(f"Unknown polytomy back-pointer mode: {mode}")

    def hard_prefix_tables(self, node, node_id: int) -> list:
        """Rebuild the prefix arrays the hard DP folds over but does not store."""
        cached = self._hard_cache.get(node_id)
        if cached is not None:
            return cached

        children = list(self.clustering.dp_children[node_id])
        max_states = len(self.table(node_id)) - 1
        prefix = self.table(self.clustering.node_to_id[children[0]])
        tables = [prefix]
        n_leaves = self.clustering.num_leaves_per_node[children[0]]

        for child in children[1:]:
            child_dp = self.table(self.clustering.node_to_id[child])
            n_leaves += self.clustering.num_leaves_per_node[child]
            n_new = min(n_leaves, max_states)
            new = np.full(n_new + 1, np.inf, dtype=prefix.dtype)
            for k in range(1, n_new + 1):
                max_i = min(k - 1, len(prefix) - 1)
                min_i = max(0, k - 1 - (len(child_dp) - 1))
                if min_i > max_i:
                    continue
                sums = (
                    prefix[min_i : max_i + 1]
                    + child_dp[k - 1 - max_i : k - min_i][::-1]
                )
                if np.isfinite(sums).any():
                    new[k] = sums.min()
            tables.append(new)
            prefix = new

        self._hard_cache[node_id] = tables
        return tables

    def enumerate_hard_prefix(
        self, node, node_id: int, fold: int, state_index: int
    ) -> list:
        """Solutions covering children[:fold + 1] at state ``state_index``."""
        children = list(self.clustering.dp_children[node_id])
        first_id = self.clustering.node_to_id[children[0]]
        if fold == 0:
            return self.enumerate_state_partitions(first_id, state_index)

        key = ("hard", node_id, fold, state_index)
        cached = self.partition_cache.get(key)
        if cached is not None:
            return cached

        tables = self.hard_prefix_tables(node, node_id)
        child_id = self.clustering.node_to_id[children[fold]]
        prefix_dp, child_dp = tables[fold - 1], self.table(child_id)

        out: list = []
        for i, j in self.find_tied_child_states(prefix_dp, child_dp, state_index):
            prefix_partitions = self.enumerate_hard_prefix(node, node_id, fold - 1, i)
            if self.append_partition_combinations(
                out, prefix_partitions, self.enumerate_state_partitions(child_id, j)
            ):
                break

        self.partition_cache[key] = out
        return out

    def soft_tables(self, node, node_id: int) -> tuple:
        """Recompute the subset DP, keeping every optimal predecessor."""
        cached = self._soft_cache.get(node_id)
        if cached is not None:
            return cached

        clustering = self.clustering
        children = list(self.clustering.dp_children[node_id])
        child_count = len(children)
        all_children_mask = (1 << child_count) - 1
        max_cluster_count = len(self.table(node_id)) - 1
        min_cluster_size = getattr(clustering, "min_cluster_size", 1)

        atom_costs = [child_subtree_cluster_cost(c, clustering) for c in children]
        child_leaves = [clustering.num_leaves_per_node[c] for c in children]

        block_leaves = [0] * (all_children_mask + 1)
        block_cost = [0.0] * (all_children_mask + 1)
        for mask in range(1, all_children_mask + 1):
            lsb = mask & -mask
            i = lsb.bit_length() - 1
            prev = mask ^ lsb
            block_leaves[mask] = block_leaves[prev] + child_leaves[i]
            block_cost[mask] = block_cost[prev] + atom_costs[i]

        transitions: dict[tuple, list] = {}
        dp: dict[int, dict[int, float]] = {0: {0: 0.0}}

        for mask in range(all_children_mask + 1):
            if mask not in dp:
                continue
            remaining_children_mask = all_children_mask ^ mask
            if remaining_children_mask == 0:
                continue

            first_remaining_child_bit = (
                remaining_children_mask & -remaining_children_mask
            )
            t = first_remaining_child_bit.bit_length() - 1
            child_dp = self.table(clustering.node_to_id[children[t]])

            for current_cluster_count, base in dp[mask].items():
                for s, cost in enumerate(child_dp):
                    if not np.isfinite(cost):
                        continue
                    next_cluster_count = current_cluster_count + s + 1
                    if next_cluster_count > max_cluster_count + 1:
                        continue
                    next_children_mask = mask | first_remaining_child_bit
                    candidate_cost = base + cost
                    dp.setdefault(next_children_mask, {})
                    if candidate_cost < dp[next_children_mask].get(
                        next_cluster_count, np.inf
                    ):
                        dp[next_children_mask][next_cluster_count] = candidate_cost
                    transitions.setdefault(
                        (next_children_mask, next_cluster_count), []
                    ).append(
                        (
                            mask,
                            current_cluster_count,
                            first_remaining_child_bit,
                            s,
                            candidate_cost,
                        )
                    )

                sub = remaining_children_mask ^ first_remaining_child_bit
                while sub:
                    block = first_remaining_child_bit | sub
                    if (
                        block.bit_count() >= 2
                        and np.isfinite(block_cost[block])
                        and block_leaves[block] >= min_cluster_size
                    ):
                        next_cluster_count = current_cluster_count + 1
                        if next_cluster_count <= max_cluster_count + 1:
                            next_children_mask = mask | block
                            candidate_cost = base + block_cost[block]
                            dp.setdefault(next_children_mask, {})
                            if candidate_cost < dp[next_children_mask].get(
                                next_cluster_count, np.inf
                            ):
                                dp[next_children_mask][
                                    next_cluster_count
                                ] = candidate_cost
                            transitions.setdefault(
                                (next_children_mask, next_cluster_count), []
                            ).append(
                                (
                                    mask,
                                    current_cluster_count,
                                    block,
                                    None,
                                    candidate_cost,
                                )
                            )
                    sub = (sub - 1) & (
                        remaining_children_mask ^ first_remaining_child_bit
                    )

        optimal: dict[tuple, list] = {}
        for key, candidates in transitions.items():
            best = dp[key[0]][key[1]]
            cost_tolerance = self.relative_cost_tolerance * abs(best)
            keep = [c[:4] for c in candidates if c[4] <= best + cost_tolerance]
            if keep:
                optimal[key] = keep

        result = (dp, optimal)
        self._soft_cache[node_id] = result
        return result

    def enumerate_soft_groups(self, node, node_id: int, mask: int, q: int) -> list:
        if mask == 0:
            return [()] if q == 0 else []

        key = ("soft", node_id, mask, q)
        cached = self.partition_cache.get(key)
        if cached is not None:
            return cached

        children = list(self.clustering.dp_children[node_id])
        _, optimal = self.soft_tables(node, node_id)

        out: list = []
        for (
            previous_children_mask,
            previous_cluster_count,
            block,
            child_state,
        ) in optimal.get((mask, q), []):
            bits = [i for i in range(len(children)) if block & (1 << i)]
            if child_state is None:
                tail = [((tuple(children[i] for i in bits)),)]
            else:
                child_id = self.clustering.node_to_id[children[bits[0]]]
                tail = self.enumerate_state_partitions(child_id, child_state)
            if self.append_partition_combinations(
                out,
                self.enumerate_soft_groups(
                    node, node_id, previous_children_mask, previous_cluster_count
                ),
                tail,
            ):
                break

        self.partition_cache[key] = out
        return out


def _as_cluster_map(clustering, solution) -> dict[Any, int]:
    """Turn a tuple of node groups into the terminal -> cluster id map."""
    name_leaves = getattr(clustering, "name_leaves_per_node", None)
    clusters: dict[Any, int] = {}
    for cluster_id, group in enumerate(solution):
        for group_node in group:
            terms = name_leaves.get(group_node) if name_leaves is not None else None
            if terms is None:
                terms = terminals(group_node)
            for terminal in terms:
                clusters[terminal] = cluster_id
    return clusters


def _signature(clusters: dict[Any, int]) -> frozenset:
    """Label-independent identity of a partition, for deduplication."""
    groups: dict[int, list] = {}
    for terminal, cluster_id in clusters.items():
        groups.setdefault(cluster_id, []).append(
            getattr(terminal, "name", None) or id(terminal)
        )
    return frozenset(frozenset(g) for g in groups.values())


def enumerate_tied_optima(
    clustering,
    k: int,
    *,
    max_solutions: Optional[int] = None,
    chosen: Optional[dict[Any, int]] = None,
) -> dict[str, Any]:
    """Return up to the configured limit of partitions tied for optimal cost.

    Parameters
    ----------
    k : int
        Number of clusters.
    max_solutions : int, optional
        Cap on returned partitions; defaults to ``clustering.max_tied_optima``. The
        enumeration stops at the max_solutions and reports ``truncated``, so a tie whose
        returned list stays bounded.
    chosen : dict, optional
        The partition ``backtrack`` returned for this k, listed first in the
        result. ``backtrack`` passes it in directly; when omitted it is looked
        up with ``clustering.get_clusters(k)``.

    Returns
    -------
    dict
        ``k``, ``score`` (the optimal objective value), ``n_solutions``,
        ``truncated``, and ``solutions`` (list of terminal -> cluster id maps).
        The first entry matches what ``backtrack`` returns for the same k.
    """
    if clustering.outlier.prefer_fewer:
        raise ConfigurationError(
            "Tied-optima enumeration is not supported with "
            "outlier.prefer_fewer=True, which minimises the outlier count "
            "before cost and so does not select on cost ties alone."
        )
    if not getattr(clustering, "_dp_ready", False):
        raise MissingDPTableError(
            "DP table not computed. Call compute_dp_table(clustering) first."
        )
    if isinstance(k, (bool, np.bool_)) or not isinstance(k, Integral) or k <= 0:
        raise InvalidKError("k must be a positive integer.")

    max_solutions = (
        max_solutions
        if max_solutions is not None
        else getattr(clustering, "max_tied_optima", 100)
    )
    if (
        isinstance(max_solutions, (bool, np.bool_))
        or not isinstance(max_solutions, Integral)
        or max_solutions < 1
    ):
        raise ConfigurationError("max_solutions must be a positive integer.")

    active_tree = (
        clustering._tree_wo_outgroup if clustering.outgroup else clustering.tree
    )
    root_id = clustering.node_to_id[active_tree.root]
    root_dp = clustering.dp_table[root_id]
    if root_dp is None:
        raise MissingDPTableError("DP table at root is missing.")
    if k - 1 >= len(root_dp) or not np.isfinite(root_dp[k - 1]):
        raise InvalidKError(f"Partition into {k} clusters is infeasible.")

    enumerator = _OptimalPartitionEnumerator(clustering, max_solutions)
    raw_solutions = enumerator.enumerate_state_partitions(root_id, k - 1)

    seen: set = set()
    solutions: list[dict[Any, int]] = []
    for solution in raw_solutions:
        clusters = _as_cluster_map(clustering, solution)
        signature = _signature(clusters)
        if signature in seen:
            continue
        seen.add(signature)
        solutions.append(clusters)

    if chosen is None:
        chosen = clustering.get_clusters(k)
    if chosen is not None:
        chosen_sig = _signature(chosen)
        for i, sol in enumerate(solutions):
            if _signature(sol) == chosen_sig:
                if i:
                    solutions.insert(0, solutions.pop(i))
                break
        else:
            solutions.insert(0, chosen.copy())
            if len(solutions) > max_solutions:
                solutions.pop()
                enumerator.truncated = True

    return {
        "k": int(k),
        "score": float(root_dp[k - 1]),
        "n_solutions": len(solutions),
        "truncated": bool(enumerator.truncated),
        "solutions": solutions,
    }
