import itertools
import random
from collections import Counter
from itertools import product

import numpy as np
import pytest
from Bio.Phylo.BaseTree import Clade, Tree

from conftest import POLYTOMY_TREE, load_tree, tree_from_newick
from phytclust import PhytClust
from phytclust.algo.dp.polytomy import compute_polytomy_dp
from phytclust.config import OutlierConfig
from phytclust.exceptions import ConfigurationError


def _cluster_groups(cmap):
    groups = {}
    for node, cluster_id in cmap.items():
        groups.setdefault(cluster_id, []).append(node.name)
    return sorted(sorted(names) for names in groups.values())


def _assert_valid_partition(groups, expected_k, expected_leaves):
    assert len(groups) == expected_k
    assert sorted(sum(groups, [])) == sorted(expected_leaves)


def test_internal_polytomy_exact_k_backtracking():
    tree = tree_from_newick("((A:1,B:1,C:1):1,(D:1,E:1):1);")
    pc = PhytClust(tree)

    result = pc.run(k=4, plot_scores=False)

    assert result["mode"] == "k"
    groups = _cluster_groups(result["clusters"][0])
    _assert_valid_partition(
        groups, expected_k=4, expected_leaves=["A", "B", "C", "D", "E"]
    )


def test_root_polytomy_exact_k_backtracking():
    # Hard mode requires each cluster to be contained within exactly one child
    # subtree of the polytomy. With 4 root children the minimum feasible k is 4
    # (one cluster per root child); k=2 would need cross-child merges (soft mode).
    tree = tree_from_newick("(A:0.1,B:0.1,C:10,(D:1,E:1):1);")
    pc = PhytClust(tree, polytomy_mode="hard")

    result = pc.run(k=4, plot_scores=False)

    assert result["mode"] == "k"
    groups = _cluster_groups(result["clusters"][0])
    _assert_valid_partition(
        groups, expected_k=4, expected_leaves=["A", "B", "C", "D", "E"]
    )


def test_soft_polytomy_mode_partition_is_valid():
    tree = tree_from_newick("(A:0.1,B:0.1,C:10,(D:1,E:1):1);")
    pc = PhytClust(
        tree,
        polytomy_mode="soft",
        soft_polytomy_max_degree=8,
    )

    result = pc.run(k=2, plot_scores=False)

    assert result["mode"] == "k"
    groups = _cluster_groups(result["clusters"][0])
    _assert_valid_partition(
        groups, expected_k=2, expected_leaves=["A", "B", "C", "D", "E"]
    )


@pytest.mark.parametrize("polytomy_mode", ["hard", "soft"])
def test_sample_polytomy_k4_recovers_four_branch_groups(polytomy_mode):
    """
    sample_polytomy.newick has a 4-way root polytomy with branches:
      P1 (A1-A5, 5 leaves, len 0.60)
      X1 (B1-B7, 7 leaves, len 0.50)
      X3 (C1-C10, 10 leaves, len 0.55)
      X4 (D1-D14, 14 leaves, len 0.58)
    k=4 must recover exactly these four groups (verified by cluster sizes).
    """
    tree = load_tree(POLYTOMY_TREE)
    pc = PhytClust(tree, polytomy_mode=polytomy_mode, soft_polytomy_max_degree=8)

    result = pc.run(k=4, plot_scores=False)

    assert result["mode"] == "k"
    cluster_map = result["clusters"][0]
    size_counts = Counter(Counter(cluster_map.values()).values())
    assert dict(size_counts) == {5: 1, 7: 1, 10: 1, 14: 1}


def test_default_polytomy_mode_is_soft():
    tree = tree_from_newick("(A:0.1,B:0.1,C:10,(D:1,E:1):1);")
    assert PhytClust(tree).polytomy_mode == "soft"


def child_subtree(prefix):
    leaves = [Clade(name=f"{prefix}{index}", branch_length=0) for index in range(4)]
    return Clade(
        branch_length=1,
        clades=[
            Clade(branch_length=0, clades=leaves[:2]),
            Clade(branch_length=0, clades=leaves[2:]),
        ],
    )


@pytest.mark.parametrize("mode", ["hard", "soft"])
@pytest.mark.parametrize("cluster_count", [6, 7, 8])
@pytest.mark.parametrize("reverse_children", [False, True])
@pytest.mark.parametrize("use_float32", [False, True])
def test_polytomy_chooses_even_child_cluster_counts(
    mode, cluster_count, reverse_children, use_float32
):
    children = [child_subtree(prefix) for prefix in ["A", "B", "C"]]
    if reverse_children:
        children.reverse()
    pc = PhytClust(
        Tree(root=Clade(clades=children)),
        polytomy_mode=mode,
        max_k=9,
        dp_float32=use_float32,
        save_tied_optima=True,
        max_tied_optima=1,
    )
    chosen = pc.get_clusters(cluster_count)
    allocations = []
    child_cluster_ids = []
    for child in children:
        cluster_ids = {chosen[leaf] for leaf in child.get_terminals()}
        allocations.append(len(cluster_ids))
        child_cluster_ids.append(cluster_ids)
    # These zero-cost optima keep clusters inside each child subtree.
    assert sum(allocations) == cluster_count
    assert not any(
        a & b
        for index, a in enumerate(child_cluster_ids)
        for b in child_cluster_ids[index + 1 :]
    )
    feasible_allocations = [
        counts
        for counts in product(range(1, 5), repeat=3)
        if sum(counts) == cluster_count
    ]
    assert sum(count**2 for count in allocations) == min(
        sum(count**2 for count in counts) for counts in feasible_allocations
    )
    assert max(allocations) - min(allocations) <= 1
    result = pc.tied_optima[cluster_count]
    assert result["n_solutions"] == 1
    assert result["solutions"][0] == chosen


def partitions(items):
    if not items:
        yield []
        return
    first, *rest = items
    for partition in partitions(rest):
        yield [[first], *partition]
        for index in range(len(partition)):
            yield partition[:index] + [[first, *partition[index]]] + partition[
                index + 1 :
            ]


@pytest.mark.parametrize(
    "minimum_size,support_enabled,prefer_fewer",
    list(product([1, 2], [False, True], [False, True])),
)
def test_soft_star_matches_exhaustive_partitions(
    minimum_size, support_enabled, prefer_fewer
):
    pc = PhytClust(
        "(a:1,b:2,c:4,d:8,e:16);",
        min_cluster_size=minimum_size,
        use_branch_support=support_enabled,
        max_k=5,
        outlier=OutlierConfig(size_threshold=3, prefer_fewer=prefer_fewer),
    )
    pc.tree.root.confidence = 50
    pc._ensure_dp(required_cap=5)
    leaves = pc.tree.get_terminals()
    root_id = pc.node_to_id[pc.tree.root]
    root_costs = pc.unpenalised_dp_table[root_id]
    for count in range(1, 6):
        candidates = []
        for partition in partitions(list(range(5))):
            if len(partition) != count or any(
                len(group) < minimum_size for group in partition
            ):
                continue
            cost = sum(
                sum(leaves[i].branch_length for i in group)
                for group in partition
                if len(group) > 1
            )
            if count == 1 and support_enabled:
                cost /= 0.5
            outliers = sum(len(group) < 3 for group in partition)
            candidates.append((cost, outliers))
        if not candidates:
            assert np.isinf(root_costs[count - 1])
            continue
        expected = min(
            candidates, key=lambda item: (item[1], item[0]) if prefer_fewer else item
        )
        assert root_costs[count - 1] == expected[0]
        assert pc._outlier_counts[root_id][count - 1] == expected[1]
        clusters = pc.get_clusters(count)
        assert set(clusters) == set(leaves)
        assert len(set(clusters.values())) == count


def test_dispatch_rejects_unknown_mode():
    pc = PhytClust("(a:1,b:2,c:3);")
    pc.polytomy_mode = "unknown"
    with pytest.raises(ConfigurationError, match="polytomy_mode"):
        compute_polytomy_dp(pc.tree.root, pc, 3, 1, None, np.float64)


def _hard_allocations(pc, k):
    """Brute force: every allocation (m_1..m_q) at the root, with its cost."""
    rows = [pc.unpenalised_dp_table[pc.node_to_id[c]] for c in pc.tree.root.clades]
    out = []
    for alloc in itertools.product(*[range(1, len(r)) for r in rows]):
        if sum(alloc) != k:
            continue
        cost = sum(float(r[m - 1]) for r, m in zip(rows, alloc))
        if np.isfinite(cost):
            out.append((cost, sum(m * m for m in alloc)))
    return out


@pytest.mark.parametrize("seed", range(12))
def test_hard_polytomy_ties_minimise_sum_of_squares(seed):
    """Among cost-tied allocations at a node, the smallest sum of m_j^2 wins."""
    rng = random.Random(seed)
    names = [f"t{i}" for i in range(60)]

    def sub(n):
        length = rng.choice([1, 1, 2])
        if n == 1:
            return f"{names.pop()}:{length}"
        left = rng.randint(1, n - 1)
        return f"({sub(left)},{sub(n - left)}):{length}"

    sizes = [rng.choice([1, 2, 3, 5, 9]) for _ in range(rng.randint(3, 5))]
    newick = "(" + ",".join(sub(s) for s in sizes) + ");"
    pc = PhytClust(
        tree=tree_from_newick(newick), polytomy_mode="hard", preserve_dp_tables=True
    )
    n = pc.num_terminals
    pc.get_clusters(n)
    leafsets = [set(pc.name_leaves_per_node[c]) for c in pc.tree.root.clades]

    for k in range(len(sizes), n + 1):
        cluster_map = pc.get_clusters(k)
        chosen = [len({cluster_map[t] for t in leaves}) for leaves in leafsets]
        candidates = _hard_allocations(pc, k)
        best = min(c for c, _ in candidates)
        tied = [sq for c, sq in candidates if abs(c - best) <= 1e-9 * max(1.0, best)]
        assert sum(m * m for m in chosen) == min(tied)
