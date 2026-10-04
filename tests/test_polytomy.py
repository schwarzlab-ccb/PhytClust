import pathlib
from io import StringIO

from Bio import Phylo

from phytclust import PhytClust

POLYTOMY_PATH = (
    pathlib.Path(__file__).parent.parent
    / "examples"
    / "sample_polytomy.newick"
)


def _load_tree(newick: str):
    return Phylo.read(StringIO(newick), "newick")


def _cluster_groups(cmap):
    groups = {}
    for node, cluster_id in cmap.items():
        groups.setdefault(cluster_id, []).append(node.name)
    return sorted(sorted(names) for names in groups.values())


def _assert_valid_partition(groups, expected_k, expected_leaves):
    assert len(groups) == expected_k
    assert sorted(sum(groups, [])) == sorted(expected_leaves)


def test_internal_polytomy_exact_k_backtracking():
    tree = _load_tree("((A:1,B:1,C:1):1,(D:1,E:1):1);")
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
    tree = _load_tree("(A:0.1,B:0.1,C:10,(D:1,E:1):1);")
    pc = PhytClust(tree, polytomy_mode="hard")

    result = pc.run(k=4, plot_scores=False)

    assert result["mode"] == "k"
    groups = _cluster_groups(result["clusters"][0])
    _assert_valid_partition(
        groups, expected_k=4, expected_leaves=["A", "B", "C", "D", "E"]
    )



def test_soft_polytomy_mode_partition_is_valid():
    tree = _load_tree("(A:0.1,B:0.1,C:10,(D:1,E:1):1);")
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


def test_sample_polytomy_hard_k4_recovers_four_branch_groups():
    """
    sample_polytomy.newick has a 4-way root polytomy with branches:
      P1 (A1-A5, 5 leaves, len 0.60)
      X1 (B1-B7, 7 leaves, len 0.50)
      X3 (C1-C10, 10 leaves, len 0.55)
      X4 (D1-D14, 14 leaves, len 0.58)
    k=4 must recover exactly these four groups (verified by cluster sizes).
    """
    tree = Phylo.read(POLYTOMY_PATH, "newick")
    pc = PhytClust(tree, polytomy_mode="hard")

    result = pc.run(k=4, plot_scores=False)

    assert result["mode"] == "k"
    cmap = result["clusters"][0]
    from collections import Counter

    size_counts = Counter(
        sum(1 for v in cmap.values() if v == cid) for cid in set(cmap.values())
    )
    assert dict(size_counts) == {5: 1, 7: 1, 10: 1, 14: 1}


def test_sample_polytomy_soft_k4_recovers_four_branch_groups():
    """Same grouping must hold under soft polytomy mode."""
    tree = Phylo.read(POLYTOMY_PATH, "newick")
    pc = PhytClust(
        tree,
        polytomy_mode="soft",
        soft_polytomy_max_degree=8,
    )

    result = pc.run(k=4, plot_scores=False)

    assert result["mode"] == "k"
    cmap = result["clusters"][0]
    from collections import Counter

    size_counts = Counter(
        sum(1 for v in cmap.values() if v == cid) for cid in set(cmap.values())
    )
    assert dict(size_counts) == {5: 1, 7: 1, 10: 1, 14: 1}


def test_default_polytomy_mode_is_soft():
    tree = _load_tree("(A:0.1,B:0.1,C:10,(D:1,E:1):1);")
    assert PhytClust(tree).polytomy_mode == "soft"
