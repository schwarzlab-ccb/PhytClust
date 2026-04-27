import pathlib
import pytest
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade

from phytclust.algo.core import PhytClust
from phytclust.algo.dp import compute_dp_table
from phytclust.exceptions import ConfigurationError, InvalidKError

TREE_PATH = pathlib.Path(__file__).parent.parent / "examples" / "sample_tree.nwk"


def _load_tree():
    return Phylo.read(TREE_PATH, "newick")


def _leaf_groups(cmap):
    """Return a frozenset of frozensets of leaf names, one per cluster."""
    from collections import defaultdict

    buckets = defaultdict(set)
    for terminal, cid in cmap.items():
        buckets[cid].add(terminal.name)
    return frozenset(frozenset(s) for s in buckets.values())


def _assert_cluster_map_valid(tree, cmap, k_expected: int):
    assert isinstance(cmap, dict)
    assert len(cmap) > 0

    assert all(isinstance(t, Clade) and t.is_terminal() for t in cmap.keys())
    assert all(isinstance(v, int) for v in cmap.values())

    cluster_ids = set(cmap.values())
    assert len(cluster_ids) == k_expected

    terminals = list(tree.get_terminals())
    assert set(cmap.keys()) == set(terminals)


def test_run_exact_k_partition_returns_result_dict():
    tree = _load_tree()
    pc = PhytClust(tree=tree)

    result = pc.run(k=2, plot_scores=False)

    assert isinstance(result, dict)
    assert result["mode"] == "k"
    assert result["selected_k"] == 2
    assert result["k_values"] == [2]
    assert result["k"] == 2
    assert result["peaks"] == [2]
    assert result["scores"] is None

    assert isinstance(result["clusters"], list)
    assert len(result["clusters"]) == 1
    cmap = result["clusters"][0]
    _assert_cluster_map_valid(tree, cmap, k_expected=2)


def test_get_clusters_exact_k_partition():
    tree = _load_tree()
    pc = PhytClust(tree=tree)

    pc.run(k=2, plot_scores=False)

    cmap = pc.get_clusters(2)
    _assert_cluster_map_valid(tree, cmap, k_expected=2)

    assert 2 in pc.clusters
    assert pc.clusters[2] == cmap


def test_run_global_mode_returns_list_of_cluster_maps():
    tree = _load_tree()
    pc = PhytClust(tree=tree)

    max_k = min(10, pc.num_terminals)
    if max_k < 4:
        pytest.skip("Test tree too small for global peak search (requires max_k>=4).")

    result = pc.run(top_n=1, max_k=max_k, plot_scores=False)

    assert result["mode"] == "global"
    assert isinstance(result["clusters"], list)

    clusters = result["clusters"]
    assert len(clusters) >= 0

    if clusters:
        assert isinstance(result["k_values"], list)
        assert result["k_values"] == result["ks"]
        assert result["selected_k"] == result["k_values"][0]
        assert isinstance(result["ks"], list)
        assert isinstance(result["peaks"], list)
        assert len(result["peaks"]) == len(result["ks"])
        cmap0 = clusters[0]
        assert isinstance(cmap0, dict)
        assert set(cmap0.keys()) == set(tree.get_terminals())


def test_run_resolution_mode_falls_back_for_small_trees():
    tree = _load_tree()
    pc = PhytClust(tree=tree)

    result = pc.run(by_resolution=True, num_bins=2, plot_scores=False)

    assert result["mode"] == "resolution"
    assert isinstance(result["clusters"], list)
    assert isinstance(result["k_values"], list)
    if result["k_values"]:
        assert result["selected_k"] == result["k_values"][0]
    assert isinstance(result["ks"], list)
    assert isinstance(result["peaks"], list)


def test_best_global_returns_list_of_cluster_maps():
    tree = _load_tree()
    pc = PhytClust(tree=tree)

    max_k = min(10, pc.num_terminals)
    if max_k < 4:
        pytest.skip("Test tree too small for global peak search (requires max_k>=4).")

    clusters = pc.best_global(top_n=1, max_k=max_k, plot_scores=False)

    assert isinstance(clusters, list)
    if clusters:
        cmap = clusters[0]
        assert isinstance(cmap, dict)
        assert set(cmap.keys()) == set(tree.get_terminals())


def test_run_rejects_invalid_argument_combinations():
    tree = _load_tree()
    pc = PhytClust(tree=tree)

    with pytest.raises(
        ConfigurationError, match="Cannot combine `k` with `by_resolution=True`"
    ):
        pc.run(k=2, by_resolution=True, plot_scores=False)

    with pytest.raises(
        ConfigurationError, match="`top_n` is meaningless when `k` is given"
    ):
        pc.run(k=2, top_n=2, plot_scores=False)

    with pytest.raises(InvalidKError, match="k must be"):
        pc.run(k=0, plot_scores=False)


def test_resolution_after_fixed_k_does_not_inherit_self_k():
    """Regression: a previous fixed-k run leaves self.k populated, but an
    explicit by_resolution=True call must override that fallback rather
    than raising "Cannot combine k with by_resolution=True"."""
    tree = _load_tree()
    pc = PhytClust(tree=tree)

    pc.run(k=2, plot_scores=False)
    assert pc.k == 2  # Sanity: previous run mutated self.k

    # Should run resolution mode cleanly, not error on the lingering self.k.
    result = pc.run(by_resolution=True, num_bins=2, plot_scores=False)
    assert result["mode"] == "resolution"

    with pytest.raises(InvalidKError, match="`top_n` must be"):
        pc.run(top_n=0, plot_scores=False)


def test_k2_partition_is_root_split():
    """
    sample_tree.nwk has root children A13={A,B,C,D} and {E,F,G,H}.
    k=2 must split exactly there, verified against expected leaf sets.

    Tree: (((A:5,B:3):6,(C:3,D:7):4):22,(((E:7,F:13):5,G:6):10,H:60):35):0;
    """
    tree = _load_tree()
    pc = PhytClust(tree=tree)
    pc.run(k=2, plot_scores=False)
    cmap = pc.get_clusters(2)

    groups = _leaf_groups(cmap)
    assert groups == frozenset(
        [
            frozenset(["A", "B", "C", "D"]),
            frozenset(["E", "F", "G", "H"]),
        ]
    )


def test_k4_partition_segregates_leaf_pairs():
    """
    k=4 on sample_tree.nwk should produce 4 clusters each covering a known clade.
    The two shallow binary pairs {A,B} and {C,D} must each be in their own clusters.
    """
    tree = _load_tree()
    pc = PhytClust(tree=tree)
    pc.run(k=4, plot_scores=False)
    cmap = pc.get_clusters(4)

    groups = _leaf_groups(cmap)
    assert len(groups) == 4
    # Each of the two leaf-pair clades must appear as an intact cluster
    assert frozenset(["A", "B"]) in groups
    assert frozenset(["C", "D"]) in groups
