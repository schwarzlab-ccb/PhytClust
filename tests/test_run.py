import pytest
from Bio.Phylo.BaseTree import Clade

from conftest import (
    POLYTOMY_TREE,
    SAMPLE_TREE,
    leaf_groups,
    load_tree,
    random_newick,
    tree_from_newick,
)
from phytclust import PeakConfig, PhytClust
from phytclust.exceptions import ConfigurationError, InvalidKError


def _assert_cluster_map_valid(tree, cluster_map, k_expected: int):
    assert isinstance(cluster_map, dict)
    assert len(cluster_map) > 0

    assert all(isinstance(t, Clade) and t.is_terminal() for t in cluster_map.keys())
    assert all(isinstance(v, int) for v in cluster_map.values())

    cluster_ids = set(cluster_map.values())
    assert len(cluster_ids) == k_expected

    leaves = list(tree.get_terminals())
    assert set(cluster_map.keys()) == set(leaves)


def test_run_exact_k_partition_returns_result_dict():
    tree = load_tree()
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
    cluster_map = result["clusters"][0]
    _assert_cluster_map_valid(tree, cluster_map, k_expected=2)


def test_get_clusters_exact_k_partition():
    tree = load_tree()
    pc = PhytClust(tree=tree)

    pc.run(k=2, plot_scores=False)

    cluster_map = pc.get_clusters(2)
    _assert_cluster_map_valid(tree, cluster_map, k_expected=2)

    assert 2 in pc.cluster_maps
    assert pc.cluster_maps[2] == cluster_map


def test_run_global_mode_returns_list_of_cluster_maps():
    tree = load_tree()
    pc = PhytClust(tree=tree)

    result = pc.run(top_n=1, max_k=8, plot_scores=False)

    assert result["mode"] == "global"
    assert result["k_values"] == [3]
    assert result["selected_k"] == 3
    assert result["ks"] == result["k_values"]
    assert result["peaks"] == result["k_values"]
    assert len(result["clusters"]) == 1
    _assert_cluster_map_valid(tree, result["clusters"][0], k_expected=3)


def test_run_resolution_mode_falls_back_for_small_trees(caplog):
    tree = load_tree()
    pc = PhytClust(tree=tree)

    with caplog.at_level("WARNING", logger="phytclust"):
        result = pc.run(by_resolution=True, num_bins=2, plot_scores=False)

    # The fallback is a global search, so the result says so, and warns.
    assert result["mode"] == "global"
    assert pc.resolution_fell_back
    assert "Resolution mode needs at least 50" in caplog.text
    assert isinstance(result["clusters"], list)
    assert isinstance(result["k_values"], list)
    if result["k_values"]:
        assert result["selected_k"] == result["k_values"][0]
    assert isinstance(result["ks"], list)
    assert isinstance(result["peaks"], list)


def test_best_global_returns_list_of_cluster_maps():
    tree = load_tree()
    pc = PhytClust(tree=tree)

    max_k = min(10, pc.num_terminals)
    if max_k < 4:
        pytest.skip("Test tree too small for global peak search (requires max_k>=4).")

    clusters = pc.best_global(top_n=1, max_k=max_k, plot_scores=False)

    assert isinstance(clusters, list)
    if clusters:
        cluster_map = clusters[0]
        assert isinstance(cluster_map, dict)
        assert set(cluster_map.keys()) == set(tree.get_terminals())


def test_run_rejects_invalid_argument_combinations():
    tree = load_tree()
    pc = PhytClust(tree=tree)

    with pytest.raises(
        ConfigurationError, match="Cannot combine `k` with `by_resolution=True`"
    ):
        pc.run(k=2, by_resolution=True, plot_scores=False)

    with pytest.raises(
        ConfigurationError, match="top_n must be 1 when k is specified"
    ):
        pc.run(k=2, top_n=2, plot_scores=False)

    with pytest.raises(InvalidKError, match="k must be"):
        pc.run(k=0, plot_scores=False)


def test_resolution_after_fixed_k_does_not_inherit_self_k():
    """run(k=...) must not affect later runs."""
    tree = load_tree()
    pc = PhytClust(tree=tree)

    fresh = PhytClust(tree=load_tree()).run(top_n=3, plot_scores=False)
    pc.run(k=2, plot_scores=False)
    assert pc.k is None
    assert pc.exact_k == 2
    assert pc.run(top_n=3, plot_scores=False)["k_values"] == fresh["k_values"]
    assert pc.exact_k is None
    assert pc.run(plot_scores=False)["mode"] == "global"

    result = pc.run(by_resolution=True, num_bins=2, plot_scores=False)
    # Too few scores for resolution mode on this tree, so it runs globally.
    assert result["mode"] == "global"

    with pytest.raises(InvalidKError, match="`top_n` must be"):
        pc.run(top_n=0, plot_scores=False)


def test_k2_partition_is_root_split():
    """
    sample_tree.nwk has root children A13={A,B,C,D} and {E,F,G,H}.
    k=2 must split exactly there, verified against expected leaf sets.

    Tree: (((A:5,B:3):6,(C:3,D:7):4):22,(((E:7,F:13):5,G:6):10,H:60):35):0;
    """
    tree = load_tree()
    pc = PhytClust(tree=tree)
    pc.run(k=2, plot_scores=False)
    cluster_map = pc.get_clusters(2)

    groups = leaf_groups(cluster_map)
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
    tree = load_tree()
    pc = PhytClust(tree=tree)
    pc.run(k=4, plot_scores=False)
    cluster_map = pc.get_clusters(4)

    groups = leaf_groups(cluster_map)
    assert len(groups) == 4
    # Each of the two leaf-pair clades must appear as an intact cluster
    assert frozenset(["A", "B"]) in groups
    assert frozenset(["C", "D"]) in groups


def test_max_k_limit_is_honoured_and_not_sticky():
    pc = PhytClust(tree=tree_from_newick(SAMPLE_TREE.read_text()))
    pc.run(top_n=1, max_k_limit=0.5, plot_scores=False)
    assert pc.max_k == 4
    pc.run(top_n=1, max_k=4, plot_scores=False)
    pc.run(top_n=1, plot_scores=False)
    assert pc.max_k == 8


def test_explicit_max_k_on_the_instance_is_kept():
    pc = PhytClust(tree=tree_from_newick(SAMPLE_TREE.read_text()), max_k=5)
    pc.run(top_n=1, plot_scores=False)
    pc.run(top_n=1, plot_scores=False)
    assert pc.max_k == 5


def test_resolution_bins_start_at_one():
    """Clade-level bins span 1..K, so k=3 and the next coarse level separate."""
    pc = PhytClust(tree=tree_from_newick(random_newick(7, 300)))
    pc.run(by_resolution=True, num_bins=5, plot_scores=False)
    bins = pc.bin_ranges_current
    assert bins[0][0] == 1
    assert bins[0][1] <= 4
    for (_, hi), (lo, _) in zip(bins, bins[1:]):
        assert lo == hi + 1


def test_resolution_fallback_threshold_is_configurable():
    pc = PhytClust(POLYTOMY_TREE)
    fallback = pc.run(by_resolution=True, plot_scores=False)
    assert fallback["mode"] == "global"
    lowered = pc.run(
        by_resolution=True,
        plot_scores=False,
        peak_config=PeakConfig(resolution_min_scores=10),
    )
    assert lowered["mode"] == "resolution"
    assert not pc.resolution_fell_back
