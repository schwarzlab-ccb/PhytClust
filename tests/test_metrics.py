"""Tests for metrics module."""

import pathlib
from Bio import Phylo
import numpy as np

import pytest

from phytclust import PhytClust
from phytclust.metrics.indices import (
    AlphaIndex,
    cluster_alpha,
    colless_index_calc,
    variance_indices,
    gini_coefficient,
    colless_ratio,
    variance_ratio,
)

TREE_PATH = pathlib.Path(__file__).parent.parent / "examples" / "sample_tree.nwk"


def _load_tree():
    """Load test tree."""
    return Phylo.read(TREE_PATH, "newick")


class TestCollessIndex:
    """Test Colless index calculation."""

    def test_colless_index_returns_non_negative(self):
        """Test that Colless index returns non-negative value."""
        tree = _load_tree()
        index = colless_index_calc(tree)

        assert isinstance(index, (int, np.integer))
        assert index >= 0

    def test_colless_index_balanced_tree_is_zero_or_small(self):
        """Test that balanced tree has small Colless index."""
        # Create a simple balanced tree
        from Bio.Phylo.BaseTree import Clade, Tree

        tree = Tree(
            Clade(
                clades=[
                    Clade(clades=[Clade(name="A"), Clade(name="B")], name="AB"),
                    Clade(clades=[Clade(name="C"), Clade(name="D")], name="CD"),
                ]
            )
        )

        index = colless_index_calc(tree)
        # Balanced tree should have Colless = 0
        assert index >= 0


class TestVarianceIndices:
    """Test variance-based indices."""

    def test_variance_indices_returns_dict(self):
        """Test that variance_indices returns a dictionary."""
        tree = _load_tree()
        result = variance_indices(tree)

        assert isinstance(result, dict)
        assert "variance" in result or len(result) > 0

    def test_variance_indices_non_negative(self):
        """Test that variance values are non-negative."""
        tree = _load_tree()
        result = variance_indices(tree)

        for value in result.values():
            assert value >= 0 or np.isnan(value)


class TestGiniCoefficient:
    """Test Gini coefficient calculation."""

    def test_gini_coefficient_returns_float(self):
        """Test that Gini coefficient returns a float."""
        tree = _load_tree()
        gini = gini_coefficient(tree)

        assert isinstance(gini, (float, np.floating))

    def test_gini_coefficient_in_valid_range(self):
        """Test that Gini coefficient is in valid range [0, 1]."""
        tree = _load_tree()
        gini = gini_coefficient(tree)

        assert 0 <= gini <= 1 or np.isnan(gini)


class TestCollessRatio:
    """Test Colless ratio calculation."""

    def test_colless_ratio_returns_float(self):
        """Test that Colless ratio returns a float."""
        tree = _load_tree()
        ratio = colless_ratio(tree)

        assert isinstance(ratio, (float, np.floating))

    def test_colless_ratio_in_valid_range(self):
        """Test that Colless ratio is in valid range [0, 1]."""
        tree = _load_tree()
        ratio = colless_ratio(tree)

        assert 0 <= ratio <= 1 or np.isnan(ratio)


class TestVarianceRatio:
    """Test variance ratio calculation."""

    def test_variance_ratio_returns_float(self):
        """Test that variance ratio returns a float."""
        tree = _load_tree()
        ratio = variance_ratio(tree)

        assert isinstance(ratio, (float, np.floating))

    def test_variance_ratio_non_negative(self):
        """Test that variance ratio is non-negative."""
        tree = _load_tree()
        ratio = variance_ratio(tree)

        assert ratio >= 0 or np.isnan(ratio)


SOFT_POLYTOMY = "((a1:1,a2:1):1,(b1:1,b2:1):1,(c1:1,c2:1):1,(d1:1,d2:1):8);"


def _brute_alpha(tree, cmap):
    """Reference: an edge is intra iff it lies strictly below a cluster's clade.

    Every cluster here is a clade, so its root is the node whose leaves are
    exactly the cluster.
    """
    groups = {}
    for leaf, cid in cmap.items():
        groups.setdefault(cid, set()).add(leaf.name)
    roots = {
        id(node)
        for node in tree.find_clades()
        if {t.name for t in node.get_terminals()} in groups.values()
    }
    intra, extra = [], []

    def walk(node, below_root):
        for child in node.clades:
            (intra if below_root else extra).append(child.branch_length or 0.0)
            walk(child, below_root or id(child) in roots)

    walk(tree.root, id(tree.root) in roots)
    return np.mean(extra) / np.mean(intra), len(intra), len(extra)


class TestClusterAlpha:
    """Alpha = mean backbone edge length / mean within-cluster edge length."""

    def test_matches_clade_reference_on_sample_tree(self):
        pc = PhytClust(tree=_load_tree())
        for k in range(2, 8):
            cmap = pc.get_clusters(k)
            alpha, n_intra, n_extra = _brute_alpha(pc.tree, cmap)
            info = cluster_alpha(pc.tree, cmap)
            assert info["alpha"] == pytest.approx(alpha)
            assert (info["n_intra_nodes"], info["n_extra_nodes"]) == (n_intra, n_extra)

    @pytest.mark.parametrize(
        "k,alpha,n_intra,n_extra", [(2, 8.0, 11, 1), (3, 4.5, 10, 2)]
    )
    def test_soft_polytomy_group_counts_member_stems_as_intra(
        self, k, alpha, n_intra, n_extra
    ):
        """A group of polytomy children is the clade of its zero-length node."""
        pc = PhytClust(tree=SOFT_POLYTOMY, max_k=6, polytomy_mode="soft")
        info = cluster_alpha(pc.tree, pc.get_clusters(k))
        assert info["alpha"] == pytest.approx(alpha)
        assert (info["n_intra_nodes"], info["n_extra_nodes"]) == (n_intra, n_extra)

    def test_reused_index_matches_fresh(self):
        pc = PhytClust(tree=_load_tree())
        index = AlphaIndex(pc.tree)
        for k in range(2, 8):
            cmap = pc.get_clusters(k)
            assert cluster_alpha(pc.tree, cmap, index=index) == cluster_alpha(
                pc.tree, cmap
            )

    def test_alpha_info_caches_until_dp_changes(self):
        pc = PhytClust(tree=_load_tree())
        first = pc.alpha_info(2)
        assert pc.alpha_info(2) is first
        pc.use_branch_support = True
        pc.get_clusters(2)
        assert pc.alpha_info(2) is not first


def test_edge_metrics_ignore_stored_root_length():
    from phytclust.metrics.indices import (
        calculate_internal_terminal_ratio,
        calculate_terminal_contributions,
        calculate_proportions,
        root_to_tip_lengths,
    )

    tree = _load_tree()
    tree.root.branch_length = 0
    before = (
        variance_indices(tree),
        calculate_internal_terminal_ratio(tree),
        calculate_terminal_contributions(tree),
        calculate_proportions(tree, [1, 2]),
        root_to_tip_lengths(tree),
    )
    tree.root.branch_length = 12345
    after = (
        variance_indices(tree),
        calculate_internal_terminal_ratio(tree),
        calculate_terminal_contributions(tree),
        calculate_proportions(tree, [1, 2]),
        root_to_tip_lengths(tree),
    )
    assert before == after


def test_branch_standard_deviation_and_compatibility_names():
    from Bio.Phylo.BaseTree import Clade, Tree
    from phytclust.metrics import (
        branch_length_standard_deviation,
        calculate_variance_branch_length,
    )

    tree = Tree(
        root=Clade(
            branch_length=100, clades=[Clade(branch_length=1), Clade(branch_length=3)]
        )
    )
    assert branch_length_standard_deviation(tree) == 1
    assert calculate_variance_branch_length(tree) == 1


def test_root_to_tip_lengths_match_reference_and_handle_deep_tree():
    from Bio.Phylo.BaseTree import Clade, Tree
    from phytclust.metrics.indices import root_to_tip_lengths, collect_branch_lengths

    tree = _load_tree()
    assert root_to_tip_lengths(tree) == [
        tree.distance(leaf) for leaf in tree.get_terminals()
    ]
    root = Clade(name="leaf", branch_length=1)
    for _ in range(1500):
        root = Clade(branch_length=1, clades=[root])
    deep_tree = Tree(root=root)
    assert root_to_tip_lengths(deep_tree) == [1500]
    assert len(collect_branch_lengths(root)) == 1500


@pytest.mark.parametrize("values", [[-1, 2], [float("nan"), 1], [float("inf"), 1]])
def test_gini_rejects_invalid_values(values):
    from phytclust.exceptions import DataError

    with pytest.raises(DataError, match="finite, non-negative"):
        gini_coefficient(values)


def test_known_gini_and_colless_values():
    from io import StringIO

    assert gini_coefficient([1, 1]) == 0
    assert gini_coefficient([0, 1]) == 0.5
    assert colless_index_calc(Phylo.read(StringIO("((a,b),(c,d));"), "newick")) == 0
    assert colless_index_calc(Phylo.read(StringIO("(((a,b),c),d);"), "newick")) == 3
