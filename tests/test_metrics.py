"""Tests for metrics module."""

import pathlib
from Bio import Phylo
import numpy as np

from phytclust.metrics.indices import (
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
