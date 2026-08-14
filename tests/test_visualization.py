"""Tests for visualization modules."""

import pathlib
import pytest
import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo

from phytclust.algo.core import PhytClust
from phytclust.exceptions import InvalidClusteringError, ConfigurationError
from phytclust.viz.scores import plot_scores
from phytclust.viz.cluster import plot_clusters

TREE_PATH = pathlib.Path(__file__).parent.parent / "examples" / "sample_tree.nwk"


def _load_tree():
    """Load test tree."""
    return Phylo.read(TREE_PATH, "newick")


class TestScoresPlot:
    """Test scores plotting function."""

    def test_plot_scores_returns_figure(self):
        """Test that plot_scores returns a matplotlib figure."""
        tree = _load_tree()
        pc = PhytClust(tree=tree)
        pc.run(top_n=1, max_k=10, plot_scores=False)

        fig = plot_scores(pc)

        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_scores_with_peaks(self):
        """Test plot_scores with peak annotations."""
        tree = _load_tree()
        pc = PhytClust(tree=tree)
        pc.run(top_n=2, max_k=10, plot_scores=False)

        peaks = [2, 3]
        fig = plot_scores(pc, peaks=peaks)

        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_scores_with_custom_range(self):
        """Test plot_scores with custom k range."""
        tree = _load_tree()
        pc = PhytClust(tree=tree)
        pc.run(top_n=1, max_k=10, plot_scores=False)

        fig = plot_scores(pc, k_start=2, k_end=8)

        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_scores_empty_scores_raises_error(self):
        """Test that plot_scores raises error with empty scores."""
        tree = _load_tree()
        pc = PhytClust(tree=tree)
        pc.scores = np.array([])

        with pytest.raises(InvalidClusteringError, match="Scores are empty"):
            plot_scores(pc)

    def test_plot_scores_invalid_axis_mode_raises_error(self):
        """Test that invalid x_axis_mode raises error."""
        tree = _load_tree()
        pc = PhytClust(tree=tree)
        pc.run(top_n=1, max_k=10, plot_scores=False)

        with pytest.raises(ConfigurationError, match="x_axis_mode must be"):
            plot_scores(pc, x_axis_mode="invalid")


class TestClusterPlot:
    """Test cluster plotting function."""

    def test_plot_clusters_returns_or_saves(self):
        """Test that plot_clusters completes successfully."""
        tree = _load_tree()
        pc = PhytClust(tree=tree)
        pc.run(k=2, plot_scores=False)

        # Should not raise
        plot_clusters(pc, save=False)

    def test_plot_clusters_with_exact_k(self):
        """Test plotting clusters with specific k."""
        tree = _load_tree()
        pc = PhytClust(tree=tree)
        pc.run(k=2, plot_scores=False)

        plot_clusters(pc, save=False)

    def test_plot_clusters_with_multiple_peaks(self):
        """Test plotting multiple cluster solutions."""
        tree = _load_tree()
        pc = PhytClust(tree=tree)
        pc.run(top_n=2, max_k=10, plot_scores=False)

        plot_clusters(pc, top_n=2, save=False)

    def test_plot_clusters_with_custom_cmap(self):
        """Test plot_clusters with custom colormap."""
        tree = _load_tree()
        pc = PhytClust(tree=tree)
        pc.run(k=2, plot_scores=False)

        # plt.get_cmap, not matplotlib.cm.get_cmap, which 3.9 removed.
        plot_clusters(pc, cmap=plt.get_cmap("viridis"), save=False)

    def test_plot_clusters_with_no_clusters_handled(self):
        """Test that plot_clusters handles missing clusters gracefully."""
        tree = _load_tree()
        pc = PhytClust(tree=tree)

        # Don't run() to leave pc without clusters
        plot_clusters(pc, save=False)


class TestVisualizationEdgeCases:
    """Test edge cases in visualization."""

    def test_plot_scores_with_log_scale(self):
        """Test plot_scores with log scale."""
        tree = _load_tree()
        pc = PhytClust(tree=tree)
        pc.run(top_n=1, max_k=10, plot_scores=False)

        fig = plot_scores(pc, log_scale_y=True)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_scores_with_resolution_bins(self):
        """Test plot_scores with resolution bins."""
        tree = _load_tree()
        pc = PhytClust(tree=tree)
        pc.run(top_n=1, max_k=10, plot_scores=False)

        fig = plot_scores(pc, resolution_on=True, num_bins=3)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_scores_linear_scale(self):
        """Test plot_scores with linear x-axis."""
        tree = _load_tree()
        pc = PhytClust(tree=tree)
        pc.run(top_n=1, max_k=10, plot_scores=False)

        fig = plot_scores(pc, x_axis_mode="linear")
        assert isinstance(fig, plt.Figure)
        plt.close(fig)
