from unittest.mock import patch

import matplotlib.pyplot as plt
import pytest

from phytclust import PhytClust
from phytclust.exceptions import ConfigurationError
from phytclust.viz.cluster import plot_clusters, plot_multiple_k


@pytest.fixture
def clustering():
    return PhytClust("((a:1,b:1):1,(c:1,d:1):1);", max_k=4)


@pytest.mark.parametrize("count", [True, 1.5, 0, -1])
def test_invalid_counts(clustering, count):
    with pytest.raises(ConfigurationError, match="positive integer"):
        plot_clusters(clustering, k=count)


def test_multiple_counts_prepare_once_and_keep_all_files(clustering, tmp_path):
    with patch.object(clustering, "_ensure_dp", wraps=clustering._ensure_dp) as prepare:
        plot_multiple_k(
            clustering, [2, 3, 2], results_dir=tmp_path, filename="clusters.png"
        )
    assert prepare.call_count == 1
    assert {path.name for path in tmp_path.iterdir()} == {
        "clusters_k2.png",
        "clusters_k3.png",
    }


def test_ranked_counts_keep_all_files(clustering, tmp_path):
    clustering.peaks_by_rank = [2, 3]
    plot_clusters(clustering, top_n=2, results_dir=tmp_path, filename="clusters.svg")
    assert {path.name for path in tmp_path.iterdir()} == {
        "clusters_k2.svg",
        "clusters_k3.svg",
    }


def test_label_spacing_and_style_do_not_change_global_settings(clustering):
    original_font_size = plt.rcParams["font.size"]
    with patch("matplotlib.pyplot.show"):
        plot_clusters(clustering, k=2, show_terminal_labels=True)
    figure = plt.gcf()
    try:
        assert figure.get_size_inches()[1] >= 2.5
        assert figure.axes[0].get_title(loc="left") == "PhytClust’s clusters at k = 2"
        assert plt.rcParams["font.size"] == original_font_size
    finally:
        plt.close(figure)


def test_clustering_failure_is_not_hidden(clustering):
    with patch.object(clustering, "_ensure_dp", side_effect=RuntimeError("failed")):
        with pytest.raises(RuntimeError, match="failed"):
            plot_multiple_k(clustering, [2])


def test_failed_save_closes_figure(clustering, tmp_path):
    existing = set(plt.get_fignums())
    with patch("matplotlib.figure.Figure.savefig", side_effect=OSError("failed")):
        with pytest.raises(OSError, match="failed"):
            plot_clusters(clustering, k=2, results_dir=tmp_path)
    assert set(plt.get_fignums()) == existing


def test_outgroup_plot_uses_partition_tree(tmp_path):
    clustering = PhytClust("((a:1,b:1):1,(c:1,d:1):1,out:3);", outgroup="out", max_k=4)
    with patch(
        "phytclust.viz.cluster.plot_cluster",
        wraps=__import__("phytclust.viz.plots", fromlist=["plot_cluster"]).plot_cluster,
    ) as draw:
        plot_clusters(
            clustering, k=2, results_dir=tmp_path, colour_branches_by_cluster=True
        )
    assert draw.call_args.kwargs["tree"] is clustering._tree_wo_outgroup
    assert set(draw.call_args.kwargs["cluster"]).issubset(
        set(clustering._tree_wo_outgroup.get_terminals())
    )


def test_supplied_axis_is_styled_without_resizing(clustering):
    figure, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].set_title("Other plot")
    with patch("matplotlib.pyplot.show"):
        plot_clusters(clustering, k=2, ax=axes[1])
    try:
        assert axes[0].get_title() == "Other plot"
        assert axes[1].get_title(loc="left") == "PhytClust’s clusters at k = 2"
        assert tuple(figure.get_size_inches()) == (10, 4)
    finally:
        plt.close(figure)


@pytest.mark.parametrize("title", ["My clusters", ""])
@pytest.mark.parametrize("plot", [plot_clusters, plot_multiple_k])
def test_custom_title(clustering, title, plot):
    options = {"k": 2} if plot is plot_clusters else {"k_values": [2]}
    with patch("matplotlib.pyplot.show"):
        plot(clustering, title=title, **options)
    figure = plt.gcf()
    try:
        assert figure.axes[0].get_title(loc="left") == title
        assert figure.axes[0].get_title(loc="center") == ""
    finally:
        plt.close(figure)
