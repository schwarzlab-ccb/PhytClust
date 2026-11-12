from unittest.mock import patch
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from Bio.Phylo.BaseTree import Clade, Tree

from conftest import POLYTOMY_TREE, SAMPLE_TREE
from phytclust import ClusterPlotConfig, PhytClust, PlotConfig, RuntimeConfig
from phytclust.exceptions import ConfigurationError, ValidationError
from phytclust.viz import plots
from phytclust.viz.cluster import plot_clusters, plot_multiple_k


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


@pytest.fixture
def clustering():
    return PhytClust("((a:0.023,b:0.031):0.08,(c:0.027,d:0.04):0.1);", max_k=4)


def test_zero_marker_size_and_line_width(clustering):
    figure = plots.plot_tree(
        clustering.tree, marker_size=0, line_width=0, show_branch_lengths=False
    )
    assert len(figure.axes[0].collections) == 2
    assert all(
        np.equal(collection.get_linewidths(), 0).all()
        for collection in figure.axes[0].collections
    )


def test_branch_labels_include_unnamed_lengths_and_zero_support():
    node = Clade(branch_length=0.023, confidence=0)
    label = plots._branch_label_formatter(None, True, True)
    assert label(node) == "0.023 / 0"
    assert plots._branch_label_formatter({}, True)(node) is None
    assert plots._branch_label_formatter(lambda _: "custom", True)(node) == "custom"


def test_root_incoming_length_is_not_labelled():
    tree = Tree(
        root=Clade(branch_length=0.123, clades=[Clade(name="a", branch_length=0.023)])
    )
    figure = plots.plot_tree(tree)
    labels = [text.get_text() for text in figure.axes[0].texts]
    assert "0.023" in labels
    assert "0.123" not in labels


def test_coordinates_computed_once_with_both_overlays(clustering):
    partition = clustering.get_clusters(2)
    with (
        patch.object(
            plots, "_get_x_positions", wraps=plots._get_x_positions
        ) as x_positions,
        patch.object(
            plots, "_get_y_positions", wraps=plots._get_y_positions
        ) as y_positions,
    ):
        plots.plot_cluster(
            partition, clustering.tree, show_cluster_bars=True, show_cluster_boxes=True
        )
    assert x_positions.call_count == y_positions.call_count == 1


def test_color_callback_error_reaches_caller(clustering):
    def fail(_):
        raise RuntimeError("colour failed")

    with pytest.raises(RuntimeError, match="colour failed"):
        plots.plot_tree(clustering.tree, branch_colour_func=fail)


def test_formatted_names_keep_cluster_colors(clustering):
    partition = clustering.get_clusters(2)
    figure = plots.plot_cluster(
        partition,
        clustering.tree,
        label_func=lambda node: node.name.upper() if node.name else None,
    )
    colors = [text.get_color() for text in figure.axes[0].texts]
    assert len(set(colors)) == 2
    assert all(color != "black" for color in colors)


def test_save_uses_supplied_figure(clustering, tmp_path):
    figure, axis = plt.subplots()
    other = plt.figure()
    with (
        patch.object(figure, "savefig") as save,
        patch.object(other, "savefig") as other_save,
    ):
        plots.plot_tree(clustering.tree, ax=axis, output_name=str(tmp_path / "tree"))
    save.assert_called_once()
    other_save.assert_not_called()


def test_plot_options_apply_to_supplied_axis(clustering):
    figure, axis = plt.subplots()
    other_axis = plt.figure().add_subplot()
    plots.plot_tree(clustering.tree, ax=axis, xlabel="Custom branch axis")
    assert axis.get_xlabel() == "Custom branch axis"
    assert other_axis.get_xlabel() == ""


def test_outgroup_uses_literal_name_and_appears_at_top():
    tree = Tree(root=Clade(clades=[Clade(name="a"), Clade(name="out.group")]))
    positions = plots._get_y_positions(tree, outgroup="out.group")
    assert positions[tree.root.clades[1]] == 0
    assert positions[tree.root.clades[0]] == 1


def test_deep_tree_geometry_and_drawing_do_not_recurse():
    root = Clade(name="tip", branch_length=1)
    for _ in range(1200):
        root = Clade(branch_length=0.001, clades=[root])
    figure = plots.plot_tree(
        Tree(root=root),
        show_terminal_labels=False,
        show_branch_lengths=False,
        height_scale=0.001,
    )
    assert len(figure.axes[0].collections[0].get_segments()) == 1201


def test_scattered_cluster_ancestor_and_invalid_members(clustering):
    leaves = clustering.tree.get_terminals()
    ancestors = plots._cluster_common_ancestors(
        clustering.tree, {0: [leaves[0], leaves[2]]}
    )
    assert ancestors[0] is clustering.tree.root
    with pytest.raises(ValidationError, match="plotted tree"):
        plots._cluster_common_ancestors(clustering.tree, {0: [Clade(name="foreign")]})


def test_empty_custom_palette_is_rejected(clustering):
    with pytest.raises(ValidationError, match="palette"):
        plots.plot_cluster(clustering.get_clusters(2), clustering.tree, palette=[])


def test_short_custom_palette_cycles_with_warning(clustering):
    with pytest.warns(UserWarning, match="colours repeat"):
        plots.plot_cluster(clustering.get_clusters(2), clustering.tree, palette=["red"])


def test_custom_lower_level_title(clustering):
    figure = plots.plot_cluster(
        clustering.get_clusters(2), clustering.tree, title="My clusters"
    )
    assert figure.axes[0].get_title(loc="left") == "My clusters"


def test_comparison_rows_match_nonuniform_tree_positions(clustering):
    table = pd.DataFrame(
        {"k = 2": [0, 0, 1, 1], "k = 3": [0, 1, 2, 2]}, index=["a", "b", "c", "d"]
    )
    figure = plots.plot_multiple_clusters(
        table,
        final_tree=clustering.tree,
        show_internal_nodes=True,
        fixed_x_range=(10000, 50000),
        cbar_width_ratio=0,
    )
    tree_axis, heat_axis = figure.axes
    positions = plots._get_y_positions(clustering.tree, adjust=True)
    expected = sorted(positions[leaf] for leaf in clustering.tree.get_terminals())
    cells = heat_axis.collections[0].get_paths()
    for path, coordinate in zip(cells[::2], expected):
        assert path.vertices[:, 1].min() < coordinate < path.vertices[:, 1].max()
    np.testing.assert_allclose(
        [text.get_position()[1] for text in tree_axis.texts], expected
    )
    np.testing.assert_allclose(heat_axis.get_xticks(), [20000, 40000])
    np.testing.assert_allclose(tree_axis.get_ylim(), heat_axis.get_ylim())


def test_comparison_accepts_long_format_and_preserves_names():
    index = pd.MultiIndex.from_tuples(
        [("a", "first"), ("b", "first"), ("a", "second"), ("b", "second")],
        names=["leaf_name", "comparison_IDs"],
    )
    table = pd.DataFrame({"cluster_ID": [0, 1, 0, 0]}, index=index)
    figure = plots.plot_multiple_clusters(table, cbar_width_ratio=0)
    assert [label.get_text() for label in figure.axes[0].get_xticklabels()] == [
        "first",
        "second",
    ]


@pytest.mark.parametrize("values", [[0, np.nan], [0, 1.5]])
def test_invalid_heatmap_assignments(values):
    with pytest.raises(ValidationError, match="integer"):
        plots.plot_multiple_clusters(pd.DataFrame({"k = 2": values}, index=["a", "b"]))


def test_missing_comparison_leaf_is_rejected(clustering):
    with pytest.raises(ValidationError, match="exactly"):
        plots.plot_multiple_clusters(
            pd.DataFrame({"k = 2": [0]}, index=["a"]), final_tree=clustering.tree
        )


def test_bar_placement_respects_formatted_label_width(clustering):
    figure = plots.plot_cluster(
        clustering.get_clusters(2),
        clustering.tree,
        show_cluster_bars=True,
        label_func=lambda node: (
            (node.name + " extra label text") if node.name else None
        ),
    )
    figure.canvas.draw()
    axis = figure.axes[0]
    renderer = figure.canvas.get_renderer()
    label_right = max(text.get_window_extent(renderer).x1 for text in axis.texts)
    bar_left = min(
        axis.transData.transform((patch.get_x(), 0))[0] for patch in axis.patches
    )
    assert bar_left > label_right


def test_unnamed_internal_nodes_can_have_markers(clustering):
    figure = plots.plot_tree(
        clustering.tree, hide_internal_nodes=False, show_branch_lengths=False
    )
    marker_collection = figure.axes[0].collections[-1]
    assert len(marker_collection.get_offsets()) == 7


def test_topology_fallback_is_labelled_honestly():
    tree = Tree(
        root=Clade(
            clades=[Clade(name="a", branch_length=0), Clade(name="b", branch_length=0)]
        )
    )
    figure = plots.plot_tree(tree)
    assert figure.axes[0].get_xlabel() == "Topology depth"


def test_cluster_wrapper_preserves_custom_axis_label(clustering):
    from phytclust.viz.cluster import plot_clusters

    with patch("matplotlib.pyplot.show"):
        plot_clusters(clustering, k=2, xlabel="Custom axis")
    assert plt.gcf().axes[0].get_xlabel() == "Custom axis"


def test_box_labels_do_not_overlap_leaf_labels(clustering):
    figure = plots.plot_cluster(
        clustering.get_clusters(2), clustering.tree, show_cluster_boxes=True
    )
    figure.canvas.draw()
    axis = figure.axes[0]
    renderer = figure.canvas.get_renderer()
    leaf_labels = [text for text in axis.texts if text.get_text().startswith(" ")]
    cluster_labels = [text for text in axis.texts if text.get_text().startswith("C")]
    for cluster_label in cluster_labels:
        cluster_bounds = cluster_label.get_window_extent(renderer)
        assert all(
            not cluster_bounds.overlaps(text.get_window_extent(renderer))
            for text in leaf_labels
        )


def test_cladogram_disables_branch_axis(clustering):
    figure = plots.plot_cluster(
        clustering.get_clusters(2),
        clustering.tree,
        layout="cladogram",
        show_cluster_bars=True,
    )
    assert not figure.axes[0].xaxis.get_visible()
    assert len(figure.axes[0].patches) == 2


def test_heatmap_default_height_keeps_large_tree_labels_apart():
    table = pd.DataFrame(
        {"k = 3": [index % 3 for index in range(177)]},
        index=[f"Taxon_{index:03d}" for index in range(177)],
    )
    figure = plots.plot_multiple_clusters(table, cbar_width_ratio=0)
    assert figure.get_size_inches()[1] == pytest.approx(43.48)
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    bounds = sorted(
        (
            label.get_window_extent(renderer)
            for label in figure.axes[0].get_yticklabels()
        ),
        key=lambda box: box.y0,
    )
    assert all(not first.overlaps(second) for first, second in zip(bounds, bounds[1:]))


def test_heatmap_respects_explicit_size():
    table = pd.DataFrame({"k = 2": [0, 1]}, index=["a", "b"])
    figure = plots.plot_multiple_clusters(table, figsize=(12, 4))
    np.testing.assert_allclose(figure.get_size_inches(), [12, 4])


def test_unlabelled_large_tree_has_compact_height():
    tree = Tree(
        root=Clade(
            clades=[
                Clade(name=f"Taxon_{index}", branch_length=1) for index in range(1103)
            ]
        )
    )
    figure = plots.plot_tree(
        tree, show_terminal_labels=False, show_branch_lengths=False, marker_size=0
    )
    assert figure.get_size_inches()[1] == 18


def test_heatmap_reserves_space_for_long_names():
    names = [
        "Long_taxon_name_with_many_characters_a",
        "Long_taxon_name_with_many_characters_b",
    ]
    tree = Tree(
        root=Clade(clades=[Clade(name=name, branch_length=1) for name in names])
    )
    figure = plots.plot_multiple_clusters(
        pd.DataFrame({"k = 2": [0, 1]}, index=names),
        final_tree=tree,
        cbar_width_ratio=0,
    )
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    tree_axis, heat_axis = figure.axes
    heat_left = heat_axis.get_window_extent(renderer).x0
    assert len(tree_axis.texts) == len(names)
    assert not heat_axis.get_yticklabels()
    assert all(
        label.get_window_extent(renderer).x1 < heat_left for label in tree_axis.texts
    )


def test_identical_assignments_form_a_continuous_color_band():
    figure = plots.plot_multiple_clusters(
        pd.DataFrame({"k = 1": [0, 0, 0]}, index=["a", "b", "c"]), cbar_width_ratio=0
    )
    figure.canvas.draw()
    axis = figure.axes[0]
    image = np.asarray(figure.canvas.buffer_rgba())
    pixels = []
    for y_coordinate in [0, 0.5, 1, 1.5, 2]:
        x_pixel, y_pixel = axis.transData.transform((0.5, y_coordinate))
        pixels.append(image[image.shape[0] - 1 - int(y_pixel), int(x_pixel)])
    assert all(np.array_equal(pixel, pixels[0]) for pixel in pixels[1:])


def test_comparison_names_stay_beside_tree_tips(clustering):
    table = pd.DataFrame({"k = 2": [0, 0, 1, 1]}, index=["a", "b", "c", "d"])
    figure = plots.plot_multiple_clusters(
        table, final_tree=clustering.tree, cbar_width_ratio=0
    )
    figure.canvas.draw()
    axis = figure.axes[0]
    leaf_by_name = {leaf.name: leaf for leaf in clustering.tree.get_terminals()}
    positions = plots._get_x_positions(clustering.tree)
    for text in axis.texts:
        leaf = leaf_by_name[text.get_text().strip()]
        tip_x = axis.transData.transform((positions[leaf], text.get_position()[1]))[0]
        text_x = axis.transData.transform(text.get_position())[0]
        assert 0 < text_x - tip_x < 15


def test_comparison_separators_are_only_between_solutions():
    table = pd.DataFrame(
        {"k = 2": [0, 1], "k = 3": [0, 1], "k = 4": [0, 1]}, index=["a", "b"]
    )
    figure = plots.plot_multiple_clusters(
        table, fixed_x_range=(10, 40), cbar_width_ratio=0
    )
    lines = figure.axes[0].lines
    assert len(lines) == 2
    np.testing.assert_allclose(
        [line.get_xdata() for line in lines], [[20, 20], [30, 30]]
    )


@pytest.fixture
def simple_clustering():
    return PhytClust("((a:1,b:1):1,(c:1,d:1):1);", max_k=4)


@pytest.mark.parametrize("count", [True, 1.5, 0, -1])
def test_invalid_counts(simple_clustering, count):
    with pytest.raises(ConfigurationError, match="positive integer"):
        plot_clusters(simple_clustering, k=count)


def test_multiple_counts_prepare_once_and_keep_all_files(simple_clustering, tmp_path):
    with patch.object(simple_clustering, "_ensure_dp", wraps=simple_clustering._ensure_dp) as prepare:
        plot_multiple_k(
            simple_clustering, [2, 3, 2], results_dir=tmp_path, filename="clusters.png"
        )
    assert prepare.call_count == 1
    assert {path.name for path in tmp_path.iterdir()} == {
        "clusters_k2.png",
        "clusters_k3.png",
    }


def test_ranked_counts_keep_all_files(simple_clustering, tmp_path):
    simple_clustering.peaks_by_rank = [2, 3]
    plot_clusters(simple_clustering, top_n=2, results_dir=tmp_path, filename="clusters.svg")
    assert {path.name for path in tmp_path.iterdir()} == {
        "clusters_k2.svg",
        "clusters_k3.svg",
    }


def test_label_spacing_and_style_do_not_change_global_settings(simple_clustering):
    original_font_size = plt.rcParams["font.size"]
    with patch("matplotlib.pyplot.show"):
        plot_clusters(simple_clustering, k=2, show_terminal_labels=True)
    figure = plt.gcf()
    try:
        assert figure.get_size_inches()[1] >= 2.5
        assert figure.axes[0].get_title(loc="left") == "PhytClust’s clusters at k = 2"
        assert plt.rcParams["font.size"] == original_font_size
    finally:
        plt.close(figure)


def test_clustering_failure_is_not_hidden(simple_clustering):
    with patch.object(simple_clustering, "_ensure_dp", side_effect=RuntimeError("failed")):
        with pytest.raises(RuntimeError, match="failed"):
            plot_multiple_k(simple_clustering, [2])


def test_failed_save_closes_figure(simple_clustering, tmp_path):
    existing = set(plt.get_fignums())
    with patch("matplotlib.figure.Figure.savefig", side_effect=OSError("failed")):
        with pytest.raises(OSError, match="failed"):
            plot_clusters(simple_clustering, k=2, results_dir=tmp_path)
    assert set(plt.get_fignums()) == existing


def test_outgroup_plot_uses_partition_tree(tmp_path):
    simple_clustering = PhytClust("((a:1,b:1):1,(c:1,d:1):1,out:3);", outgroup="out", max_k=4)
    with patch(
        "phytclust.viz.cluster.plot_cluster",
        wraps=__import__("phytclust.viz.plots", fromlist=["plot_cluster"]).plot_cluster,
    ) as draw:
        plot_clusters(
            simple_clustering, k=2, results_dir=tmp_path, colour_branches_by_cluster=True
        )
    assert draw.call_args.kwargs["tree"] is simple_clustering._tree_wo_outgroup
    assert set(draw.call_args.kwargs["cluster"]).issubset(
        set(simple_clustering._tree_wo_outgroup.get_terminals())
    )


def test_supplied_axis_is_styled_without_resizing(simple_clustering):
    figure, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].set_title("Other plot")
    with patch("matplotlib.pyplot.show"):
        plot_clusters(simple_clustering, k=2, ax=axes[1])
    try:
        assert axes[0].get_title() == "Other plot"
        assert axes[1].get_title(loc="left") == "PhytClust’s clusters at k = 2"
        assert tuple(figure.get_size_inches()) == (10, 4)
    finally:
        plt.close(figure)


@pytest.mark.parametrize("title", ["My clusters", ""])
@pytest.mark.parametrize("plot", [plot_clusters, plot_multiple_k])
def test_custom_title(simple_clustering, title, plot):
    options = {"k": 2} if plot is plot_clusters else {"k_values": [2]}
    with patch("matplotlib.pyplot.show"):
        plot(simple_clustering, title=title, **options)
    figure = plt.gcf()
    try:
        assert figure.axes[0].get_title(loc="left") == title
        assert figure.axes[0].get_title(loc="center") == ""
    finally:
        plt.close(figure)


def test_configured_height_scale_is_honoured(tmp_path):
    from PIL import Image

    pc_default = PhytClust(POLYTOMY_TREE)
    pc_tall = PhytClust(
        POLYTOMY_TREE,
        runtime_config=RuntimeConfig(
            plot=PlotConfig(cluster=ClusterPlotConfig(height_scale=0.8))
        ),
    )
    sizes = []
    for pc in (pc_default, pc_tall):
        pc.get_clusters(4)
        plot_clusters(pc, k=4, results_dir=str(tmp_path), save=True)
        with Image.open(tmp_path / "tree_k4.png") as image:
            sizes.append(image.size[1])
    assert sizes[1] != sizes[0]


def test_continuous_colormap_gives_distinct_colours():
    pc = PhytClust(SAMPLE_TREE)
    figure = plots.plot_cluster(pc.get_clusters(3), pc.tree, cmap="viridis")
    colours = {
        mcolors.to_hex(colour)
        for collection in figure.axes[0].collections
        for colour in collection.get_facecolors()
    }
    assert len(colours) >= 3
    plt.close(figure)
