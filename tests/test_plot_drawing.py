from unittest.mock import patch

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from Bio.Phylo.BaseTree import Clade, Tree

from phytclust import PhytClust
from phytclust.exceptions import ValidationError
from phytclust.viz import plots


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
        plots.plot_tree(clustering.tree, branch_color_func=fail)


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


@pytest.mark.parametrize("palette", [[], ["red"]])
def test_short_custom_palette_is_rejected(clustering, palette):
    with pytest.raises(ValidationError, match="palette"):
        plots.plot_cluster(clustering.get_clusters(2), clustering.tree, palette=palette)


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


def test_peak_labels_and_coordinates_use_cluster_counts():
    figure = plots.plot_peaks([0.2, 0.8, 0.3], [1], k_start=2, show_plot=False)
    np.testing.assert_allclose(figure.axes[0].lines[0].get_xdata(), [2, 3, 4])
    assert figure.axes[0].texts[0].get_text() == "3"


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
