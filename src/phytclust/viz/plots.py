from __future__ import annotations

import os
import logging
from dataclasses import dataclass
from numbers import Integral
from typing import Any, Callable, Optional, Union

import matplotlib as mpl
import matplotlib.collections as mpcollections

from ..exceptions import ValidationError
from ..utils.traversal import iter_clades
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.textpath import TextToPath
from matplotlib.font_manager import FontProperties
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

COLORS = {
    "allele_a": mpl.colors.to_rgba("orange"),
    "allele_b": mpl.colors.to_rgba("teal"),
    "clonal": mpl.colors.to_rgba("lightgrey"),
    "normal": mpl.colors.to_rgba("dimgray"),
    "gain": mpl.colors.to_rgba("red"),
    "wgd": mpl.colors.to_rgba("green"),
    "loss": mpl.colors.to_rgba("blue"),
    "chr_label": mpl.colors.to_rgba("grey"),
    "vlines": "#1f77b4",
    "marker_internal": "#1f77b4",
    "marker_terminal": "black",
    "marker_normal": "green",
    "summary_label": "grey",
    "branch_label": mpl.colors.to_rgba("dimgray"),
    "background": "white",
    "background_hatch": "lightgray",
    "patch_background": "white",
}

LINEWIDTHS = {"copy_numbers": 2, "chr_boundary": 1, "segment_boundary": 0.5}
ALPHAS = {"patches": 0.15, "patches_wgd": 0.3, "clonal": 0.3}
SIZES = {
    "tree_marker": 40,
    "ylabel_font": 8,
    "ylabel_tick": 6,
    "xlabel_font": 10,
    "xlabel_tick": 8,
    "chr_label": 8,
    "branch_label": 6,
}


class PlotError(Exception):
    pass


def _format_branch_value(
    value: Optional[float], *, show_zero: bool = False
) -> Optional[str]:
    """Format finite values with three significant digits."""
    if value is None or not np.isfinite(value) or (value == 0 and not show_zero):
        return None
    return f"{value:.3g}"


def _branch_label_formatter(
    branch_labels, show_branch_lengths, show_branch_support=False
):
    """Format lengths and support, or use the supplied labels."""
    if branch_labels is not None:
        if isinstance(branch_labels, dict):
            return branch_labels.get
        if not callable(branch_labels):
            raise ValidationError("branch_labels must be a dict or a callable.")

        def custom_label(clade):
            value = branch_labels(clade)
            return (
                value
                if isinstance(value, str)
                else _format_branch_value(value, show_zero=True)
            )

        return custom_label

    def format_label(clade):
        parts = []
        if show_branch_lengths:
            parts.append(_format_branch_value(getattr(clade, "branch_length", None)))
        if show_branch_support:
            parts.append(
                _format_branch_value(getattr(clade, "confidence", None), show_zero=True)
            )
        return " / ".join(part for part in parts if part is not None) or None

    return format_label


@dataclass
class _TreeGeometry:
    """Coordinates and traversal order shared by tree drawing layers."""

    nodes: list
    leaves: list
    x_positions: dict
    y_positions: dict


def _tree_geometry(
    tree, *, layout="rectangular", hide_internal_nodes=True, outgroup=None
):
    nodes = list(iter_clades(tree.root))
    return _TreeGeometry(
        nodes=nodes,
        leaves=[node for node in nodes if not node.clades],
        x_positions=_get_x_positions(tree, layout=layout, nodes=nodes),
        y_positions=_get_y_positions(
            tree, adjust=not hide_internal_nodes, outgroup=outgroup, nodes=nodes
        ),
    )


def plot_tree(
    input_tree: Any,
    label_func: Optional[Callable[[Any], str]] = None,
    title: str = "",
    ax: Optional[plt.Axes] = None,
    output_name: Optional[str] = None,
    outgroup: Optional[str] = None,
    width_scale: float = 1.0,
    height_scale: Optional[float] = None,
    show_terminal_labels: bool = True,
    show_branch_lengths: bool = True,
    show_branch_support: bool = False,
    branch_labels: Optional[Union[dict[Any, str], Callable[[Any], str]]] = None,
    label_colors: Optional[Union[dict[str, str], Callable[[str], str]]] = None,
    hide_internal_nodes: bool = True,
    marker_size: Optional[int] = None,
    line_width: Optional[float] = None,
    layout: str = "rectangular",
    branch_color_func: Optional[Callable[[Any], Any]] = None,
    show_branch_axis: bool = True,
    _geometry: Optional[_TreeGeometry] = None,
    **kwargs: Any,
) -> plt.Figure:
    """Draw a tree with batched branches and optional leaf and branch labels.

    Rectangular layout uses branch lengths; cladogram layout aligns leaves
    and hides the branch axis. Zero marker size or line width disables that
    element. Callbacks receive clades, and callback errors reach the caller.
    """
    marker_size = SIZES["tree_marker"] if marker_size is None else marker_size
    line_width = LINEWIDTHS["segment_boundary"] if line_width is None else line_width
    for name, value in (("marker_size", marker_size), ("line_width", line_width)):
        if not np.isfinite(value) or value < 0:
            raise ValidationError(f"{name} must be finite and zero or greater.")
    label_func = label_func or (lambda clade: getattr(clade, "name", None))
    geometry = _geometry or _tree_geometry(
        input_tree,
        layout=layout,
        hide_internal_nodes=hide_internal_nodes,
        outgroup=outgroup,
    )

    horizontal_lines: list[list[tuple[float, float]]] = []
    vertical_lines: list[list[tuple[float, float]]] = []
    horizontal_colors: list[Any] = []
    vertical_colors: list[Any] = []
    horizontal_line_widths: list[float] = []
    vertical_line_widths: list[float] = []

    marker_x: list[float] = []
    marker_y: list[float] = []
    marker_sizes: list[float] = []
    marker_colors: list[Any] = []

    text_x: list[float] = []
    text_y: list[float] = []
    texts: list[str] = []
    text_colors: list[Any] = []

    branch_text_x: list[float] = []
    branch_text_y: list[float] = []
    branch_texts: list[str] = []

    if ax is None:
        node_count = len(geometry.nodes)
        plot_height = (
            max(
                2.5,
                (
                    len(geometry.leaves) * 0.24 + 1
                    if show_terminal_labels
                    else min(18, len(geometry.leaves) * 0.035 + 1)
                ),
            )
            if height_scale is None
            else max(1.5, height_scale * node_count * 0.25)
        )
        plot_width = max(8, 5 + width_scale)
        fig, ax = plt.subplots(figsize=(min(250, plot_width), min(250, plot_height)))
    else:
        fig = ax.figure

    if label_colors is None:
        clade_colors: dict[str, Any] = {}
        for clade in geometry.nodes:
            name = getattr(clade, "name", None)
            if not name:
                continue
            is_term = clade.is_terminal()
            clade_colors[name] = (
                COLORS["marker_terminal"] if is_term else COLORS["marker_internal"]
            )
            if outgroup is not None and name == outgroup:
                clade_colors[name] = COLORS["marker_normal"]

        def get_label_color(label):
            return clade_colors.get(label, "black")

    else:
        get_label_color = (
            label_colors
            if callable(label_colors)
            else (lambda label: label_colors.get(label, "black"))
        )

    def marker_func(node):
        name = getattr(node, "name", None)
        if label_colors is None and name is None:
            color = (
                COLORS["marker_internal"] if node.clades else COLORS["marker_terminal"]
            )
        else:
            color = get_label_color(name)
        return marker_size, color

    ax.axes.get_yaxis().set_visible(False)
    for spine in ("right", "left", "top"):
        ax.spines[spine].set_visible(False)
    is_cladogram = layout == "cladogram"
    if is_cladogram or not show_branch_axis:
        ax.spines["bottom"].set_visible(False)
        ax.axes.get_xaxis().set_visible(False)
    else:
        ax.xaxis.set_major_locator(mpl.ticker.AutoLocator())
        ax.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())
        ax.xaxis.set_tick_params(labelsize=SIZES["xlabel_tick"])
        ax.xaxis.label.set_size(SIZES["xlabel_font"])
    ax.set_title(
        title,
        loc="left",
        pad=14,
        fontweight="semibold",
        fontsize=12,
        zorder=10,
    )

    x_positions = geometry.x_positions
    y_positions = geometry.y_positions

    xmax = max(x_positions.values(), default=0.0) or 1.0
    ax.set_xlim(-0.05 * xmax, 1.05 * xmax)
    top_margin = 0.5
    ymax = (max(y_positions.values()) if y_positions else 0.0) + top_margin
    ax.set_ylim(ymax, -0.5)
    ax_scale = ax.get_xlim()[1] - ax.get_xlim()[0]

    format_branch_label = _branch_label_formatter(
        branch_labels, show_branch_lengths, show_branch_support
    )

    def draw_clade_lines(
        *,
        use_linecollection: bool,
        orientation: str,
        y_here: float = 0.0,
        x_start: float = 0.0,
        x_here: float = 0.0,
        y_bot: float = 0.0,
        y_top: float = 0.0,
        color: Any = "black",
        branch_width: float = 0.1,
    ) -> None:
        if use_linecollection and orientation == "horizontal":
            horizontal_lines.append([(x_start, y_here), (x_here, y_here)])
            horizontal_colors.append(color)
            horizontal_line_widths.append(branch_width)
        elif use_linecollection and orientation == "vertical":
            vertical_lines.append([(x_here, y_bot), (x_here, y_top)])
            vertical_colors.append(color)
            vertical_line_widths.append(branch_width)

    def draw_clade(clade: Any, x_start: float, color: Any, branch_width: float) -> None:
        stack: list[tuple[Any, float, Any, float]] = [
            (clade, x_start, color, branch_width)
        ]
        while stack:
            clade, x_start, color, branch_width = stack.pop()
            x_here = x_positions.get(clade, 0.0)
            y_here = y_positions.get(clade, 0.0)

            if hasattr(clade, "color") and clade.color is not None:
                try:
                    color = clade.color.to_hex()
                except Exception:
                    color = clade.color

            if line_width > 0 and hasattr(clade, "width") and clade.width is not None:
                branch_width = float(clade.width) * float(
                    plt.rcParams["lines.linewidth"]
                )

            edge_color = color
            if branch_color_func is not None:
                branch_color = branch_color_func(clade)
                if branch_color is not None:
                    edge_color = branch_color

            draw_clade_lines(
                use_linecollection=True,
                orientation="horizontal",
                y_here=y_here,
                x_start=x_start,
                x_here=x_here,
                color=edge_color,
                branch_width=branch_width,
            )

            branch_label = (
                None if clade is input_tree.root else format_branch_label(clade)
            )
            if branch_label:
                branch_text_x.append((x_start + x_here) / 2.0)
                branch_text_y.append(y_here)
                branch_texts.append(str(branch_label))

            if marker_size > 0 and not (
                hide_internal_nodes and not clade.is_terminal()
            ):
                marker = marker_func(clade)
                if marker is not None:
                    node_marker_size, node_marker_color = marker
                    marker_x.append(x_here)
                    marker_y.append(y_here)
                    marker_sizes.append(node_marker_size)
                    marker_colors.append(node_marker_color)

            label = label_func(clade)
            is_terminal = clade.is_terminal()
            if (
                label not in (None, clade.__class__.__name__)
                and not (hide_internal_nodes and not is_terminal)
                and (show_terminal_labels or not is_terminal)
            ):
                text_x.append(x_here + min(0.02 * ax_scale, 1.0))
                text_y.append(y_here)
                texts.append(f" {label}")
                text_colors.append(get_label_color(getattr(clade, "name", "")))

            if clade.clades:
                child_positions = [y_positions[child] for child in clade.clades]
                y_top = min(child_positions)
                y_bot = max(child_positions)
                draw_clade_lines(
                    use_linecollection=True,
                    orientation="vertical",
                    x_here=x_here,
                    y_bot=y_bot,
                    y_top=y_top,
                    color=edge_color,
                    branch_width=branch_width,
                )
                for child in reversed(clade.clades):
                    stack.append((child, x_here, edge_color, branch_width))

    line_width = float(
        line_width if line_width is not None else plt.rcParams["lines.linewidth"]
    )
    draw_clade(input_tree.root, 0.0, "k", line_width)

    if horizontal_lines:
        h = mpcollections.LineCollection(
            horizontal_lines,
            colors=horizontal_colors,
            linewidths=horizontal_line_widths,
        )
        ax.add_collection(h)
    if vertical_lines:
        v = mpcollections.LineCollection(
            vertical_lines, colors=vertical_colors, linewidths=vertical_line_widths
        )
        ax.add_collection(v)

    if marker_x:
        ax.scatter(marker_x, marker_y, s=marker_sizes, c=marker_colors, zorder=3)

    for x, y, text, color in zip(text_x, text_y, texts, text_colors):
        ax.text(x, y, text, va="center", color=color, fontsize=9)

    for x, y, text in zip(branch_text_x, branch_text_y, branch_texts):
        ax.text(
            x,
            y,
            text,
            ha="center",
            va="bottom",
            fontsize=SIZES["branch_label"],
            color=COLORS["branch_label"],
        )

    if not is_cladogram and show_branch_axis:
        axis_label = (
            "Branch length"
            if any(node.branch_length for node in geometry.nodes)
            else "Topology depth"
        )
        ax.set_xlabel(axis_label, fontsize=10, labelpad=8)
        ax.spines["bottom"].set_color("#999999")
        ax.spines["bottom"].set_linewidth(0.7)
        ax.tick_params(axis="x", labelsize=9, colors="#555555")
    ax.set_ylabel("taxa")

    for key, value in kwargs.items():
        method = getattr(ax, str(key), None) or getattr(ax, f"set_{key}", None)
        if method is None:
            raise ValidationError(f"Unknown plot option: {key}.")
        if isinstance(value, dict):
            method(**value)
        elif (
            isinstance(value, tuple)
            and len(value) == 2
            and isinstance(value[0], tuple)
            and isinstance(value[1], dict)
        ):
            method(*value[0], **value[1])
        else:
            method(*(value if isinstance(value, (list, tuple)) else (value,)))

    if output_name is not None:
        fig.savefig(str(output_name) + ".png", bbox_inches="tight")

    return fig


def _get_x_positions(tree, layout="rectangular", *, nodes=None) -> dict[Any, float]:
    """Calculate branch-length coordinates or aligned cladogram coordinates."""
    if layout not in ("rectangular", "cladogram"):
        raise ValidationError("layout must be 'rectangular' or 'cladogram'.")
    nodes = list(iter_clades(tree.root)) if nodes is None else nodes
    depths = {tree.root: 0}
    lengths = {}
    for node in nodes:
        length = node.branch_length or 0.0
        if not np.isfinite(length) or length < 0:
            raise ValidationError("Branch lengths must be finite and zero or greater.")
        lengths[node] = length
        for child in node.clades:
            depths[child] = depths[node] + 1
    if layout == "cladogram":
        max_depth = max(max(depths.values()), 1)
        return {
            node: float(depths[node]) / max_depth if node.clades else 1.0
            for node in nodes
        }
    positions = {tree.root: float(lengths[tree.root])}
    for node in nodes:
        for child in node.clades:
            positions[child] = positions[node] + lengths[child]
            if not np.isfinite(positions[child]):
                raise ValidationError("Total branch lengths are too large to plot.")
    # With no branch lengths, use topology so the tree remains visible.
    if not max(positions.values()):
        return {node: float(depth) for node, depth in depths.items()}
    return positions


def _get_y_positions(
    tree, adjust=False, outgroup=None, *, nodes=None
) -> dict[Any, float]:
    """Place leaves in tree order and internal nodes between their children.

    An outgroup clade's leaves appear first, at the top. With adjust=True,
    internal nodes also receive space for markers and labels.
    """
    nodes = list(iter_clades(tree.root)) if nodes is None else nodes
    leaves = [node for node in nodes if not node.clades]
    if outgroup is not None:
        matches = [node for node in nodes if node.name == outgroup]
        if len(matches) != 1:
            raise PlotError(f"Outgroup {outgroup!r} must match exactly one clade.")
        outgroup_leaves = {node for node in iter_clades(matches[0]) if not node.clades}
        leaves = [leaf for leaf in leaves if leaf in outgroup_leaves] + [
            leaf for leaf in leaves if leaf not in outgroup_leaves
        ]
    positions = {leaf: float(index) for index, leaf in enumerate(leaves)}
    for node in reversed(nodes):
        if node.clades:
            positions[node] = (
                min(positions[child] for child in node.clades)
                + max(positions[child] for child in node.clades)
            ) / 2
    if adjust:
        ordered_nodes = sorted(positions, key=positions.get)
        positions = {node: float(index) for index, node in enumerate(ordered_nodes)}
    return positions


_MIXED_CLUSTER = object()


def _cluster_common_ancestors(
    tree: Any, clusters_by_id: dict[int, list], *, nodes=None
) -> dict[int, Any]:
    """Find cluster ancestors, including clusters spread across subtrees."""
    nodes = list(iter_clades(tree.root)) if nodes is None else nodes
    tree_leaves = {node for node in nodes if not node.clades}
    if any(
        not members or any(member not in tree_leaves for member in members)
        for members in clusters_by_id.values()
    ):
        raise ValidationError("Cluster members must be leaves of the plotted tree.")
    sizes = {cluster_id: len(members) for cluster_id, members in clusters_by_id.items()}
    leaf_cluster_ids: dict[Any, int] = {}
    for cluster_id, members in clusters_by_id.items():
        for member in members:
            leaf_cluster_ids[member] = cluster_id

    uniform_cluster_ids: dict[Any, Any] = {}
    assigned_leaf_counts: dict[Any, int] = {}
    cluster_ancestors: dict[int, Any] = {}

    for node in reversed(nodes):
        if node.is_terminal():
            cluster_id = leaf_cluster_ids.get(node)
            uniform_cluster_ids[node] = cluster_id
            assigned_leaf_counts[node] = 1 if cluster_id is not None else 0
        else:
            uniform_cluster_id: Any = None
            assigned_leaf_count = 0
            for child_index, child in enumerate(node.clades):
                child_cluster_id = uniform_cluster_ids.get(child)
                assigned_leaf_count += assigned_leaf_counts.get(child, 0)
                if child_index == 0:
                    uniform_cluster_id = child_cluster_id
                elif (
                    uniform_cluster_id is _MIXED_CLUSTER
                    or child_cluster_id is _MIXED_CLUSTER
                    or child_cluster_id != uniform_cluster_id
                ):
                    uniform_cluster_id = _MIXED_CLUSTER
            uniform_cluster_ids[node] = uniform_cluster_id
            assigned_leaf_counts[node] = assigned_leaf_count

        uniform_cluster_id = uniform_cluster_ids[node]
        if (
            uniform_cluster_id is not None
            and uniform_cluster_id is not _MIXED_CLUSTER
            and uniform_cluster_id in sizes
            and assigned_leaf_counts[node] == sizes[uniform_cluster_id]
            and uniform_cluster_id not in cluster_ancestors
        ):
            cluster_ancestors[uniform_cluster_id] = node

    if len(cluster_ancestors) == len(clusters_by_id):
        return cluster_ancestors
    parents = {tree.root: None}
    depths = {tree.root: 0}
    for node in nodes:
        for child in node.clades:
            parents[child] = node
            depths[child] = depths[node] + 1
    for cluster_id, members in clusters_by_id.items():
        if cluster_id in cluster_ancestors:
            continue
        ancestor = members[0]
        for member in members[1:]:
            left, right = ancestor, member
            if left not in parents or right not in parents:
                raise ValidationError(
                    "Cluster members must belong to the plotted tree."
                )
            while depths[left] > depths[right]:
                left = parents[left]
            while depths[right] > depths[left]:
                right = parents[right]
            while left is not right:
                left, right = parents[left], parents[right]
            ancestor = left
        cluster_ancestors[cluster_id] = ancestor
    return cluster_ancestors


def _text_right_edges(axis):
    """Measure the right edges of visible labels in tree coordinates."""
    if not axis.texts:
        return []
    canvas = axis.figure.canvas
    if not hasattr(canvas, "get_renderer"):
        canvas.draw()
    renderer = canvas.get_renderer() if hasattr(canvas, "get_renderer") else None
    inverse_transform = axis.transData.inverted()
    return [
        (
            text.get_position()[1],
            inverse_transform.transform(
                (text.get_window_extent(renderer=renderer).x1, 0)
            )[0],
        )
        for text in axis.texts
        if text.get_visible()
    ]


def _draw_cluster_boxes(
    ax: plt.Axes,
    tree: Any,
    cluster: dict[Any, int],
    cluster_to_color: dict[int, Any],
    *,
    outgroup: Optional[str],
    hide_internal_nodes: bool,
    layout: str = "rectangular",
    mark_outliers: bool = True,
    box_alpha: float = 0.18,
    box_pad_y: float = 0.4,
    box_pad_x_frac: float = 0.01,
    show_labels: bool = True,
    geometry: Optional[_TreeGeometry] = None,
) -> None:
    """Draw cluster boxes from their common ancestors to their leaf tips."""
    geometry = geometry or _tree_geometry(
        tree, layout=layout, hide_internal_nodes=hide_internal_nodes, outgroup=outgroup
    )
    leaves = [
        leaf for leaf in geometry.leaves if leaf.name != outgroup or outgroup is None
    ]
    if not leaves:
        return

    x_positions = geometry.x_positions
    y_positions = geometry.y_positions
    xmax = max(x_positions.values()) if x_positions else 1.0

    clusters_by_id: dict[int, list] = {}
    for leaf in leaves:
        cluster_id = cluster.get(leaf)
        if cluster_id is None:
            continue
        clusters_by_id.setdefault(int(cluster_id), []).append(leaf)

    pad_x = xmax * box_pad_x_frac
    label_edges = _text_right_edges(ax)

    ancestors_by_cluster = _cluster_common_ancestors(
        tree, clusters_by_id, nodes=geometry.nodes
    )

    for cluster_id, members in clusters_by_id.items():
        if cluster_id < 0 and mark_outliers:
            colour = cluster_to_color.get(cluster_id, "#888888")
            outlier = True
        else:
            colour = cluster_to_color.get(cluster_id, "#888888")
            outlier = False
        if colour is None:
            continue

        ys = [y_positions.get(m, 0.0) for m in members]
        if not ys:
            continue
        y_min = min(ys) - box_pad_y
        y_max = max(ys) + box_pad_y

        mrca = ancestors_by_cluster.get(cluster_id, members[0])
        x_left = x_positions.get(mrca, 0.0) - pad_x
        member_y_positions = {y_positions[member] for member in members}
        label_right = max(
            (right for y, right in label_edges if y in member_y_positions), default=0.0
        )
        x_right = (
            max(max(x_positions[member] for member in members), label_right) + pad_x
        )

        rect = plt.Rectangle(
            (x_left, y_min),
            x_right - x_left,
            y_max - y_min,
            facecolor="none" if outlier else colour,
            edgecolor=colour,
            linewidth=1.2,
            linestyle="--" if outlier else "-",
            alpha=box_alpha if not outlier else max(box_alpha * 2, 0.4),
            zorder=1,
            clip_on=False,
        )
        ax.add_patch(rect)

        if show_labels:
            ax.text(
                x_right + xmax * 0.005,
                (y_min + y_max) / 2,
                "outlier" if outlier else f"C{cluster_id}",
                fontsize=7,
                fontweight="600",
                color=colour,
                style="italic" if outlier else "normal",
                va="center",
                ha="left",
                clip_on=False,
            )


def _draw_cluster_bars(
    ax: plt.Axes,
    tree: Any,
    cluster: dict[Any, int],
    cluster_to_color: dict[int, Any],
    *,
    outgroup: Optional[str],
    hide_internal_nodes: bool,
    bar_width_frac: float = 0.04,
    bar_alpha: float = 0.85,
    layout: str = "rectangular",
    geometry: Optional[_TreeGeometry] = None,
) -> None:
    """Draw a side bar, joining adjacent leaves assigned to the same cluster."""
    geometry = geometry or _tree_geometry(
        tree, layout=layout, hide_internal_nodes=hide_internal_nodes, outgroup=outgroup
    )
    leaves = [
        leaf for leaf in geometry.leaves if leaf.name != outgroup or outgroup is None
    ]
    if not leaves:
        return

    y_positions = geometry.y_positions
    leaves.sort(key=lambda c: y_positions.get(c, 0))

    xmax = ax.get_xlim()[1]
    if not np.isfinite(xmax) or xmax <= 0:
        return

    label_pad = xmax * 0.02
    label_right = max((right for _, right in _text_right_edges(ax)), default=xmax)
    label_right = max(label_right, xmax)
    bar_left = label_right + label_pad
    bar_width = max(xmax * bar_width_frac, xmax * 0.02)

    leaf_y = np.array([y_positions[leaf] for leaf in leaves], dtype=float)
    leaf_count = len(leaf_y)
    bounds = _row_boundaries(leaf_y)

    cluster_ids = [cluster.get(c) for c in leaves]
    run_start = 0
    run_cluster_id = cluster_ids[0]
    for i in range(1, leaf_count + 1):
        cluster_id = cluster_ids[i] if i < leaf_count else object()
        if cluster_id != run_cluster_id:
            if run_cluster_id is not None:
                colour = cluster_to_color.get(int(run_cluster_id))
                if colour is not None:
                    ax.add_patch(
                        plt.Rectangle(
                            (bar_left, bounds[run_start]),
                            bar_width,
                            bounds[i] - bounds[run_start],
                            facecolor=colour,
                            edgecolor="none",
                            alpha=bar_alpha,
                            zorder=2,
                            clip_on=False,
                        )
                    )
            run_start = i
            run_cluster_id = cluster_id


def plot_peaks(
    scores_subset: list[float],
    peaks: list[int],
    k_start: int,
    k_end: Optional[int] = None,
    fig_width: int = 10,
    fig_height: int = 10,
    log_scale_x: bool = False,
    log_scale_y: bool = False,
    show_plot: bool = True,
) -> plt.Figure:
    """Plot scores starting at k_start and mark peaks supplied as array indices.

    k_end limits the displayed range. Peak labels show cluster counts, not
    array indices. Return the figure and optionally display it.
    """
    scores = np.asarray(scores_subset, dtype=float)
    if scores.ndim != 1 or not scores.size or not np.isfinite(scores).all():
        raise ValidationError(
            "scores_subset must be a nonempty vector of finite scores."
        )
    if isinstance(k_start, bool) or not isinstance(k_start, Integral) or k_start < 1:
        raise ValidationError("k_start must be a positive integer.")
    if any(
        isinstance(peak, bool)
        or not isinstance(peak, Integral)
        or peak < 0
        or peak >= len(scores)
        for peak in peaks
    ):
        raise ValidationError("Peak indices must refer to entries in scores_subset.")
    if k_end is not None and (
        isinstance(k_end, bool) or not isinstance(k_end, Integral) or k_end < k_start
    ):
        raise ValidationError("k_end must be an integer at least k_start.")
    if log_scale_y and not (scores > 0).all():
        raise ValidationError("Logarithmic score axes require positive scores.")
    title_fontsize, label_fontsize = 16, 14
    tick_labelsize, legend_fontsize, peak_labelsize = 12, 12, 12

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    x_vals = np.arange(len(scores)) + k_start
    ax.plot(x_vals, scores_subset, label="Scores")

    peak_x_vals = np.array(peaks) + k_start
    peak_y_vals = np.array([scores_subset[p] for p in peaks])
    ax.plot(
        peak_x_vals, peak_y_vals, "x", markersize=10, color="red", label="Top Peaks"
    )

    ax.set_title("Cluster score peaks", fontsize=title_fontsize)
    ax.set_xlabel("k", fontsize=label_fontsize)
    ax.set_ylabel("Score", fontsize=label_fontsize)
    ax.tick_params(axis="both", labelsize=tick_labelsize)

    if log_scale_x:
        ax.set_xscale("log")
    if log_scale_y:
        ax.set_yscale("log")

    if k_end is not None:
        ax.set_xlim(k_start, k_end + 1)
        ax.set_title(
            f"Cluster score peaks from k = {k_start} to {k_end}",
            fontsize=title_fontsize,
        )

    data_min, data_max = min(scores_subset), max(scores_subset)
    offset_y = 0.02 if data_max == data_min else 0.02 * (data_max - data_min)

    for px, py in zip(peak_x_vals, peak_y_vals):
        ax.text(
            px,
            py + offset_y,
            str(int(px)),
            fontsize=peak_labelsize,
            ha="center",
            va="bottom",
        )

    ax.legend(fontsize=legend_fontsize)
    if show_plot:
        plt.show()
    return fig


def plot_cluster(
    cluster: dict[Any, int],
    tree: Any,
    *,
    cmap: str | mcolors.Colormap = "phytclust",
    save: bool = False,
    filename: str | None = None,
    results_dir: str | None = None,
    outlier: bool = True,
    hide_internal_nodes: bool = True,
    show_terminal_labels: bool = True,
    width_scale: float = 2.0,
    height_scale: Optional[float] = None,
    label_func: Callable[[Any], str] | None = None,
    show_branch_lengths: bool = False,
    marker_size: int = 40,
    outgroup: str | None = None,
    scores: list[float] | np.ndarray | None = None,
    show_cluster_bars: bool = False,
    cluster_bar_width: float = 0.04,
    cluster_bar_alpha: float = 0.85,
    show_cluster_boxes: bool = False,
    cluster_box_alpha: float = 0.18,
    show_cluster_box_labels: bool = True,
    colour_branches_by_cluster: bool = False,
    layout: str = "rectangular",
    palette: list | None = None,
    show_branch_axis: bool = True,
    title: Optional[str] = None,
    **kwargs: Any,
) -> plt.Figure:
    """Draw one partition with optional cluster bars, boxes, and coloured branches.

    Cluster keys must be leaves of the plotted tree. Negative IDs represent
    outliers; their boxes are dashed when outlier=True. Use title to replace
    the default heading, or an empty string to omit it. A custom palette must
    contain a colour for every cluster. This function returns the figure,
    including when save=True, so the caller can close or edit it.
    """
    geometry = _tree_geometry(
        tree, layout=layout, hide_internal_nodes=hide_internal_nodes, outgroup=outgroup
    )
    tree_leaves = set(geometry.leaves)
    if not cluster or any(leaf not in tree_leaves for leaf in cluster):
        raise ValidationError("Cluster members must be leaves of the plotted tree.")
    if any(
        isinstance(cluster_id, bool) or not isinstance(cluster_id, Integral)
        for cluster_id in cluster.values()
    ):
        raise ValidationError("Cluster IDs must be integers.")
    cluster_count = len(set(cluster.values()))
    if palette is not None:
        resolved_palette = list(palette)
    elif isinstance(cmap, str) and cmap == "phytclust":
        from .palette import expand_palette

        resolved_palette = expand_palette(max(cluster_count, 8))
    elif isinstance(cmap, str):
        cmap_obj = plt.get_cmap(cmap)
        resolved_palette = list(
            cmap_obj.colors
            if hasattr(cmap_obj, "colors")
            else cmap_obj(np.linspace(0, 1, getattr(cmap_obj, "N", 20)))
        )
    else:
        resolved_palette = list(
            cmap.colors
            if hasattr(cmap, "colors")
            else cmap(np.linspace(0, 1, getattr(cmap, "N", 20)))
        )

    if not resolved_palette or (
        palette is not None and cluster_count > len(resolved_palette)
    ):
        raise ValidationError("palette must contain a colour for every cluster.")

    palette = [mcolors.to_rgba(col) for col in resolved_palette]

    ids = np.fromiter(cluster.values(), dtype=int)
    unique_ids = sorted(set(int(v) for v in ids.tolist()))
    color_index_by_cluster = {cluster_id: i for i, cluster_id in enumerate(unique_ids)}
    colours: list[Any] = [
        palette[color_index_by_cluster[int(cluster_id)] % len(palette)]
        for cluster_id in ids
    ]

    if outgroup is not None:
        for i, leaf in enumerate(cluster.keys()):
            if getattr(leaf, "name", None) == outgroup:
                colours[i] = "grey"

    leaf_names = [getattr(leaf, "name", None) for leaf in cluster.keys()]
    if any(not isinstance(name, str) or not name for name in leaf_names) or len(
        set(leaf_names)
    ) != len(leaf_names):
        raise ValidationError("Clustered leaves must have unique, nonempty names.")
    leaf_colors = dict(zip(leaf_names, colours))

    if title is None:
        title = f"PhytClust’s clusters at k = {cluster_count}"
        if scores is not None and cluster_count <= len(scores):
            score = float(scores[cluster_count - 1])
            if np.isfinite(score):
                title += f" (score = {score:.3g})"

    cluster_to_color = {
        int(cluster_id): palette[color_index_by_cluster[int(cluster_id)] % len(palette)]
        for cluster_id in unique_ids
    }

    branch_color_func: Optional[Callable[[Any], Any]] = None
    if colour_branches_by_cluster:
        subtree_cluster_ids: dict[Any, Optional[int]] = {}
        for clade in reversed(geometry.nodes):
            if not clade.clades:
                subtree_cluster_ids[clade] = cluster.get(clade)
                continue
            child_reps = [subtree_cluster_ids.get(c) for c in clade.clades]
            if any(subtree_cluster_id is None for subtree_cluster_id in child_reps):
                subtree_cluster_ids[clade] = None
            elif all(
                subtree_cluster_id == child_reps[0] for subtree_cluster_id in child_reps
            ):
                subtree_cluster_ids[clade] = child_reps[0]
            else:
                subtree_cluster_ids[clade] = None
        default_branch_color = "#888888"

        def _resolve_branch_color(clade: Any) -> Any:
            subtree_cluster_id = subtree_cluster_ids.get(clade)
            if subtree_cluster_id is None or int(subtree_cluster_id) < 0:
                return default_branch_color
            return cluster_to_color.get(int(subtree_cluster_id), default_branch_color)

        branch_color_func = _resolve_branch_color

    effective_label_colors = None if show_cluster_bars else leaf_colors

    fig = plot_tree(
        tree,
        title=title,
        label_colors=effective_label_colors,
        hide_internal_nodes=hide_internal_nodes,
        show_terminal_labels=show_terminal_labels,
        width_scale=width_scale,
        height_scale=height_scale,
        label_func=label_func,
        show_branch_lengths=show_branch_lengths,
        marker_size=marker_size,
        outgroup=outgroup,
        layout=layout,
        branch_color_func=branch_color_func,
        show_branch_axis=show_branch_axis,
        _geometry=geometry,
        **kwargs,
    )

    target_ax = kwargs.get("ax") or fig.axes[0]

    if show_cluster_boxes:
        _draw_cluster_boxes(
            target_ax,
            tree,
            cluster,
            cluster_to_color,
            outgroup=outgroup,
            hide_internal_nodes=hide_internal_nodes,
            layout=layout,
            box_alpha=cluster_box_alpha,
            show_labels=show_cluster_box_labels,
            mark_outliers=outlier,
            geometry=geometry,
        )

    if show_cluster_bars:
        _draw_cluster_bars(
            target_ax,
            tree,
            cluster,
            cluster_to_color,
            outgroup=outgroup,
            hide_internal_nodes=hide_internal_nodes,
            bar_width_frac=cluster_bar_width,
            bar_alpha=cluster_bar_alpha,
            layout=layout,
            geometry=geometry,
        )

    if save:
        results_dir = results_dir or "."
        os.makedirs(results_dir, exist_ok=True)
        filename = filename or f"tree_k{cluster_count}.png"
        fig.savefig(os.path.join(results_dir, filename), bbox_inches="tight")

    return fig


def _comparison_table(input_df: pd.DataFrame) -> pd.DataFrame:
    """Accept a leaf-by-solution table or the older long table format."""
    if isinstance(input_df.index, pd.MultiIndex):
        required_levels = {"leaf_name", "comparison_IDs"}
        if set(input_df.index.names) != required_levels or "cluster_ID" not in input_df:
            raise ValidationError(
                "Long comparison tables need leaf_name and comparison_IDs index levels and a cluster_ID column."
            )
        if input_df.index.has_duplicates:
            raise ValidationError("Each leaf must have one assignment per solution.")
        solution_order = input_df.index.get_level_values("comparison_IDs").unique()
        table = (
            input_df["cluster_ID"]
            .unstack("comparison_IDs")
            .reindex(columns=solution_order)
        )
    else:
        table = input_df.copy()
    if table.empty or table.index.has_duplicates or table.columns.has_duplicates:
        raise ValidationError(
            "Comparison tables need unique leaves and solutions and at least one assignment."
        )
    if any(not isinstance(name, str) or not name for name in table.index):
        raise ValidationError("Leaf names must be nonempty strings.")
    try:
        values = table.to_numpy(dtype=float)
    except (TypeError, ValueError) as error:
        raise ValidationError("Cluster assignments must be finite integers.") from error
    if not np.isfinite(values).all() or not np.equal(values, np.floor(values)).all():
        raise ValidationError(
            "Every leaf needs a finite integer assignment for every solution."
        )
    return table


def _row_boundaries(positions):
    """Place row edges halfway between neighbouring leaf coordinates."""
    positions = np.asarray(positions, dtype=float)
    if len(positions) == 1:
        return np.array([positions[0] - 0.5, positions[0] + 0.5])
    boundaries = np.empty(len(positions) + 1)
    boundaries[1:-1] = (positions[:-1] + positions[1:]) / 2
    boundaries[0] = positions[0] - (boundaries[1] - positions[0])
    boundaries[-1] = positions[-1] + (positions[-1] - boundaries[-2])
    return boundaries


def _fit_tree_labels(axis):
    """Keep tip labels inside the tree panel without moving them from their tips."""
    if not axis.texts:
        return
    canvas = axis.figure.canvas
    if not hasattr(canvas, "get_renderer"):
        canvas.draw()
    renderer = canvas.get_renderer() if hasattr(canvas, "get_renderer") else None
    axis_width = axis.get_window_extent(renderer).width
    left, right = axis.get_xlim()
    for text in axis.texts:
        label_width = text.get_window_extent(renderer).width
        available_fraction = max(0.1, 1 - (label_width + 8) / axis_width)
        right = max(right, left + (text.get_position()[0] - left) / available_fraction)
    axis.set_xlim(left, right)


def plot_multiple_clusters(
    input_df: pd.DataFrame,
    final_tree: Optional[Any] = None,
    y_posns: Optional[dict[str, float]] = None,
    cmax: Optional[int] = None,
    tree_width_ratio: float = 1.0,
    cbar_width_ratio: float = 0.05,
    figsize: Optional[tuple[float, float]] = None,
    tree_marker_size: int = 0,
    show_internal_nodes: bool = False,
    title: str = "PhytClust’s cluster assignments",
    tree_label_func: Optional[Callable[[Any], str]] = None,
    cmap: str = "phytclust",
    outgroup: Optional[str] = None,
    fixed_x_range: Optional[tuple[float, float]] = None,
) -> plt.Figure:
    """Compare assignments across solutions, with an optional aligned tree.

    Rows are leaf names and columns are solution names. The older long format
    with leaf_name and comparison_IDs index levels and a cluster_ID column is
    also accepted. Each leaf needs one integer assignment in every solution.
    Leaf labels sit beside tree tips when a tree is supplied. Assignment
    bands meet halfway between neighbouring leaves. Vertical separators
    distinguish solutions; there are no horizontal cell borders.
    IDs are local to each partition; matching colours across columns do not
    imply that the clusters have identical members. Negative IDs are grey.
    y_posns supplies leaf coordinates. fixed_x_range changes the heatmap's
    horizontal span. cbar_width_ratio controls the colour key's width.
    cmax, when supplied, is an upper limit for assignment IDs. Without figsize,
    the figure grows to keep leaf labels readable. Supply figsize for a fixed size.
    """
    custom_y_positions = y_posns is not None
    table = _comparison_table(input_df)
    if not np.isfinite(tree_width_ratio) or tree_width_ratio <= 0:
        raise ValidationError("tree_width_ratio must be finite and greater than zero.")
    if not np.isfinite(cbar_width_ratio) or cbar_width_ratio < 0:
        raise ValidationError("cbar_width_ratio must be finite and zero or greater.")
    values = table.to_numpy(dtype=float)
    if cmax is not None and (not np.isfinite(cmax) or cmax < values.max()):
        raise ValidationError(
            "cmax must be finite and at least the largest cluster ID."
        )
    geometry = None
    if final_tree is not None:
        geometry = _tree_geometry(
            final_tree, hide_internal_nodes=not show_internal_nodes, outgroup=outgroup
        )
        leaf_names = [leaf.name for leaf in geometry.leaves]
        if len(set(leaf_names)) != len(leaf_names) or set(leaf_names) != set(
            table.index
        ):
            raise ValidationError(
                "The comparison table must contain exactly the tree's unique leaf names."
            )
        if y_posns is None:
            y_posns = {
                leaf.name: geometry.y_positions[leaf] for leaf in geometry.leaves
            }
    if y_posns is None:
        y_posns = {name: float(index) for index, name in enumerate(table.index)}
    if any(
        name not in y_posns or not np.isfinite(y_posns[name]) for name in table.index
    ):
        raise ValidationError("y_posns must give a finite coordinate for every leaf.")
    leaf_order = sorted(table.index, key=y_posns.__getitem__)
    row_centres = np.array([y_posns[name] for name in leaf_order], dtype=float)
    if len(np.unique(row_centres)) != len(row_centres):
        raise ValidationError("Leaf coordinates must be distinct.")
    table = table.loc[leaf_order]
    if geometry is not None and custom_y_positions:
        for leaf in geometry.leaves:
            geometry.y_positions[leaf] = y_posns[leaf.name]
        for node in reversed(geometry.nodes):
            if node.clades:
                geometry.y_positions[node] = (
                    min(geometry.y_positions[child] for child in node.clades)
                    + max(geometry.y_positions[child] for child in node.clades)
                ) / 2
    column_count = len(table.columns)
    start, end = (0.0, float(column_count)) if fixed_x_range is None else fixed_x_range
    if not np.isfinite([start, end]).all() or end <= start:
        raise ValidationError("fixed_x_range must have two finite, increasing values.")
    column_edges = np.linspace(start, end, column_count + 1)
    row_edges = _row_boundaries(row_centres)
    cluster_ids = sorted(set(table.to_numpy(dtype=float).ravel()))
    color_indices = {cluster_id: index for index, cluster_id in enumerate(cluster_ids)}
    data = np.array(
        [[color_indices[value] for value in row] for row in table.to_numpy(dtype=float)]
    )
    if cmap == "phytclust":
        from .palette import expand_palette

        colors = expand_palette(max(len(cluster_ids), 8))[: len(cluster_ids)]
    else:
        colors = list(plt.get_cmap(cmap)(np.linspace(0, 1, len(cluster_ids))))
    colors = [
        "#aaaaaa" if cluster_id < 0 else color
        for cluster_id, color in zip(cluster_ids, colors)
    ]
    color_map = mcolors.ListedColormap(colors)
    color_norm = mcolors.BoundaryNorm(
        np.arange(len(cluster_ids) + 1) - 0.5, len(cluster_ids)
    )
    text_metrics = TextToPath()
    label_font = FontProperties(size=9)
    label_width_inches = (
        max(
            text_metrics.get_text_width_height_descent(name, label_font, False)[0]
            for name in leaf_order
        )
        / 72
    )
    if figsize is None:
        minimum_gap = np.diff(row_centres).min() if len(row_centres) > 1 else 1.0
        row_span = (row_edges[-1] - row_edges[0]) / minimum_gap
        figsize = (
            max(
                10,
                column_count * 0.65
                + label_width_inches
                + (4 if geometry is not None else 1),
            ),
            min(250, max(6, row_span * 0.24 + 1)),
        )
    if geometry is None:
        fig, heat_axis = plt.subplots(figsize=figsize)
    else:
        fig, (tree_axis, heat_axis) = plt.subplots(
            ncols=2,
            figsize=figsize,
            gridspec_kw={"width_ratios": [tree_width_ratio, 1], "wspace": 0.04},
        )
        plot_tree(
            final_tree,
            ax=tree_axis,
            outgroup=outgroup,
            label_func=tree_label_func
            or (lambda node: node.name if not node.clades else None),
            hide_internal_nodes=not show_internal_nodes,
            show_terminal_labels=True,
            show_branch_lengths=False,
            line_width=0.8,
            marker_size=tree_marker_size,
            _geometry=geometry,
        )
        tree_axis.set_axis_off()
        tree_axis.set_ylim(row_edges[-1], row_edges[0])
    left, top = np.meshgrid(column_edges[:-1], row_edges[:-1])
    right, bottom = np.meshgrid(column_edges[1:], row_edges[1:])
    vertices = np.stack(
        (
            np.stack((left, top), axis=-1),
            np.stack((right, top), axis=-1),
            np.stack((right, bottom), axis=-1),
            np.stack((left, bottom), axis=-1),
        ),
        axis=2,
    ).reshape(-1, 4, 2)
    mesh = mpcollections.PolyCollection(
        vertices,
        array=data.ravel(),
        cmap=color_map,
        norm=color_norm,
        edgecolors="none",
        linewidths=0,
        antialiased=False,
        rasterized=True,
    )
    heat_axis.add_collection(mesh)
    for boundary in column_edges[1:-1]:
        heat_axis.axvline(boundary, color="#333333", linewidth=0.7, zorder=3)
    heat_axis.set_xlim(start, end)
    heat_axis.set_xticks((column_edges[:-1] + column_edges[1:]) / 2)
    heat_axis.set_xticklabels([str(name) for name in table.columns], fontsize=9)
    heat_axis.tick_params(
        axis="x", labeltop=True, labelbottom=False, top=False, bottom=False, pad=8
    )
    heat_axis.set_yticks(row_centres if geometry is None else [])
    if geometry is None:
        heat_axis.set_yticklabels(leaf_order, fontsize=9)
    heat_axis.tick_params(axis="y", length=0)
    heat_axis.set_ylim(row_edges[-1], row_edges[0])
    for spine in heat_axis.spines.values():
        spine.set_visible(False)
    if cbar_width_ratio:
        colorbar = fig.colorbar(
            mesh,
            ax=heat_axis,
            fraction=cbar_width_ratio,
            pad=0.04,
            ticks=np.arange(len(cluster_ids)),
        )
        colorbar.ax.set_yticklabels(
            [str(int(cluster_id)) for cluster_id in cluster_ids]
        )
        colorbar.ax.tick_params(labelsize=8, length=0)
        colorbar.set_label("Cluster ID", fontsize=9)
        colorbar.outline.set_visible(False)
    figure_height = figsize[1]
    fig.subplots_adjust(
        top=1 - min(0.25, 0.75 / figure_height), bottom=min(0.2, 0.45 / figure_height)
    )
    if geometry is not None:
        _fit_tree_labels(tree_axis)
    fig.suptitle(
        title,
        fontsize=12,
        fontweight="semibold",
        x=0.125,
        y=1 - min(0.1, 0.2 / figure_height),
        ha="left",
    )
    fig.set_facecolor("white")
    return fig
