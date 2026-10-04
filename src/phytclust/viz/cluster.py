import os
import logging
import warnings
from typing import Any, Optional, Callable, List
from numbers import Integral
from pathlib import Path
import matplotlib.pyplot as plt

from ..exceptions import ConfigurationError
from .plots import plot_cluster

logger = logging.getLogger(__name__)


def _positive_count(value, name: str) -> int:
    """Check that a plot count is a positive integer."""
    if isinstance(value, bool) or not isinstance(value, Integral) or value < 1:
        raise ConfigurationError(f"{name} must be a positive integer.")
    return int(value)


def _partition_maps(pc, counts) -> list[tuple[int, dict]]:
    """Prepare the tables once and retrieve each requested partition."""
    counts = list(dict.fromkeys(_positive_count(count, "k") for count in counts))
    if not counts:
        return []
    pc._ensure_dp(required_cap=max(counts))
    return [(count, pc._clusters(count)) for count in counts]


def _resolve_plot_targets(pc, *, top_n, k, n) -> list[tuple[int, dict]]:
    """Use explicit k, then the old n argument, stored k, or ranked peaks."""
    top_n = _positive_count(top_n, "top_n")
    if n is not None and k is None:
        warnings.warn(
            "Parameter 'n' is deprecated; use 'k' instead.",
            DeprecationWarning,
            stacklevel=3,
        )
    selected_count = k if k is not None else n if n is not None else pc.k
    counts = (
        [selected_count]
        if selected_count is not None
        else (pc.peaks_by_rank or [])[:top_n]
    )
    return _partition_maps(pc, counts)


def _style_cluster_figure(
    figure, cluster_count, leaf_count, *, labels, resize, axis=None, title=None
):
    """Space leaf labels and apply a consistent title and branch axis."""
    if resize:
        width = figure.get_size_inches()[0]
        height = leaf_count * 0.24 + 1 if labels else min(18, leaf_count * 0.035 + 1)
        figure.set_size_inches(max(width, 8), min(250, max(2.5, height)))
    axis = figure.axes[0] if axis is None else axis
    axis.set_title(
        f"PhytClust’s clusters at k = {cluster_count}" if title is None else title,
        loc="left",
        fontsize=12,
        fontweight="semibold",
        pad=14,
    )
    # The renderer uses a left-aligned centre title; remove it before adding ours.
    axis.set_title("", loc="center")
    axis.set_ylabel("")
    if axis.xaxis.get_visible():
        axis.set_xlabel(axis.get_xlabel() or "Branch length", fontsize=10, labelpad=8)
        axis.tick_params(axis="x", labelsize=9, colors="#555555")
        axis.spines["bottom"].set_color("#999999")
        axis.spines["bottom"].set_linewidth(0.7)
    figure.set_facecolor("white")


def _save_or_show(figure, *, results_dir, save, filename, count, multiple, dpi=180):
    """Save and close a figure, or display it."""
    if not (save or results_dir):
        plt.show()
        return
    output_name = Path(filename or f"tree_k{count}.png")
    if filename and multiple:
        output_name = output_name.with_name(
            f"{output_name.stem}_k{count}{output_name.suffix}"
        )
    output_path = Path(results_dir or ".") / output_name
    try:
        figure.savefig(output_path, bbox_inches="tight", dpi=dpi, facecolor="white")
    finally:
        plt.close(figure)


def plot_clusters(
    pc,
    results_dir: Optional[str] = None,
    top_n: int = 1,
    n: Optional[int] = None,
    k: Optional[int] = None,
    cmap=None,
    show_terminal_labels: bool = False,
    outlier: bool = False,
    save: bool = False,
    filename: Optional[str] = None,
    hide_internal_nodes: Optional[bool] = None,
    width_scale: Optional[float] = None,
    height_scale: Optional[float] = None,
    label_func: Optional[Callable[[Any], str]] = None,
    show_branch_lengths: Optional[bool] = None,
    marker_size: Optional[int] = None,
    show_cluster_bars: bool = False,
    show_cluster_boxes: bool = False,
    colour_branches_by_cluster: bool = False,
    layout: str = "rectangular",
    palette: Optional[List] = None,
    show_branch_axis: bool = True,
    title: Optional[str] = None,
    dpi: int = 180,
    **kwargs,
) -> None:
    """Draw selected partitions and optionally save each figure.

    Explicit k takes precedence over stored k and ranked peaks. When several
    partitions share a filename, append their k values. Labelled trees receive
    extra vertical space unless an axis or height_scale is supplied. Use title
    to replace the default heading; an empty string removes it.
    """
    dpi = _positive_count(dpi, "dpi")
    automatic_height = height_scale is None and kwargs.get("ax") is None
    cluster_cfg = getattr(getattr(pc, "plot_config", None), "cluster", None)

    def _cfg(value, field, fallback):
        if value is not None:
            return value
        return getattr(cluster_cfg, field, fallback)

    cmap = _cfg(cmap, "cmap", "phytclust")
    hide_internal_nodes = _cfg(hide_internal_nodes, "hide_internal_nodes", True)
    width_scale = _cfg(width_scale, "width_scale", 2.0)
    height_scale = _cfg(height_scale, "height_scale", 0.1)
    show_branch_lengths = _cfg(show_branch_lengths, "show_branch_lengths", False)
    marker_size = _cfg(marker_size, "marker_size", 40)

    clusters_to_plot = _resolve_plot_targets(pc, top_n=top_n, k=k, n=n)

    if not clusters_to_plot:
        logger.warning("No partitions selected. Provide k or run peak selection first.")
        return

    if (save or results_dir) and results_dir is not None:
        os.makedirs(results_dir, exist_ok=True)

    plot_tree = (
        pc._tree_wo_outgroup
        if pc.outgroup and pc._tree_wo_outgroup is not None
        else pc.tree
    )
    kwargs.setdefault("line_width", 0.9)
    for k_val, clmap in clusters_to_plot:
        fig = plot_cluster(
            cluster=clmap,
            tree=plot_tree,
            cmap=cmap,
            outlier=outlier,
            hide_internal_nodes=hide_internal_nodes,
            show_terminal_labels=show_terminal_labels,
            width_scale=width_scale,
            height_scale=height_scale,
            label_func=label_func,
            show_branch_lengths=show_branch_lengths,
            marker_size=marker_size,
            outgroup=None,
            results_dir=None,
            show_cluster_bars=show_cluster_bars,
            show_cluster_boxes=show_cluster_boxes,
            colour_branches_by_cluster=colour_branches_by_cluster,
            layout=layout,
            palette=palette,
            show_branch_axis=show_branch_axis,
            **kwargs,
        )

        _style_cluster_figure(
            fig,
            k_val,
            len(clmap),
            labels=show_terminal_labels,
            resize=automatic_height,
            axis=kwargs.get("ax"),
            title=title,
        )
        _save_or_show(
            fig,
            results_dir=results_dir,
            save=save,
            filename=filename,
            count=k_val,
            multiple=len(clusters_to_plot) > 1,
            dpi=dpi,
        )


def plot_multiple_k(
    pc,
    k_values: Optional[List[int]] = None,
    results_dir: Optional[str] = None,
    top_n: int = 3,
    cmap="phytclust",
    show_terminal_labels: bool = False,
    hide_internal_nodes: bool = True,
    width_scale: float = 1.5,
    height_scale: Optional[float] = None,
    label_func: Optional[Callable[[Any], str]] = None,
    show_branch_lengths: bool = False,
    marker_size: int = 30,
    save: bool = False,
    title: Optional[str] = None,
    dpi: int = 180,
    **kwargs,
) -> None:
    """Draw the supplied cluster counts, or the highest-ranked peaks.

    Prepare the clustering tables once for the largest requested count.
    Plotting and clustering errors are reported to the caller. Use title to
    give every figure the same heading; an empty string removes it.
    """
    dpi = _positive_count(dpi, "dpi")
    top_n = _positive_count(top_n, "top_n")
    counts = k_values if k_values is not None else (pc.peaks_by_rank or [])[:top_n]
    targets = _partition_maps(pc, counts)
    if not targets:
        logger.warning(
            "No partitions selected. Provide k_values or run peak selection first."
        )
        return
    cluster_cfg = dict(
        cmap=cmap,
        show_terminal_labels=show_terminal_labels,
        hide_internal_nodes=hide_internal_nodes,
        width_scale=width_scale,
        height_scale=0.08 if height_scale is None else height_scale,
        label_func=label_func,
        show_branch_lengths=show_branch_lengths,
        marker_size=marker_size,
    )
    filename = kwargs.pop("filename", None)
    if results_dir is not None:
        os.makedirs(results_dir, exist_ok=True)
    plot_tree = (
        pc._tree_wo_outgroup
        if pc.outgroup and pc._tree_wo_outgroup is not None
        else pc.tree
    )
    kwargs.setdefault("line_width", 0.9)
    for count, partition in targets:
        figure = plot_cluster(
            cluster=partition, tree=plot_tree, outgroup=None, **cluster_cfg, **kwargs
        )
        _style_cluster_figure(
            figure,
            count,
            len(partition),
            labels=show_terminal_labels,
            resize=height_scale is None and kwargs.get("ax") is None,
            axis=kwargs.get("ax"),
            title=title,
        )
        _save_or_show(
            figure,
            results_dir=results_dir,
            save=save,
            filename=filename,
            count=count,
            multiple=len(targets) > 1,
            dpi=dpi,
        )
