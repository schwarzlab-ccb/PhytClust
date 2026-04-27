import os
import logging
import warnings
import numpy as np
from typing import Any, Optional, Callable, Tuple, List
import matplotlib.pyplot as plt

from ..algo.dp import cluster_map
from .draw import plot_cluster

logger = logging.getLogger(__name__)


def _get_map(pc, k: int):
    """Retrieve or compute a cluster map for a given k."""
    try:
        get_clusters = getattr(pc, "get_clusters", None)
        if callable(get_clusters):
            return get_clusters(int(k))
    except Exception:
        pass
    return cluster_map(pc, int(k))


def _resolve_plot_targets(
    pc,
    *,
    top_n: int,
    k: Optional[int],
    n: Optional[int],
) -> list[tuple[int, dict]]:
    """Resolve which k-values to plot.

    Precedence is preserved for backward compatibility:
    explicit ``k`` > compatibility alias ``n`` > ``pc.k`` > ranked peaks.
    """
    clusters_to_plot: list[tuple[int, dict]] = []
    if n is not None and k is None:
        warnings.warn(
            "Parameter 'n' is deprecated; use 'k' instead.",
            DeprecationWarning,
            stacklevel=3,
        )

    selected_k = k if k is not None else n if n is not None else pc.k

    if selected_k is not None:
        k_val = int(selected_k)
        clmap = _get_map(pc, k_val)
        if clmap is not None:
            clusters_to_plot.append((k_val, clmap))
        return clusters_to_plot

    for k_val in (pc.peaks_by_rank or [])[:top_n]:
        kv = int(k_val)
        clmap = _get_map(pc, kv)
        if clmap is not None:
            clusters_to_plot.append((kv, clmap))
    return clusters_to_plot


def plot_clusters(
    pc,
    results_dir: Optional[str] = None,
    top_n: int = 1,
    n: Optional[int] = None,
    k: Optional[int] = None,
    cmap="phytclust",
    show_terminal_labels: bool = False,
    outlier: bool = False,
    save: bool = False,
    filename: Optional[str] = None,
    hide_internal_nodes: bool = True,
    width_scale: float = 2,
    height_scale: float = 0.1,
    label_func: Optional[Callable[[Any], Tuple[float, str]]] = None,
    show_branch_lengths: bool = False,
    marker_size: int = 40,
    show_cluster_bars: bool = False,
    show_cluster_boxes: bool = False,
    colour_branches_by_cluster: bool = False,
    layout: str = "rectangular",
    palette: Optional[List] = None,
    show_branch_axis: bool = True,
    **kwargs,
) -> None:
    if pc.clusters is None:
        pc.clusters = {}

    if pc.k is None and (pc.scores is None):
        logger.info("Scores not available - continuing without a score plot.")

    clusters_to_plot = _resolve_plot_targets(pc, top_n=top_n, k=k, n=n)

    if not clusters_to_plot:
        logger.warning("No clusters to plot - check your arguments.")
        return

    if (save or results_dir) and results_dir is not None:
        os.makedirs(results_dir, exist_ok=True)

    for k_val, clmap in clusters_to_plot:
        fig = plot_cluster(
            cluster=clmap,
            tree=pc.tree,
            cmap=cmap,
            outlier=outlier,
            hide_internal_nodes=hide_internal_nodes,
            show_terminal_labels=show_terminal_labels,
            width_scale=width_scale,
            height_scale=height_scale,
            label_func=label_func,
            show_branch_lengths=show_branch_lengths,
            marker_size=marker_size,
            outgroup=pc.outgroup,
            results_dir=None,
            show_cluster_bars=show_cluster_bars,
            show_cluster_boxes=show_cluster_boxes,
            colour_branches_by_cluster=colour_branches_by_cluster,
            layout=layout,
            palette=palette,
            show_branch_axis=show_branch_axis,
            **kwargs,
        )

        if save or results_dir:
            out_dir = results_dir or "."
            out_name = filename or f"tree_k{k_val}.png"
            out_path = os.path.join(out_dir, out_name)
            fig.savefig(out_path, bbox_inches="tight")
            plt.close(fig)
        else:
            plt.show()


def plot_multiple_k(
    pc,
    k_values: Optional[List[int]] = None,
    results_dir: Optional[str] = None,
    top_n: int = 3,
    cmap="phytclust",
    show_terminal_labels: bool = False,
    hide_internal_nodes: bool = True,
    width_scale: float = 1.5,
    height_scale: float = 0.08,
    label_func: Optional[Callable[[Any], Tuple[float, str]]] = None,
    show_branch_lengths: bool = False,
    marker_size: int = 30,
    save: bool = False,
    **kwargs,
) -> None:
    """
    Plot multiple cluster solutions (different k values).

    Parameters
    ----------
    pc : PhytClust
        The PhytClust object with computed clusters.
    k_values : list[int], optional
        Specific k values to plot. If None, uses top_n peaks.
    results_dir : str, optional
        Directory for saving output.
    top_n : int, optional
        Number of peaks to use if k_values not specified.
    cmap : colormap, optional
        Matplotlib colormap for clusters.
    show_terminal_labels : bool, optional
        Show leaf names.
    hide_internal_nodes : bool, optional
        Hide internal node markers.
    width_scale : float, optional
        Scale factor for subplot width.
    height_scale : float, optional
        Scale factor for subplot height.
    label_func : callable, optional
        Function to format labels.
    show_branch_lengths : bool, optional
        Show branch length labels.
    marker_size : int, optional
        Size of leaf markers.
    save : bool, optional
        Save figures to files named tree_k{k}.png.
    **kwargs : dict
        Additional arguments passed to plot_cluster.
    """
    if pc.clusters is None:
        pc.clusters = {}

    # Determine which k values to plot
    clusters_to_plot: List[Tuple[int, dict]] = []

    source = k_values if k_values is not None else (pc.peaks_by_rank or [])[:top_n]
    for k_val in source:
        try:
            clmap = _get_map(pc, int(k_val))
            if clmap is not None:
                clusters_to_plot.append((int(k_val), clmap))
        except Exception:
            logger.warning("Could not compute clusters for k=%d", k_val)

    if not clusters_to_plot:
        logger.warning("No clusters to plot - check your arguments.")
        return

    if (save or results_dir) and results_dir is not None:
        os.makedirs(results_dir, exist_ok=True)

    for k_val, clmap in clusters_to_plot:
        try:
            fig_single = plot_cluster(
                cluster=clmap,
                tree=pc.tree,
                cmap=cmap,
                outlier=False,
                hide_internal_nodes=hide_internal_nodes,
                show_terminal_labels=show_terminal_labels,
                width_scale=width_scale,
                height_scale=height_scale,
                label_func=label_func,
                show_branch_lengths=show_branch_lengths,
                marker_size=marker_size,
                outgroup=pc.outgroup,
                **kwargs,
            )
            if save or results_dir:
                out_dir = results_dir or "."
                out_name = f"tree_k{k_val}.png"
                out_path = os.path.join(out_dir, out_name)
                fig_single.savefig(out_path, bbox_inches="tight", dpi=150)
                plt.close(fig_single)
            else:
                plt.show()
        except Exception as e:
            logger.warning("Error plotting k=%d: %s", k_val, e)
