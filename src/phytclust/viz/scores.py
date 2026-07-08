from typing import Optional, List
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.ticker import LogLocator, LogFormatter, MaxNLocator

from ..exceptions import InvalidClusteringError, ConfigurationError
from .palette import ACCENT_HEX


def plot_scores(
    pc,
    scores_subset: Optional[np.ndarray] = None,
    k_start: int = 1,
    k_end: Optional[int] = None,
    peaks: Optional[List[int]] = None,
    resolution_on: bool = False,
    num_bins: int = 3,
    fig_width: Optional[int] = None,
    fig_height: Optional[int] = None,
    log_scale_y: Optional[bool] = None,
    x_axis_mode: Optional[str] = None,
    log_base: Optional[float] = None,
) -> plt.Figure:
    """Split-out of `_plot_scores` with optional runtime-config defaults."""
    use_full_scores = scores_subset is None
    if scores_subset is None:
        scores_subset = pc.scores
    if len(scores_subset) == 0:
        raise InvalidClusteringError("Scores are empty, please calculate scores first.")

    scores_cfg = getattr(getattr(pc, "plot_config", None), "scores", None)

    title_fontsize = getattr(scores_cfg, "title_fontsize", 50)
    axis_label_fontsize = getattr(scores_cfg, "axis_label_fontsize", 35)
    tick_labelsize = getattr(scores_cfg, "tick_labelsize", 30)
    peak_labelsize = getattr(scores_cfg, "peak_labelsize", 50)

    if fig_width is None:
        fig_width = int(getattr(scores_cfg, "fig_width", 18))
    if fig_height is None:
        fig_height = int(getattr(scores_cfg, "fig_height", 10))
    if log_scale_y is None:
        log_scale_y = bool(getattr(scores_cfg, "log_scale_y", False))
    if x_axis_mode is None:
        x_axis_mode = str(getattr(scores_cfg, "x_axis_mode", "log"))
    if log_base is None:
        log_base = getattr(scores_cfg, "log_base", None)

    clamp_negative_to_zero = bool(getattr(scores_cfg, "clamp_negative_to_zero", True))
    if use_full_scores:
        start_idx = k_start - 1
        end_idx = k_end if k_end is not None else None
        scores_slice = np.asarray(scores_subset[start_idx:end_idx], dtype=float)
    else:
        scores_slice = np.asarray(scores_subset, dtype=float)
        if k_end is not None:
            expected_len = max(0, int(k_end - k_start + 1))
            scores_slice = scores_slice[:expected_len]

    if clamp_negative_to_zero:
        scores_slice = np.maximum(scores_slice, 0.0)

    if scores_slice.size == 0:
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        ax.text(0.5, 0.5, "No scores to plot", ha="center", va="center")
        ax.set_axis_off()
        plt.tight_layout()
        return fig
    x_indices = np.arange(k_start, k_start + len(scores_slice))

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    ax.plot(x_indices, scores_slice, "o-", label="Scores", markersize=2, linewidth=7)

    mid_pts, xtick_labels = [], []
    peak_points = []
    if peaks:
        valid_peaks, valid_scores = [], []
        for p in peaks:
            idx = p - k_start
            if 0 <= idx < len(scores_slice):
                valid_peaks.append(p)
                valid_scores.append(scores_slice[idx])

        if valid_peaks:
            ax.plot(valid_peaks, valid_scores, "rx", label="Peaks", markersize=8)
            peak_points = list(zip(valid_peaks, valid_scores))

        if resolution_on and num_bins > 0:
            bin_ranges = getattr(pc, "bin_ranges_current", None)
            if bin_ranges is None:
                from ..algo.bins import define_bins

                bin_ranges = define_bins(pc, num_bins, k_lo=k_start, k_hi=k_end)

            boundaries = [
                (bin_ranges[i][1] + bin_ranges[i + 1][0]) / 2
                for i in range(len(bin_ranges) - 1)
            ]
            cont_bins = []
            cont_bins.append((bin_ranges[0][0], boundaries[0]))
            for i in range(1, len(bin_ranges) - 1):
                cont_bins.append((boundaries[i - 1], boundaries[i]))
            cont_bins.append((boundaries[-1], bin_ranges[-1][1]))

            custom = getattr(scores_cfg, "colorblind_palette", None)
            cb_palette = list(custom or ACCENT_HEX)
            colours = (cb_palette * (num_bins // len(cb_palette) + 1))[:num_bins]

            mid_pts, xtick_labels = [], []
            use_log_mid = x_axis_mode == "log"

            for i, ((lo_nom, hi_nom), (lo, hi), colour) in enumerate(
                zip(bin_ranges, cont_bins, colours), 1
            ):
                ax.axvspan(lo, hi, color=colour, alpha=0.20, zorder=-1)
                mid = (lo_nom * hi_nom) ** 0.5 if use_log_mid else (lo_nom + hi_nom) / 2
                mid_pts.append(mid)
                xtick_labels.append(f"[{lo_nom}-{hi_nom}]")
                ax.text(
                    mid,
                    1.03,
                    f"CL{i}",
                    ha="center",
                    va="bottom",
                    fontsize=20,
                    weight="bold",
                    color=colour,
                    transform=ax.get_xaxis_transform(),
                    clip_on=False,
                )

            ax.set_xticks(mid_pts)
            ax.xaxis.set_major_locator(mticker.FixedLocator(mid_pts))
            ax.xaxis.set_minor_locator(mticker.NullLocator())

    x_min_plot = max(1, x_indices[0])
    x_max_plot = x_indices[-1] + 0.05 * (x_indices[-1] - x_min_plot)
    ax.set_xlim(x_min_plot, x_max_plot)

    if x_axis_mode == "log":
        ax.set_xscale("log")
        if resolution_on:
            ax.xaxis.set_major_locator(mticker.FixedLocator(mid_pts))
            ax.xaxis.set_minor_locator(mticker.NullLocator())
        else:
            ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=10))
            ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs="auto", numticks=10))
            ax.xaxis.set_major_formatter(LogFormatter(base=10.0, labelOnlyBase=False))
            for label in ax.get_xticklabels():
                label.set_rotation(45)
                label.set_ha("right")
        x_label_str = "No. of Clusters (log)"

    elif x_axis_mode == "linear":
        ax.set_xscale("linear")
        if resolution_on:
            ax.set_xticks(mid_pts)
        else:
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        x_label_str = "No. of Clusters"

    else:
        raise ConfigurationError("x_axis_mode must be either 'log' or 'linear'.")

    if log_scale_y:
        has_nonpositive = bool(np.any(scores_slice <= 0))
        if has_nonpositive and not clamp_negative_to_zero:
            ax.set_yscale("symlog", linthresh=1e-3)
            y_label_str = "Scores (symlog)"
        else:
            ax.set_yscale("log")
            y_label_str = "Scores (log)"
    else:
        y_label_str = "Scores"

    ax.tick_params(axis="both", which="both", labelsize=tick_labelsize)
    ax.set_xlabel(x_label_str, fontsize=axis_label_fontsize)
    ax.set_ylabel(y_label_str, fontsize=axis_label_fontsize)
    title_pad = 62 if resolution_on else 40
    ax.set_title("Scores", fontsize=title_fontsize, pad=title_pad)
    ax.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
    y_min = np.min(scores_slice)
    y_max = np.max(scores_slice)
    padding = (y_max - y_min) * 0.2

    yscale = ax.get_yscale()
    if yscale == "log":
        # Keep lower bound strictly positive for log axis.
        positive_vals = scores_slice[scores_slice > 0]
        if positive_vals.size == 0:
            ax.set_yscale("linear")
            lower = 0.0
            upper = 1.0
            ax.set_ylabel("Scores", fontsize=axis_label_fontsize)
            ax.set_ylim(lower, upper)
        else:
            ymin_pos = max(float(np.min(positive_vals)), 1e-12)
            lower = max(ymin_pos * 0.8, 1e-12)
            upper_base = max(y_max + padding, ymin_pos * 1.2)
            upper = upper_base * 1.35
            ax.set_ylim(lower, upper)
    else:
        lower = 0.0 if clamp_negative_to_zero else (y_min - padding)
        label_headroom = max(padding, 0.15 * max(y_max, 1.0))
        upper = max(y_max + padding + label_headroom, lower + 1.0)
        ax.set_ylim(lower, upper)

    # Add peak labels after limits are finalized to improve placement.
    if peak_points:
        fig.canvas.draw()
        label_fontsize = int(peak_labelsize)
        top_px = ax.bbox.y1

        peak_points_sorted = sorted(peak_points, key=lambda t: t[0])
        last_x_px = None
        near_group_rank = 0
        for px, py in peak_points_sorted:
            x_px, y_px = ax.transData.transform((px, py))

            if last_x_px is None or abs(x_px - last_x_px) > 35:
                near_group_rank = 0
            else:
                near_group_rank += 1
            last_x_px = x_px

            # Keep labels away from CL headers by flipping below when near top.
            if top_px - y_px < 55:
                y_offset = -(10 + 10 * (near_group_rank % 3))
                va = "top"
            else:
                y_offset = 10 + 10 * (near_group_rank % 4)
                va = "bottom"

            ax.annotate(
                str(px),
                xy=(px, py),
                xytext=(0, y_offset),
                textcoords="offset points",
                fontsize=label_fontsize,
                color="red",
                ha="center",
                va=va,
                bbox={
                    "boxstyle": "round,pad=0.15",
                    "fc": "white",
                    "ec": "none",
                    "alpha": 0.7,
                },
            )

    if resolution_on:
        plt.tight_layout(rect=(0.0, 0.0, 1.0, 0.88))
    else:
        plt.tight_layout()
    return fig
