"""Plot partition scores and the selected cluster counts."""

from numbers import Integral
import warnings

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

from ..exceptions import ConfigurationError, InvalidClusteringError
from .palette import ACCENT_HEX


def _cluster_count(value, name):
    """Check an integer cluster count supplied to the plot."""
    if isinstance(value, bool) or not isinstance(value, Integral) or value < 1:
        raise ConfigurationError(f"{name} must be a positive integer.")
    return int(value)


def _resolution_bands(pc, axis, config, first, last, num_bins, logarithmic):
    """Draw the portions of resolution classes inside the plotted range."""
    from ..algo.bins import define_bins

    ranges = getattr(pc, "bin_ranges_current", None)
    if ranges is None:
        ranges = (
            define_bins(pc, num_bins, k_lo=first, k_hi=last)
            if last > first
            else [(first, last)]
        )
    palette = getattr(config, "colorblind_palette", None)
    palette = list(ACCENT_HEX if palette is None else palette)
    if not palette or any(not mcolors.is_color_like(color) for color in palette):
        raise ConfigurationError("colorblind_palette must contain valid colours.")
    centres, labels = [], []
    for index, (start, end) in enumerate(ranges):
        start, end = max(first, start), min(last, end)
        if start > end:
            continue
        colour = palette[index % len(palette)]
        axis.axvspan(
            max(first, start - 0.5),
            min(last, end + 0.5),
            color=colour,
            alpha=0.12,
            zorder=0,
        )
        centre = (
            np.exp((np.log(start) + np.log(end)) / 2)
            if logarithmic
            else (start + end) / 2
        )
        centres.append(centre)
        labels.append(f"[{start}–{end}]")
        axis.text(
            centre,
            1.025,
            f"CL{index + 1}",
            transform=axis.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=getattr(config, "bin_labelsize", 9),
            weight="semibold",
            color=colour,
        )
    if centres:
        axis.xaxis.set_major_locator(mticker.FixedLocator(centres))
        axis.xaxis.set_major_formatter(mticker.FixedFormatter(labels))
        axis.xaxis.set_minor_locator(mticker.NullLocator())


def _label_peaks(axis, points, fontsize, marker_size):
    """Place cluster-count labels without covering earlier labels."""
    renderer = axis.figure.canvas.get_renderer()
    occupied = []
    for count, score in sorted(points):
        annotation = axis.annotate(
            str(count),
            (count, score),
            xytext=(0, 0),
            textcoords="offset points",
            fontsize=fontsize,
            color="#b84b4b",
            ha="center",
            bbox={
                "boxstyle": "round,pad=0.15",
                "fc": "white",
                "ec": "none",
                "alpha": 0.9,
            },
        )
        for row in range(len(points) + 2):
            for direction in (1, -1):
                annotation.set_position(
                    (0, direction * (marker_size / 2 + 5 + row * (fontsize + 3)))
                )
                annotation.set_verticalalignment("bottom" if direction > 0 else "top")
                bounds = annotation.get_window_extent(renderer).expanded(1.15, 1.15)
                if (
                    axis.bbox.contains(bounds.x0, bounds.y0)
                    and axis.bbox.contains(bounds.x1, bounds.y1)
                    and not any(bounds.overlaps(previous) for previous in occupied)
                ):
                    break
            else:
                continue
            break
        occupied.append(bounds)


def plot_scores(
    pc,
    scores_subset: np.ndarray | None = None,
    k_start: int = 1,
    k_end: int | None = None,
    peaks: list[int] | None = None,
    resolution_on: bool = False,
    num_bins: int = 3,
    fig_width: float | None = None,
    fig_height: float | None = None,
    log_scale_y: bool | None = None,
    x_axis_mode: str | None = None,
    log_base: float | None = None,
    title: str | None = None,
) -> plt.Figure:
    """Draw scores against cluster count, with optional peak labels and resolution bands.

    An explicit score subset starts at ``k_start``; stored scores start at one.
    Peaks are cluster counts. Zero scores are omitted on a logarithmic axis.
    When negative scores are retained, a symmetric log axis also shows zero.
    ``log_base=None`` uses base ten. Set ``title=""`` to leave out the title.
    """
    k_start = _cluster_count(k_start, "k_start")
    if k_end is not None:
        k_end = _cluster_count(k_end, "k_end")
        if k_end < k_start:
            raise ConfigurationError("k_end must be at least k_start.")
    config = getattr(getattr(pc, "plot_config", None), "scores", None)
    width = getattr(config, "fig_width", 9) if fig_width is None else fig_width
    height = getattr(config, "fig_height", 5) if fig_height is None else fig_height
    if not all(np.isfinite(value) and value > 0 for value in (width, height)):
        raise ConfigurationError("Figure width and height must be finite and positive.")
    mode = getattr(config, "x_axis_mode", "log") if x_axis_mode is None else x_axis_mode
    if mode not in ("log", "linear"):
        raise ConfigurationError("x_axis_mode must be either 'log' or 'linear'.")
    base = getattr(config, "log_base", None) if log_base is None else log_base
    base = 10.0 if base is None else base
    if not np.isfinite(base) or base <= 1:
        raise ConfigurationError("log_base must be finite and greater than one.")
    if resolution_on:
        num_bins = _cluster_count(num_bins, "num_bins")
    use_stored = scores_subset is None
    values = getattr(pc, "scores", None) if use_stored else scores_subset
    if values is None or (use_stored and np.size(values) == 0):
        raise InvalidClusteringError("Scores are empty, please calculate scores first.")
    try:
        values = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as error:
        raise InvalidClusteringError(
            "Scores must be a finite one-dimensional array."
        ) from error
    if values.ndim != 1 or not np.all(np.isfinite(values)):
        raise InvalidClusteringError("Scores must be a finite one-dimensional array.")
    values = (
        values[k_start - 1 : k_end]
        if use_stored
        else values[: None if k_end is None else k_end - k_start + 1]
    )
    counts = np.arange(k_start, k_start + len(values))
    zero_negative_scores = getattr(config, "clamp_negative_to_zero", True)
    if zero_negative_scores:
        values = np.maximum(values, 0)
    peak_counts = (
        []
        if peaks is None
        else sorted({_cluster_count(value, "Peak cluster count") for value in peaks})
    )
    log_y = getattr(config, "log_scale_y", True) if log_scale_y is None else log_scale_y
    plotted_values = values.copy()
    scale = "linear"
    if log_y and values.size:
        if not zero_negative_scores and np.any(values <= 0):
            scale = "symlog"
        elif np.any(values > 0):
            scale = "log"
            if np.any(values <= 0):
                warnings.warn(
                    "Zero or negative scores are omitted from the logarithmic plot.",
                    UserWarning,
                    stacklevel=2,
                )
                plotted_values[values <= 0] = np.nan
    figure, axis = plt.subplots(figsize=(width, height))
    if not values.size:
        axis.text(0.5, 0.5, "No scores to plot", ha="center", va="center")
        axis.set_axis_off()
        figure.tight_layout()
        return figure
    axis.plot(
        counts,
        plotted_values,
        "o-",
        color="#3f408a",
        markersize=2,
        linewidth=1.6,
        zorder=2,
    )
    if mode == "log":
        axis.set_xscale("log", base=base)
        axis.xaxis.set_major_locator(mticker.LogLocator(base=base, numticks=8))
        axis.xaxis.set_major_formatter(
            mticker.LogFormatter(base=base, labelOnlyBase=False)
        )
    else:
        axis.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    axis.set_xlim(max(0.5, counts[0] - 0.5), counts[-1] + 0.5)
    if scale == "log":
        axis.set_yscale("log", base=base)
        positive = values[values > 0]
        axis.set_ylim(
            max(np.nextafter(0.0, 1.0), positive.min() * 0.8), positive.max() * 1.4
        )
    else:
        if scale == "symlog":
            nonzero = np.abs(values[values != 0])
            axis.set_yscale(
                "symlog",
                base=base,
                linthresh=float(nonzero.min()) if nonzero.size else 1.0,
            )
        span = float(np.ptp(values)) or float(np.max(np.abs(values))) or 1.0
        axis.set_ylim(
            0 if zero_negative_scores else values.min() - span * 0.15,
            values.max() + span * 0.25,
        )
    marker_size = getattr(config, "peak_markersize", 7)
    points = [
        (count, plotted_values[count - k_start])
        for count in peak_counts
        if k_start <= count < k_start + len(values)
        and np.isfinite(plotted_values[count - k_start])
    ]
    if points:
        axis.plot(
            *zip(*points),
            linestyle="none",
            marker=getattr(config, "peak_marker", "o"),
            markersize=marker_size,
            markerfacecolor="#b84b4b",
            markeredgecolor="white",
            markeredgewidth=1,
            zorder=3,
        )
    if resolution_on:
        try:
            _resolution_bands(
                pc, axis, config, k_start, int(counts[-1]), num_bins, mode == "log"
            )
        except Exception:
            plt.close(figure)
            raise
    for name in ("top", "right"):
        axis.spines[name].set_visible(False)
    for name in ("bottom", "left"):
        axis.spines[name].set_color("#b4b4b4")
    axis.grid(axis="y", which="major", color="#dedede", linewidth=0.6)
    axis.set_axisbelow(True)
    axis.tick_params(
        axis="both",
        which="both",
        labelsize=getattr(config, "tick_labelsize", 9),
        color="#999999",
    )
    axis.set_xlabel(
        "Cluster count" + (" (log scale)" if mode == "log" else ""),
        fontsize=getattr(config, "axis_label_fontsize", 10),
    )
    axis.set_ylabel(
        "Score" + (f" ({scale})" if scale != "linear" else ""),
        fontsize=getattr(config, "axis_label_fontsize", 10),
    )
    axis.set_title(
        "PhytClust’s scores" if title is None else title,
        fontsize=getattr(config, "title_fontsize", 12),
        fontfamily=[
            getattr(config, "title_fontfamily", "Liberation Sans"),
            "DejaVu Sans",
        ],
        weight="normal",
        loc="left",
        pad=30 if resolution_on else 14,
    )
    figure.tight_layout()
    _label_peaks(axis, points, getattr(config, "peak_labelsize", 9), marker_size)
    return figure
