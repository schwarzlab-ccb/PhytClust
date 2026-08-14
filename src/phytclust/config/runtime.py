"""Runtime settings schema for plotting and saving.

Algorithm tuning lives in `phytclust.config.peak.PeakConfig`.
This module intentionally contains only plotting/save settings.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class ScorePlotConfig:
    """Score-plot specific defaults."""

    title_fontsize: int = 50
    axis_label_fontsize: int = 35
    tick_labelsize: int = 30
    peak_labelsize: int = 50
    fig_width: int = 18
    fig_height: int = 10
    # Clamp negative scores to 0 on the plot by default.
    clamp_negative_to_zero: bool = True
    log_scale_y: bool = True
    x_axis_mode: str = "log"
    log_base: float | None = None
    # If True, save/display k>=3 score curve as primary `scores.png` output.
    prefer_unsmoothed_primary: bool = True
    # If True, generate the secondary companion score figure.
    show_secondary_score_plot: bool = False
    colorblind_palette: list[str] | None = None


@dataclass
class ClusterPlotConfig:
    """Cluster/tree rendering defaults.

    These are the single source of truth for cluster-plot defaults: both the CLI
    and ``plot_clusters`` resolve unset options from here, so a
    ``RuntimeConfig`` set on the ``PhytClust`` object applies to ``pc.plot()``
    as well. Explicitly passed keyword arguments always win.
    """

    cmap: str = "phytclust"
    width_scale: float = 2.0
    height_scale: float = 0.1
    marker_size: int = 40
    show_branch_lengths: bool = False
    hide_internal_nodes: bool = True


@dataclass
class SaveConfig:
    """Output defaults used by save methods / CLI.

    The results table is always tab-separated, so the field is ``tsv_name``
    and the matching CLI flag is ``--tsv-name``.
    """

    tsv_name: str = "phytclust_results.tsv"
    outlier: bool = True


@dataclass
class PlotConfig:
    """Top-level plotting config grouped by score vs tree/cluster plots."""

    cluster: ClusterPlotConfig = field(default_factory=ClusterPlotConfig)
    scores: ScorePlotConfig = field(default_factory=ScorePlotConfig)


@dataclass
class RuntimeConfig:
    """Runtime defaults attached to `PhytClust` for plotting/save behavior."""

    plot: PlotConfig = field(default_factory=PlotConfig)
    save: SaveConfig = field(default_factory=SaveConfig)


def _deep_update_dataclass(dc_obj: Any, updates: dict[str, Any]) -> None:
    """Recursively apply dictionary overrides into nested dataclasses."""
    for key, value in updates.items():
        if not hasattr(dc_obj, key):
            continue
        cur = getattr(dc_obj, key)
        if hasattr(cur, "__dataclass_fields__") and isinstance(value, dict):
            _deep_update_dataclass(cur, value)
        else:
            setattr(dc_obj, key, value)


def build_runtime_config(
    overrides: dict[str, Any] | None = None,
) -> RuntimeConfig:
    """Construct a typed runtime config, optionally applying overrides."""
    cfg = RuntimeConfig()
    if overrides:
        _deep_update_dataclass(cfg, overrides)
    return cfg
