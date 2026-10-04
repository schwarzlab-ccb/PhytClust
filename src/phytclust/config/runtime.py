"""Default plotting and save settings."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field, fields, is_dataclass
import warnings
from typing import Any

from ..exceptions import ConfigurationError


@dataclass
class ScorePlotConfig:
    """Default sizes, colours, and axes for score plots."""

    title_fontsize: int = 20
    title_fontfamily: str = "Liberation Sans"
    axis_label_fontsize: int = 18
    tick_labelsize: int = 14
    peak_labelsize: int = 16
    bin_labelsize: int = 14
    peak_marker: str = "o"
    peak_markersize: int = 7
    fig_width: float = 9
    fig_height: float = 5
    clamp_negative_to_zero: bool = True
    log_scale_y: bool = True
    x_axis_mode: str = "log"
    log_base: float | None = None
    prefer_unsmoothed_primary: bool = True
    show_secondary_score_plot: bool = False
    colorblind_palette: list[str] | None = None


@dataclass
class ClusterPlotConfig:
    """Default settings for cluster plots."""

    cmap: str = "phytclust"
    width_scale: float = 2.0
    height_scale: float = 0.1
    marker_size: int = 40
    show_branch_lengths: bool = False
    hide_internal_nodes: bool = True


@dataclass
class SaveConfig:
    """Default filename and outlier labels for saved results."""

    tsv_name: str = "phytclust_results.tsv"
    outlier: bool = True


@dataclass
class PlotConfig:
    """Default settings for cluster and score plots."""

    cluster: ClusterPlotConfig = field(default_factory=ClusterPlotConfig)
    scores: ScorePlotConfig = field(default_factory=ScorePlotConfig)


@dataclass
class RuntimeConfig:
    """Plotting and save settings used by PhytClust."""

    plot: PlotConfig = field(default_factory=PlotConfig)
    save: SaveConfig = field(default_factory=SaveConfig)


def _apply_runtime_overrides(
    config: Any, overrides: Mapping[str, Any], setting_path: str = ""
) -> None:
    """Apply nested overrides and warn about unknown settings."""
    field_names = {config_field.name for config_field in fields(config)}
    for name, value in overrides.items():
        full_name = f"{setting_path}.{name}" if setting_path else str(name)
        if name not in field_names:
            warnings.warn(
                f"Unknown runtime setting {full_name!r}; ignored.",
                UserWarning,
                stacklevel=3,
            )
            continue
        current_value = getattr(config, name)
        if is_dataclass(current_value):
            if not isinstance(value, Mapping):
                raise ConfigurationError(f"{full_name} must be a mapping of settings.")
            _apply_runtime_overrides(current_value, value, full_name)
        elif isinstance(value, Mapping):
            raise ConfigurationError(f"{full_name} must be a value, not a mapping.")
        else:
            setattr(config, name, value)


def build_runtime_config(
    overrides: Mapping[str, Any] | None = None,
) -> RuntimeConfig:
    """Create plotting and save settings with optional nested overrides."""
    config = RuntimeConfig()
    if overrides is not None:
        if not isinstance(overrides, Mapping):
            raise ConfigurationError("Runtime overrides must be a mapping of settings.")
        _apply_runtime_overrides(config, overrides)
    return config
