"""Cluster colours and accent colours used in PhytClust plots."""

from __future__ import annotations

import math
import warnings
from numbers import Integral
from collections.abc import Sequence

from ..exceptions import ConfigurationError
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap

BASE_HEX: list[str] = [
    "#b84b4b",
    "#4f8f4a",
    "#da63aa",
    "#ceb94b",
    "#3f408a",
    "#c06f2e",
    "#5b6bb3",
    "#849060",
    "#6e3f8a",
    "#2f8a85",
    "#8f4b7f",
    "#3d7c74",
    "#ad5c7a",
    "#2f6f93",
    "#7a5d3b",
    "#3f648a",
]

ACCENT_HEX: list[str] = [
    "#0072B2",
    "#009E73",
    "#D55E00",
    "#56B4E9",
    "#E69F00",
    "#CC79A7",
    "#F0E442",
    "#000000",
]

BASE_RGBA: list[tuple[float, ...]] = [mcolors.to_rgba(h) for h in BASE_HEX]


def _validate_color_count(value, *, allow_zero=True) -> int:
    """Check the requested number of colours."""
    minimum = 0 if allow_zero else 1
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        requirement = "zero or greater" if allow_zero else "greater than zero"
        raise ConfigurationError(f"Colour count must be an integer {requirement}.")
    return int(value)


def expand_palette(
    n: int, base: Sequence[str | tuple[float, ...]] | None = None
) -> list[tuple[float, ...]]:
    """Return n colours, using lighter shades after the base colours.

    Each round multiplies the original opacity by the round's opacity.
    Lightening stops at 84% toward white; colours repeat after seven rounds.
    Base colours retain their original order.
    """
    color_count = _validate_color_count(n)
    base_colors = BASE_HEX if base is None else list(base)
    if not base_colors:
        raise ConfigurationError("The base palette must contain at least one colour.")
    try:
        base_rgba = [mcolors.to_rgba(color) for color in base_colors]
    except (ValueError, TypeError) as error:
        raise ConfigurationError(f"Invalid palette colour: {error}") from error
    if any(not math.isfinite(channel) for color in base_rgba for channel in color):
        raise ConfigurationError("Palette colour channels must be finite.")
    if color_count > len(base_rgba) * 7:
        warnings.warn(
            f"Colours repeat beyond {len(base_rgba) * 7} entries. Use cluster labels for larger partitions.",
            UserWarning,
            stacklevel=2,
        )
    colors = []
    for color_index in range(color_count):
        round_index, base_index = divmod(color_index, len(base_rgba))
        lightening_fraction = min(round_index * 0.14, 0.84)
        round_opacity = max(1 - round_index * 0.12, 0.58)
        color = _blend_toward_white(base_rgba[base_index], lightening_fraction)
        colors.append((*color[:3], color[3] * round_opacity))
    return colors


def get_cmap(n: int = 20) -> ListedColormap:
    """Return a PhytClust colormap with n colours; n must be positive."""
    color_count = _validate_color_count(n, allow_zero=False)
    return ListedColormap(expand_palette(color_count), name="phytclust")


def _blend_toward_white(rgba: tuple[float, ...], fraction: float) -> tuple[float, ...]:
    """Blend RGB channels toward white, retaining the original opacity."""
    return (*(channel + (1 - channel) * fraction for channel in rgba[:3]), rgba[3])
