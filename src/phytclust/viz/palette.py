"""PhytClust brand palette — single source of truth for all colors.

The base palette starts from logo-derived tones and extends with harmonized
accents. When more colors are needed, lighter/less opaque variants are
generated automatically while preserving a readable alpha floor.
"""

from __future__ import annotations

import numpy as np
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap

# Base colors for cluster fills. Same brand tones, led by red as before, but
# ordered so that consecutive entries are always perceptually far apart (every
# adjacent pair is >=~50 Lab units): with few clusters, neighbouring clusters
# never get near-identical colours. This fixes the original order's collisions
# (e.g. steel blue directly followed by indigo, olive next to teal).
BASE_HEX: list[str] = [
    "#b84b4b",  # red
    "#4f8f4a",  # green
    "#da63aa",  # pink
    "#ceb94b",  # gold
    "#3f408a",  # indigo
    "#c06f2e",  # amber brown
    "#5b6bb3",  # slate blue
    "#849060",  # olive
    "#6e3f8a",  # purple
    "#2f8a85",  # cyan teal
    "#8f4b7f",  # magenta plum
    "#3d7c74",  # teal
    "#ad5c7a",  # rose
    "#2f6f93",  # ocean blue
    "#7a5d3b",  # earth
    "#3f648a",  # steel blue
]

# Colorblind-safe accent palette (for score-plot peaks / bin labels)
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


def expand_palette(
    n: int,
    base: list[str] | None = None,
) -> list[tuple[float, ...]]:
    """Return *n* RGBA colors by cycling base colors with decreasing alpha.

    For n <= len(base), returns the first n base colors.
    For larger n, lightens and reduces alpha in successive rounds,
    matching the GUI's ``generateClusterColors`` behavior.
    """
    base = base or BASE_HEX
    base_rgba = [mcolors.to_rgba(h) for h in base]
    if n <= len(base_rgba):
        return base_rgba[:n]

    colors: list[tuple[float, ...]] = []
    min_alpha = 0.58
    alpha_step = 0.12
    light_step = 0.14
    repeats = int(np.ceil(n / len(base_rgba)))
    for r in range(repeats):
        factor = r * light_step
        alpha = max(1.0 - r * alpha_step, min_alpha)
        for rgba in base_rgba:
            adjusted = _adjust_lightness(rgba, factor)
            colors.append((*adjusted[:3], alpha))
    return colors[:n]


def get_cmap(n: int = 20) -> ListedColormap:
    """Return a PhytClust ``ListedColormap`` with *n* colors."""
    return ListedColormap(expand_palette(n), name="phytclust")


def _adjust_lightness(rgba: tuple[float, ...], factor: float) -> tuple[float, ...]:
    """Shift lightness of an RGBA color by *factor* (-1..1)."""
    r, g, b = rgba[:3]
    if factor > 0:
        r = r + (1 - r) * factor
        g = g + (1 - g) * factor
        b = b + (1 - b) * factor
    elif factor < 0:
        r = r * (1 + factor)
        g = g * (1 + factor)
        b = b * (1 + factor)
    return (
        max(0.0, min(1.0, r)),
        max(0.0, min(1.0, g)),
        max(0.0, min(1.0, b)),
        rgba[3],
    )
