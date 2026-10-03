from numbers import Integral
from typing import Optional

import numpy as np

from ..exceptions import ConfigurationError, InvalidKError


def define_bins(
    clustering,
    num_bins: int = 3,
    *,
    k_lo: int = 1,
    k_hi: Optional[int] = None,
) -> list[tuple[int, int]]:
    """Return inclusive, non-overlapping logarithmic ranges of cluster counts.

    Cover k_lo through k_hi. Default k_hi to the number of leaves.
    Rounding may produce fewer than num_bins ranges.

    Resolution mode starts at 1 to preserve detail at small cluster counts,
    although 1 is not selected as a score peak.
    """
    if k_hi is None:
        k_hi = clustering.num_terminals
    for name, value, error in (
        ("num_bins", num_bins, ConfigurationError),
        ("k_lo", k_lo, InvalidKError),
        ("k_hi", k_hi, InvalidKError),
    ):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value < 1:
            raise error(f"{name} must be a positive integer.")
    if k_hi <= k_lo:
        raise InvalidKError("k_hi must be > k_lo")

    logarithmic_edges = np.geomspace(k_lo, k_hi, num_bins + 1)
    integer_edges = np.unique(np.round(logarithmic_edges).astype(int))
    integer_edges[0], integer_edges[-1] = k_lo, k_hi

    bins = []
    for index in range(len(integer_edges) - 1):
        lower_bound = int(integer_edges[index]) if index == 0 else bins[-1][1] + 1
        upper_bound = int(integer_edges[index + 1])
        if lower_bound <= upper_bound:
            bins.append((lower_bound, upper_bound))
    return bins
