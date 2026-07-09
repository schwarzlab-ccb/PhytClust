from dataclasses import dataclass
from typing import Optional


@dataclass
class OutlierConfig:
    """Configuration for outlier detection and handling in DP and output.

    Controls which clusters are treated as outliers, how they affect
    the DP cost function, and whether they are marked as -1 in output.

    Usage::

        from phytclust.config import OutlierConfig
        cfg = OutlierConfig(size_threshold=5, prefer_fewer=True)
        pc = PhytClust(tree=tree, outlier=cfg)
    """

    # Clusters with fewer leaves than this are outliers.
    # None disables threshold-based outlier detection.
    size_threshold: Optional[int] = 2

    # If True, DP minimises outlier count first, then breaks ties by cost.
    # If False (default), cost is minimised first, outlier count is tie-breaker.
    prefer_fewer: bool = False

    # Additive penalty applied, during the DP, to every cluster smaller than
    # ``size_threshold`` — steering the optimiser away from tiny/outlier
    # clusters. It is a *penalised* objective: ``raw_dp_table`` keeps the pure
    # within-cluster dispersion (used for scoring), while ``dp_table`` carries
    # raw + penalty (what the DP minimises). Opt-in and OFF by default so
    # existing results are unchanged; enable with ``penalty_enabled=True``.
    #
    # The penalty is applied per formed cluster of size ``s < size_threshold``
    # as ``ratio_weight * shape(s)``, which is additive over clusters and so
    # decomposes exactly over the tree DP. ``ratio_mode`` selects ``shape``:
    #   "power"   -> (size_threshold - s)      # linear in the size deficit
    #   "inverse" -> 1.0 / s                   # size-based, mild
    #   "exp"     -> exp(size_threshold - s) - 1.0   # steep
    # NOTE: a non-linear penalty on the *total outlier count* (rather than
    # per-cluster) would NOT decompose over the tree DP and cannot be optimised
    # exactly; that is why the penalty is defined per-cluster. Penalty is
    # currently supported for fully-bifurcating trees only.
    penalty_enabled: bool = False
    ratio_weight: float = 10.0
    ratio_mode: str = "exp"  # "exp", "inverse", or "power"
