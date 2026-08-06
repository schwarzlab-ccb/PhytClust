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

    # If True, the outlier-ratio penalty shapes the DP cost function.
    penalty_enabled: bool = False

    # Weight and mode for the outlier-ratio penalty in the DP cost function.
    ratio_weight: float = 10.0
    ratio_mode: str = "exp"  # "exp", "inverse", or "power"
