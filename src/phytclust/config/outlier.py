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

    size_threshold: Optional[int] = 2  # None disables outlier detection

    # If True, minimise outlier count first and break ties on cost.
    prefer_fewer: bool = False

    # Penalty is per formed cluster, not on the total outlier count: a
    # non-linear penalty on the total does not decompose over the tree DP and
    # cannot be optimised exactly. Binary trees only.
    penalty_enabled: bool = False
    ratio_weight: float = 10.0
    # shape(s) for a cluster of size s < size_threshold:
    #   power -> size_threshold - s;  inverse -> 1/s;  exp -> e**(size_threshold - s) - 1
    ratio_mode: str = "exp"
