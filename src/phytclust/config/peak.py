from dataclasses import dataclass
from typing import Optional

from ..exceptions import ConfigurationError

RANKING_MODES = ("raw", "adjusted")
RESOLUTION_FALLBACK_MODES = ("none", "max_score")


@dataclass
class PeakConfig:
    """Configuration for peak detection and ranking in find_score_peaks.

    Usage::

        cfg = PeakConfig(prominence_weight=0.5)
        result = pc.run(top_n=3, peak_config=cfg)
    """

    # --- Ranking ---
    # 1.0 ranks on prominence alone, 0.0 on absolute score alone. Ignored
    # when ranking_mode="raw".
    prominence_weight: float = 0.7
    ranking_mode: str = "adjusted"

    # --- Boundary candidate (k=2) ---
    boundary_window_size: int = 5
    boundary_ratio_threshold: float = 1.5

    # --- Detection ---
    min_prominence: Optional[float] = None  # None = 1% of score range
    use_log_peak_input: bool = False
    log_peak_offset: float = 1e-12
    use_relative_prominence: bool = False
    min_relative_prominence: Optional[float] = None
    # Keep peak at k only if prominence >= min_prominence * k**power.
    prominence_k_power: float = 0.0
    min_k: int = 2
    resolution_fallback_mode: str = "none"  # resolution mode only
    exclude_k2: bool = True

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Check enumerated and bounded fields.

        Called again by find_score_peaks, since config files and CLI flags
        overlay values onto an already-built instance.
        """
        if self.ranking_mode not in RANKING_MODES:
            raise ConfigurationError(
                f"ranking_mode must be one of {RANKING_MODES}, "
                f"got {self.ranking_mode!r}."
            )
        if self.resolution_fallback_mode not in RESOLUTION_FALLBACK_MODES:
            raise ConfigurationError(
                f"resolution_fallback_mode must be one of "
                f"{RESOLUTION_FALLBACK_MODES}, "
                f"got {self.resolution_fallback_mode!r}."
            )
        if not 0.0 <= float(self.prominence_weight) <= 1.0:
            raise ConfigurationError(
                f"prominence_weight must be between 0 and 1, "
                f"got {self.prominence_weight!r}."
            )
