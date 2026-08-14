from dataclasses import dataclass
from typing import Optional

from ..exceptions import ConfigurationError

#: Accepted values for :attr:`PeakConfig.ranking_mode`.
RANKING_MODES = ("raw", "adjusted")

#: Accepted values for :attr:`PeakConfig.resolution_fallback_mode`.
RESOLUTION_FALLBACK_MODES = ("none", "max_score")


@dataclass
class PeakConfig:
    """Configuration for peak detection and ranking in find_score_peaks.

    Usage::

        cfg = PeakConfig(prominence_weight=0.5)
        result = pc.run(top_n=3, peak_config=cfg)
    """

    # --- Ranking ---
    # Blend between the two ranking signals when ranking_mode="adjusted":
    #   metric = prominence_weight * norm(prominence)
    #          + (1 - prominence_weight) * norm(score)
    # 1.0 ranks purely by peak prominence (how much a peak stands out from its
    # neighbours), 0.0 purely by absolute score height. Ignored when
    # ranking_mode="raw", which ranks by prominence alone.
    prominence_weight: float = 0.7
    ranking_mode: str = "adjusted"  # one of RANKING_MODES

    # --- Boundary candidate (k=2) ---
    boundary_window_size: int = 5  # right-window size for k=2 comparison
    boundary_ratio_threshold: float = 1.5  # min ratio vs right-window median

    # --- Detection ---
    min_prominence: Optional[float] = None  # None = auto (1% of score range)
    # Optional: detect peaks on log-transformed scores instead of raw scores.
    # Useful when large-k tails dominate absolute score scale.
    use_log_peak_input: bool = False
    log_peak_offset: float = 1e-12
    # Optional: keep detection on linear scores but rank/filter peaks using
    # relative prominence (fold-change over local baseline) instead of
    # absolute prominence.
    use_relative_prominence: bool = False
    min_relative_prominence: Optional[float] = None
    # Optional k-scaled prominence threshold: keep peak at k only if
    # prominence >= min_prominence * (k ** prominence_k_power).
    # Set to 0.0 to disable scaling (default behavior).
    prominence_k_power: float = 0.0
    min_k: int = 2
    # Resolution mode only: what to do with a bin that contains no detected
    # peak. "none" leaves the bin empty; "max_score" falls back to the
    # highest-scoring k in that bin so every bin returns a k.
    resolution_fallback_mode: str = "none"  # one of RESOLUTION_FALLBACK_MODES
    # If True, k=2 is never returned as a peak. In top_n / global mode the
    # next-ranked peak is chosen instead; in resolution mode the bin that
    # would otherwise pick k=2 advances to its next-ranked candidate.
    exclude_k2: bool = True

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Check the enumerated and bounded fields; raise on a bad value.

        Called on construction and again by ``find_score_peaks``, since config
        files and CLI flags overlay values onto an already-built instance.
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
