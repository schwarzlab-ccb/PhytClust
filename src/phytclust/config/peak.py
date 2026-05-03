from dataclasses import dataclass
from typing import Optional


@dataclass
class PeakConfig:
    """Configuration for peak detection and ranking in find_score_peaks.

    Usage::

        cfg = PeakConfig(lambda_weight=0.5)
        result = pc.run(top_n=3, peak_config=cfg)
    """

    # --- Ranking ---
    lambda_weight: float = 0.7
    ranking_mode: str = "adjusted"  # "raw" or "adjusted"

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
    # Resolution mode: fallback for bins with no detected peaks.
    resolution_fallback_mode: str = "none"  # "none" or "max_score"
    # If True, k=2 is never returned as a peak. In top_n / global mode the
    # next-ranked peak is chosen instead; in resolution mode the bin that
    # would otherwise pick k=2 advances to its next-ranked candidate.
    exclude_k2: bool = True
