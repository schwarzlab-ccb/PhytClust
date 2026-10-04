"""Settings for selecting peaks in the score curve."""

from dataclasses import dataclass
import math
from numbers import Integral, Real

from ..exceptions import ConfigurationError

RANKING_MODES = ("raw", "adjusted")
RESOLUTION_FALLBACK_MODES = ("none", "max_score")


def _finite_nonnegative(value, name):
    """Check a numeric threshold or weight."""
    if (
        isinstance(value, bool)
        or not isinstance(value, Real)
        or not math.isfinite(value)
        or value < 0
    ):
        raise ConfigurationError(f"{name} must be finite and zero or greater.")


def _positive_integer(value, name):
    """Check a window size or cluster count."""
    if isinstance(value, bool) or not isinstance(value, Integral) or value < 1:
        raise ConfigurationError(f"{name} must be a positive integer.")


@dataclass
class PeakConfig:
    """Settings for finding and ranking score peaks."""

    prominence_weight: float = 0.7
    ranking_mode: str = "adjusted"
    partition_preference: str = "none"
    partition_weight: float = 0.5
    boundary_window_size: int = 5
    boundary_ratio_threshold: float = 1.5
    min_prominence: float | None = None
    use_log_peak_input: bool = False
    log_peak_offset: float = 1e-12
    use_relative_prominence: bool = False
    min_relative_prominence: float | None = None
    # Raise the prominence threshold by cluster_count ** prominence_k_power.
    prominence_k_power: float = 0.0
    # Smallest cluster count considered for automatic peak selection.
    min_k: int = 2
    resolution_fallback_mode: str = "none"
    exclude_k2: bool = False

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Check settings before peak detection, including changes made after creation."""
        if self.ranking_mode not in RANKING_MODES:
            raise ConfigurationError(
                f"ranking_mode must be one of {RANKING_MODES}, got {self.ranking_mode!r}."
            )
        if self.resolution_fallback_mode not in RESOLUTION_FALLBACK_MODES:
            raise ConfigurationError(
                f"resolution_fallback_mode must be one of {RESOLUTION_FALLBACK_MODES}, got {self.resolution_fallback_mode!r}."
            )
        if self.partition_preference not in {"none", "fewer_outliers", "balanced"}:
            raise ConfigurationError(
                "partition_preference must be none, fewer_outliers, or balanced."
            )
        if self.partition_preference != "none":
            _finite_nonnegative(self.partition_weight, "partition_weight")
            if self.partition_weight > 1:
                raise ConfigurationError("partition_weight must be between 0 and 1.")
        for name in ("use_log_peak_input", "use_relative_prominence", "exclude_k2"):
            if not isinstance(getattr(self, name), bool):
                raise ConfigurationError(f"{name} must be True or False.")
        _positive_integer(self.min_k, "min_k")
        if self.ranking_mode == "adjusted":
            _finite_nonnegative(self.prominence_weight, "prominence_weight")
            if self.prominence_weight > 1:
                raise ConfigurationError("prominence_weight must be between 0 and 1.")
        if not self.exclude_k2 and self.min_k <= 2:
            _positive_integer(self.boundary_window_size, "boundary_window_size")
            _finite_nonnegative(
                self.boundary_ratio_threshold, "boundary_ratio_threshold"
            )
        if self.min_prominence is not None:
            _finite_nonnegative(self.min_prominence, "min_prominence")
        if self.use_log_peak_input:
            _finite_nonnegative(self.log_peak_offset, "log_peak_offset")
            if self.log_peak_offset == 0:
                raise ConfigurationError("log_peak_offset must be greater than zero.")
        if self.use_relative_prominence:
            if self.min_relative_prominence is not None:
                _finite_nonnegative(
                    self.min_relative_prominence, "min_relative_prominence"
                )
        elif self.min_prominence is not None:
            _finite_nonnegative(self.prominence_k_power, "prominence_k_power")
