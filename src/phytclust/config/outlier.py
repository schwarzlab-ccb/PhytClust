"""Settings for identifying and penalizing small clusters."""

from dataclasses import dataclass
import math
from numbers import Integral, Real

from ..exceptions import ConfigurationError


@dataclass
class OutlierConfig:
    """Settings for small-cluster counts and penalties.

    Clusters smaller than size_threshold count as outliers. Setting
    prefer_fewer prioritizes their count before partition cost.
    With no threshold, prefer_fewer counts singleton clusters as outliers.
    Penalties require both penalty_enabled and an explicit size_threshold.
    Saved outlier labels are controlled separately by the save options.
    """

    size_threshold: int | None = None
    prefer_fewer: bool = False
    penalty_enabled: bool = False
    ratio_weight: float = 10.0
    ratio_mode: str = "exp"

    @property
    def counting_threshold(self) -> int | None:
        """Use the chosen threshold, or singletons when fewer outliers are preferred."""
        if self.size_threshold is not None:
            return self.size_threshold
        return 2 if self.prefer_fewer else None

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Check active outlier settings, including changes made after creation."""
        for name in ("prefer_fewer", "penalty_enabled"):
            if not isinstance(getattr(self, name), bool):
                raise ConfigurationError(f"outlier.{name} must be True or False.")
        if self.size_threshold is not None:
            if (
                isinstance(self.size_threshold, bool)
                or not isinstance(self.size_threshold, Integral)
                or self.size_threshold < 1
            ):
                raise ConfigurationError(
                    "outlier.size_threshold must be a positive integer."
                )
        if not self.penalty_enabled or self.size_threshold is None:
            return
        if (
            isinstance(self.ratio_weight, bool)
            or not isinstance(self.ratio_weight, Real)
            or not math.isfinite(self.ratio_weight)
            or self.ratio_weight < 0
        ):
            raise ConfigurationError(
                "outlier.ratio_weight must be finite and zero or greater."
            )
        if self.ratio_weight > 0 and self.ratio_mode not in ("exp", "inverse", "power"):
            raise ConfigurationError(
                "outlier.ratio_mode must be 'exp', 'inverse', or 'power'."
            )
