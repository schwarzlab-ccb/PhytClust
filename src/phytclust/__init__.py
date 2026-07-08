from .algo.core import PhytClust
from .config import (
    OutlierConfig,
    PeakConfig,
    RuntimeConfig,
    PlotConfig,
    ScorePlotConfig,
    ClusterPlotConfig,
)
from .exceptions import (
    PhytClustError,
    ValidationError,
    ConfigurationError,
    ComputationError,
    DataError,
    InvalidKError,
    InvalidTreeError,
    MissingDPTableError,
    InvalidClusteringError,
)
from importlib.metadata import version, PackageNotFoundError

__all__ = [
    "PhytClust",
    "OutlierConfig",
    "PeakConfig",
    "RuntimeConfig",
    "PlotConfig",
    "ScorePlotConfig",
    "ClusterPlotConfig",
    "PhytClustError",
    "ValidationError",
    "ConfigurationError",
    "ComputationError",
    "DataError",
    "InvalidKError",
    "InvalidTreeError",
    "MissingDPTableError",
    "InvalidClusteringError",
]

try:
    __version__ = version("phytclust")
except PackageNotFoundError:
    __version__ = "0.0.0"
