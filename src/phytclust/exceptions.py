"""Custom exceptions for PhytClust."""


class PhytClustError(Exception):
    """Base exception for all PhytClust errors."""
    pass


class ValidationError(PhytClustError):
    """Raised when input validation fails."""
    pass


class ConfigurationError(PhytClustError):
    """Raised when configuration parameters are invalid."""
    pass


class ComputationError(PhytClustError):
    """Raised when an internal computation fails."""
    pass


class DataError(PhytClustError):
    """Raised when data is malformed or missing."""
    pass


class InvalidKError(ConfigurationError):
    """Raised when k value is invalid."""
    pass


class InvalidTreeError(DataError):
    """Raised when a tree or a requested node is invalid or unavailable."""
    pass


class MissingDPTableError(ComputationError):
    """Raised when required dynamic programming state is missing or invalid."""
    pass

class InvalidClusteringError(ComputationError):
    """Raised when a requested clustering or required result is unavailable."""
    pass
