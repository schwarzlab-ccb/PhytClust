from typing import Optional
import math
from numbers import Integral

import numpy as np

from ...exceptions import ConfigurationError, InvalidKError, InvalidTreeError
from ...utils.traversal import iter_clades


_EPS_CACHE: dict = {}


def floating_point_precision(dtype) -> float:
    """Return the smallest relative spacing near 1 for a floating-point type."""
    precision = _EPS_CACHE.get(dtype)
    if precision is None:
        precision = float(np.finfo(dtype).eps)
        _EPS_CACHE[dtype] = precision
    return precision


def validate_clustering_parameters(clustering) -> None:
    """Check counts, cost settings, and support values used by clustering."""
    def check_count(name, value, minimum=1, error=ConfigurationError):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value < minimum:
            raise error(f"{name} must be a whole number of at least {minimum}.")

    for name in ("k", "max_k"):
        value = getattr(clustering, name, None)
        if value is not None:
            check_count(name, value, error=InvalidKError)
    check_count("min_cluster_size", clustering.min_cluster_size)
    check_count("max_tied_optima", getattr(clustering, "max_tied_optima", 100))
    if clustering.outlier.size_threshold is not None:
        check_count("outlier.size_threshold", clustering.outlier.size_threshold)
    if not _finite_number(clustering.max_k_limit):
        raise ConfigurationError("max_k_limit must be a finite number.")
    if not _finite_number(clustering.outlier.ratio_weight):
        raise ConfigurationError("outlier.ratio_weight must be a finite number.")
    if clustering.use_branch_support:
        if not _finite_number(clustering.support_weight) or clustering.support_weight < 0:
            raise ConfigurationError("support_weight must be a finite number of zero or greater.")
        if not _finite_number(clustering.min_support) or not 0 < clustering.min_support <= 1:
            raise ConfigurationError("min_support must be greater than zero and at most 1.")
        for node in iter_clades(clustering.tree.root):
            branch_support_fraction(node, clustering)
    if clustering.k is not None and clustering.k < 1:
        raise InvalidKError("k must be ≥ 1 if provided.")
    if not 0 < clustering.max_k_limit <= 1:
        raise ConfigurationError("max_k_limit must be between 0 and 1")

    if clustering.outlier.ratio_weight < 0:
        raise ConfigurationError("outlier.ratio_weight must be ≥ 0")

    if clustering.outlier.ratio_mode not in {"exp", "inverse", "power"}:
        raise ConfigurationError(
            "outlier.ratio_mode must be one of: 'exp', 'inverse', 'power'"
        )

    if clustering.outlier.prefer_fewer and clustering.outlier.size_threshold is None:
        raise ConfigurationError(
            "outlier.prefer_fewer=True requires outlier.size_threshold to be set."
        )

    polytomy_mode = getattr(clustering, "polytomy_mode", "soft")
    if polytomy_mode not in {"hard", "soft"}:
        raise ConfigurationError("polytomy_mode must be 'hard' or 'soft'.")

    if polytomy_mode == "soft":
        max_deg = getattr(clustering, "soft_polytomy_max_degree", 12)
        check_count("soft_polytomy_max_degree", max_deg, minimum=2)


def _finite_number(value) -> bool:
    """Return whether a value is a finite number."""
    try:
        return math.isfinite(value)
    except (TypeError, ValueError, OverflowError):
        return False


def small_cluster_penalty_enabled(clustering) -> bool:
    """Return whether small-cluster penalties are enabled and have nonzero weight."""
    outlier_config = clustering.outlier
    return bool(
        getattr(outlier_config, "penalty_enabled", False)
        and outlier_config.size_threshold is not None
        and outlier_config.ratio_weight > 0
    )


def small_cluster_penalty(cluster_size: int, clustering) -> float:
    """Return the weighted penalty for a cluster below the size threshold."""
    outlier_config = clustering.outlier
    size_threshold = outlier_config.size_threshold
    if not small_cluster_penalty_enabled(clustering) or cluster_size >= size_threshold:
        return 0.0
    cluster_size = max(1, int(cluster_size))
    mode = outlier_config.ratio_mode
    if mode == "inverse":
        shape = 1.0 / cluster_size
    elif mode == "exp":
        with np.errstate(over="raise", invalid="raise"):
            try:
                shape = float(np.expm1(size_threshold - cluster_size))
            except (FloatingPointError, OverflowError) as exc:
                raise ConfigurationError(
                    "The exponential outlier penalty is too large. "
                    "Reduce size_threshold or use inverse or power mode."
                ) from exc
    else:
        shape = float(size_threshold - cluster_size)
    penalty = float(outlier_config.ratio_weight) * shape
    if not math.isfinite(penalty):
        raise ConfigurationError("The outlier penalty is too large. Reduce ratio_weight or size_threshold.")
    return penalty


def branch_support_fraction(node, clustering) -> float:
    """Convert confidence percentages to fractions; missing confidence means full support."""
    if not getattr(clustering, "use_branch_support", False):
        return 1.0
    confidence = node.confidence if node.confidence is not None else 100.0
    if not _finite_number(confidence) or not 0 <= confidence <= 100:
        raise InvalidTreeError(f"Confidence for node {node.name!r} must be between 0 and 100.")
    return max(confidence / 100.0, clustering.min_support)


def is_better_partition(
    candidate_cost: float,
    candidate_outlier_count: int,
    best_cost: float,
    best_outlier_count: int,
    *,
    outlier_detection_enabled: bool,
    prioritize_fewer_outliers: bool,
    relative_cost_tolerance: float = 0.0,
    candidate_cluster_count_sum_of_squares: Optional[int] = None,
    best_cluster_count_sum_of_squares: Optional[int] = None,
) -> bool:
    """Return whether the candidate partition should replace the current best.

    When requested, compare outlier counts before costs. Otherwise compare
    costs first and use outlier counts to break ties. For remaining ties,
    prefer a more even distribution of clusters among child groups.

    This comparator is used for polytomies, where small-cluster penalties
    are disabled. Binary merges with penalties compare penalized costs first.
    Treat costs as tied within the relative cost tolerance.
    """
    if not np.isfinite(candidate_cost):
        return False
    if not np.isfinite(best_cost):
        return True

    cost_tolerance = relative_cost_tolerance * max(abs(candidate_cost), abs(best_cost))
    tied = abs(candidate_cost - best_cost) <= cost_tolerance

    if outlier_detection_enabled and prioritize_fewer_outliers:
        if candidate_outlier_count != best_outlier_count:
            return candidate_outlier_count < best_outlier_count
        if candidate_cost < best_cost - cost_tolerance:
            return True
        if tied and candidate_cluster_count_sum_of_squares is not None and best_cluster_count_sum_of_squares is not None:
            return candidate_cluster_count_sum_of_squares < best_cluster_count_sum_of_squares
        return False

    if candidate_cost < best_cost - cost_tolerance:
        return True

    if tied:
        if outlier_detection_enabled and candidate_outlier_count != best_outlier_count:
            return candidate_outlier_count < best_outlier_count
        if candidate_cluster_count_sum_of_squares is not None and best_cluster_count_sum_of_squares is not None:
            return candidate_cluster_count_sum_of_squares < best_cluster_count_sum_of_squares

    return False


def child_subtree_cluster_cost(child, clustering) -> float:
    """
    Return the cost of one cluster containing the child subtree at its parent.
    """
    child_cost = clustering.cluster_cost[child]
    if not np.isfinite(child_cost):
        return np.inf
    return child_cost + clustering.num_leaves_per_node[child] * effective_branch_length(child, clustering)


def effective_branch_length(node, clustering):
    """
    Return branch length plus the support penalty when support is enabled.

    The penalty is ``support_weight * -log(support_fraction)``.
    """
    base = node.branch_length or 0.0
    if not clustering.use_branch_support:
        return base

    support_fraction = branch_support_fraction(node, clustering)
    penalty = -np.log(support_fraction)
    return base + clustering.support_weight * penalty


def subtree_has_only_zero_length_branches(node, clustering, zero_length_tolerance: float = 1e-12):
    """Return whether every branch below the node is effectively zero length."""
    stack = list(node.clades)
    while stack:
        child = stack.pop()
        if effective_branch_length(child, clustering) > zero_length_tolerance:
            return False
        if not child.is_terminal():
            stack.extend(child.clades)
    return True
