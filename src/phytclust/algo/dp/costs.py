from typing import Optional
import logging
import math

import numpy as np

from ...validation import positive_int
from ...exceptions import ConfigurationError, InvalidKError, InvalidTreeError
from ...utils.traversal import iter_clades

logger = logging.getLogger("phytclust")


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
    for name, value in (("k", clustering.k), ("max_k", clustering.max_k)):
        if value is not None:
            positive_int(value, name, error=InvalidKError)
    positive_int(clustering.min_cluster_size, "min_cluster_size")
    positive_int(clustering.max_tied_optima, "max_tied_optima")
    clustering.outlier.validate()
    if not _finite_number(clustering.max_k_limit):
        raise ConfigurationError("max_k_limit must be a finite number.")
    if clustering.use_branch_support:
        if not _finite_number(clustering.support_weight) or clustering.support_weight < 0:
            raise ConfigurationError("support_weight must be a finite number of zero or greater.")
        if not _finite_number(clustering.min_support) or not 0 < clustering.min_support <= 1:
            raise ConfigurationError("min_support must be greater than zero and at most 1.")
        for node in iter_clades(clustering.tree.root):
            branch_support_fraction(node, clustering)
    if not 0 < clustering.max_k_limit <= 1:
        raise ConfigurationError("max_k_limit must be between 0 and 1")

    polytomy_mode = clustering.polytomy_mode
    if polytomy_mode not in {"hard", "soft"}:
        raise ConfigurationError("polytomy_mode must be 'hard' or 'soft'.")

    if polytomy_mode == "soft":
        max_deg = clustering.soft_polytomy_max_degree
        positive_int(max_deg, "soft_polytomy_max_degree", minimum=2)


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
        outlier_config.penalty_enabled
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


def detect_support_scale(tree) -> float:
    """Return 1.0 when every branch support in ``tree`` is at most 1, else 100.0.

    Bootstrap percentages (0-100) and posterior probabilities (0-1) are both
    common; a tree whose supports never exceed 1 is read as 0-1.
    """
    supports = [
        clade.confidence
        for clade in iter_clades(tree.root)
        if clade.confidence is not None and _finite_number(clade.confidence)
    ]
    if supports and max(supports) <= 1.0:
        return 1.0
    return 100.0


def branch_support_fraction(node, clustering) -> float:
    """Convert a clade's support to a fraction; missing support means full support.

    The scale comes from detect_support_scale, run once per DP build.
    """
    if not clustering.use_branch_support:
        return 1.0
    scale = clustering._support_scale or 100.0
    confidence = node.confidence if node.confidence is not None else scale
    if not _finite_number(confidence) or not 0 <= confidence <= scale:
        raise InvalidTreeError(
            f"Confidence for node {node.name!r} must be between 0 and {scale:g}."
        )
    return max(confidence / scale, clustering.min_support)


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

    With prefer_fewer, outlier counts are compared before cost; otherwise cost
    first, then outlier count, then the more even spread of clusters across
    children. Costs within the relative tolerance are tied. Used at polytomies,
    where small-cluster penalties do not apply.
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
