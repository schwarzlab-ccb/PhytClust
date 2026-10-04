"""Checks for clustering costs and parameter validation."""

import numpy as np
import pytest
from phytclust import PhytClust
from phytclust.config import OutlierConfig
from phytclust.exceptions import ConfigurationError, InvalidTreeError
from phytclust.algo.dp.costs import (
    floating_point_precision, is_better_partition, small_cluster_penalty,
)


def test_precision_is_relative_spacing():
    assert floating_point_precision(np.float64) == np.finfo(np.float64).eps


@pytest.mark.parametrize('weight', [-1, float('nan'), float('inf')])
def test_support_settings_checked_only_when_enabled(weight):
    pc = PhytClust('(a:1,b:2);', support_weight=weight)
    pc.tree.root.confidence = 200
    pc.get_clusters(1)
    pc.use_branch_support = True
    with pytest.raises(ConfigurationError, match='support_weight'):
        pc.get_clusters(1)


@pytest.mark.parametrize('confidence', [-1, 101, float('nan'), float('inf')])
def test_invalid_confidence_rejected_when_support_enabled(confidence):
    pc = PhytClust('(a:1,b:2);', use_branch_support=True)
    pc.tree.root.confidence = confidence
    with pytest.raises(InvalidTreeError, match='Confidence'):
        pc.get_clusters(1)


@pytest.mark.parametrize('name,value', [('min_cluster_size', 1.5), ('soft_polytomy_max_degree', True), ('min_support', 0)])
def test_invalid_parameters(name, value):
    pc = PhytClust('(a:1,b:2);', use_branch_support=True, **{name: value})
    with pytest.raises(ConfigurationError):
        pc.get_clusters(1)


def test_overflowing_penalty_has_clear_error():
    pc = PhytClust('(a:1,b:2);', outlier=OutlierConfig(size_threshold=1000, penalty_enabled=True))
    with pytest.raises(ConfigurationError, match='exponential outlier penalty'):
        small_cluster_penalty(1, pc)
    pc.outlier.penalty_enabled = False
    assert small_cluster_penalty(1, pc) == 0


def test_comparison_priorities_and_relative_ties():
    assert is_better_partition(10, 0, 1, 1, outlier_detection_enabled=True, prioritize_fewer_outliers=True)
    assert not is_better_partition(10, 0, 1, 1, outlier_detection_enabled=True, prioritize_fewer_outliers=False)
    for scale in [1, 1000]:
        assert is_better_partition(1.000001 * scale, 0, 1 * scale, 1,
                                  outlier_detection_enabled=True, prioritize_fewer_outliers=False,
                                  relative_cost_tolerance=1e-5)
