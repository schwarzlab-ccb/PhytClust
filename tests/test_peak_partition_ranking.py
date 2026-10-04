"""Partition preferences reorder peaks without changing scores or partitions."""

from types import SimpleNamespace

import numpy as np
import pytest

from phytclust import PhytClust
from phytclust.config import PeakConfig
from phytclust.exceptions import ConfigurationError
from phytclust.algo.scoring import _rank_peak_partitions, find_score_peaks


def test_balanced_ranking_penalizes_uneven_sizes_and_all_singletons():
    partitions = {
        2: dict(enumerate([0] + [1] * 99)),
        3: dict(enumerate([0] * 34 + [1] * 33 + [2] * 33)),
        100: dict(enumerate(range(100))),
    }
    clustering = SimpleNamespace(_clusters=lambda count: partitions[count])
    candidates = [{"k": count, "combined_metric": 1.0} for count in partitions]
    ranked = _rank_peak_partitions(
        clustering,
        candidates,
        PeakConfig(partition_preference="balanced", partition_weight=1),
    )
    assert [item["k"] for item in ranked] == [3, 2, 100]
    assert ranked[-1]["partition_quality"] == 0
    assert ranked[0]["size_balance"] > 0.99


def test_fewer_outliers_uses_singleton_fraction():
    clustering = SimpleNamespace(
        _clusters=lambda count: dict(
            enumerate([0] * (100 - count + 1) + list(range(1, count)))
        )
    )
    ranked = _rank_peak_partitions(
        clustering,
        [{"k": 2, "combined_metric": 1.0}, {"k": 4, "combined_metric": 1.0}],
        PeakConfig(partition_preference="fewer_outliers", partition_weight=1),
    )
    assert [item["singleton_count"] for item in ranked] == [1, 3]
    assert ranked[0]["partition_quality"] == pytest.approx(0.99)


@pytest.mark.parametrize(
    "settings",
    [PeakConfig(), PeakConfig(partition_preference="balanced", partition_weight=0)],
)
def test_disabled_preference_does_not_reconstruct_partitions(settings):
    def fail(count):
        pytest.fail("Disabled preference reconstructed a partition.")

    candidates = [{"k": 2, "combined_metric": 0.7}]
    assert (
        _rank_peak_partitions(SimpleNamespace(_clusters=fail), candidates, settings)
        is candidates
    )


@pytest.mark.parametrize(
    "settings",
    [
        dict(partition_preference="bad"),
        dict(partition_preference="balanced", partition_weight=-1),
        dict(partition_preference="balanced", partition_weight=1.1),
        dict(partition_preference="balanced", partition_weight=np.nan),
    ],
)
def test_invalid_partition_preferences(settings):
    with pytest.raises(ConfigurationError, match="partition_"):
        PeakConfig(**settings)


@pytest.mark.parametrize("resolution", [False, True])
def test_ranking_preserves_detected_candidates_scores_and_partitions(resolution):
    clustering = PhytClust("((a:1,b:1):1,(c:1,d:1):1,(e:1,f:1):1,g:1);", max_k=7)
    clustering.get_clusters(1)
    clustering.scores = np.array([0.0, 2.0, 1.0, 3.0, 1.0, 4.0, 1.0])
    original = clustering.scores.copy()
    find_score_peaks(
        clustering,
        plot=False,
        global_peaks=5,
        resolution_on=resolution,
        peak_config=PeakConfig(exclude_k2=True),
    )
    candidates = {item["k"] for item in clustering.peak_ranking_details}
    partitions = {
        count: {
            leaf.name: cluster
            for leaf, cluster in clustering.get_clusters(count).items()
        }
        for count in candidates
    }
    find_score_peaks(
        clustering,
        plot=False,
        global_peaks=5,
        resolution_on=resolution,
        peak_config=PeakConfig(exclude_k2=True, partition_preference="balanced"),
    )
    assert {item["k"] for item in clustering.peak_ranking_details} == candidates
    np.testing.assert_array_equal(clustering.scores, original)
    for count in candidates:
        assert {
            leaf.name: cluster
            for leaf, cluster in clustering.get_clusters(count).items()
        } == partitions[count]


@pytest.mark.parametrize("values", [np.zeros(5), np.array([])])
def test_no_peaks_clear_previous_ranking_details(values):
    clustering = SimpleNamespace(
        scores=values, peak_ranking_details=[{"k": 3}], no_split_zero_length=True
    )
    assert (
        find_score_peaks(
            clustering, plot=False, peak_config=PeakConfig(exclude_k2=True)
        )
        == []
    )
    assert clustering.peak_ranking_details == []
