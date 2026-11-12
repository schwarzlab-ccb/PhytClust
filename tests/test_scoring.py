import logging
from types import SimpleNamespace

import numpy as np
import pytest
from Bio.Phylo.BaseTree import Tree

from conftest import SAMPLE_TREE, tree_from_newick
from phytclust import PhytClust
from phytclust.algo import scoring
from phytclust.algo.bins import define_bins
from phytclust.algo.scoring import _rank_peak_partitions, find_score_peaks
from phytclust.config import PeakConfig
from phytclust.exceptions import ConfigurationError, InvalidKError


def scoring_state(costs):
    tree = Tree()
    return SimpleNamespace(
        tree=tree,
        outgroup=None,
        node_to_id={tree.root: 0},
        unpenalised_dp_table=[np.array(costs, dtype=float)],
        dp_table=[np.array(costs, dtype=float)],
        max_k=len(costs),
        num_terminals=len(costs),
        k=None,
        _dp_cache_sig=("same",),
        _dp_cap=len(costs),
        _score_arrays_cache=None,
        _score_arrays_sig=None,
        _score_arrays_cap=None,
        use_penalized_beta_for_scoring=False,
    )


@pytest.mark.parametrize(
    "costs",
    [[10, 5, 0, float("inf")], [float("inf"), 5, 2, 0], [10, float("nan"), 3, 1]],
)
def test_scalar_and_vector_scores_agree(costs):
    pc = scoring_state(costs)
    vector = scoring._score_all_cluster_counts(pc)
    scalar = np.array(
        [scoring._score_cluster_count(pc, k=k) for k in range(1, len(costs) + 1)]
    )
    for index in range(3):
        np.testing.assert_allclose(vector[index], scalar[:, index], equal_nan=True)


def test_impossible_costs_do_not_create_elbows_or_shift_k():
    pc = scoring_state([10, float("inf"), 5, 2, 1])
    scoring.calculate_scores(pc)
    assert len(pc.scores) >= 4
    assert pc.scores[1] == 0
    assert pc.scores[2] == 0
    assert pc.scores[3] > 0
    assert np.isinf(pc.beta_values[1])


@pytest.mark.parametrize(
    "scores,start,end,min_k",
    [
        ([0, 10, 0, 4, 0, 0], 4, 5, 2),
        ([0, 10, 0, 0, 0, 0], 3, 6, 2),
        ([0, 10, 0, 0, 0, 0], 1, 6, 3),
        ([0, 10, 0, 0, 0, 0], 1, 1, 2),
    ],
)
def test_peaks_obey_range_and_minimum(monkeypatch, scores, start, end, min_k):
    pc = SimpleNamespace(scores=np.array(scores, dtype=float))
    peaks = scoring.find_score_peaks(
        pc, plot=False, k_start=start, k_end=end, peak_config=PeakConfig(min_k=min_k)
    )
    assert all(max(start, min_k) <= k <= end for k in peaks)
    if start == 4:
        assert peaks == [4]
    else:
        assert peaks == []


@pytest.mark.parametrize("resolution", [False, True])
@pytest.mark.parametrize("log_input", [False, True])
@pytest.mark.parametrize("relative", [False, True])
@pytest.mark.parametrize(
    "curve,expected",
    [
        ([420.6, 269.3, 141.3], []),
        ([82.2, 622.6, 186.6], [3]),
        ([10, 5], []),
        ([10], []),
        ([10, 10, 5], []),
        ([1, 10, 10, 1], [3]),
        ([1, 2, 8, 2, 1], [4]),
        ([1, 2, 3], []),
    ],
)
def test_interior_peaks_use_observed_k2_context(
    resolution, log_input, relative, curve, expected
):
    pc = SimpleNamespace(scores=np.array([0, *curve]), no_split_zero_length=True)
    peaks = scoring.find_score_peaks(
        pc,
        plot=False,
        resolution_on=resolution,
        num_bins=1,
        peak_config=PeakConfig(
            exclude_k2=True,
            use_log_peak_input=log_input,
            use_relative_prominence=relative,
        ),
    )
    assert peaks == expected


@pytest.mark.parametrize("resolution", [False, True])
@pytest.mark.parametrize("log_input", [False, True])
@pytest.mark.parametrize("relative", [False, True])
def test_eligibility_preserves_peak_prominence(resolution, log_input, relative):
    pc = SimpleNamespace(
        scores=np.array([0, 8, 10, 1, 6, 2]), no_split_zero_length=True
    )
    common = dict(
        use_log_peak_input=log_input,
        use_relative_prominence=relative,
        min_prominence=0,
    )

    def detect(**settings):
        scoring.find_score_peaks(
            pc,
            plot=False,
            resolution_on=resolution,
            num_bins=1,
            global_peaks=10,
            peak_config=PeakConfig(**common, **settings),
        )
        return {item["k"]: item["prominence"] for item in pc.peak_ranking_details}

    original = detect()
    assert set(original) == {3, 5}
    assert detect(exclude_k2=True) == original
    assert detect(min_k=4) == {5: original[5]}
    expected = 10 / 8 if relative else np.log(10 / 8) if log_input else 2
    assert original[3] == pytest.approx(expected)


@pytest.mark.parametrize("resolution", [False, True])
@pytest.mark.parametrize("restriction", [dict(min_k=4), dict(k_start=4)])
def test_restricted_boundary_does_not_create_peak(resolution, restriction):
    pc = SimpleNamespace(
        scores=np.array([0, 1, 10, 7, 2, 5, 1]), no_split_zero_length=True
    )
    settings = PeakConfig(min_k=restriction.get("min_k", 2))
    assert scoring.find_score_peaks(
        pc,
        plot=False,
        resolution_on=resolution,
        num_bins=1,
        k_start=restriction.get("k_start"),
        peak_config=settings,
    ) == [6]


@pytest.mark.parametrize("resolution", [False, True])
def test_k2_boundary_selection_is_preserved(resolution):
    pc = SimpleNamespace(
        scores=np.array([0, 420.6, 269.3, 141.3]), no_split_zero_length=True
    )
    assert scoring.find_score_peaks(pc, plot=False, resolution_on=resolution) == [2]


@pytest.mark.parametrize("resolution", [False, True])
@pytest.mark.parametrize(
    "settings",
    [
        {"boundary_ratio_threshold": 1000},
        {"use_relative_prominence": True, "min_relative_prominence": 100},
        {
            "boundary_ratio_threshold": 1000,
            "use_relative_prominence": True,
            "min_relative_prominence": 100,
        },
    ],
)
def test_rejected_k2_is_not_reintroduced_by_fallback(resolution, settings):
    pc = SimpleNamespace(scores=np.array([0, 10, 5, 1]), no_split_zero_length=True)
    assert (
        scoring.find_score_peaks(
            pc, plot=False, resolution_on=resolution, peak_config=PeakConfig(**settings)
        )
        == []
    )
    assert pc.peak_ranking_details == []


def test_only_nonzero_k2_still_obeys_relative_threshold():
    pc = SimpleNamespace(scores=np.array([0, 10, 0, 0]), no_split_zero_length=True)
    assert (
        scoring.find_score_peaks(
            pc,
            plot=False,
            peak_config=PeakConfig(
                use_relative_prominence=True, min_relative_prominence=1e20
            ),
        )
        == []
    )


@pytest.mark.parametrize("scores", [[0, 10], [0, 10, 0, 0]])
def test_only_nonzero_k2_remains_selectable(scores):
    pc = SimpleNamespace(scores=np.array(scores), no_split_zero_length=True)
    assert scoring.find_score_peaks(pc, plot=False) == [2]
    assert pc.peak_ranking_details[0]["k"] == 2


@pytest.mark.parametrize("min_k,expected", [(3, [8, 60]), (20, [60])])
def test_resolution_fallback_runs_when_entire_curve_has_no_peaks(min_k, expected):
    pc = SimpleNamespace(scores=np.arange(60, dtype=float), no_split_zero_length=True)
    assert (
        scoring.find_score_peaks(
            pc,
            plot=False,
            resolution_on=True,
            num_bins=2,
            peak_config=PeakConfig(min_k=min_k, resolution_fallback_mode="max_score"),
        )
        == expected
    )
    selected = [item for items in pc.resolution_info.values() for item in items]
    assert [item["k"] for item in selected] == expected
    assert all(item["selection_note"] == "fallback_max_score" for item in selected)


def test_resolution_fallback_disabled_leaves_peakless_curve_empty():
    pc = SimpleNamespace(scores=np.arange(60, dtype=float), no_split_zero_length=True)
    assert (
        scoring.find_score_peaks(
            pc,
            plot=False,
            resolution_on=True,
            peak_config=PeakConfig(resolution_fallback_mode="none"),
        )
        == []
    )


def test_resolution_fallback_respects_requested_range():
    pc = SimpleNamespace(scores=np.arange(60, dtype=float), no_split_zero_length=True)
    assert scoring.find_score_peaks(
        pc,
        plot=False,
        resolution_on=True,
        num_bins=2,
        k_start=10,
        k_end=30,
        peak_config=PeakConfig(resolution_fallback_mode="max_score"),
    ) == [30]


def test_resolution_fallback_reaches_public_run_and_backtracks_partitions(monkeypatch):
    from phytclust import PhytClust
    from phytclust.algo import core

    subtrees = [f"leaf{index}:1" for index in range(64)]
    while len(subtrees) > 1:
        subtrees = [
            f"({subtrees[index]},{subtrees[index + 1]}):1"
            for index in range(0, len(subtrees), 2)
        ]
    clustering = PhytClust(subtrees[0] + ";", max_k=60)

    def increasing_scores(state):
        state.scores = np.arange(60, dtype=float)

    monkeypatch.setattr(core, "calculate_scores", increasing_scores)
    result = clustering.run(
        by_resolution=True,
        num_bins=2,
        plot_scores=False,
        peak_config=PeakConfig(resolution_fallback_mode="max_score"),
    )
    assert result["mode"] == "resolution"
    assert len(result["k_values"]) == len(result["clusters"]) == 2
    for count, partition in zip(result["k_values"], result["clusters"]):
        assert len(partition) == 64
        assert len(set(partition.values())) == count
    assert all(
        item["selection_note"] == "fallback_max_score"
        for choices in clustering.resolution_info.values()
        for item in choices
    )


def test_relative_prominence_uses_original_score_baseline(monkeypatch):
    clustering = SimpleNamespace(scores=np.array([0.0, 0.0, 1.0, 10.0, 1.0, 1.0]))
    peaks = scoring.find_score_peaks(
        clustering,
        plot=False,
        peak_config=PeakConfig(
            use_log_peak_input=True,
            use_relative_prominence=True,
            min_relative_prominence=5,
        ),
    )
    assert peaks == [4]
    assert clustering.peak_ranking_details[0]["prominence"] == pytest.approx(10)


@pytest.mark.parametrize(
    "settings, field_name",
    [
        ({"ranking_mode": "other"}, "ranking_mode"),
        ({"resolution_fallback_mode": "other"}, "resolution_fallback_mode"),
        ({"prominence_weight": np.nan}, "prominence_weight"),
        ({"prominence_weight": 1.1}, "prominence_weight"),
        ({"prominence_weight": "0.7"}, "prominence_weight"),
        ({"boundary_window_size": 2.5}, "boundary_window_size"),
        ({"boundary_window_size": True}, "boundary_window_size"),
        ({"boundary_window_size": 0}, "boundary_window_size"),
        ({"boundary_ratio_threshold": np.inf}, "boundary_ratio_threshold"),
        ({"min_k": 0}, "min_k"),
        ({"min_k": 2.5}, "min_k"),
        ({"min_prominence": -1}, "min_prominence"),
        ({"use_log_peak_input": True, "log_peak_offset": 0}, "log_peak_offset"),
        ({"use_log_peak_input": True, "log_peak_offset": np.nan}, "log_peak_offset"),
        (
            {"use_relative_prominence": True, "min_relative_prominence": -1},
            "min_relative_prominence",
        ),
        ({"min_prominence": 1.0, "prominence_k_power": np.inf}, "prominence_k_power"),
        ({"use_relative_prominence": "false"}, "use_relative_prominence"),
    ],
)
def test_invalid_active_settings(settings, field_name):
    with pytest.raises(ConfigurationError, match=field_name):
        PeakConfig(**settings)


@pytest.mark.parametrize(
    "settings",
    [
        {"ranking_mode": "raw", "prominence_weight": np.nan},
        {
            "exclude_k2": True,
            "boundary_window_size": 0,
            "boundary_ratio_threshold": np.nan,
        },
        {"min_k": 3, "boundary_window_size": 0},
        {"use_log_peak_input": False, "log_peak_offset": np.nan},
        {"use_log_peak_input": False, "log_peak_offset": "unused"},
        {"use_relative_prominence": False, "min_relative_prominence": np.nan},
        {"min_prominence": None, "prominence_k_power": np.nan},
        {
            "use_relative_prominence": True,
            "min_prominence": 0.0,
            "prominence_k_power": np.nan,
        },
    ],
)
def test_unused_settings_are_ignored(settings, monkeypatch):
    config = PeakConfig(**settings)
    clustering = SimpleNamespace(scores=np.array([0.0, 0.0, 1.0, 10.0, 1.0, 1.0]))
    assert scoring.find_score_peaks(clustering, peak_config=config, plot=False) == [4]


def test_changes_after_creation_are_validated(monkeypatch):
    config = PeakConfig()
    config.min_k = 0
    with pytest.raises(ConfigurationError, match="min_k"):
        scoring.find_score_peaks(
            SimpleNamespace(scores=[]), peak_config=config, plot=False
        )


def test_numpy_numbers_and_zero_thresholds():
    PeakConfig(
        min_k=np.int64(2),
        boundary_window_size=np.int64(3),
        min_prominence=0.0,
        prominence_weight=np.float64(0.5),
        boundary_ratio_threshold=0.0,
    )


@pytest.mark.parametrize("log_input", [False, True])
def test_relative_ratio_uses_higher_of_the_two_bases(monkeypatch, log_input):
    clustering = SimpleNamespace(scores=np.array([0.0, 0.0, 2.0, 10.0, 4.0, 4.0]))
    scoring.find_score_peaks(
        clustering,
        plot=False,
        peak_config=PeakConfig(
            use_relative_prominence=True, use_log_peak_input=log_input
        ),
    )
    assert clustering.peak_ranking_details[0]["prominence"] == pytest.approx(2.5)


def test_unused_relative_threshold_does_not_filter_boundary_peak(monkeypatch):
    clustering = SimpleNamespace(scores=np.array([0.0, 10.0, 1.0, 5.0, 1.0, 1.0]))
    peaks = scoring.find_score_peaks(
        clustering,
        plot=False,
        global_peaks=2,
        peak_config=PeakConfig(min_relative_prominence=np.nan),
    )
    assert 2 in peaks


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


class _PC:
    def __init__(self, n):
        self.num_terminals = n


def test_define_bins_basic_properties():
    pc = _PC(100)
    bins = define_bins(pc, num_bins=3, k_lo=1, k_hi=100)

    assert len(bins) >= 1
    assert bins[0][0] == 1
    assert bins[-1][1] == 100

    for (lo1, hi1), (lo2, hi2) in zip(bins, bins[1:]):
        assert lo1 <= hi1
        assert lo2 <= hi2
        assert hi1 + 1 == lo2

    assert all(1 <= lo <= hi <= 100 for lo, hi in bins)


def test_define_bins_rejects_invalid_range():
    pc = _PC(10)
    with pytest.raises(InvalidKError, match="k_hi must be > k_lo"):
        define_bins(pc, num_bins=3, k_lo=5, k_hi=5)


@pytest.mark.parametrize("parameter", ["num_bins", "k_lo", "k_hi"])
@pytest.mark.parametrize("value", [0, -1, 1.5, True, float("nan"), float("inf"), "3"])
def test_bin_parameters_require_positive_integers(parameter, value):
    from phytclust.exceptions import ConfigurationError
    arguments = {"num_bins": 3, "k_lo": 1, "k_hi": 10}
    arguments[parameter] = value
    error = ConfigurationError if parameter == "num_bins" else InvalidKError
    with pytest.raises(error, match=parameter):
        define_bins(_PC(10), **arguments)


def test_rounding_reduces_bin_count_without_gaps():
    bins = define_bins(_PC(3), num_bins=10)
    assert len(bins) < 10
    assert [k for lower, upper in bins for k in range(lower, upper + 1)] == [1, 2, 3]
    assert all(type(boundary) is int for pair in bins for boundary in pair)


def test_custom_range_is_covered_exactly():
    bins = define_bins(_PC(100), num_bins=4, k_lo=3, k_hi=20)
    assert [k for lower, upper in bins for k in range(lower, upper + 1)] == list(range(3, 21))


def test_resolution_min_scores_is_validated():
    with pytest.raises(Exception, match="resolution_min_scores"):
        PeakConfig(resolution_min_scores=0)


def test_zero_distance_warning_ignores_separated_leaves(caplog):
    tree = "(((A:5,B:5):0,(C:1,D:1):3):1,((E:1,F:1):3,(G:1,H:1):3):1);"
    pc = PhytClust(tree)
    with caplog.at_level(logging.WARNING, logger="phytclust"):
        pc.run(top_n=3, plot_scores=False)
    assert "zero distance" not in caplog.text
    assert pc.zero_length_splits == {}


def test_zero_distance_warning_reports_split_identical_leaves(caplog):
    tree = "(((A:0,B:0):2,(C:1,D:1):3):1,((E:1,F:1):3,(G:1,H:1):3):1);"
    pc = PhytClust(tree)
    with caplog.at_level(logging.WARNING, logger="phytclust"):
        pc.run(k=8)
    assert pc.zero_length_splits == {8: 1}
    assert "'A' and 'B'" in caplog.text
    pc.run(k=3)
    assert pc.zero_length_splits == {}


def test_no_peaks_is_a_warning(caplog):
    with caplog.at_level(logging.WARNING, logger="phytclust"):
        result = PhytClust(SAMPLE_TREE, min_cluster_size=3).run(top_n=3, plot_scores=False)
    assert result["k_values"] == []
    assert "no partition was selected" in caplog.text


def test_k2_is_a_candidate_by_default():
    """Scoring selects k* from k in {2, ..., N-1}."""
    assert PeakConfig().exclude_k2 is False
    newick = "(((A:1,B:1):1,(C:1,D:1):1):40,((E:1,F:1):1,(G:1,H:1):1):40);"

    result = PhytClust(tree=tree_from_newick(newick)).run(top_n=1, plot_scores=False)
    assert result["k_values"] == [2]

    excluded = PhytClust(tree=tree_from_newick(newick)).run(
        top_n=1, plot_scores=False, peak_config=PeakConfig(exclude_k2=True)
    )
    assert 2 not in excluded["k_values"]
