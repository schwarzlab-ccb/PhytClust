"""Regression checks for scoring caches, invalid costs, and peak ranges."""

from types import SimpleNamespace
import numpy as np
import pytest
from Bio.Phylo.BaseTree import Tree
from phytclust.algo import scoring
from phytclust.config import PeakConfig


def scoring_state(costs):
    tree = Tree()
    return SimpleNamespace(
        tree=tree,
        outgroup=None,
        node_to_id={tree.root: 0},
        raw_dp_table=[np.array(costs, dtype=float)],
        dp_table=[np.array(costs, dtype=float)],
        max_k=len(costs),
        num_terminals=len(costs),
        k=None,
        _dp_cache_sig=("same",),
        _dp_cap=len(costs),
        _score_raw_arrays=None,
        _score_raw_base_sig=None,
        _score_raw_cap=None,
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


@pytest.mark.parametrize(
    "floor_name", ["score_beta_floor_abs", "score_beta_floor_frac"]
)
def test_cost_floor_changes_invalidate_raw_score_cache(floor_name):
    pc = scoring_state([10, 5, 2, 1])
    before = scoring._cached_cluster_count_scores(pc)[2]
    setattr(pc, floor_name, 100)
    after = scoring._cached_cluster_count_scores(pc)[2]
    assert not np.array_equal(before, after)
    np.testing.assert_array_equal(after, scoring._score_all_cluster_counts(pc)[2])


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
    monkeypatch.setattr(scoring, "_first_zero_length_pair_split", lambda *args: None)
    pc = SimpleNamespace(scores=np.array(scores, dtype=float))
    peaks = scoring.find_score_peaks(
        pc, plot=False, k_start=start, k_end=end, peak_config=PeakConfig(min_k=min_k)
    )
    assert all(max(start, min_k) <= k <= end for k in peaks)
    if start == 4:
        assert peaks == [4]
    else:
        assert peaks == []
