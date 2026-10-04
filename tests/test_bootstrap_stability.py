"""Checks for bootstrap pair consistency and independent replicate clustering."""

from io import StringIO
import numpy as np
import pytest
from Bio import Phylo
from phytclust.algo.bootstrap_stability import (
    stability_for_k,
    compute_coassoc_for_k,
    choose_k_by_stability,
    _pairwise_coassociation,
)
from phytclust.exceptions import DataError, ConfigurationError


def tree():
    return Phylo.read(StringIO("((a:0.123456789,b:2)X:3,(c:4,d:5)Y:6);"), "newick")


@pytest.mark.parametrize(
    "frequency,expected", [(0, 1), (0.1, 0.8), (0.5, 0), (0.9, 0.8), (1, 1)]
)
def test_consistency_rewards_grouping_and_separation(frequency, expected):
    assert stability_for_k(np.array([[1, frequency], [frequency, 1]])) == pytest.approx(
        expected
    )


def test_unobserved_pairs_are_excluded():
    matrix = np.array([[1, 0.5, 0], [0.5, 1, 0], [0, 0, 1]])
    counts = np.array([[2, 2, 0], [2, 2, 0], [0, 0, 0]])
    assert stability_for_k(matrix, counts) == 0
    assert stability_for_k(np.ones((1, 1))) == 0


def test_coassociation_matches_direct_pair_counts():
    labels = np.array([[1_000_000_000, 1_000_000_000, -1], [0, 1, 1], [5, 5, 7]])
    result = _pairwise_coassociation(labels)
    for first in range(3):
        for second in range(3):
            valid = (labels[:, first] >= 0) & (labels[:, second] >= 0)
            expected = (
                np.mean(labels[valid, first] == labels[valid, second])
                if valid.any()
                else 0
            )
            assert result[first, second] == expected


@pytest.mark.parametrize("jobs", [0, -2, True, 1.5])
def test_invalid_worker_counts(jobs):
    with pytest.raises(ConfigurationError, match="n_jobs"):
        compute_coassoc_for_k([tree()], 2, n_jobs=jobs)


def test_empty_and_duplicate_taxa_rejected():
    with pytest.raises(DataError, match="No bootstrap trees"):
        compute_coassoc_for_k([], 2)
    duplicate = tree()
    duplicate.get_terminals()[0].name = "b"
    with pytest.raises(DataError, match="duplicate"):
        compute_coassoc_for_k([duplicate], 2)


def test_serial_clustering_preserves_input_tree():
    original = tree()
    original.root.name = None
    original.root.extra = {"values": [1]}
    taxa, labels, matrix = compute_coassoc_for_k([original], 2)
    assert original.root.name is None
    assert original.root.extra == {"values": [1]}
    assert original.find_any(name="a").branch_length == 0.123456789
    assert labels.shape == (1, 4)
    assert matrix.shape == (4, 4)
    assert taxa == ["a", "b", "c", "d"]


def test_serial_and_parallel_results_agree():
    trees = [tree(), tree()]
    serial = compute_coassoc_for_k(trees, 2, n_jobs=1)
    parallel = compute_coassoc_for_k(trees, 2, n_jobs=2)
    assert serial[0] == parallel[0]
    np.testing.assert_array_equal(serial[1], parallel[1])
    np.testing.assert_array_equal(serial[2], parallel[2])


def test_candidates_reuse_dp_tables(monkeypatch):
    from phytclust.algo import core

    real = core.compute_dp_table
    calls = []

    def counted(clustering):
        calls.append(clustering)
        return real(clustering)

    monkeypatch.setattr(core, "compute_dp_table", counted)
    result = choose_k_by_stability([tree(), tree()], [2, 3], n_jobs=1)
    assert len(calls) == 2
    assert set(result["scores"]) == {2, 3}
    assert result["best_k"] == 2


@pytest.mark.parametrize(
    "matrix",
    [
        np.ones(3),
        np.ones((2, 3)),
        np.array([[1, 2], [2, 1]]),
        np.array([[1, float("nan")], [float("nan"), 1]]),
    ],
)
def test_invalid_coassociation_rejected(matrix):
    with pytest.raises(DataError):
        stability_for_k(matrix)
