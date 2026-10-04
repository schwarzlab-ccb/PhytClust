"""Regression checks for export selection and output metadata."""

import pandas as pd
import pytest
from Bio.Phylo.BaseTree import Clade
from phytclust import PhytClust
from phytclust.exceptions import (
    ConfigurationError,
    InvalidClusteringError,
    InvalidKError,
)


def clustering():
    return PhytClust("((a:1,b:2):3,(c:4,d:5):6);", max_k=4)


def test_explicit_count_overrides_previous_run(tmp_path):
    pc = clustering()
    pc.run(k=2, plot_scores=False)
    path = pc.save(str(tmp_path), k=3, outlier=False)
    assert list(pd.read_csv(path, sep="\t").columns) == ["Node Name", "clusters_k3"]
    assert pd.read_csv(tmp_path / "alphas.tsv", sep="\t")["k"].tolist() == [3]
    assert (tmp_path / "peaks_by_rank.txt").read_text() == ""


@pytest.mark.parametrize(
    "filename",
    ["alphas.tsv", "peaks_by_rank.txt", "../results.tsv", "a/b.tsv", "a\\b.tsv", ""],
)
def test_invalid_filenames_rejected_without_creating_directory(tmp_path, filename):
    destination = tmp_path / "new"
    with pytest.raises(ValueError):
        clustering().save(str(destination), filename=filename)
    assert not destination.exists()


@pytest.mark.parametrize("count", [0, -1, 1.5, True])
def test_invalid_explicit_counts_rejected(tmp_path, count):
    with pytest.raises(InvalidKError):
        clustering().save(str(tmp_path), k=count)


@pytest.mark.parametrize("count", [0, -1, 1.5, True])
def test_invalid_peak_limits_rejected(tmp_path, count):
    with pytest.raises(ConfigurationError):
        clustering().save(str(tmp_path), top_n=count)


def test_no_selection_does_not_create_directory(tmp_path):
    destination = tmp_path / "new"
    assert clustering().save(str(destination)) is None
    assert not destination.exists()


def test_peak_file_only_lists_exported_peaks(tmp_path):
    pc = clustering()
    pc.peaks_by_rank = [3, 2]
    pc.save(str(tmp_path), top_n=1)
    assert (tmp_path / "peaks_by_rank.txt").read_text() == "Rank 1: 3 clusters\n"
    pc.save(str(tmp_path), k=2)
    assert (tmp_path / "peaks_by_rank.txt").read_text() == ""


def test_duplicate_leaf_names_are_not_silently_merged(tmp_path, monkeypatch):
    pc = clustering()
    monkeypatch.setattr(
        pc, "_clusters", lambda count: {Clade(name="same"): 0, Clade(name="same"): 1}
    )
    monkeypatch.setattr(pc, "alpha_info", lambda *args: {"k": 2, "alpha": 1})
    with pytest.raises(InvalidClusteringError, match="duplicate"):
        pc.save(str(tmp_path), k=2)
