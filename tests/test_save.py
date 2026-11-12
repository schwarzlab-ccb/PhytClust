import pandas as pd
import pytest
from Bio.Phylo.BaseTree import Clade

from conftest import POLYTOMY_TREE, SAMPLE_TREE
from phytclust import PhytClust, RuntimeConfig
from phytclust.config.runtime import SaveConfig
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


def test_save_all_skips_infeasible_k_but_keeps_feasible_ones(tmp_path):
    pc = PhytClust("(a:1,b:1,c:1,d:1);", polytomy_mode="hard", max_k=4)

    path = pc.save(str(tmp_path), output_all=True)

    assert path is not None
    header = (tmp_path / "phytclust_results.tsv").read_text().splitlines()[0]
    assert header == "Node Name\tclusters_k1\tclusters_k4"


def test_save_all_counts_keeps_peak_ranks(tmp_path):
    pc = PhytClust(POLYTOMY_TREE)
    selected = pc.run(top_n=2, plot_scores=False)["k_values"]
    pc.save(str(tmp_path), top_n=2, output_all=True)
    ranks = (tmp_path / "peaks_by_rank.txt").read_text().splitlines()
    assert ranks == [f"Rank {i}: {k} clusters" for i, k in enumerate(selected, 1)]


def test_save_uses_save_config_defaults(tmp_path):
    pc = PhytClust(
        SAMPLE_TREE,
        runtime_config=RuntimeConfig(save=SaveConfig(tsv_name="custom.tsv", outlier=False)),
    )
    pc.run(k=5)
    pc.save(str(tmp_path))
    rows = (tmp_path / "custom.tsv").read_text().splitlines()[1:]
    assert all(not row.endswith("\t-1") for row in rows)
    assert not (tmp_path / "phytclust_results.tsv").exists()
