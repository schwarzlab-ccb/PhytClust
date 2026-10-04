"""Regression checks for outgroups, all-k output, and local GUI exports."""

from collections import OrderedDict
from io import StringIO

import pytest
from fastapi import HTTPException

from phytclust import PhytClust
from phytclust.exceptions import InvalidTreeError
from phytclust.gui import api


def test_internal_clade_outgroup_counts_all_removed_leaves():
    pc = PhytClust("((a:1,b:1)X:1,(c:1,d:1):1);", outgroup="X")

    assert pc.num_terminals == 2
    assert pc.max_k == 2
    assert {leaf.name for leaf in pc.get_clusters(2)} == {"c", "d"}


def test_entire_tree_cannot_be_pruned_as_outgroup():
    with pytest.raises(InvalidTreeError, match="entire tree"):
        PhytClust("((a:1,b:1):1,c:1)ROOT;", outgroup="ROOT")


def test_save_all_skips_infeasible_k_but_keeps_feasible_ones(tmp_path):
    pc = PhytClust("(a:1,b:1,c:1,d:1);", polytomy_mode="hard", max_k=4)

    path = pc.save(str(tmp_path), output_all=True)

    assert path is not None
    header = (tmp_path / "phytclust_results.tsv").read_text().splitlines()[0]
    assert header == "Node Name\tclusters_k1\tclusters_k4"


def test_historical_export_uses_requested_run_snapshot(monkeypatch):
    monkeypatch.setattr(api, "_CACHE", OrderedDict())
    monkeypatch.setattr(api, "LAST_PC", None)
    monkeypatch.setattr(api, "LAST_RESULT", None)
    monkeypatch.setattr(api, "LAST_CONSTRUCTION_KEY", None)

    newick = "((a:1,b:1):1,(c:1,d:1):1);"
    old = api._run_phytclust(api.PhytclustRequest(newick=newick, mode="k", k=2))
    api._run_phytclust(api.PhytclustRequest(newick=newick, mode="k", k=3))

    response = api.export_tsv(
        api.ExportTSVRequest(run_id=old["run_id"], outlier=False)
    )
    assert StringIO(response.body.decode()).readline().strip() == "Node Name\tclusters_k2"
    with pytest.raises(HTTPException) as exc:
        api.export_tsv(api.ExportTSVRequest(run_id="missing"))
    assert exc.value.status_code == 404


def test_local_gui_save_writes_displayed_tsv(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(api, "LAST_PC", object())
    monkeypatch.setattr(api, "LAST_RESULT", {"run_id": "current"})
    displayed = "Node Name\tclusters_k3\na\t0\nb\t-1\n"

    api.save_results(
        api.SaveRequest(
            results_dir="results",
            filename="phytclust_k3.tsv",
            tsv=displayed,
        )
    )

    assert (tmp_path / "results" / "phytclust_k3.tsv").read_text() == displayed
