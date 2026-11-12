"""GUI request errors, stored exports, and JSON values."""

from collections import OrderedDict
import threading
import asyncio
import json
from pathlib import Path
import re
from types import SimpleNamespace

import numpy as np
import pytest
from fastapi import HTTPException

from phytclust.gui import api

TREE = "((a:1,b:1):1,(c:1,d:1):1);"


@pytest.fixture(autouse=True)
def isolated_runs(monkeypatch):
    for name in ("LAST_PC", "LAST_RESULT", "LAST_NEWICK", "LAST_CONSTRUCTION_KEY"):
        monkeypatch.setattr(api, name, None)
    monkeypatch.setattr(api, "LAST_CONSTRUCTION_NOTES", [])
    monkeypatch.setattr(api, "_CACHE", OrderedDict())
    monkeypatch.setattr(api, "_RUN_LOCK", threading.Lock())


@pytest.fixture
def client():
    def request(method, path, json_body=None, host="127.0.0.1"):
        async def send_request():
            body = json.dumps(json_body).encode() if json_body is not None else b""
            received = False
            messages = []

            async def receive():
                nonlocal received
                if not received:
                    received = True
                    return {"type": "http.request", "body": body, "more_body": False}
                await asyncio.Event().wait()

            async def send(message):
                messages.append(message)

            scope = {
                "type": "http",
                "asgi": {"version": "3.0"},
                "http_version": "1.1",
                "method": method,
                "scheme": "http",
                "path": path,
                "raw_path": path.encode(),
                "root_path": "",
                "query_string": b"",
                "headers": [
                    (b"content-type", b"application/json"),
                    (b"host", host.encode()),
                ],
                "server": ("localhost", 80),
                "client": ("127.0.0.1", 1234),
            }

            async def keep_loop_active():
                while True:
                    await asyncio.sleep(0.01)

            heartbeat = asyncio.create_task(keep_loop_active())
            try:
                await asyncio.wait_for(api.app(scope, receive, send), timeout=5)
            finally:
                heartbeat.cancel()
                await asyncio.gather(heartbeat, return_exceptions=True)
            start = next(
                message
                for message in messages
                if message["type"] == "http.response.start"
            )
            content = b"".join(
                message.get("body", b"")
                for message in messages
                if message["type"] == "http.response.body"
            )
            return SimpleNamespace(
                status_code=start["status"],
                content=content,
                headers={
                    name.decode(): value.decode() for name, value in start["headers"]
                },
                json=lambda: json.loads(content),
            )

        return asyncio.run(send_request())

    return SimpleNamespace(
        post=lambda path, json, host="127.0.0.1": request("POST", path, json, host),
        get=lambda path, host="127.0.0.1": request("GET", path, host=host),
    )


def test_non_loopback_host_header_is_rejected(client):
    """A DNS-rebinding page sends its own hostname; the API must refuse it."""
    response = client.get("/", host="evil.example.com")
    assert response.status_code == 400
    assert client.get("/", host="localhost:8000").status_code == 200


def test_save_filename_must_be_tabular_text(client, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    run = client.post("/api/run", json={"newick": "((a:1,b:1):1,(c:1,d:1):1);", "mode": "k", "k": 2})
    assert run.status_code == 200
    for filename in ("sitecustomize.py", ".bashrc", "notes.html"):
        response = client.post(
            "/api/save",
            json={"results_dir": "out", "filename": filename, "tsv": "a\tb\n"},
        )
        assert response.status_code == 400, filename
        assert "must end in" in response.json()["detail"]
    assert not (tmp_path / "out").exists() or not list((tmp_path / "out").iterdir())
    ok = client.post(
        "/api/save", json={"results_dir": "out", "filename": "ok.tsv", "tsv": "a\tb\n"}
    )
    assert ok.status_code == 200
    assert (tmp_path / "out" / "ok.tsv").read_text() == "a\tb\n"


@pytest.mark.parametrize(
    "settings",
    [
        {"top_n": 0},
        {"k": True},
        {"max_k": 1.5},
        {"num_bins": -1},
        {"mode": "other"},
        {"polytomy_mode": "other"},
    ],
)
def test_invalid_request_fields(client, settings):
    response = client.post("/api/run", json={"newick": TREE, **settings})
    assert response.status_code == 422
    assert "detail" in response.json()


def test_bad_newick_reports_parse_error(client):
    response = client.post("/api/run", json={"newick": "((a:1,b:1);"})
    assert response.status_code == 400
    assert "Could not parse Newick" in response.json()["detail"]


def test_new_runs_parse_once(monkeypatch):
    original = api.Phylo.read
    parsed = []

    def read(*args, **kwargs):
        result = original(*args, **kwargs)
        parsed.append(result)
        return result

    monkeypatch.setattr(api.Phylo, "read", read)
    result = api.run_phytclust(api.PhytclustRequest(newick=TREE, mode="k", k=2))
    assert result["k_values"] == [2]
    assert len(parsed) == 1


def test_exact_k_ignores_unused_peak_settings(client):
    response = client.post(
        "/api/run",
        json={
            "newick": TREE,
            "mode": "k",
            "k": 2,
            "prominence_weight": -1,
            "boundary_window_size": 0,
        },
    )
    assert response.status_code == 200
    assert response.json()["k_values"] == [2]


def test_active_peak_errors_reported(client):
    response = client.post("/api/run", json={"newick": TREE, "prominence_weight": -1})
    assert response.status_code == 400
    assert "prominence_weight" in response.json()["detail"]


def test_nonfinite_scores_are_json_null(client, monkeypatch):
    original = api.PhytClust.run

    def run(self, **kwargs):
        result = original(self, **kwargs)
        result["scores"] = np.array([np.nan, np.inf, -np.inf, 2.0])
        result["alpha_details"] = [
            {"score": np.float64(np.nan), "count": np.int64(3), "label": "example"}
        ]
        return result

    monkeypatch.setattr(api.PhytClust, "run", run)
    response = client.post("/api/run", json={"newick": TREE, "mode": "k", "k": 2})
    assert response.status_code == 200
    assert response.json()["scores"] == [None, None, None, 2.0]
    assert response.json()["alpha_details"][0] == {
        "score": None,
        "count": 3.0,
        "label": "example",
    }


def test_latest_export_uses_snapshot_without_live_save(monkeypatch):
    api.run_phytclust(api.PhytclustRequest(newick=TREE, mode="k", k=2))
    monkeypatch.setattr(
        api.LAST_PC, "save", lambda **kwargs: pytest.fail("export must use snapshot")
    )
    response = api.export_tsv(api.ExportTSVRequest(outlier=False))
    assert response.body.decode().splitlines()[0] == "Node Name\tclusters_k2"


def test_save_requested_run_uses_snapshot(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    old = api.run_phytclust(api.PhytclustRequest(newick=TREE, mode="k", k=2))
    api.run_phytclust(api.PhytclustRequest(newick=TREE, mode="k", k=3))
    api.save_results(api.SaveRequest(results_dir="results", run_id=old["run_id"]))
    assert (tmp_path / "results/phytclust_results.tsv").read_text().splitlines()[
        0
    ] == "Node Name\tclusters_k2"


def test_save_waits_for_run_lock(monkeypatch):
    saved = threading.Event()
    monkeypatch.setattr(api, "_save_results", lambda request: saved.set())
    with api._RUN_LOCK:
        thread = threading.Thread(
            target=api.save_results,
            args=(api.SaveRequest(results_dir="results"),),
            daemon=True,
        )
        thread.start()
        assert not saved.wait(0.1)
    assert saved.wait(1)
    thread.join(1)
    assert not thread.is_alive()


def test_export_serializes_snapshot_after_releasing_lock(monkeypatch):
    api.run_phytclust(api.PhytclustRequest(newick=TREE, mode="k", k=2))
    original = api._run_tsv

    def serialize(result, request):
        assert not api._RUN_LOCK.locked()
        return original(result, request)

    monkeypatch.setattr(api, "_run_tsv", serialize)
    assert api.export_tsv(api.ExportTSVRequest()).status_code == 200


def test_invalid_filename_does_not_create_directory(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    api.run_phytclust(api.PhytclustRequest(newick=TREE, mode="k", k=2))
    with pytest.raises(HTTPException) as error:
        api.save_results(api.SaveRequest(results_dir="new", filename="../outside.tsv"))
    assert error.value.status_code == 400
    assert not (tmp_path / "new").exists()


def test_favicon_uses_existing_logo(client):
    response = client.get("/favicon.ico")
    assert response.status_code == 200
    assert response.headers["content-type"] == "image/png"
    assert response.content.startswith(b"\x89PNG")


@pytest.mark.parametrize(
    "newick, names",
    [
        ('("A,B":1,c:1);', {"A,B", "c"}),
        ("('O''Brien':1,'A,(B): C':1);", {"O'Brien", "A,(B): C"}),
        ("[&R](a:1[notes, here],b:1);", {"a", "b"}),
        ("(Taxon one:1,Taxon two:1);", {"Taxon one", "Taxon two"}),
    ],
)
def test_browser_and_server_keep_quoted_names_and_comments(client, newick, names):
    response = client.post("/api/run", json={"newick": newick, "mode": "k", "k": 1})
    assert response.status_code == 200
    assert set(response.json()["clusters"][0]) == names


def test_unclosed_quoted_name_reports_request_error(client):
    response = client.post("/api/run", json={"newick": "('a,b);", "mode": "k", "k": 1})
    assert response.status_code == 400
    assert "not closed" in response.json()["detail"]


REPAIRED_TREE = (
    "((((A:1,A:1):1,(C:1,D:1):1):5,((E:1,F:1):1,(G:1,H:1):1):5):3,"
    "(((I:1,J:1):1,(K:1,:1):1):5,((M:1,N:-1):1,(O:1,P:1)95:1):5):3);"
)


def _drawn_leaves(newick):
    """Leaf names of the returned tree: labels that follow '(' or ','."""
    names = []
    for match in re.finditer(r"(?:^|[(,])('(?:[^']|'')*'|[^(),:;']+)(?=[:,)])", newick):
        names.append(match.group(1).strip("'").replace("''", "'"))
    return names


def test_input_repairs_are_returned_as_notes(client):
    response = client.post("/api/run", json={"newick": REPAIRED_TREE, "mode": "k", "k": 3})
    assert response.status_code == 200
    notes = response.json()["notes"]
    assert any("renaming to 'A_1'" in note for note in notes)
    assert any("negative length" in note for note in notes)
    # A reused clustering instance still reports how its tree was repaired.
    again = client.post("/api/run", json={"newick": REPAIRED_TREE, "mode": "k", "k": 4})
    assert again.json()["notes"] == notes


def test_clean_tree_has_no_notes(client):
    response = client.post("/api/run", json={"newick": TREE, "mode": "k", "k": 2})
    assert response.json()["notes"] == []


@pytest.mark.parametrize(
    "newick",
    [REPAIRED_TREE, "((A:1,B:1)A:1,(C:1,D:1):1);", "('O''Brien':1,'A,(B): C':1);"],
)
def test_every_drawn_leaf_has_a_cluster(client, newick):
    """The returned tree is the clustered one, so leaf names are cluster keys."""
    response = client.post("/api/run", json={"newick": newick, "mode": "k", "k": 1})
    assert response.status_code == 200
    data = response.json()
    leaves = _drawn_leaves(data["newick"])
    assert sorted(leaves) == sorted(data["clusters"][0])
    assert len(set(leaves)) == len(leaves)


def test_returned_tree_keeps_support_labels_and_repaired_lengths(client):
    data = client.post("/api/run", json={"newick": REPAIRED_TREE, "mode": "k", "k": 1}).json()
    assert ")95:1" in data["newick"]
    assert "N:0)" in data["newick"]
    assert "internal_node_" not in data["newick"].replace("internal_node_P0", "")


@pytest.mark.parametrize("length", ["nan", "inf", "-inf", "abc", "1e999", "0x10"])
def test_nonfinite_or_nonnumeric_branch_lengths_are_rejected(client, length):
    response = client.post(
        "/api/run", json={"newick": f"((A:{length},B:1):1,(C:1,D:1):1);", "mode": "k", "k": 2}
    )
    assert response.status_code == 400
    assert f"Branch length '{length}' must be a finite number" in response.json()["detail"]


def test_numeric_branch_length_forms_are_accepted(client):
    response = client.post(
        "/api/run",
        json={"newick": "((A: 1 ,B:1e-3):.5,(C:+1.,D:1E2)[c]:1);", "mode": "k", "k": 2},
    )
    assert response.status_code == 200


def test_resolution_fallback_reports_global_mode(client):
    response = client.post("/api/run", json={"newick": REPAIRED_TREE, "mode": "resolution"})
    assert response.status_code == 200
    data = response.json()
    assert data["requested_mode"] == "resolution"
    if data["mode"] == "global":
        assert data["bin_ranges"] is None
        assert any("resolution" in note.lower() for note in data["notes"])


def test_save_expands_home_inside_working_directory(client, tmp_path, monkeypatch):
    home = tmp_path / "home"
    work = home / "project"
    work.mkdir(parents=True)
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.chdir(work)
    assert client.post("/api/run", json={"newick": TREE, "mode": "k", "k": 2}).status_code == 200
    response = client.post(
        "/api/save",
        json={"results_dir": "~/project/gui_results", "filename": "ok.tsv", "tsv": "a\tb\n"},
    )
    assert response.status_code == 200
    assert response.json()["results_dir"] == str(work / "gui_results")
    assert (work / "gui_results" / "ok.tsv").is_file()
    assert not (work / "~").exists()

    outside = client.post(
        "/api/save", json={"results_dir": "~/gui_results", "filename": "ok.tsv", "tsv": "x"}
    )
    assert outside.status_code == 403
    assert "outside the server's working directory" in outside.json()["detail"]
    assert not (home / "gui_results").exists()
    assert not (work / "~").exists()


GUI_DIR = Path(api.__file__).parent


def test_help_example_clusters_with_default_settings(client):
    help_html = (GUI_DIR / "templates/partials/_help_sidebar.html").read_text()
    sidebar_html = (GUI_DIR / "templates/partials/_sidebar.html").read_text()
    example = re.search(r"Example: <code>([^<]+)</code>", help_html).group(1)
    assert f'placeholder="{example}"' in sidebar_html
    default_request = {
        "newick": example,
        "mode": "global",
        "ranking_mode": "adjusted",
        "use_relative_prominence": False,
        "exclude_k2": False,
        "polytomy_mode": "soft",
    }
    response = client.post("/api/run", json=default_request)
    assert response.status_code == 200
    assert response.json()["clusters"], "the help example must produce a clustering"


def test_dark_theme_branches_are_light_enough_to_see():
    dark = (GUI_DIR / "static/css/base.css").read_text().split('[data-theme="dark"]', 1)[1]
    branch = re.search(r"--pc-tree-branch:\s*(#[0-9a-fA-F]{6})", dark).group(1)
    assert int(branch[1:3], 16) > 0x80
