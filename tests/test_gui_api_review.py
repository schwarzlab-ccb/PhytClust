"""GUI request errors, stored exports, and JSON values."""

from collections import OrderedDict
import threading
import asyncio
import json
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
    monkeypatch.setattr(api, "_CACHE", OrderedDict())
    monkeypatch.setattr(api, "PUBLIC_MODE", False)
    monkeypatch.setattr(api, "_RUN_LOCK", threading.Lock())


@pytest.fixture
def client():
    def request(method, path, json_body=None):
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
                "headers": [(b"content-type", b"application/json")],
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
        post=lambda path, json: request("POST", path, json),
        get=lambda path: request("GET", path),
    )


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
    original = api._historical_tsv

    def serialize(result, request):
        assert not api._RUN_LOCK.locked()
        return original(result, request)

    monkeypatch.setattr(api, "_historical_tsv", serialize)
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


def test_public_mode_limit_and_export(client, monkeypatch):
    monkeypatch.setattr(api, "PUBLIC_MODE", True)
    monkeypatch.setattr(api, "PUBLIC_MAX_TIPS", 3)
    assert (
        client.post("/api/run", json={"newick": TREE, "mode": "k", "k": 2}).status_code
        == 400
    )
    assert client.post("/api/save", json={"results_dir": "results"}).status_code == 403


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
