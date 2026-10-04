"""GUI startup validation, browser timing, and environment cleanup."""

import os
import threading
from types import SimpleNamespace

import pytest

from phytclust.gui import launch


@pytest.mark.parametrize("port", ["0", "-1", "65536", "abc", "1.5"])
def test_invalid_ports(port):
    with pytest.raises(SystemExit) as error:
        launch.build_gui_parser().parse_args(["--port", port])
    assert error.value.code == 2


@pytest.mark.parametrize(
    "host, expected",
    [
        ("0.0.0.0", "127.0.0.1"),
        ("::", "[::1]"),
        ("::1", "[::1]"),
        ("[::1]", "[::1]"),
        ("localhost", "localhost"),
    ],
)
def test_browser_addresses(host, expected):
    assert launch._browser_url(host, 8000) == f"http://{expected}:8000"


@pytest.fixture
def fake_uvicorn(monkeypatch):
    state = SimpleNamespace(
        started=False, should_exit=False, configs=[], reload_calls=[]
    )

    def config(app, **options):
        state.configs.append((app, options))
        return SimpleNamespace(**options)

    state.run = lambda: setattr(state, "started", True)
    fake = SimpleNamespace(
        Config=config,
        Server=lambda config: state,
        run=lambda *args, **kwargs: state.reload_calls.append((args, kwargs)),
    )
    monkeypatch.setattr(
        launch.importlib,
        "import_module",
        lambda name: fake if name == "uvicorn" else SimpleNamespace(),
    )
    return state


def test_no_browser_and_ipv6_host(fake_uvicorn, monkeypatch):
    monkeypatch.setattr(
        launch.webbrowser, "open", lambda url: pytest.fail("browser must not open")
    )
    assert launch.main(["--host", "[::1]", "--no-browser"]) == 0
    assert fake_uvicorn.configs == [
        ("phytclust.gui.api:app", {"host": "::1", "port": 8000})
    ]


def test_reload_does_not_open_browser(fake_uvicorn, monkeypatch):
    monkeypatch.setattr(
        launch.webbrowser, "open", lambda url: pytest.fail("browser must not open")
    )
    assert launch.main(["--reload"]) == 0
    assert fake_uvicorn.reload_calls == [
        (
            ("phytclust.gui.api:app",),
            {"host": "127.0.0.1", "port": 8000, "reload": True},
        )
    ]


@pytest.mark.parametrize("previous", [None, "0", "1"])
@pytest.mark.parametrize("failure", [False, True])
def test_public_mode_restored_after_shutdown(
    fake_uvicorn, monkeypatch, previous, failure
):
    if previous is None:
        monkeypatch.delenv("PHYTCLUST_PUBLIC_MODE", raising=False)
    else:
        monkeypatch.setenv("PHYTCLUST_PUBLIC_MODE", previous)

    def run():
        assert os.environ["PHYTCLUST_PUBLIC_MODE"] == "1"
        if failure:
            raise OSError("port occupied")
        fake_uvicorn.started = True

    fake_uvicorn.run = run
    assert launch.main(["--public", "--no-browser"]) == (1 if failure else 0)
    assert os.environ.get("PHYTCLUST_PUBLIC_MODE") == previous


def test_existing_public_mode_kept_without_flag(fake_uvicorn, monkeypatch):
    monkeypatch.setenv("PHYTCLUST_PUBLIC_MODE", "1")
    assert launch.main(["--no-browser"]) == 0
    assert os.environ["PHYTCLUST_PUBLIC_MODE"] == "1"


def test_dependency_error_uses_stderr(monkeypatch, capsys):
    def missing(name):
        raise ImportError("missing dependency")

    monkeypatch.setattr(launch.importlib, "import_module", missing)
    assert launch.main([]) == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "pip install 'phytclust[gui]'" in captured.err


def test_startup_failure_exit_status(fake_uvicorn):
    fake_uvicorn.run = lambda: None
    assert launch.main(["--no-browser"]) == 1

    def fail():
        raise SystemExit(3)

    fake_uvicorn.run = fail
    assert launch.main(["--no-browser"]) == 3


def test_browser_waits_for_started_server(monkeypatch):
    server = SimpleNamespace(started=False, should_exit=False)
    opened = threading.Event()
    stopped = threading.Event()
    monkeypatch.setattr(launch.webbrowser, "open", lambda url: opened.set() or True)
    thread = launch._open_browser_when_ready(server, "http://localhost:8000", stopped)
    try:
        assert thread.daemon
        assert not opened.wait(0.1)
        server.started = True
        assert opened.wait(1)
    finally:
        stopped.set()
        thread.join(1)
    assert not thread.is_alive()


def test_browser_cancelled_before_startup(monkeypatch):
    stopped = threading.Event()
    monkeypatch.setattr(
        launch.webbrowser, "open", lambda url: pytest.fail("browser must not open")
    )
    thread = launch._open_browser_when_ready(
        SimpleNamespace(started=False, should_exit=False),
        "http://localhost:8000",
        stopped,
    )
    stopped.set()
    thread.join(1)
    assert not thread.is_alive()


def test_browser_failure_warns(monkeypatch, caplog):
    monkeypatch.setattr(launch.webbrowser, "open", lambda url: False)
    thread = launch._open_browser_when_ready(
        SimpleNamespace(started=True, should_exit=False),
        "http://localhost:8000",
        threading.Event(),
    )
    thread.join(1)
    assert "Open http://localhost:8000 manually" in caplog.text
