"""Start the PhytClust web GUI."""

from __future__ import annotations

import argparse
import importlib
import logging
import os
import sys
import threading
import webbrowser


logger = logging.getLogger("uvicorn.error")


def _port_number(value: str) -> int:
    """Check a TCP port supplied on the command line."""
    try:
        port = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "port must be an integer from 1 to 65535."
        ) from error
    if not 1 <= port <= 65535:
        raise argparse.ArgumentTypeError("port must be an integer from 1 to 65535.")
    return port


def build_gui_parser() -> argparse.ArgumentParser:
    """Define the GUI startup options."""
    parser = argparse.ArgumentParser(
        prog="phytclust gui",
        description="Start the PhytClust web GUI.",
    )
    parser.add_argument(
        "--host", default="127.0.0.1", help="Address to listen on (default: 127.0.0.1)."
    )
    parser.add_argument(
        "--port",
        type=_port_number,
        default=8000,
        help="Port to listen on (default: 8000).",
    )
    parser.add_argument(
        "--reload",
        action="store_true",
        help="Restart the server when code changes. Browser opening is disabled.",
    )
    parser.add_argument(
        "--public",
        action="store_true",
        help="Limit tree size and disable saving files on the server.",
    )
    parser.add_argument(
        "--no-browser", action="store_true", help="Do not open a browser automatically."
    )
    return parser


def _browser_url(host: str, port: int) -> str:
    """Use a local address for wildcard listeners and bracket IPv6 addresses."""
    host = host.strip("[]")
    host = {"0.0.0.0": "127.0.0.1", "::": "::1"}.get(host, host)
    if ":" in host:
        host = f"[{host}]"
    return f"http://{host}:{port}"


def _open_browser_when_ready(
    server, url: str, stopped: threading.Event
) -> threading.Thread:
    """Open the browser after startup unless the server has stopped."""

    def open_browser():
        while not stopped.is_set():
            if server.should_exit:
                return
            if server.started:
                try:
                    if not webbrowser.open(url):
                        logger.warning(
                            "Could not open a browser. Open %s manually.", url
                        )
                except (webbrowser.Error, OSError) as error:
                    logger.warning(
                        "Could not open a browser: %s. Open %s manually.", error, url
                    )
                return
            stopped.wait(0.05)

    thread = threading.Thread(
        target=open_browser, name="phytclust-browser", daemon=True
    )
    thread.start()
    return thread


def main(argv=None) -> int:
    """Start the GUI and return its exit status."""
    args = build_gui_parser().parse_args(argv)
    try:
        uvicorn = importlib.import_module("uvicorn")
        for module_name in ("fastapi", "jinja2"):
            importlib.import_module(module_name)
    except ImportError as error:
        print(
            f"ERROR: GUI dependencies are unavailable: {error}\nInstall them with: pip install 'phytclust[gui]'",
            file=sys.stderr,
        )
        return 1

    previous_public_mode = os.environ.get("PHYTCLUST_PUBLIC_MODE")
    if args.public:
        os.environ["PHYTCLUST_PUBLIC_MODE"] = "1"
    stopped = threading.Event()
    browser_thread = None
    url = _browser_url(args.host, args.port)
    print(f"Starting PhytClust GUI at {url} (Ctrl-C to stop)", flush=True)
    try:
        if args.reload:
            uvicorn.run(
                "phytclust.gui.api:app",
                host=args.host.strip("[]"),
                port=args.port,
                reload=True,
            )
        else:
            config = uvicorn.Config(
                "phytclust.gui.api:app", host=args.host.strip("[]"), port=args.port
            )
            server = uvicorn.Server(config)
            if not args.no_browser:
                browser_thread = _open_browser_when_ready(server, url, stopped)
            server.run()
            if not server.started:
                print("ERROR: The GUI server did not start.", file=sys.stderr)
                return 1
    except KeyboardInterrupt:
        return 0
    except SystemExit as error:
        return error.code if isinstance(error.code, int) else 1
    except OSError as error:
        print(f"ERROR: Could not start the GUI: {error}", file=sys.stderr)
        return 1
    finally:
        stopped.set()
        if browser_thread is not None:
            browser_thread.join(timeout=0.2)
        if args.public:
            if previous_public_mode is None:
                os.environ.pop("PHYTCLUST_PUBLIC_MODE", None)
            else:
                os.environ["PHYTCLUST_PUBLIC_MODE"] = previous_public_mode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
