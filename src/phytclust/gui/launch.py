"""Launch the PhytClust web GUI (a thin uvicorn wrapper).

Exposed as the ``phytclust gui`` subcommand; see ``phytclust gui --help``.
"""
from __future__ import annotations

import argparse
import os


def build_gui_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="phytclust gui",
        description="Launch the PhytClust web GUI (FastAPI + uvicorn).",
    )
    p.add_argument(
        "--host", default="127.0.0.1", help="Host to bind (default: 127.0.0.1)."
    )
    p.add_argument(
        "--port", type=int, default=8000, help="Port to bind (default: 8000)."
    )
    p.add_argument(
        "--reload", action="store_true", help="Auto-reload on code changes (dev)."
    )
    p.add_argument(
        "--public",
        action="store_true",
        help="Enable public-safe mode: cap tip count and disable server-side save.",
    )
    p.add_argument(
        "--no-browser",
        action="store_true",
        help="Do not open a browser window automatically.",
    )
    return p


def _open_browser_soon(url: str, delay: float = 1.5) -> None:
    """Best-effort: open the default browser a moment after the server binds."""
    import threading
    import webbrowser

    def _open() -> None:
        try:
            webbrowser.open(url)
        except Exception:
            pass

    threading.Timer(delay, _open).start()


def main(argv=None) -> int:
    args = build_gui_parser().parse_args(argv)

    try:
        import uvicorn
    except ImportError:
        print(
            "The GUI requires the optional 'gui' dependencies. Install with:\n"
            "    pip install 'phytclust[gui]'\n"
            "(or: pip install 'uvicorn[standard]' fastapi jinja2 python-multipart)"
        )
        return 1

    if args.public:
        os.environ["PHYTCLUST_PUBLIC_MODE"] = "1"

    url = f"http://{args.host}:{args.port}"
    # Under --reload uvicorn spawns a reloader subprocess, so this parent would
    # open the browser before the server is actually ready — skip it there.
    if not args.no_browser and not args.reload:
        _open_browser_soon(url)

    print(f"PhytClust GUI starting at {url}  (press Ctrl-C to stop)")
    uvicorn.run(
        "phytclust.gui.api:app",
        host=args.host,
        port=args.port,
        reload=args.reload,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
