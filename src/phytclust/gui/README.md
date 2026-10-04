# PhytClust Web Frontend

A FastAPI + D3.js web GUI for PhytClust. Paste or upload a Newick tree to
visualise it, run clustering, and explore/compare/export the results
interactively.

## Installation

Install PhytClust with the optional GUI dependencies:

```bash
pip install -e '.[gui]'
```

(equivalently: `pip install 'uvicorn[standard]' fastapi jinja2 python-multipart`)

## Running

The simplest way is the bundled subcommand, which starts the server and opens a
browser:

```bash
phytclust gui                 # http://127.0.0.1:8000
phytclust gui --port 8080     # custom port
phytclust gui --public        # public-safe mode (tip cap, no server-side save)
phytclust gui --reload        # auto-reload for development
```

Ports must be between 1 and 65535. The browser opens after server startup;
`--no-browser` disables it. With `--reload`, open the URL manually. Wildcard
listeners (`--host 0.0.0.0` or `--host ::`) use a loopback URL for the local
browser.

See `phytclust gui --help` for all options.

Or run uvicorn directly (useful for custom deployment):

```bash
uvicorn phytclust.gui.api:app --host 127.0.0.1 --port 8000 --reload
```

Then open `http://127.0.0.1:8000` in a browser.

## Public-safe mode

Set `PHYTCLUST_PUBLIC_MODE=1` (or pass `--public`) when exposing the server
beyond localhost. It caps the input tip count (`PHYTCLUST_MAX_TIPS`, default
10000) and disables the server-side `/api/save` endpoint; use **Export TSV** in
the UI instead.

Run requests require positive integer counts and one of the modes `k`, `global`, or `resolution`. Exact-k requests ignore peak settings. Nonfinite scores and alpha values are returned as JSON null.

TSV exports use stored assignments. An optional `run_id` selects an earlier run; omitting it selects the latest result. Server-side saves also accept `run_id` for a stored TSV, or `tsv` for the table currently displayed in the GUI. Saves that use the live clustering instance wait for any active run to finish.
