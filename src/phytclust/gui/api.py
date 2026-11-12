"""Run clustering and serve GUI results and exports."""

import csv
import logging
import math
import os
import re
import threading
import uuid
from collections import OrderedDict
from contextlib import contextmanager
from io import StringIO
from pathlib import Path
from typing import Annotated, Any, Literal, Optional
from numbers import Real

from Bio import Phylo
from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import FileResponse, HTMLResponse, PlainTextResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from pydantic import BaseModel, Field, field_validator
from starlette.middleware.trustedhost import TrustedHostMiddleware

from phytclust import __version__
from phytclust.algo.core import PhytClust
from phytclust.config import OutlierConfig, PeakConfig
from phytclust.io.save import _output_filename
from phytclust.exceptions import (
    PhytClustError,
)
from phytclust.utils.traversal import match_child_order, tree_fingerprint

logger = logging.getLogger("phytclust.gui")

# Same grammar as DECIMAL_NUMBER in static/js/utils.js, so the server accepts
# exactly the branch lengths the browser can draw.
_DECIMAL_NUMBER = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?")


def _is_finite_decimal(text: str) -> bool:
    """Return True for a plain decimal number that is finite as a float."""
    return bool(_DECIMAL_NUMBER.fullmatch(text)) and math.isfinite(float(text))


def _normalize_newick(newick: str) -> str:
    """Use single quotes for quoted names and names containing spaces."""
    result: list[str] = []
    position = 0
    while position < len(newick):
        character = newick[position]
        if character in "(),;":
            result.append(character)
            position += 1
        elif character == "[":
            start = position
            depth = 1
            position += 1
            while position < len(newick) and depth:
                if newick[position] == "[":
                    depth += 1
                elif newick[position] == "]":
                    depth -= 1
                position += 1
            if depth:
                raise ValueError("The Newick comment is not closed.")
            result.append(newick[start:position])
        elif character in "'\"":
            quote = character
            label: list[str] = []
            position += 1
            while position < len(newick):
                character = newick[position]
                position += 1
                if character != quote:
                    label.append(character)
                elif position < len(newick) and newick[position] == quote:
                    label.append(quote)
                    position += 1
                else:
                    break
            else:
                raise ValueError("The quoted Newick name is not closed.")
            result.append("'" + "".join(label).replace("'", "''") + "'")
        elif character == ":":
            start = position
            position += 1
            while position < len(newick) and newick[position] not in "(),;[":
                position += 1
            length = newick[start + 1 : position].strip()
            if length and not _is_finite_decimal(length):
                raise ValueError(
                    f"Branch length {length!r} must be a finite number."
                )
            result.append(newick[start:position])
        else:
            start = position
            while position < len(newick) and newick[position] not in "(),:;[]'\"":
                position += 1
            if start == position:
                raise ValueError("Unexpected character in Newick tree.")
            token = newick[start:position]
            label = token.strip()
            if label and any(character.isspace() for character in label):
                result.append("'" + label.replace("'", "''") + "'")
            else:
                result.append(token)
    return "".join(result)


_CACHE: OrderedDict[str, tuple[PhytClust, dict[str, Any]]] = OrderedDict()
_CACHE_MAX = 20

_RUN_LOCK = threading.Lock()

app = FastAPI(title="PhytClust Web API", version=__version__)

LOOPBACK_HOSTS = ["127.0.0.1", "localhost", "[::1]"]


def _allowed_hosts() -> list[str]:
    """Accept loopback Host headers unless PHYTCLUST_ALLOWED_HOSTS widens the list.

    Other Host values are rejected to block DNS rebinding against a loopback
    server. The launcher sets the variable when binding a non-loopback address.
    """
    configured = os.environ.get("PHYTCLUST_ALLOWED_HOSTS", "").strip()
    if configured:
        return [host.strip() for host in configured.split(",") if host.strip()]
    return list(LOOPBACK_HOSTS)


app.add_middleware(TrustedHostMiddleware, allowed_hosts=_allowed_hosts())

SAVE_EXTENSIONS = {".tsv", ".txt", ".csv"}


def _export_filename(filename: str) -> str:
    """Require a plain filename with a tabular text extension."""
    name = _output_filename(filename)
    if os.path.splitext(name)[1].casefold() not in SAVE_EXTENSIONS:
        raise ValueError(
            f"Output filename {filename!r} must end in one of: "
            + ", ".join(sorted(SAVE_EXTENSIONS))
            + "."
        )
    return name

STATIC_DIR = Path(__file__).parent / "static"
TEMPLATES_DIR = Path(__file__).parent / "templates"
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")
templates = Jinja2Templates(directory=TEMPLATES_DIR)

LAST_PC: Optional[PhytClust] = None
LAST_RESULT: Optional[dict[str, Any]] = None
LAST_NEWICK: Optional[str] = None
LAST_CONSTRUCTION_KEY: Optional[tuple] = None
LAST_CONSTRUCTION_NOTES: list[str] = []


def _cache_put(run_id: str, pc: PhytClust, result: dict[str, Any]) -> None:
    """Store a run and discard the oldest entries above the cache limit."""
    _CACHE[run_id] = (pc, result)
    _CACHE.move_to_end(run_id)
    while len(_CACHE) > _CACHE_MAX:
        _CACHE.popitem(last=False)


def _cache_get(run_id: Optional[str]) -> tuple[PhytClust, dict[str, Any]]:
    """Retrieve a stored run or the latest result."""
    if run_id is not None:
        if run_id in _CACHE:
            return _CACHE[run_id]
        raise HTTPException(
            status_code=404, detail="Requested run is no longer available."
        )
    if LAST_PC is not None and LAST_RESULT is not None:
        return LAST_PC, LAST_RESULT
    raise HTTPException(
        status_code=400,
        detail="No PhytClust result available. Run PhytClust before exporting results.",
    )


@app.get("/", response_class=HTMLResponse)
def serve_frontend(request: Request):
    """Serve the GUI page."""
    return templates.TemplateResponse(
        request, "index.html", {"phytclust_version": __version__}
    )


@app.get("/favicon.ico", include_in_schema=False)
def serve_favicon():
    """Use the existing logo as the browser tab icon."""
    return FileResponse(
        STATIC_DIR / "phytclust_logo_colour.png",
        media_type="image/png",
        headers={"Cache-Control": "public, max-age=86400"},
    )


PositiveCount = Annotated[int, Field(strict=True, ge=1)]


class PhytclustRequest(BaseModel):
    """Tree and settings submitted by the GUI."""

    newick: str
    mode: Literal["k", "global", "resolution"] = "global"
    k: PositiveCount | None = None
    top_n: PositiveCount = 1
    by_resolution: bool = False
    num_bins: PositiveCount | None = None
    max_k: PositiveCount | None = None
    max_k_limit: Optional[float] = None
    outgroup: Optional[str] = None
    root_taxon: Optional[str] = None
    min_cluster_size: PositiveCount | None = None
    compute_all_clusters: bool = False

    use_branch_support: bool = False
    min_support: Optional[float] = None
    support_weight: Optional[float] = None

    outlier_size_threshold: PositiveCount | None = None
    outlier_prefer_fewer: bool = False
    outlier_ratio_weight: Optional[float] = None
    outlier_ratio_mode: Optional[str] = None

    no_split_zero_length: bool = False
    polytomy_mode: Literal["hard", "soft"] = "soft"

    prominence_weight: float = 0.7
    ranking_mode: str = "adjusted"
    min_prominence: Optional[float] = None
    use_relative_prominence: bool = False
    boundary_window_size: Optional[int] = None
    boundary_ratio_threshold: Optional[float] = None
    exclude_k2: bool = False

    @field_validator("mode", mode="before")
    @classmethod
    def normalize_mode(cls, value):
        """Accept selection mode names without case differences."""
        return value.lower() if isinstance(value, str) else value


def _leaf_key_to_name(leaf: Any) -> str:
    name = getattr(leaf, "name", None)
    if name is not None:
        return str(name)
    return str(leaf)


def _json_float(value) -> float | None:
    """Return a finite float, or None for missing or nonfinite values."""
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return number if math.isfinite(number) else None


def _serialize_clusters(raw_clusters: Any) -> list[dict[str, int]]:
    """Convert leaf objects to names in each cluster map."""
    if raw_clusters is None:
        return []
    cluster_maps = [raw_clusters] if isinstance(raw_clusters, dict) else raw_clusters
    return [
        {
            _leaf_key_to_name(leaf): int(cluster_id)
            for leaf, cluster_id in cluster_map.items()
        }
        for cluster_map in cluster_maps
    ]


class _WarningCollector(logging.Handler):
    """Keep the distinct warning messages the library logs on this thread."""

    def __init__(self) -> None:
        """Collect WARNING and higher records from the creating thread only."""
        super().__init__(level=logging.WARNING)
        self.thread = threading.get_ident()
        self.messages: list[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        """Store the formatted message unless it is a GUI record or a repeat."""
        if record.thread != self.thread or record.name.startswith("phytclust.gui"):
            return
        try:
            message = record.getMessage()
        except Exception:
            return
        if message not in self.messages:
            self.messages.append(message)


@contextmanager
def _library_warnings():
    """Yield a list that fills with warnings logged by PhytClust in the block.

    PhytClust reports input repairs (renamed nodes, negative or missing
    lengths) through its logger; they are returned as run notes.
    """
    collector = _WarningCollector()
    library_logger = logging.getLogger("phytclust")
    library_logger.addHandler(collector)
    try:
        yield collector.messages
    finally:
        library_logger.removeHandler(collector)


_NEWICK_QUOTED_CHARACTERS = set("()[]':;,\"")


def _newick_label(name: str) -> str:
    """Quote a node name when the browser's Newick reader would split it."""
    if name and not any(
        character in _NEWICK_QUOTED_CHARACTERS or character.isspace()
        for character in name
    ):
        return name
    return "'" + name.replace("'", "''") + "'"


def _newick_number(value: Any) -> Optional[str]:
    """Format a finite number for Newick, or return None."""
    number = _json_float(value)
    if number is None:
        return None
    if number.is_integer() and abs(number) < 1e15:
        return str(int(number))
    return repr(number)


def _clustered_newick(tree: Any, input_names: set[str]) -> str:
    """Write the tree PhytClust clustered, so drawn leaves match cluster keys.

    PhytClust roots the tree, collapses single-child nodes, zeroes negative
    lengths and renames duplicate or unnamed nodes before clustering, so the
    drawn tree must be the prepared one. Generated internal names are replaced
    by the node's support value.
    """

    def suffix(clade: Any) -> str:
        label = ""
        name = None if clade.name is None else str(clade.name)
        if not clade.clades:
            label = _newick_label(name or "")
        elif name and (name in input_names or not name.startswith("internal_node_")):
            label = _newick_label(name)
        else:
            label = _newick_number(clade.confidence) or ""
        length = _newick_number(clade.branch_length)
        return label + (":" + length if length is not None else "")

    parts: list[str] = []
    stack: list[Any] = [tree.root]
    while stack:
        item = stack.pop()
        if isinstance(item, str):
            parts.append(item)
            continue
        if not item.clades:
            parts.append(suffix(item))
            continue
        parts.append("(")
        stack.append(")" + suffix(item))
        for index, child in enumerate(reversed(item.clades)):
            if index:
                stack.append(",")
            stack.append(child)
    return "".join(parts) + ";"


def _construction_key(request: PhytclustRequest, tree_key: Any) -> tuple:
    """Identify the tree and settings needed to reuse a clustering instance."""
    return (
        tree_key,
        request.outgroup,
        request.root_taxon,
        request.min_cluster_size,
        request.use_branch_support,
        request.min_support,
        request.support_weight,
        request.outlier_size_threshold,
        request.outlier_prefer_fewer,
        request.outlier_ratio_weight,
        request.outlier_ratio_mode,
        request.no_split_zero_length,
        request.polytomy_mode,
    )


@app.post("/api/run")
def run_phytclust(request: PhytclustRequest):
    """Run clustering while protecting the shared instance and cache."""
    with _RUN_LOCK:
        return _run_phytclust(request)


def _run_phytclust(request: PhytclustRequest):
    """Parse the tree, select partitions, and store the response."""
    if not request.newick.strip():
        raise HTTPException(status_code=400, detail="Empty Newick string.")

    try:
        newick = _normalize_newick(request.newick)
        parsed = Phylo.read(StringIO(newick.strip()), "newick")
    except (ValueError, Phylo.NewickIO.NewickError) as error:
        raise HTTPException(
            status_code=400, detail=f"Could not parse Newick: {error}"
        ) from error
    tree_key = tree_fingerprint(parsed.root)
    input_names = {
        str(clade.name) for clade in parsed.find_clades() if clade.name is not None
    }

    mode = request.mode.lower()
    if mode not in {"k", "global", "resolution"}:
        raise HTTPException(
            status_code=400, detail="mode must be 'k', 'global', or 'resolution'."
        )

    current_key = _construction_key(request, tree_key)
    if (
        LAST_CONSTRUCTION_KEY is not None
        and LAST_PC is not None
        and LAST_CONSTRUCTION_KEY == current_key
    ):
        pc = LAST_PC
        if newick != LAST_NEWICK and parsed is not None:
            match_child_order(pc.tree.root, parsed.root)
        pc.compute_all_clusters = request.compute_all_clusters
        if mode != "k":
            pc.exact_k = None
        construction_notes = list(LAST_CONSTRUCTION_NOTES)
        logger.debug("Reusing the previous clustering instance.")
    else:
        try:
            outlier_kwargs: dict[str, Any] = {}
            if request.outlier_size_threshold is not None:
                outlier_kwargs["size_threshold"] = request.outlier_size_threshold
            if request.outlier_prefer_fewer:
                outlier_kwargs["prefer_fewer"] = True
            if request.outlier_ratio_weight is not None:
                outlier_kwargs["ratio_weight"] = request.outlier_ratio_weight
            if request.outlier_ratio_mode is not None:
                outlier_kwargs["ratio_mode"] = request.outlier_ratio_mode

            kwargs: dict[str, Any] = dict(
                tree=parsed,
                outgroup=request.outgroup,
                compute_all_clusters=request.compute_all_clusters,
                use_branch_support=request.use_branch_support,
                no_split_zero_length=request.no_split_zero_length,
                polytomy_mode=request.polytomy_mode,
            )
            if request.root_taxon:
                kwargs["root_taxon"] = request.root_taxon
            if request.min_cluster_size is not None:
                kwargs["min_cluster_size"] = request.min_cluster_size
            if request.min_support is not None:
                kwargs["min_support"] = request.min_support
            if request.support_weight is not None:
                kwargs["support_weight"] = request.support_weight
            if outlier_kwargs:
                kwargs["outlier"] = OutlierConfig(**outlier_kwargs)

            with _library_warnings() as construction_notes:
                pc = PhytClust(**kwargs)
        except PhytClustError as e:
            raise HTTPException(status_code=400, detail=str(e))
        except (ValueError, TypeError) as error:
            raise HTTPException(status_code=400, detail=str(error)) from error
        except Exception:
            logger.exception("Could not initialize clustering")
            raise HTTPException(
                status_code=500, detail="Could not initialize clustering."
            )

    try:
        with _library_warnings() as run_notes:
            result = _run_selection(pc, request, mode)
    except HTTPException:
        raise
    except PhytClustError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception:
        logger.exception("PhytClust run failed")
        raise HTTPException(status_code=500, detail="PhytClust internal error.")
    notes = list(construction_notes)
    notes.extend(note for note in run_notes if note not in notes)

    return _store_run(
        pc,
        result,
        request,
        notes,
        newick=newick,
        clustered_newick=_clustered_newick(pc.tree, input_names),
        construction_key=current_key,
        construction_notes=construction_notes,
    )


def _run_selection(pc: PhytClust, request: PhytclustRequest, mode: str) -> dict:
    """Run the selection mode requested by the GUI and return the result."""
    peak_settings = None
    if mode != "k":
        peak_kwargs: dict[str, Any] = {"prominence_weight": request.prominence_weight}
        if request.ranking_mode:
            peak_kwargs["ranking_mode"] = request.ranking_mode
        if request.min_prominence is not None:
            peak_kwargs["min_prominence"] = request.min_prominence
        if request.use_relative_prominence:
            peak_kwargs["use_relative_prominence"] = True
        if request.boundary_window_size is not None:
            peak_kwargs["boundary_window_size"] = request.boundary_window_size
        if request.boundary_ratio_threshold is not None:
            peak_kwargs["boundary_ratio_threshold"] = request.boundary_ratio_threshold
        peak_kwargs["exclude_k2"] = bool(request.exclude_k2)

        peak_settings = PeakConfig(**peak_kwargs)

    if mode == "k":
        if request.k is None:
            raise HTTPException(status_code=400, detail="k must be provided in k-mode.")
        return pc.run(
            k=request.k,
            top_n=1,
            by_resolution=False,
            plot_scores=False,
            peak_config=peak_settings,
        )
    if mode == "resolution":
        return pc.run(
            k=None,
            top_n=request.top_n,
            by_resolution=True,
            num_bins=request.num_bins,
            max_k=request.max_k,
            max_k_limit=request.max_k_limit,
            plot_scores=False,
            peak_config=peak_settings,
        )
    if mode == "global":
        return pc.run(
            k=None,
            top_n=request.top_n,
            by_resolution=False,
            max_k=request.max_k,
            max_k_limit=request.max_k_limit,
            plot_scores=False,
            peak_config=peak_settings,
        )
    raise HTTPException(status_code=400, detail="Invalid mode.")


def _store_run(
    pc: PhytClust,
    result: dict[str, Any],
    request: PhytclustRequest,
    notes: list[str],
    *,
    newick: str,
    clustered_newick: str,
    construction_key: tuple,
    construction_notes: list[str],
) -> dict[str, Any]:
    """Build the JSON response, cache it, and remember the instance for reuse.

    ``newick`` is the normalized request tree used to detect child-order
    changes; ``clustered_newick`` is the tree returned for drawing.
    """
    global LAST_PC, LAST_RESULT, LAST_NEWICK, LAST_CONSTRUCTION_KEY
    global LAST_CONSTRUCTION_NOTES

    clusters_json = _serialize_clusters(result.get("clusters"))
    scores = result.get("scores")
    if scores is not None:
        scores = [_json_float(value) for value in scores]

    all_clusters_json = None
    all_ks = None
    all_alphas = None
    if request.compute_all_clusters and pc.cluster_maps:
        sorted_ks = sorted(pc.cluster_maps.keys())
        all_ks = sorted_ks
        all_clusters_json = _serialize_clusters([pc.cluster_maps[k] for k in sorted_ks])
        all_alphas = [
            _json_float(pc.alpha_info(k, pc.cluster_maps[k])["alpha"]) for k in sorted_ks
        ]

    alpha_details = result.get("alpha_details") or []
    alphas = result.get("alphas") or []

    bin_ranges = None
    peaks_by_resolution = None
    if result.get("mode") == "resolution":
        raw_ranges = pc.bin_ranges_current
        if raw_ranges:
            bin_ranges = [[int(lo), int(hi)] for lo, hi in raw_ranges]
        raw_pbr = pc.peaks_by_resolution
        if isinstance(raw_pbr, dict):
            peaks_by_resolution = {
                str(label): [int(count) for count in (cluster_counts or [])]
                for label, cluster_counts in raw_pbr.items()
            }

    payload = {
        "mode": result.get("mode"),
        "requested_mode": request.mode,
        "selected_k": result.get("selected_k"),
        "k_values": result.get("k_values"),
        "k": result.get("k"),
        "ks": result.get("ks"),
        "peaks": result.get("peaks"),
        "bin_ranges": bin_ranges,
        "peaks_by_resolution": peaks_by_resolution,
        "alphas": [_json_float(a) for a in alphas],
        "alpha_details": [
            {
                name: (
                    _json_float(value)
                    if isinstance(value, Real) and not isinstance(value, bool)
                    else value
                )
                for name, value in detail.items()
            }
            for detail in alpha_details
        ],
        "notes": [str(n) for n in notes if n is not None],
        "newick": clustered_newick,
        "scores": scores,
        "clusters": clusters_json,
        "all_clusters": all_clusters_json,
        "all_ks": all_ks,
        "all_alphas": all_alphas,
        "outlier_size_threshold": pc.outlier.size_threshold,
    }

    run_id = str(uuid.uuid4())
    payload["run_id"] = run_id
    _cache_put(run_id, pc, payload)
    LAST_PC = pc
    LAST_NEWICK = newick
    LAST_CONSTRUCTION_KEY = construction_key
    LAST_CONSTRUCTION_NOTES = list(construction_notes)
    LAST_RESULT = payload
    return payload


class SaveRequest(BaseModel):
    """Settings for saving displayed or calculated assignments."""

    results_dir: str
    filename: str = "phytclust_results.tsv"
    top_n: PositiveCount = 1
    outlier: bool = True
    output_all: bool = False
    tsv: Optional[str] = None
    run_id: Optional[str] = None


def _results_dir(directory: str) -> Path:
    """Require an output directory inside the server working directory.

    A leading ``~`` is expanded to the home directory before the check, so
    ``~/results`` is accepted only when it lies inside the working directory.
    """
    root = Path.cwd().resolve()
    try:
        candidate = Path(directory).expanduser()
        candidate = (
            candidate if candidate.is_absolute() else root / candidate
        ).resolve()
    except (OSError, RuntimeError, ValueError) as e:
        raise HTTPException(
            status_code=400, detail=f"Invalid output directory {directory!r}: {e}"
        )
    if candidate != root and not candidate.is_relative_to(root):
        raise HTTPException(
            status_code=403,
            detail=(
                f"Output directory {directory!r} is outside the server's working "
                f"directory ({root}). Use a path inside it, or restart the GUI "
                "from where you want the results written."
            ),
        )
    return candidate


@app.post("/api/save")
def save_results(request: SaveRequest):
    """Save assignments without overlapping a clustering run."""
    with _RUN_LOCK:
        return _save_results(request)


def _save_results(request: SaveRequest):
    """Write a displayed table, a stored run, or the current clustering result."""
    pc, result = _cache_get(request.run_id)
    try:
        filename = _export_filename(request.filename)
    except ValueError as error:
        raise HTTPException(status_code=400, detail=str(error)) from error

    out_dir = _results_dir(request.results_dir)
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        raise HTTPException(
            status_code=400,
            detail=f"Cannot create output directory '{request.results_dir}': {e}",
        )
    target = out_dir / filename
    if target.is_symlink() or (target.exists() and not target.is_file()):
        raise HTTPException(
            status_code=400,
            detail=f"Output path {filename!r} exists and is not a regular file.",
        )

    try:
        if request.tsv is not None:
            target.write_text(request.tsv, encoding="utf-8")
        elif request.run_id is not None:
            export_settings = ExportTSVRequest(
                top_n=request.top_n,
                outlier=request.outlier,
                output_all=request.output_all,
            )
            target.write_text(
                _run_tsv(result, export_settings), encoding="utf-8"
            )
        else:
            saved_path = pc.save(
                results_dir=str(out_dir),
                top_n=request.top_n,
                filename=filename,
                outlier=request.outlier,
                output_all=request.output_all,
            )
            if saved_path is None:
                raise HTTPException(
                    status_code=400, detail="No cluster assignments to save."
                )
        return {
            "status": "ok",
            "results_dir": str(out_dir),
            "filename": filename,
        }
    except HTTPException:
        raise
    except (ValueError, PhytClustError, OSError) as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception:
        logger.exception("Failed to save results")
        raise HTTPException(status_code=500, detail="Failed to save results.")


class ExportTSVRequest(BaseModel):
    """Settings for downloading assignments from a stored run."""

    filename: str = "phytclust_results.tsv"
    top_n: PositiveCount = 1
    outlier: bool = True
    output_all: bool = False
    run_id: Optional[str] = None


def _run_tsv(result: dict[str, Any], request: ExportTSVRequest) -> str:
    """Build a TSV from saved assignments without rerunning clustering."""
    if request.output_all:
        cluster_counts = result.get("all_ks") or []
        cluster_maps = result.get("all_clusters") or []
        if not cluster_maps:
            raise ValueError("All-k clusters were not stored for this run.")
    else:
        cluster_counts = (result.get("k_values") or [])[: request.top_n]
        cluster_maps = (result.get("clusters") or [])[: request.top_n]
    if not cluster_counts or len(cluster_counts) != len(cluster_maps):
        raise ValueError("No clusters are available for the requested run.")

    threshold = result.get("outlier_size_threshold")
    by_name: dict[str, dict[int, int]] = {}
    for k, cluster_map in zip(cluster_counts, cluster_maps):
        counts: dict[int, int] = {}
        for cluster_id in cluster_map.values():
            counts[cluster_id] = counts.get(cluster_id, 0) + 1
        for name, cluster_id in cluster_map.items():
            output_cluster_id = cluster_id
            if request.outlier and counts[cluster_id] < (
                threshold if threshold is not None else 2
            ):
                output_cluster_id = -1
            by_name.setdefault(name, {})[int(k)] = output_cluster_id

    output = StringIO()
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    sorted_cluster_counts = sorted(int(k) for k in cluster_counts)
    writer.writerow(["Node Name", *(f"clusters_k{k}" for k in sorted_cluster_counts)])
    for name in sorted(by_name):
        writer.writerow(
            [name, *(by_name[name].get(k, "") for k in sorted_cluster_counts)]
        )
    return output.getvalue()


@app.post("/api/export_tsv")
def export_tsv(request: ExportTSVRequest):
    """Download the requested run, or the latest run when no ID is supplied."""
    with _RUN_LOCK:
        _, result = _cache_get(request.run_id)
    try:
        _export_filename(request.filename)
        return PlainTextResponse(
            content=_run_tsv(result, request),
            media_type="text/tab-separated-values",
        )
    except ValueError as error:
        raise HTTPException(status_code=400, detail=str(error)) from error
