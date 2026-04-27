# src/phytclust/gui/api.py

import logging
import math
import os
import tempfile
import uuid
from collections import OrderedDict
from io import StringIO
from pathlib import Path
from typing import Any, Optional

from Bio import Phylo
from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import HTMLResponse, PlainTextResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from pydantic import BaseModel

from phytclust.algo.core import PhytClust
from phytclust.config import OutlierConfig, PeakConfig

logger = logging.getLogger("phytclust.gui")


def _normalize_newick(newick: str) -> str:
    """Quote any unquoted names containing spaces so JS and BioPython parse the same names.

    The JS newick tokenizer splits only on ( ) , : ; so spaces are part of names.
    BioPython splits on whitespace too, so it only keeps the last word. This
    pre-pass wraps such names in single quotes before BioPython ever sees the string.
    """
    result: list[str] = []
    i = 0
    n = len(newick)

    while i < n:
        ch = newick[i]

        if ch in "(),;":
            result.append(ch)
            i += 1

        elif ch in "'\"":
            # Already-quoted string — pass through unchanged
            quote = ch
            j = i + 1
            while j < n and newick[j] != quote:
                j += 1
            result.append(newick[i : j + 1])
            i = j + 1

        elif ch == ":":
            result.append(":")
            i += 1
            # Read the branch-length value verbatim (digits, dot, sign, exponent)
            while i < n and newick[i] not in "(),;":
                result.append(newick[i])
                i += 1

        else:
            # Unquoted name token (leaf or internal node label)
            j = i
            while j < n and newick[j] not in "(),;:'\"\\":
                j += 1
            token = newick[i:j]
            stripped = token.strip()
            if stripped and " " in stripped:
                leading = token[: len(token) - len(token.lstrip())]
                result.append(leading + "'" + stripped + "'")
            else:
                result.append(token)
            i = j

    return "".join(result)

# Set PHYTCLUST_PUBLIC_MODE=1 to enable public-safe restrictions:
#   - leaf count capped at PUBLIC_MAX_TIPS
#   - /api/save disabled (arbitrary server-path writes)
#   - per-run result cache so concurrent users don't clobber each other
PUBLIC_MODE: bool = os.getenv("PHYTCLUST_PUBLIC_MODE", "0") == "1"
PUBLIC_MAX_TIPS: int = int(os.getenv("PHYTCLUST_MAX_TIPS", "10000"))

# LRU result cache keyed by run_id (UUID). Holds last _CACHE_MAX results.
_CACHE: OrderedDict[str, tuple[PhytClust, dict[str, Any]]] = OrderedDict()
_CACHE_MAX = 20

# ------------------------------------------------------
# Create the FastAPI app
# ------------------------------------------------------
app = FastAPI(title="PhytClust Web API")

STATIC_DIR = Path(__file__).parent / "static"
TEMPLATES_DIR = Path(__file__).parent / "templates"
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")
templates = Jinja2Templates(directory=TEMPLATES_DIR)

# Single-user fallback (local use). In PUBLIC_MODE the per-run cache is used
# instead, so concurrent users don't clobber each other's results.
LAST_PC: Optional[PhytClust] = None
LAST_RESULT: Optional[dict[str, Any]] = None
LAST_NEWICK: Optional[str] = None
LAST_CONSTRUCTION_KEY: Optional[tuple] = None


def _cache_put(run_id: str, pc: PhytClust, result: dict[str, Any]) -> None:
    _CACHE[run_id] = (pc, result)
    _CACHE.move_to_end(run_id)
    while len(_CACHE) > _CACHE_MAX:
        _CACHE.popitem(last=False)


def _cache_get(run_id: Optional[str]) -> tuple[PhytClust, dict[str, Any]]:
    if run_id and run_id in _CACHE:
        return _CACHE[run_id]
    if LAST_PC is not None and LAST_RESULT is not None:
        return LAST_PC, LAST_RESULT
    raise HTTPException(
        status_code=400,
        detail="No PhytClust result available. Please run PhytClust first.",
    )


# ------------------------------------------------------
# Serve index.html at root URL
# ------------------------------------------------------
@app.get("/", response_class=HTMLResponse)
def serve_frontend(request: Request):
    return templates.TemplateResponse(request, "index.html")


# ------------------------------------------------------
# Request model
# ------------------------------------------------------
class PhytclustRequest(BaseModel):
    newick: str
    mode: str = "global"
    k: Optional[int] = None
    top_n: int = 1
    by_resolution: bool = False
    num_bins: Optional[int] = None
    max_k: Optional[int] = None
    max_k_limit: Optional[float] = None
    outgroup: Optional[str] = None
    root_taxon: Optional[str] = None
    min_cluster_size: Optional[int] = None
    compute_all_clusters: bool = False

    # Branch-support options
    use_branch_support: bool = False
    min_support: Optional[float] = None
    support_weight: Optional[float] = None

    # Outlier config
    outlier_size_threshold: Optional[int] = None
    outlier_prefer_fewer: bool = False
    outlier_ratio_weight: Optional[float] = None
    outlier_ratio_mode: Optional[str] = None

    # Zero-length / polytomy handling
    no_split_zero_length: bool = False
    optimize_polytomies: bool = True

    # Peak config
    lambda_weight: float = 0.7
    ranking_mode: str = "adjusted"
    min_prominence: Optional[float] = None
    use_relative_prominence: bool = False
    boundary_window_size: Optional[int] = None
    boundary_ratio_threshold: Optional[float] = None


# ------------------------------------------------------
# Helpers
# ------------------------------------------------------
def _require_last_result(run_id: Optional[str] = None) -> tuple[PhytClust, dict[str, Any]]:
    return _cache_get(run_id)


def _leaf_key_to_name(x: Any) -> str:
    name = getattr(x, "name", None)
    if name is not None:
        return str(name)
    return str(x)


def _safe_float(v) -> float | None:
    """Return a JSON-safe float, converting nan/inf → None."""
    if v is None:
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def _serialize_clusters(raw_clusters: Any) -> list[dict[str, int]]:
    out: list[dict[str, int]] = []

    if raw_clusters is None:
        return out

    if isinstance(raw_clusters, dict):
        out.append({_leaf_key_to_name(k): int(v) for k, v in raw_clusters.items()})
        return out

    for cmap in raw_clusters:
        out.append({_leaf_key_to_name(k): int(v) for k, v in cmap.items()})

    return out


# ------------------------------------------------------
# API endpoint
# ------------------------------------------------------
def _construction_key(req: PhytclustRequest, newick: str) -> tuple:
    """Hashable key covering all PhytClust constructor parameters."""
    return (
        newick,
        req.outgroup,
        req.root_taxon,
        req.min_cluster_size,
        req.use_branch_support,
        req.min_support,
        req.support_weight,
        req.outlier_size_threshold,
        req.outlier_prefer_fewer,
        req.outlier_ratio_weight,
        req.outlier_ratio_mode,
        req.no_split_zero_length,
        req.optimize_polytomies,
    )


@app.post("/api/run")
def run_phytclust(req: PhytclustRequest):
    global LAST_PC, LAST_RESULT, LAST_NEWICK, LAST_CONSTRUCTION_KEY
    notes: list[str] = []

    if not req.newick.strip():
        raise HTTPException(status_code=400, detail="Empty Newick string.")

    newick = _normalize_newick(req.newick)

    if PUBLIC_MODE:
        try:
            _t = Phylo.read(StringIO(newick.strip()), "newick")
            n_tips = _t.count_terminals()
        except Exception:
            n_tips = 0
        if n_tips > PUBLIC_MAX_TIPS:
            raise HTTPException(
                status_code=400,
                detail=f"Tree has {n_tips} tips; the public demo is limited to {PUBLIC_MAX_TIPS}.",
            )

    mode = req.mode.lower()
    if mode not in {"k", "global", "resolution"}:
        raise HTTPException(
            status_code=400, detail="mode must be 'k', 'global', or 'resolution'."
        )

    current_key = _construction_key(req, newick)
    if (
        LAST_CONSTRUCTION_KEY is not None
        and LAST_PC is not None
        and LAST_CONSTRUCTION_KEY == current_key
    ):
        pc = LAST_PC
        pc.compute_all_clusters = req.compute_all_clusters
        # A previous fixed-k run leaves pc.k populated; PhytClust._run_inner
        # falls back to self.k when k=None is passed, which would then collide
        # with by_resolution=True on the next run. Reset run-mode state on
        # reuse so this run's mode is what actually takes effect.
        if mode != "k":
            pc.k = None
        notes.append("Reusing previous PhytClust instance.")
    else:
        try:
            # Build outlier config
            outlier_kwargs: dict[str, Any] = {}
            if req.outlier_size_threshold is not None:
                outlier_kwargs["size_threshold"] = req.outlier_size_threshold
            if req.outlier_prefer_fewer:
                outlier_kwargs["prefer_fewer"] = True
            if req.outlier_ratio_weight is not None:
                outlier_kwargs["ratio_weight"] = req.outlier_ratio_weight
            if req.outlier_ratio_mode is not None:
                outlier_kwargs["ratio_mode"] = req.outlier_ratio_mode

            kwargs: dict[str, Any] = dict(
                tree=newick,
                outgroup=req.outgroup,
                compute_all_clusters=req.compute_all_clusters,
                use_branch_support=req.use_branch_support,
                no_split_zero_length=req.no_split_zero_length,
                optimize_polytomies=req.optimize_polytomies,
            )
            if req.root_taxon:
                kwargs["root_taxon"] = req.root_taxon
            if req.min_cluster_size is not None:
                kwargs["min_cluster_size"] = req.min_cluster_size
            if req.min_support is not None:
                kwargs["min_support"] = req.min_support
            if req.support_weight is not None:
                kwargs["support_weight"] = req.support_weight
            if outlier_kwargs:
                kwargs["outlier"] = OutlierConfig(**outlier_kwargs)

            pc = PhytClust(**kwargs)
        except Exception as e:
            raise HTTPException(status_code=400, detail=f"Could not parse Newick: {e}")

    try:
        # Build peak config
        peak_kwargs: dict[str, Any] = {"lambda_weight": req.lambda_weight}
        if req.ranking_mode:
            peak_kwargs["ranking_mode"] = req.ranking_mode
        if req.min_prominence is not None:
            peak_kwargs["min_prominence"] = req.min_prominence
        if req.use_relative_prominence:
            peak_kwargs["use_relative_prominence"] = True
        if req.boundary_window_size is not None:
            peak_kwargs["boundary_window_size"] = req.boundary_window_size
        if req.boundary_ratio_threshold is not None:
            peak_kwargs["boundary_ratio_threshold"] = req.boundary_ratio_threshold

        cfg = PeakConfig(**peak_kwargs)

        if mode == "k":
            if req.k is None:
                raise HTTPException(
                    status_code=400, detail="k must be provided in k-mode."
                )
            result = pc.run(
                k=req.k,
                top_n=1,
                by_resolution=False,
                plot_scores=False,
                peak_config=cfg,
            )

        elif mode == "resolution":
            result = pc.run(
                k=None,
                top_n=req.top_n,
                by_resolution=True,
                num_bins=req.num_bins,
                max_k=req.max_k,
                max_k_limit=req.max_k_limit,
                plot_scores=False,
                peak_config=cfg,
            )

        elif mode == "global":
            result = pc.run(
                k=None,
                top_n=req.top_n,
                by_resolution=False,
                max_k=req.max_k,
                max_k_limit=req.max_k_limit,
                plot_scores=False,
                peak_config=cfg,
            )
        else:
            raise HTTPException(status_code=400, detail="Invalid mode.")

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"PhytClust error: {e}")

    clusters_json = _serialize_clusters(result.get("clusters"))
    scores = result.get("scores")
    if scores is not None:
        try:
            scores = list(scores)
        except Exception:
            pass

    # When compute_all_clusters is on, include every cached k->cluster mapping
    all_clusters_json = None
    all_ks = None
    all_alphas = None
    if req.compute_all_clusters and pc.clusters:
        sorted_ks = sorted(pc.clusters.keys())
        all_ks = sorted_ks
        all_clusters_json = _serialize_clusters([pc.clusters[kv] for kv in sorted_ks])
        from ..metrics.indices import cluster_alpha

        active_tree = (
            pc._tree_wo_outgroup
            if (pc.outgroup and pc._tree_wo_outgroup is not None)
            else pc.tree
        )
        all_alphas = []
        for kv in sorted_ks:
            cached = pc.alpha_by_k.get(int(kv))
            if cached is None:
                info = cluster_alpha(active_tree, pc.clusters[kv])
                cached = {"k": int(kv), **info}
                pc.alpha_by_k[int(kv)] = cached
            all_alphas.append(_safe_float(cached["alpha"]))

    alpha_details = result.get("alpha_details") or []
    alphas = result.get("alphas") or []

    payload = {
        "mode": result.get("mode"),
        "selected_k": result.get("selected_k"),
        "k_values": result.get("k_values"),
        "k": result.get("k"),
        "ks": result.get("ks"),
        "peaks": result.get("peaks"),
        "alphas": [_safe_float(a) for a in alphas],
        "alpha_details": [
            {k: (_safe_float(v) if isinstance(v, (int, float)) else v) for k, v in d.items()}
            for d in alpha_details
        ],
        "outgroup": result.get("outgroup"),
        "notes": [str(n) for n in notes if n is not None],
        "newick": newick,
        "scores": scores,
        "clusters": clusters_json,
        "all_clusters": all_clusters_json,
        "all_ks": all_ks,
        "all_alphas": all_alphas,
    }

    run_id = str(uuid.uuid4())
    payload["run_id"] = run_id
    _cache_put(run_id, pc, payload)
    LAST_PC = pc
    LAST_NEWICK = newick
    LAST_CONSTRUCTION_KEY = current_key
    LAST_RESULT = payload
    return payload


# ------------------------------------------------------
# Save endpoint
# ------------------------------------------------------
class SaveRequest(BaseModel):
    results_dir: str
    filename: str = "phytclust_results.tsv"
    top_n: int = 1
    outlier: bool = True
    output_all: bool = False


@app.post("/api/save")
def save_results(req: SaveRequest):
    if PUBLIC_MODE:
        raise HTTPException(
            status_code=403,
            detail="Save to server is disabled in public mode. Use Export TSV instead.",
        )
    pc, _ = _require_last_result()

    try:
        out_dir = Path(req.results_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        raise HTTPException(
            status_code=400,
            detail=f"Cannot create output directory '{req.results_dir}': {e}",
        )

    try:
        pc.save(
            results_dir=req.results_dir,
            top_n=req.top_n,
            filename=req.filename,
            outlier=req.outlier,
            output_all=req.output_all,
        )
        return {
            "status": "ok",
            "results_dir": req.results_dir,
            "filename": req.filename,
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to save results: {e}")


class ExportTSVRequest(BaseModel):
    filename: str = "phytclust_results.tsv"
    top_n: int = 1
    outlier: bool = True
    output_all: bool = False
    run_id: Optional[str] = None


@app.post("/api/export_tsv")
def export_tsv(req: ExportTSVRequest):
    pc, _ = _require_last_result(req.run_id)

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            pc.save(
                results_dir=tmpdir,
                top_n=req.top_n,
                filename=req.filename,
                outlier=req.outlier,
                output_all=req.output_all,
            )
            tsv_path = Path(tmpdir) / req.filename
            if not tsv_path.exists():
                raise HTTPException(status_code=500, detail="TSV file not generated.")
            text = tsv_path.read_text(encoding="utf-8")
            return PlainTextResponse(
                content=text, media_type="text/tab-separated-values"
            )
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to export TSV: {e}")
