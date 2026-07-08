import os
import logging
from collections import Counter
from typing import Any, Optional

import pandas as pd

from ..algo.dp import cluster_map
from ..metrics.indices import cluster_alpha

logger = logging.getLogger(__name__)


def _safe_filename(filename: str) -> str:
    """Reduce a user-supplied output filename to a safe basename.

    Rejects path separators, parent-dir references, and absolute paths so a
    caller can never write outside ``results_dir`` (path-traversal guard).
    """
    name = os.path.basename(filename)
    if name != filename or "/" in filename or "\\" in filename or name in ("", ".", ".."):
        raise ValueError(
            f"Invalid output filename {filename!r}: must be a plain filename "
            "with no path separators or '..'."
        )
    return name


def save_clusters(
    pc,
    results_dir: str,
    top_n: int = 1,
    filename: str = "phytclust_results.tsv",
    outlier: bool = True,
    n: Optional[int] = None,
    output_all: bool = False,
) -> Optional[str]:
    """Write cluster assignments to a TSV file.

    Returns the path to the written file, or None if nothing was saved.
    """
    filename = _safe_filename(filename)
    os.makedirs(results_dir, exist_ok=True)

    if pc.k is not None:
        ks = [pc.k]
    elif output_all:
        ks = list(range(1, pc.max_k + 1))
    elif n is not None:
        ks = [n]
    else:
        ks = (pc.peaks_by_rank or [])[:top_n]

    outlier_thresh = pc.outlier.size_threshold
    records: list[dict[str, Any]] = []
    alpha_rows: list[dict[str, Any]] = []
    active_tree = (
        pc._tree_wo_outgroup
        if (pc.outgroup and pc._tree_wo_outgroup is not None)
        else pc.tree
    )

    for k_val in ks:
        cmap = cluster_map(pc, k_val)
        if cmap is None:
            continue
        counts = Counter(cmap.values())

        cached_alpha = getattr(pc, "alpha_by_k", {}).get(int(k_val))
        if cached_alpha is None:
            alpha_info = cluster_alpha(active_tree, cmap)
            alpha_info = {"k": int(k_val), **alpha_info}
            if hasattr(pc, "alpha_by_k"):
                pc.alpha_by_k[int(k_val)] = alpha_info
        else:
            alpha_info = cached_alpha
        alpha_rows.append(alpha_info)

        for node, cid in cmap.items():
            if not outlier:
                # `outlier=False` (e.g. --no-outlier) is a master off switch:
                # never mark, regardless of any configured size threshold.
                mark = cid
            elif outlier_thresh is not None:
                mark = -1 if counts[cid] < outlier_thresh else cid
            else:
                mark = -1 if counts[cid] == 1 else cid
            records.append({
                "Node Name": node.name,
                "k": k_val,
                "Cluster ID": mark,
            })

    if not records:
        logger.warning("No clusters to save.")
        return None

    df = pd.DataFrame(records)
    pivot = df.pivot_table(
        index="Node Name", columns="k",
        values="Cluster ID", aggfunc="first",
    )
    pivot.columns = [f"clusters_k{col}" for col in pivot.columns]
    pivot.reset_index(inplace=True)

    out_path = os.path.join(results_dir, filename)
    pivot.to_csv(out_path, index=False, sep="\t")
    logger.info("Wrote clusters to %s", out_path)

    if getattr(pc, "peaks_by_rank", None):
        peaks_path = os.path.join(results_dir, "peaks_by_rank.txt")
        with open(peaks_path, "w") as fh:
            for rank, k_val in enumerate(pc.peaks_by_rank, 1):
                fh.write(f"Rank {rank}: {k_val} clusters\n")
        logger.info("Wrote peaks_by_rank.txt")

    if alpha_rows:
        alpha_path = os.path.join(results_dir, "alphas.tsv")
        pd.DataFrame(alpha_rows).to_csv(alpha_path, sep="\t", index=False)
        logger.info("Wrote alphas.tsv")

    return out_path
