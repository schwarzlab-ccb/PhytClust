import os
import logging
from collections import Counter
from numbers import Integral
from typing import Any, Optional

import numpy as np
import pandas as pd

from ..exceptions import ConfigurationError, InvalidClusteringError, InvalidKError

logger = logging.getLogger(__name__)


def _safe_filename(filename: str) -> str:
    """Require a plain filename without path separators or reserved output names."""
    if not isinstance(filename, str) or not filename or "\x00" in filename:
        raise ValueError(
            "Output filename must be a nonempty string without null characters."
        )
    name = os.path.basename(filename)
    if name != filename or "/" in filename or "\\" in filename or name in (".", ".."):
        raise ValueError(
            f"Invalid output filename {filename!r}: use a plain filename without path separators."
        )
    if name.casefold() in {"alphas.tsv", "peaks_by_rank.txt"}:
        raise ValueError(
            f"Output filename {filename!r} is reserved for export metadata."
        )
    return name


def _positive_count(value, name, error):
    if (
        isinstance(value, (bool, np.bool_))
        or not isinstance(value, Integral)
        or value < 1
    ):
        raise error(f"{name} must be a positive integer.")


def save_clusters(
    pc,
    results_dir: str,
    top_n: int = 1,
    filename: str = "phytclust_results.tsv",
    outlier: bool = True,
    n: Optional[int] = None,
    output_all: bool = False,
) -> Optional[str]:
    """Write selected cluster assignments and alpha statistics as TSV files.

    An explicit n takes priority over the stored fixed count, all-count export,
    and selected peaks, in that order. The public save method calls this n
    argument k. Skip impossible partitions only when exporting all counts.

    With outlier=True, mark clusters below the configured threshold as -1.
    Without a threshold, mark singleton clusters as -1.

    Write peak ranks only for exported peaks; the ranking file is empty for
    fixed-count exports. Return the assignment path, or None for no results.
    """
    clustering = pc
    filename = _safe_filename(filename)
    _positive_count(top_n, "top_n", ConfigurationError)
    if n is not None:
        _positive_count(n, "k", InvalidKError)
    saving_all = n is None and clustering.k is None and output_all
    exporting_peaks = n is None and clustering.k is None and not saving_all
    if n is not None:
        cluster_counts = [n]
    elif clustering.k is not None:
        _positive_count(clustering.k, "k", InvalidKError)
        cluster_counts = [clustering.k]
    elif saving_all:
        maximum = clustering._effective_max_k()
        cluster_counts = list(range(1, min(maximum, clustering.num_terminals) + 1))
    else:
        cluster_counts = list(dict.fromkeys((clustering.peaks_by_rank or [])[:top_n]))

    outlier_size_threshold = clustering.outlier.size_threshold
    records: list[dict[str, Any]] = []
    alpha_rows: list[dict[str, Any]] = []
    exported_counts = []
    if cluster_counts:
        for cluster_count in cluster_counts:
            _positive_count(cluster_count, "k", InvalidKError)
        clustering._ensure_dp(required_cap=max(cluster_counts))
    for cluster_count in cluster_counts:
        try:
            cluster_map = clustering._clusters(cluster_count)
        except InvalidClusteringError:
            if not saving_all:
                raise
            logger.info(
                "Skipping impossible partition with k=%d during all-count export.",
                cluster_count,
            )
            continue
        if not cluster_map:
            continue
        cluster_sizes = Counter(cluster_map.values())
        alpha_rows.append(clustering.alpha_info(cluster_count, cluster_map))
        exported_counts.append(cluster_count)
        for node, cluster_id in cluster_map.items():
            output_cluster_id = cluster_id
            if outlier:
                threshold = (
                    outlier_size_threshold if outlier_size_threshold is not None else 2
                )
                if cluster_sizes[cluster_id] < threshold:
                    output_cluster_id = -1
            records.append(
                {
                    "Node Name": node.name,
                    "k": cluster_count,
                    "Cluster ID": output_cluster_id,
                }
            )

    if not records:
        logger.warning("No cluster assignments to save.")
        return None
    frame = pd.DataFrame(records)
    if frame["Node Name"].isna().any() or (frame["Node Name"] == "").any():
        raise InvalidClusteringError("Export requires named leaves.")
    if frame.duplicated(["Node Name", "k"]).any():
        raise InvalidClusteringError(
            "Export contains duplicate leaf assignments for the same cluster count."
        )
    assignments = frame.pivot(index="Node Name", columns="k", values="Cluster ID")
    assignments = assignments.sort_index()
    assignments.columns = [f"clusters_k{count}" for count in assignments.columns]
    assignments.reset_index(inplace=True)

    os.makedirs(results_dir, exist_ok=True)
    output_path = os.path.join(results_dir, filename)
    assignments.to_csv(output_path, index=False, sep="\t")
    logger.info("Wrote cluster assignments to %s", output_path)

    selected_peaks = (
        [count for count in cluster_counts if count in exported_counts]
        if exporting_peaks
        else []
    )
    peaks_path = os.path.join(results_dir, "peaks_by_rank.txt")
    with open(peaks_path, "w", encoding="utf-8") as stream:
        for rank, cluster_count in enumerate(selected_peaks, 1):
            stream.write(f"Rank {rank}: {cluster_count} clusters\n")
    pd.DataFrame(alpha_rows).to_csv(
        os.path.join(results_dir, "alphas.tsv"), sep="\t", index=False
    )
    return output_path
