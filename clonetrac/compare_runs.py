"""Compare selected partitions between two Clonetrac runs."""

import json
from pathlib import Path

import pandas as pd


def partition_statistics(folder: Path, rank: int, cluster_count: int, threshold: int):
    """Measure small clusters and the largest cluster in one saved partition."""
    table = pd.read_csv(
        folder / f"rank_{rank}_k{cluster_count}" / "cell_assignments.tsv", sep="\t"
    )
    sizes = table.groupby("cluster").size()
    return {
        "k": cluster_count,
        "small_clusters": int((sizes < threshold).sum()),
        "cells_in_small_clusters": int(sizes[sizes < threshold].sum()),
        "largest_cluster_cells": int(sizes.max()),
        "largest_cluster_percent": round(float(sizes.max() / sizes.sum() * 100), 1),
        "median_cluster_size": float(sizes.median()),
        "relative_size_variation": round(float(sizes.std(ddof=0) / sizes.mean()), 3),
        "cluster_sizes": sorted(sizes.tolist(), reverse=True),
    }


def compare_runs(baseline: Path, trial: Path, threshold: int = 3):
    """Compare every saved rank and summarize the first ranked partition."""
    primary_rows = []
    partition_rows = []
    for folder in sorted(trial.glob("NB*")):
        if not (folder / "run_summary.json").is_file():
            continue
        primary = {"sample": folder.name}
        for label, run_folder in [
            ("before", baseline / folder.name),
            ("after", folder),
        ]:
            summary = json.loads((run_folder / "run_summary.json").read_text())
            primary[f"{label}_selected_k"] = summary["selected_k"]
            primary[f"{label}_fallback_used"] = summary["fallback_used"]
            for rank, count in enumerate(summary["selected_k"], start=1):
                stats = partition_statistics(run_folder, rank, count, threshold)
                partition_rows.append(
                    {"sample": folder.name, "run": label, "rank": rank, **stats}
                )
                if rank == 1:
                    primary.update(
                        {f"{label}_{name}": value for name, value in stats.items()}
                    )
        primary_rows.append(primary)
    comparison = pd.DataFrame(primary_rows)
    comparison.to_csv(trial / "comparison.tsv", sep="\t", index=False)
    pd.DataFrame(partition_rows).to_csv(
        trial / "comparison_all_ranks.tsv", sep="\t", index=False
    )
    return comparison
