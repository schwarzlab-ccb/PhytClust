"""Run PhytClust on the original Clonetrac trees and plot their copy-number profiles."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from Bio import Phylo

from phytclust import PhytClust
from phytclust.config.outlier import OutlierConfig
from phytclust.config.peak import PeakConfig
from phytclust.viz.palette import expand_palette
from phytclust.viz.plots import plot_cluster
from phytclust.viz.scores import plot_scores

ROOT = Path(__file__).resolve().parent


def normalize_cell_name(name: str) -> str:
    """Match tree and profile names to the normal-cell list."""
    return re.sub(r"\.final\.bam$", "", name.removeprefix("CTR_"))


def discover_samples(root: Path = ROOT) -> list[tuple[str, Path, Path]]:
    """Find original NB tree/profile pairs directly inside the project folder."""
    samples = []
    for tree_path in sorted(root.glob("NB*_final_tree.new")):
        prefix = tree_path.name.removesuffix("_final_tree.new")
        candidates = [
            root / f"{prefix}_final_cn_profiles.tsv",
            root / f"{prefix}_final_cn_profiles_without_ecdna.tsv",
        ]
        profile_paths = [path for path in candidates if path.is_file()]
        if len(profile_paths) != 1:
            raise ValueError(
                f"Expected one profile file for {tree_path.name}; found {len(profile_paths)}."
            )
        samples.append((tree_path.name.split("_", 1)[0], tree_path, profile_paths[0]))
    if not samples:
        raise ValueError(f"No original NB trees found in {root}.")
    return samples


def load_normal_cells(root: Path = ROOT) -> set[str]:
    """Read normal-cell identifiers and check the two supplied lists agree."""
    table_names = set(
        pd.read_csv(root / "clonetrac_normal_cells.tsv", sep="\t")["sample_id"]
    )
    text_names = set(
        re.findall(r'"([^\"]+)"', (root / "clonetrac_normal_cells.txt").read_text())
    )
    if table_names != text_names:
        raise ValueError("The normal-cell TSV and text lists disagree.")
    return {normalize_cell_name(name) for name in table_names}


def prepare_inputs(
    tree_path: Path, profile_path: Path, normal_cells: set[str], output_dir: Path
):
    """Remove normal cells from copies of the tree and profile table."""
    tree = Phylo.read(tree_path, "newick")
    original_leaves = tree.get_terminals()
    names = [leaf.name for leaf in original_leaves]
    if any(not name for name in names) or len(names) != len(set(names)):
        raise ValueError(
            f"Tree leaf names must be unique and nonempty: {tree_path.name}."
        )
    removed_names = [
        name for name in names if normalize_cell_name(name) in normal_cells
    ]
    if len(names) - len(removed_names) < 2:
        raise ValueError(
            "Fewer than two leaves would remain after removing normal cells."
        )
    for leaf in original_leaves:
        if leaf.name in removed_names:
            tree.prune(leaf)
    profiles = pd.read_csv(
        profile_path, sep="\t", dtype={"sample_id": str, "chrom": str}
    )
    required_columns = {"sample_id", "chrom", "start", "end", "major"}
    if not required_columns.issubset(profiles.columns):
        raise ValueError(f"Profile table requires columns {sorted(required_columns)}.")
    normal_rows = profiles["sample_id"].map(normalize_cell_name).isin(normal_cells)
    profiles = profiles.loc[~normal_rows].copy()
    profiles.to_csv(output_dir / "cn_profiles.no_normals.tsv", sep="\t", index=False)
    Phylo.write(
        tree,
        output_dir / "tree.no_normals.newick",
        "newick",
        format_branch_length="%.15g",
    )
    pd.DataFrame(
        {
            "cell": removed_names,
            "normal_cell_id": [normalize_cell_name(name) for name in removed_names],
        }
    ).to_csv(output_dir / "removed_normal_cells.tsv", sep="\t", index=False)
    return tree, profiles, len(names), removed_names


def chromosome_order(name: str):
    """Sort numbered chromosomes before X, Y, and other names."""
    label = name.removeprefix("chr")
    if label.isdigit():
        return (0, int(label))
    return (1, {"X": "0", "Y": "1", "M": "2", "MT": "2"}.get(label, label))


def profile_matrix(profiles: pd.DataFrame, leaf_names: list[str]):
    """Arrange complete copy-number profiles in the tree's leaf order."""
    selected = profiles.loc[profiles["sample_id"].isin(leaf_names)].copy()
    keys = ["sample_id", "chrom", "start", "end"]
    if selected.duplicated(keys).any():
        raise ValueError("Duplicate copy-number segments for a cell.")
    matrix = selected.pivot(
        index="sample_id", columns=["chrom", "start", "end"], values="major"
    )
    missing = set(leaf_names) - set(matrix.index)
    if missing:
        raise ValueError(f"Copy-number profiles are missing for {sorted(missing)[:5]}.")
    bins = sorted(
        matrix.columns,
        key=lambda segment: (chromosome_order(segment[0]), segment[1], segment[2]),
    )
    matrix = matrix.reindex(index=leaf_names, columns=pd.MultiIndex.from_tuples(bins))
    values = matrix.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError(
            "Copy-number profiles must have finite, nonnegative values for every cell and segment."
        )
    # Each segment occupies its genomic length; unmeasured gaps are omitted.
    widths = np.array([end - start + 1 for _, start, end in bins], dtype=float)
    if (widths <= 0).any():
        raise ValueError("Copy-number segment ends must be at or after their starts.")
    for previous, current in zip(bins, bins[1:]):
        if previous[0] == current[0] and current[1] <= previous[2]:
            raise ValueError("Copy-number segments overlap.")
    return values, bins, np.r_[0, np.cumsum(widths)]


def plot_tree_and_profiles(
    tree, clusters, profiles, sample_id: str, cluster_count: int
):
    """Draw cluster colours and copy-number rows at the same leaf coordinates."""
    leaves = tree.get_terminals()
    leaf_names = [leaf.name for leaf in leaves]
    values, bins, x_edges = profile_matrix(profiles, leaf_names)
    palette = expand_palette(max(len(set(clusters.values())), 8))
    cluster_colors = dict(zip(sorted(set(clusters.values())), palette))
    height = max(7, min(18, len(leaves) * 0.04))
    figure = plt.figure(figsize=(18, height))
    grid = figure.add_gridspec(1, 4, width_ratios=[3.3, 0.14, 8.5, 0.22], wspace=0.045)
    tree_axis, cluster_axis, profile_axis, colorbar_axis = [
        figure.add_subplot(grid[0, index]) for index in range(4)
    ]
    plot_cluster(
        cluster=clusters,
        tree=tree,
        ax=tree_axis,
        palette=palette,
        show_terminal_labels=False,
        show_cluster_bars=False,
        colour_branches_by_cluster=True,
        marker_size=3,
        show_branch_lengths=False,
        title="",
        line_width=0.65,
    )
    y_edges = np.arange(len(leaves) + 1) - 0.5
    row_colors = np.array([cluster_colors[clusters[leaf]] for leaf in leaves])
    cluster_axis.imshow(
        row_colors[:, None, :],
        extent=(0, 1, len(leaves) - 0.5, -0.5),
        aspect="auto",
        interpolation="nearest",
    )
    cluster_axis.set_axis_off()
    maximum = max(3, int(np.ceil(values.max())))
    norm = TwoSlopeNorm(vmin=0, vcenter=2, vmax=maximum)
    mesh = profile_axis.pcolormesh(
        x_edges,
        y_edges,
        values,
        cmap="coolwarm",
        norm=norm,
        shading="flat",
        linewidth=0,
        antialiased=False,
        rasterized=True,
    )
    chromosome_starts = [0] + [
        index for index in range(1, len(bins)) if bins[index][0] != bins[index - 1][0]
    ]
    chromosome_ends = chromosome_starts[1:] + [len(bins)]
    for end in chromosome_ends[:-1]:
        profile_axis.axvline(x_edges[end], color="#444444", linewidth=0.55)
    profile_axis.set_xticks(
        [
            (x_edges[start] + x_edges[end]) / 2
            for start, end in zip(chromosome_starts, chromosome_ends)
        ]
    )
    profile_axis.set_xticklabels(
        [bins[start][0].removeprefix("chr") for start in chromosome_starts],
        fontsize=9,
        rotation=60,
        ha="left",
    )
    profile_axis.tick_params(
        axis="x", length=0, labeltop=True, labelbottom=False, pad=5
    )
    profile_axis.set_yticks([])
    for spine in profile_axis.spines.values():
        spine.set_visible(False)
    for axis in (tree_axis, cluster_axis, profile_axis):
        axis.set_ylim(len(leaves) - 0.5, -0.5)
    tree_axis.set_xlabel("Branch length", fontsize=9)
    tree_axis.set_title("PhytClust clusters", fontsize=11, fontfamily="Liberation Sans")
    profile_axis.set_title(
        "Copy-number profiles", fontsize=11, fontfamily="Liberation Sans", pad=25
    )
    colorbar = figure.colorbar(mesh, cax=colorbar_axis)
    colorbar.set_label("Total copy number (major column)", fontsize=9)
    colorbar.set_ticks(np.arange(maximum + 1))
    colorbar.ax.tick_params(labelsize=8, length=0)
    colorbar.outline.set_visible(False)
    figure.suptitle(
        f"{sample_id}: PhytClust’s clusters at k = {cluster_count}",
        fontsize=14,
        fontfamily="Liberation Sans",
        y=0.985,
    )
    cluster_sizes = pd.Series(list(clusters.values())).value_counts()
    if len(cluster_colors) <= 16:
        figure.legend(
            handles=[
                Patch(
                    facecolor=color,
                    label=f"C{cluster_id} (n={cluster_sizes[cluster_id]})",
                )
                for cluster_id, color in cluster_colors.items()
            ],
            loc="lower center",
            ncol=min(8, len(cluster_colors)),
            frameon=False,
            fontsize=9,
        )
    else:
        figure.text(
            0.5,
            0.015,
            f"{len(cluster_colors)} clusters, {len(leaves)} cells",
            ha="center",
            fontsize=9,
            color="#555555",
        )
    figure.subplots_adjust(
        top=0.91,
        bottom=0.12 if len(cluster_colors) <= 16 else 0.08,
        left=0.04,
        right=0.94,
    )
    return figure, leaf_names


def process_sample(
    sample_id,
    tree_path,
    profile_path,
    normal_cells,
    *,
    output_root=ROOT / "phytclust_outputs",
    top_n=3,
    max_k_fraction=0.5,
    prefer_fewer_outliers=False,
    outlier_size_threshold=None,
    max_k_limit=None,
    partition_preference="none",
    partition_weight=0.5,
):
    """Save filtered inputs, scores, assignments, and plots for one sample."""
    outlier_settings = OutlierConfig(
        size_threshold=outlier_size_threshold, prefer_fewer=prefer_fewer_outliers
    )
    peak_settings = PeakConfig(
        partition_preference=partition_preference, partition_weight=partition_weight
    )
    output_dir = Path(output_root) / sample_id
    output_dir.mkdir(parents=True, exist_ok=True)
    tree, profiles, original_count, removed_names = prepare_inputs(
        tree_path, profile_path, normal_cells, output_dir
    )
    outgroup = (
        "diploid"
        if any(leaf.name == "diploid" for leaf in tree.get_terminals())
        else None
    )
    tumor_count = tree.count_terminals() - int(outgroup is not None)
    maximum_k = max(2, min(tumor_count, round(tumor_count * max_k_fraction)))
    if max_k_limit is not None:
        if not math.isfinite(max_k_limit) or not 0 < max_k_limit <= 1:
            raise ValueError("max_k_limit must be between 0 and 1.")
        maximum_k = max(2, math.ceil(tumor_count * max_k_limit))
    print(
        f"{sample_id}: removed {len(removed_names)} normal cells; {tumor_count} tumour cells; max k = {maximum_k}",
        flush=True,
    )
    clustering = PhytClust(
        tree,
        outgroup=outgroup,
        no_split_zero_length=True,
        compute_all_clusters=False,
        max_k=None if max_k_limit is not None else maximum_k,
        max_k_limit=0.9 if max_k_limit is None else max_k_limit,
        outlier=outlier_settings,
    )
    result = clustering.run(top_n=top_n, plot_scores=False, peak_config=peak_settings)
    pd.DataFrame(getattr(clustering, "peak_ranking_details", None) or []).to_csv(
        output_dir / "peak_ranking.tsv", sep="\t", index=False
    )
    selected_counts = list(result["peaks"])
    fallback_used = False
    if not selected_counts:
        scores = np.asarray(clustering.scores, dtype=float)
        available_indices = np.flatnonzero(np.isfinite(scores))
        if not len(available_indices):
            raise ValueError(
                f"No finite scores or selected partitions for {sample_id}."
            )
        best_index = available_indices[np.argmax(scores[available_indices])]
        selected_counts = [int(best_index + 1)]
        fallback_used = True
        print(
            f"WARNING: {sample_id} has no detected peaks; saving the highest finite score at k = {selected_counts[0]}.",
            flush=True,
        )
    # Keep the original peak ranking intact when retrieving specific partitions.
    from phytclust.algo.dp import cluster_map

    partitions = [cluster_map(clustering, count) for count in selected_counts]
    plot_tree = clustering._tree_wo_outgroup if outgroup else clustering.tree
    figure = plot_scores(
        clustering,
        peaks=selected_counts,
        title=f"{sample_id}: clustering validity index",
    )
    figure.savefig(output_dir / "scores.png", dpi=180, bbox_inches="tight")
    figure.savefig(output_dir / "scores.pdf", bbox_inches="tight")
    plt.close(figure)
    pd.DataFrame(
        {"k": np.arange(1, len(clustering.scores) + 1), "score": clustering.scores}
    ).to_csv(output_dir / "scores.tsv", sep="\t", index=False)
    all_assignments = []
    for rank, (count, partition) in enumerate(
        zip(selected_counts, partitions), start=1
    ):
        partition_dir = output_dir / f"rank_{rank}_k{count}"
        partition_dir.mkdir(exist_ok=True)
        clustering.save(str(partition_dir), k=count, outlier=False)
        assignments = pd.DataFrame(
            {
                "cell": [leaf.name for leaf in partition],
                "cluster": list(partition.values()),
            }
        )
        assignments["normalized_cell"] = assignments["cell"].map(normalize_cell_name)
        assignments.insert(0, "k", count)
        assignments.insert(0, "rank", rank)
        assignments.to_csv(
            partition_dir / "cell_assignments.tsv", sep="\t", index=False
        )
        all_assignments.append(assignments)
        figure = plot_cluster(
            partition,
            plot_tree,
            show_terminal_labels=False,
            colour_branches_by_cluster=True,
            show_cluster_bars=True,
            marker_size=5,
            height_scale=0.04,
            title=f"{sample_id}: PhytClust’s clusters at k = {count}",
        )
        figure.savefig(
            partition_dir / "tree_clusters.png", dpi=180, bbox_inches="tight"
        )
        plt.close(figure)
        figure, leaf_order = plot_tree_and_profiles(
            plot_tree, partition, profiles, sample_id, count
        )
        figure.savefig(partition_dir / "tree_and_cnp.png", dpi=180, bbox_inches="tight")
        figure.savefig(partition_dir / "tree_and_cnp.pdf", bbox_inches="tight")
        plt.close(figure)
        pd.DataFrame(
            {
                "plot_row": np.arange(len(leaf_order)),
                "cell": leaf_order,
                "cluster": [partition[leaf] for leaf in plot_tree.get_terminals()],
            }
        ).to_csv(partition_dir / "plot_row_order.tsv", sep="\t", index=False)
        print(
            f"  rank {rank}: k = {count}; {partition_dir.relative_to(output_root)}",
            flush=True,
        )
    pd.concat(all_assignments, ignore_index=True).to_csv(
        output_dir / "all_selected_assignments.tsv", sep="\t", index=False
    )
    summary = {
        "sample": sample_id,
        "tree_file": tree_path.name,
        "profile_file": profile_path.name,
        "original_leaves": original_count,
        "normal_cells_removed": len(removed_names),
        "tumour_cells": tumor_count,
        "outgroup": outgroup,
        "max_k": maximum_k,
        "max_k_limit": max_k_limit,
        "partition_preference": partition_preference,
        "partition_weight": partition_weight,
        "prefer_fewer_outliers": prefer_fewer_outliers,
        "outlier_size_threshold": outlier_size_threshold,
        "selected_k": selected_counts,
        "fallback_used": fallback_used,
        "normal_cells_remaining": sum(
            normalize_cell_name(leaf.name) in normal_cells
            for leaf in plot_tree.get_terminals()
        ),
        "output_dir": str(output_dir),
    }
    (output_dir / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sample",
        action="append",
        help="Sample ID to process; repeat for multiple samples.",
    )
    parser.add_argument("--top-n", type=int, default=3)
    parser.add_argument("--max-k-fraction", type=float, default=0.5)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "phytclust_outputs")
    parser.add_argument(
        "--prefer-fewer-outliers",
        action="store_true",
        help="Prioritize fewer clusters below the size threshold at each k.",
    )
    parser.add_argument(
        "--outlier-size-threshold",
        type=int,
        help="Clusters smaller than this count as outliers.",
    )
    parser.add_argument(
        "--max-k-limit",
        type=float,
        help="Cap k as a fraction of tumour cells; overrides max-k-fraction.",
    )
    parser.add_argument(
        "--peak-partition-preference",
        choices=("none", "fewer_outliers", "balanced"),
        default="none",
    )
    parser.add_argument("--peak-partition-weight", type=float, default=0.5)
    args = parser.parse_args()
    if args.top_n < 1 or not 0 < args.max_k_fraction <= 1:
        parser.error("top-n must be positive and max-k-fraction must be in (0, 1].")
    if args.outlier_size_threshold is not None and args.outlier_size_threshold < 1:
        parser.error("outlier-size-threshold must be positive.")
    if args.max_k_limit is not None and not 0 < args.max_k_limit <= 1:
        parser.error("max-k-limit must be in (0, 1].")
    normals = load_normal_cells()
    samples = discover_samples()
    if args.sample:
        unknown = set(args.sample) - {sample[0] for sample in samples}
        if unknown:
            parser.error(f"Unknown samples: {sorted(unknown)}")
        samples = [sample for sample in samples if sample[0] in args.sample]
    summaries = [
        process_sample(
            *sample,
            normals,
            output_root=args.output_dir,
            top_n=args.top_n,
            max_k_fraction=args.max_k_fraction,
            prefer_fewer_outliers=args.prefer_fewer_outliers,
            outlier_size_threshold=args.outlier_size_threshold,
            max_k_limit=args.max_k_limit,
            partition_preference=args.peak_partition_preference,
            partition_weight=args.peak_partition_weight,
        )
        for sample in samples
    ]
    pd.DataFrame(summaries).to_csv(
        args.output_dir / "summary.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
