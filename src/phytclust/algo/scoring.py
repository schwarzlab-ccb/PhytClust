import logging
from collections import Counter
from typing import Optional, Any

import numpy as np

from ..exceptions import (
    ConfigurationError,
    InvalidClusteringError,
    MissingDPTableError,
)
from ..viz.scores import plot_scores as _plot_scores
from .bins import define_bins as _define_bins
from ..config import RANKING_MODES, PeakConfig
from ..utils.traversal import nonterminals, terminals

logger = logging.getLogger("phytclust")


def _first_zero_length_pair_split(
    clustering, max_k: int, eps: float = 1e-12
) -> Optional[int]:
    """Find the first observed split of a two-leaf subtree on a zero-length edge.

    Check adjacent feasible cluster counts up to max_k or the diagnostic
    search limit, whichever is smaller. This does not detect every zero-length split."""
    if getattr(clustering, "no_split_zero_length", False):
        return None

    active_tree = (
        clustering._tree_wo_outgroup if clustering.outgroup else clustering.tree
    )
    zero_length_pairs: list[tuple] = []
    for node in nonterminals(active_tree.root):
        for child in node.clades:
            branch_length = child.branch_length or 0.0
            if branch_length <= eps:
                cached = getattr(clustering, "name_leaves_per_node", None)
                terms = cached.get(child) if cached is not None else None
                if terms is None:
                    terms = list(terminals(child))
                if len(terms) == 2:
                    zero_length_pairs.append((terms[0], terms[1]))

    if not zero_length_pairs:
        return None

    search_cap = int(getattr(clustering, "zero_length_split_max_k", 100))
    effective_max_k = min(max_k, search_cap)
    clustering._ensure_dp(required_cap=effective_max_k)
    try:
        prev = clustering._clusters(1)
    except (ValueError, RuntimeError, InvalidClusteringError, MissingDPTableError):
        prev = None

    for k in range(2, effective_max_k + 1):
        try:
            cur = clustering._clusters(k)
        except (ValueError, RuntimeError, InvalidClusteringError, MissingDPTableError):
            prev = None
            continue

        if prev is not None:
            for t1, t2 in zero_length_pairs:
                pc1 = prev.get(t1)
                pc2 = prev.get(t2)
                if pc1 is None or pc2 is None or pc1 != pc2:
                    continue
                cc1 = cur.get(t1)
                cc2 = cur.get(t2)
                if cc1 is None or cc2 is None:
                    continue
                if cc1 != cc2:
                    return k

        prev = cur

    return None


def _score_cluster_count(
    clustering, clusters: Optional[dict[Any, Any]] = None, k: Optional[int] = None
):
    """Return cost, relative cost reduction, and score for one cluster count.

    A supplied map selects its number of clusters; its membership is not scored."""
    if not clustering.max_k or clustering.max_k <= 0:
        raise ConfigurationError(
            "max_k must be set and positive to compute cluster scores."
        )

    active_tree = (
        clustering._tree_wo_outgroup if clustering.outgroup else clustering.tree
    )
    root = active_tree.root
    root_id = clustering.node_to_id[root]

    use_penalized_beta = getattr(clustering, "use_penalized_beta_for_scoring", False)
    dp_row = (
        clustering.dp_table[root_id]
        if use_penalized_beta
        else clustering.raw_dp_table[root_id]
    )

    if dp_row is None:
        raise MissingDPTableError("Root DP row missing.")

    dp_row = np.asarray(dp_row, dtype=float)
    clustering.beta_1 = dp_row[0]
    num_terminals = clustering.num_terminals

    if clusters is not None:
        num_clusters = len(set(clusters.values()))
        if num_clusters < 1 or num_clusters > clustering.max_k:
            return (float("inf"), float("inf"), float("inf"))
        if num_clusters - 1 >= len(dp_row):
            return (float("inf"), float("inf"), float("inf"))
        beta = dp_row[num_clusters - 1]

    elif k is not None:
        if k < 1 or k > clustering.max_k or k - 1 >= len(dp_row):
            return (float("inf"), float("inf"), float("inf"))
        num_clusters = k
        beta = dp_row[k - 1]

    else:
        raise ConfigurationError(
            "Either 'clusters' or 'k' must be provided to compute the score."
        )

    if not np.isfinite(beta):
        return (beta, float("inf"), 0.0)

    if beta == 0:
        return (beta, float("inf"), 0.0)

    relative_cost_floor = float(
        getattr(clustering, "score_beta_floor_frac", 0.0) or 0.0
    )
    absolute_cost_floor = float(getattr(clustering, "score_beta_floor_abs", 0.0) or 0.0)
    beta_floor = max(
        absolute_cost_floor, relative_cost_floor * float(clustering.beta_1)
    )
    beta_denom = max(float(beta), beta_floor)

    beta_ratios = (clustering.beta_1 - beta) / beta_denom
    norm_ratios = (num_terminals - num_clusters) / float(num_clusters)

    if not np.isfinite(beta_ratios) or not np.isfinite(norm_ratios):
        score = 0.0
    else:
        score = beta_ratios * norm_ratios

    return (beta, beta_ratios, score)


def _score_all_cluster_counts(clustering) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return costs, cost-reduction ratios, and scores for k = 1 through max_k.

    Zero or non-finite costs have zero score. Counts beyond the table have
    infinite cost and score."""
    active_tree = (
        clustering._tree_wo_outgroup if clustering.outgroup else clustering.tree
    )
    root_id = clustering.node_to_id[active_tree.root]

    use_penalized_cost = getattr(clustering, "use_penalized_beta_for_scoring", False)
    dp_row = (
        clustering.dp_table[root_id]
        if use_penalized_cost
        else clustering.raw_dp_table[root_id]
    )
    if dp_row is None:
        raise MissingDPTableError("Root DP row missing.")

    dp_row = np.asarray(dp_row, dtype=float)
    clustering.beta_1 = float(dp_row[0])
    num_terminals = clustering.num_terminals
    max_k = int(clustering.max_k)

    relative_cost_floor = float(
        getattr(clustering, "score_beta_floor_frac", 0.0) or 0.0
    )
    absolute_cost_floor = float(getattr(clustering, "score_beta_floor_abs", 0.0) or 0.0)
    beta_floor = max(absolute_cost_floor, relative_cost_floor * clustering.beta_1)

    available_cost_count = min(max_k, len(dp_row))
    betas = np.full(max_k, np.inf, dtype=float)
    betas[:available_cost_count] = dp_row[:available_cost_count]

    ks = np.arange(1, max_k + 1, dtype=float)
    norm_ratios = (num_terminals - ks) / ks

    edge = ~np.isfinite(betas) | (betas == 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        beta_denom = np.maximum(betas, beta_floor)
        beta_ratios_full = (clustering.beta_1 - betas) / beta_denom
        scores_full = beta_ratios_full * norm_ratios

    beta_ratios = np.where(edge, np.inf, beta_ratios_full)
    scores = np.where(
        edge | ~np.isfinite(scores_full),
        0.0,
        scores_full,
    )

    if available_cost_count < max_k:
        scores[available_cost_count:] = np.inf

    return betas, beta_ratios, scores


def _cached_cluster_count_scores(
    clustering,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Reuse score arrays for the current DP tables and cost-floor settings.

    A cached larger range can serve a smaller request by copying its prefix."""
    requested = int(clustering.max_k)
    base_sig = (
        clustering._dp_cache_sig,
        clustering._dp_cap,
        bool(getattr(clustering, "use_penalized_beta_for_scoring", False)),
        float(getattr(clustering, "score_beta_floor_frac", 0.0) or 0.0),
        float(getattr(clustering, "score_beta_floor_abs", 0.0) or 0.0),
    )

    cached = clustering._score_raw_arrays
    cached_sig = clustering._score_raw_base_sig
    cached_cap = clustering._score_raw_cap

    if (
        cached is not None
        and cached_sig == base_sig
        and cached_cap is not None
        and requested <= cached_cap
    ):
        betas, ratios, scores = cached
        return (
            betas[:requested].copy(),
            ratios[:requested].copy(),
            scores[:requested].copy(),
        )

    betas, ratios, scores = _score_all_cluster_counts(clustering)

    clustering._score_raw_arrays = (betas.copy(), ratios.copy(), scores.copy())
    clustering._score_raw_cap = requested
    clustering._score_raw_base_sig = base_sig
    return betas, ratios, scores


def calculate_scores(clustering) -> None:
    """Calculate elbow-weighted scores while preserving cluster-count positions."""
    if clustering.k is not None:
        from .dp import cluster_map

        cmap = cluster_map(clustering, clustering.k)
        beta, cost_reduction_ratio, score = _score_cluster_count(
            clustering, clusters=cmap
        )
        beta_values = np.array([beta], dtype=float)
        cost_reduction_ratios = np.array([cost_reduction_ratio], dtype=float)
        scores = np.array([score], dtype=float)
    else:
        if not clustering.max_k or clustering.max_k <= 0:
            raise ConfigurationError(
                "max_k must be set and positive to compute DP-based scores."
            )
        beta_values, cost_reduction_ratios, scores = _cached_cluster_count_scores(
            clustering
        )

    scores[scores < 0] = 0
    beta_values[beta_values < 0] = 0
    finite_costs = np.isfinite(beta_values)

    n = len(beta_values)
    elbow_scores = np.zeros(n, dtype=float)
    if n >= 3:
        drops = np.zeros(n - 1, dtype=float)
        adjacent_finite = finite_costs[:-1] & finite_costs[1:]
        drops[adjacent_finite] = np.maximum(
            0.0, beta_values[:-1][adjacent_finite] - beta_values[1:][adjacent_finite]
        )
        prev_drop = drops[:-1]
        next_drop = drops[1:]

        drop_eps = 1e-12 * (max(float(beta_values[0]), 0.0) if finite_costs[0] else 0.0)
        stable = (prev_drop > drop_eps) & (next_drop > drop_eps)
        ratios = np.zeros_like(prev_drop)
        ratios[stable] = prev_drop[stable] / next_drop[stable]

        ratios = np.clip(ratios, 0.0, 50.0)
        elbow_scores[1 : n - 1] = ratios

    invalid_mask = (
        np.isnan(scores)
        | np.isinf(scores)
        | np.isnan(elbow_scores)
        | np.isinf(elbow_scores)
    )
    # Keep array positions aligned with k; impossible partitions get zero score.
    scores[invalid_mask | ~finite_costs] = 0.0
    valid_mask = np.ones(n, dtype=bool)

    scores_valid = scores[valid_mask]
    beta_valid = beta_values[valid_mask]
    den_valid = cost_reduction_ratios[valid_mask]
    elbow_valid = elbow_scores[valid_mask]

    if len(scores_valid) == 0:
        clustering.scores = np.array([], dtype=float)
        clustering.beta_values = np.array([], dtype=float)
        clustering.norm_ratios = np.array([], dtype=float)
        return

    combined_scores = np.nan_to_num(
        elbow_valid * scores_valid, nan=0.0, posinf=0.0, neginf=0.0
    )

    eps = 1e-12
    nonzero_idx = np.where(np.abs(combined_scores) > eps)[0]
    if nonzero_idx.size > 0:
        last_useful = nonzero_idx[-1] + 1
    else:
        last_useful = len(combined_scores)

    combined_scores = combined_scores[:last_useful]
    beta_valid = beta_valid[:last_useful]
    den_valid = den_valid[:last_useful]

    clustering.scores = combined_scores
    clustering.beta_values = beta_valid
    clustering.norm_ratios = den_valid


def _rank_peaks(
    peak_data: list,
    ranking_mode: str,
    prominence_weight: float,
) -> list:
    """Rank peaks by prominence and score; returns list of dicts sorted best-first."""
    if ranking_mode not in RANKING_MODES:
        raise ConfigurationError(f"ranking_mode must be one of {RANKING_MODES}.")

    all_prom = [x[1] for x in peak_data]
    all_sc = [x[2] for x in peak_data]
    prom_min, prom_max = min(all_prom), max(all_prom)
    score_min, score_max = min(all_sc), max(all_sc)

    ranked_data = []
    for cluster_count, prominence, score in peak_data:
        if ranking_mode == "raw":
            prom_norm = score_norm = None
            combined_metric = prominence
        else:
            prom_norm = (
                (prominence - prom_min) / (prom_max - prom_min)
                if prom_max > prom_min
                else 1.0
            )
            score_norm = (
                (score - score_min) / (score_max - score_min)
                if score_max > score_min
                else 1.0
            )
            combined_metric = (
                prominence_weight * prom_norm + (1 - prominence_weight) * score_norm
            )

        ranked_data.append(
            {
                "k": cluster_count,
                "prominence": prominence,
                "score": score,
                "prom_norm": prom_norm,
                "score_norm": score_norm,
                "combined_metric": combined_metric,
            }
        )

    ranked_data.sort(key=lambda x: x["combined_metric"], reverse=True)
    return ranked_data


def _rank_peak_partitions(clustering, ranked_data, settings):
    """Blend peak strength with singleton avoidance or cluster-size balance."""
    if settings.partition_preference == "none" or settings.partition_weight == 0:
        return ranked_data
    maximum_strength = max(item["combined_metric"] for item in ranked_data)
    weight = settings.partition_weight
    for item in ranked_data:
        partition = clustering._clusters(int(item["k"]))
        sizes = list(Counter(partition.values()).values())
        cell_count = sum(sizes)
        singleton_count = sum(size == 1 for size in sizes)
        singleton_fraction = singleton_count / cell_count
        # Effective cluster count divided by k: one for equally sized clusters.
        balance = cell_count**2 / (len(sizes) * sum(size**2 for size in sizes))
        quality = 1 - singleton_fraction
        if settings.partition_preference == "balanced":
            quality *= balance
        strength = (
            item["combined_metric"] / maximum_strength if maximum_strength > 0 else 1.0
        )
        item.update(
            singleton_count=singleton_count,
            singleton_fraction=singleton_fraction,
            cluster_sizes=sorted(sizes, reverse=True),
            size_balance=balance,
            partition_quality=quality,
            ranking_metric=(1 - weight) * strength + weight * quality,
        )
    return sorted(ranked_data, key=lambda item: item["ranking_metric"], reverse=True)


def _plot_raw(
    peaks_to_plot: list,
    *,
    clustering,
    plot: bool,
    scores: np.ndarray,
    k_end: int,
    resolution_on: bool,
    num_bins: int,
) -> None:
    """Build and assign score plots onto clustering.plot_of_scores / clustering.plot_of_raw_scores."""
    if not plot:
        return

    scores_cfg = getattr(getattr(clustering, "plot_config", None), "scores", None)
    prefer_unsmoothed_primary = bool(
        getattr(scores_cfg, "prefer_unsmoothed_primary", True)
    )
    show_secondary_score_plot = bool(
        getattr(scores_cfg, "show_secondary_score_plot", False)
    )

    if resolution_on:
        primary_arr = scores[1:k_end].copy()
        clustering.plot_of_scores = _plot_scores(
            clustering,
            scores_subset=primary_arr,
            peaks=peaks_to_plot,
            k_start=2,
            k_end=k_end,
            resolution_on=True,
            num_bins=num_bins,
        )
        if show_secondary_score_plot:
            secondary_arr = scores[2:k_end].copy()
            clustering.plot_of_raw_scores = _plot_scores(
                clustering,
                scores_subset=secondary_arr,
                peaks=peaks_to_plot,
                k_start=3,
                k_end=k_end,
                resolution_on=False,
                num_bins=num_bins,
            )
        else:
            clustering.plot_of_raw_scores = None
        return

    if prefer_unsmoothed_primary:
        primary_arr = scores[2:k_end].copy()
        clustering.plot_of_scores = _plot_scores(
            clustering,
            scores_subset=primary_arr,
            peaks=peaks_to_plot,
            k_start=3,
            k_end=k_end,
            resolution_on=False,
            num_bins=num_bins,
        )
        if show_secondary_score_plot:
            secondary_arr = scores[1:k_end].copy()
            clustering.plot_of_raw_scores = _plot_scores(
                clustering,
                scores_subset=secondary_arr,
                peaks=peaks_to_plot,
                k_start=2,
                k_end=k_end,
                resolution_on=resolution_on,
                num_bins=num_bins,
            )
        else:
            clustering.plot_of_raw_scores = None
    else:
        primary_arr = scores[1:k_end].copy()
        clustering.plot_of_scores = _plot_scores(
            clustering,
            scores_subset=primary_arr,
            peaks=peaks_to_plot,
            k_start=2,
            k_end=k_end,
            resolution_on=resolution_on,
            num_bins=num_bins,
        )
        if show_secondary_score_plot:
            secondary_arr = scores[2:k_end].copy()
            clustering.plot_of_raw_scores = _plot_scores(
                clustering,
                scores_subset=secondary_arr,
                peaks=peaks_to_plot,
                k_start=3,
                k_end=k_end,
                resolution_on=False,
                num_bins=num_bins,
            )
        else:
            clustering.plot_of_raw_scores = None


def find_score_peaks(
    clustering,
    scores: Optional[np.ndarray] = None,
    global_peaks: int = 3,
    peaks_per_bin: int = 1,
    resolution_on: bool = False,
    num_bins: int = 3,
    k_start: Optional[int] = None,
    k_end: Optional[int] = None,
    plot: bool = True,
    peak_config: Optional[PeakConfig] = None,
) -> list[int]:
    """Find and rank score peaks within the inclusive k_start to k_end range."""
    from scipy.signal import find_peaks

    peak_settings = peak_config or PeakConfig()
    peak_settings.validate()
    clustering.peak_ranking_details = []

    min_k = peak_settings.min_k
    min_prominence = peak_settings.min_prominence
    use_log_peak_input = bool(getattr(peak_settings, "use_log_peak_input", False))
    log_peak_offset = getattr(peak_settings, "log_peak_offset", 1e-12)
    use_relative_prominence = bool(
        getattr(peak_settings, "use_relative_prominence", False)
    )
    min_relative_prominence = getattr(peak_settings, "min_relative_prominence", None)
    prominence_k_power = getattr(peak_settings, "prominence_k_power", 0.0)
    ranking_mode = peak_settings.ranking_mode
    prominence_weight = peak_settings.prominence_weight
    boundary_window_size = peak_settings.boundary_window_size
    boundary_ratio_threshold = peak_settings.boundary_ratio_threshold
    resolution_fallback_mode = peak_settings.resolution_fallback_mode
    exclude_k2 = bool(getattr(peak_settings, "exclude_k2", False))
    if exclude_k2:
        min_k = max(min_k, 3)

    if scores is None:
        scores = clustering.scores

    if scores is None or len(scores) == 0:
        clustering.peaks_by_rank = []
        clustering.resolution_info = None
        clustering.peaks_by_resolution = None

        if plot:
            clustering.plot_of_scores = _plot_scores(
                clustering,
                scores_subset=np.array([], dtype=float),
                peaks=[],
                k_start=1,
                resolution_on=resolution_on,
                num_bins=num_bins,
            )

        logger.info("No score peaks found.")
        return []

    scores = np.asarray(scores, dtype=float)
    scores = np.nan_to_num(scores, nan=0.0, posinf=0.0, neginf=0.0)

    zero_split_k = _first_zero_length_pair_split(clustering, len(scores))
    clustering.zero_length_split_k = zero_split_k
    if zero_split_k is not None:
        logger.warning(
            "At k=%d, leaves in a two-leaf subtree on a zero-length edge "
            "are assigned to different clusters. Check peaks at k >= %d.",
            zero_split_k,
            zero_split_k,
        )

    eps = 1e-12
    nonzero_idx = np.flatnonzero(np.abs(scores) > eps)

    k_start = k_start if k_start is not None else 1
    k_end = k_end if k_end is not None else len(scores)

    if len(scores) < 2:
        raise InvalidClusteringError(
            f"At least two scores are required to find peaks. scores = {scores}"
        )
    k_start = max(1, int(k_start))
    k_end = min(int(k_end), len(scores))
    min_k = max(min_k, k_start)
    allow_k2 = not exclude_k2 and min_k <= 2 <= k_end

    def _emit_plot(peaks_to_plot):
        _plot_raw(
            peaks_to_plot,
            clustering=clustering,
            plot=plot,
            scores=scores,
            k_end=k_end,
            resolution_on=resolution_on,
            num_bins=num_bins,
        )

    if allow_k2 and len(scores) > 1 and nonzero_idx.size == 1 and nonzero_idx[0] == 1:
        clustering.peaks_by_rank = [2]
        if not resolution_on:
            clustering.resolution_info = None
            clustering.peaks_by_resolution = None
        else:
            clustering.resolution_info = {
                "special_case": [(2, float(scores[1]), float(scores[1]), 1.0)]
            }
            clustering.peaks_by_resolution = {"special_case": [2]}
        _emit_plot([2])
        return clustering.peaks_by_rank

    peak_data = []

    if len(scores) > 2 and allow_k2:
        score_k2 = scores[1]
        w = min(boundary_window_size, len(scores) - 2)
        right_window = scores[2 : 2 + w]

        if len(right_window) > 0 and score_k2 > right_window[0]:
            window_baseline = float(np.median(right_window))
            window_ratio = score_k2 / (window_baseline + eps)

            if window_ratio > boundary_ratio_threshold:
                boundary_prom = score_k2 - window_baseline
                if use_relative_prominence:
                    boundary_prom = (score_k2 + eps) / (window_baseline + eps)

                if (
                    not use_relative_prominence
                    or min_relative_prominence is None
                    or boundary_prom >= float(min_relative_prominence)
                ):
                    peak_data.append((2, boundary_prom, score_k2))

    interior_start = 2
    interior_scores = scores[interior_start:k_end].astype(float)
    interior_scores = np.nan_to_num(interior_scores, nan=0.0, posinf=0.0, neginf=0.0)

    # A trailing zero supplies the descending side of the preceding peak.
    peak_input = interior_scores

    if use_log_peak_input:
        peak_input = np.log(np.maximum(peak_input, 0.0) + max(log_peak_offset, 1e-15))

    if len(peak_input) >= 2:
        score_range = float(np.max(peak_input) - np.min(peak_input))
        if min_prominence is None:
            auto_prom = 0.01 * score_range
        else:
            auto_prom = min_prominence

        sentinel = float(np.min(peak_input)) - 1.0
        padded_input = np.concatenate([[sentinel], peak_input])
        peaks_idx_padded, props = find_peaks(padded_input, prominence=auto_prom)
        peaks_idx = peaks_idx_padded - 1
        prominences = props["prominences"]

        for i, pidx in enumerate(peaks_idx):
            cluster_count = pidx + 3
            if cluster_count >= min_k:
                prom_raw = float(prominences[i])
                score_at_k = float(scores[cluster_count - 1])
                baseline = max(score_at_k - prom_raw, 0.0)
                prom_used = prom_raw

                if use_relative_prominence:
                    if use_log_peak_input:
                        # SciPy's bases are positions; read their original scores.
                        left_base = int(props["left_bases"][i]) - 1
                        right_base = int(props["right_bases"][i]) - 1
                        baseline = max(
                            float(interior_scores[left_base])
                            if left_base >= 0
                            else 0.0,
                            float(interior_scores[right_base])
                            if right_base >= 0
                            else 0.0,
                            0.0,
                        )
                    prom_used = (score_at_k + eps) / (baseline + eps)

                if use_relative_prominence:
                    if min_relative_prominence is None or prom_used >= float(
                        min_relative_prominence
                    ):
                        peak_data.append((cluster_count, prom_used, score_at_k))
                elif min_prominence is None or prominence_k_power <= 0:
                    peak_data.append((cluster_count, prom_used, score_at_k))
                else:
                    min_prom_k = float(min_prominence) * (
                        float(cluster_count) ** float(prominence_k_power)
                    )
                    if prom_used >= min_prom_k:
                        peak_data.append((cluster_count, prom_used, score_at_k))

    dedup = {}
    for cluster_count, prominence, score in peak_data:
        if cluster_count not in dedup or prominence > dedup[cluster_count][0]:
            dedup[cluster_count] = (prominence, score)
    peak_data = [
        (cluster_count, prominence, score)
        for cluster_count, (prominence, score) in dedup.items()
    ]

    if len(peak_data) == 0:
        if allow_k2 and len(scores) > 2 and scores[1] > scores[2]:
            clustering.peaks_by_rank = [2]
            if not resolution_on:
                clustering.resolution_info = None
                clustering.peaks_by_resolution = None
            else:
                clustering.resolution_info = {
                    "fallback_k2": [(2, float(scores[1]), float(scores[1]), 1.0)]
                }
                clustering.peaks_by_resolution = {"fallback_k2": [2]}
            _emit_plot([2])
            return clustering.peaks_by_rank

        logger.info("No score peaks found.")
        clustering.peaks_by_rank = []
        clustering.resolution_info = None
        clustering.peaks_by_resolution = None
        _emit_plot([])
        return clustering.peaks_by_rank

    if not resolution_on:
        clustering.resolution_info = None
        clustering.peaks_by_resolution = None

        ranked_data = _rank_peaks(peak_data, ranking_mode, prominence_weight)
        ranked_data = _rank_peak_partitions(clustering, ranked_data, peak_settings)

        chosen = ranked_data[:global_peaks]
        final_peaks = [int(x["k"]) for x in chosen]
        clustering.peaks_by_rank = final_peaks
        clustering.peak_ranking_details = ranked_data

    else:
        bin_ranges = _define_bins(clustering, num_bins=num_bins, k_lo=1, k_hi=k_end)
        clustering.bin_ranges_current = bin_ranges
        clustering.resolution_info = {}
        clustering.peaks_by_resolution = {}

        ranked_data = _rank_peaks(peak_data, ranking_mode, prominence_weight)
        ranked_data = _rank_peak_partitions(clustering, ranked_data, peak_settings)
        clustering.peak_ranking_details = ranked_data

        all_picked_peaks = []
        for i, (start_k, end_k) in enumerate(bin_ranges, start=1):
            bin_label = f"Bin {i}: {start_k}-{end_k}"
            candidates = [item for item in ranked_data if start_k <= item["k"] <= end_k]
            chosen = candidates[:peaks_per_bin]

            bin_k_lo = max(int(start_k), int(min_k), 2)
            bin_k_hi = min(int(end_k), len(scores), k_end)
            bin_best_k = None
            bin_best_score = None
            if bin_k_lo <= bin_k_hi:
                bin_ks = np.arange(bin_k_lo, bin_k_hi + 1, dtype=int)
                bin_vals = scores[bin_ks - 1]
                if bin_vals.size > 0:
                    best_idx = int(np.argmax(bin_vals))
                    bin_best_k = int(bin_ks[best_idx])
                    bin_best_score = float(bin_vals[best_idx])

            if len(chosen) == 0 and resolution_fallback_mode == "max_score":
                if bin_best_k is not None and bin_best_score is not None:
                    chosen = [
                        {
                            "k": bin_best_k,
                            "prominence": 0.0,
                            "score": bin_best_score,
                            "combined_metric": bin_best_score,
                            "selection_note": "fallback_max_score",
                        }
                    ]

            chosen_kvals = [x["k"] for x in chosen]
            clustering.resolution_info[bin_label] = chosen
            clustering.peaks_by_resolution[bin_label] = chosen_kvals
            all_picked_peaks.extend(chosen_kvals)

        final_peaks = sorted({int(cluster_count) for cluster_count in all_picked_peaks})
        clustering.peaks_by_rank = final_peaks

    _emit_plot(clustering.peaks_by_rank)
    return clustering.peaks_by_rank
