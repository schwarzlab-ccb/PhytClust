import logging
from typing import Optional, Any

import numpy as np

from ..exceptions import (
    ConfigurationError,
    InvalidClusteringError,
    MissingDPTableError,
)
from ..viz.scores import plot_scores as _plot_scores
from .bins import define_bins as _define_bins
from ..config import PeakConfig

logger = logging.getLogger("phytclust")


def _find_zero_length_split_k(pc, max_k: int, eps: float = 1e-12) -> Optional[int]:
    """
    Find the smallest k where the optimal partition splits terminals
    connected by a zero-length branch into different clusters.

    Returns that k value, or None if no zero-length split occurs up to max_k.

    Notes
    -----
    This is a diagnostic / warning helper only. It must not dominate runtime,
    so:

    1. If the user already set ``pc.no_split_zero_length=True``, the DP has
       already forbidden any such split; return ``None`` immediately.
    2. Otherwise walk k in order and return as soon as the first split is
       observed (early exit). The common case is a small k, so the loop
       almost never runs to ``max_k``.
    3. Cluster lookups go through ``pc.get_clusters(k)``, which hits the
       per-k backtrack cache on ``pc.clusters`` and cooperates with the
       cache invalidation in ``_ensure_dp``.
    """
    # 1. Short-circuit: if splits are structurally forbidden, there's nothing
    #    to find. This alone turns an O(max_k) loop into O(1) for configs
    #    where the user already opted into `no_split_zero_length=True`.
    if getattr(pc, "no_split_zero_length", False):
        return None

    # 2. Collect zero-length sibling pairs (clades with a zero-length branch
    #    whose subtree contains exactly two terminals).
    active_tree = pc._tree_wo_outgroup if pc.outgroup else pc.tree
    zero_length_pairs: list[tuple] = []
    for node in active_tree.get_nonterminals():
        for child in node.clades:
            bl = child.branch_length or 0.0
            if bl <= eps:
                # Prefer the pre-cached terminal list over a fresh DFS.
                cached = getattr(pc, "name_leaves_per_node", None)
                terms = cached.get(child) if cached is not None else None
                if terms is None:
                    terms = list(child.get_terminals())
                if len(terms) == 2:
                    zero_length_pairs.append((terms[0], terms[1]))

    if not zero_length_pairs:
        return None

    # 3. Walk k in order, using the cached backtrack. Early-exit on the first
    #    split — in practice this terminates at very small k (often k=2).
    #    Cap the search at `zero_length_split_max_k` (default 100). The peak
    #    detector's warning is most actionable at small k; scanning into the
    #    thousands provides little extra information and dominates runtime
    #    on large trees with many feasible k values.
    search_cap = int(getattr(pc, "zero_length_split_max_k", 100))
    effective_max_k = min(max_k, search_cap)
    try:
        prev = pc.get_clusters(1)
    except (ValueError, RuntimeError, InvalidClusteringError, MissingDPTableError):
        prev = None

    for k in range(2, effective_max_k + 1):
        try:
            cur = pc.get_clusters(k)
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


def _single_cluster_score(
    pc, clusters: Optional[dict[Any, Any]] = None, k: Optional[int] = None
):
    if not pc.max_k or pc.max_k <= 0:
        raise ConfigurationError(
            "max_k must be set and positive to compute cluster scores."
        )

    active_tree = pc._tree_wo_outgroup if pc.outgroup else pc.tree
    root = active_tree.root
    root_id = pc.node_to_id[root]

    use_penalized_beta = getattr(pc, "use_penalized_beta_for_scoring", False)
    dp_row = pc.dp_table[root_id] if use_penalized_beta else pc.raw_dp_table[root_id]

    if dp_row is None:
        raise MissingDPTableError("Root DP row missing.")

    dp_row = np.asarray(dp_row, dtype=float)
    pc.beta_1 = dp_row[0]
    num_terminals = pc.num_terminals

    if clusters is not None:
        num_clusters = len(set(clusters.values()))
        if num_clusters < 1 or num_clusters > pc.max_k:
            return (float("inf"), float("inf"), float("inf"))
        if num_clusters - 1 >= len(dp_row):
            return (float("inf"), float("inf"), float("inf"))
        beta = dp_row[num_clusters - 1]

    elif k is not None:
        if k < 1 or k > pc.max_k or k - 1 >= len(dp_row):
            return (float("inf"), float("inf"), float("inf"))
        num_clusters = k
        beta = dp_row[k - 1]

    else:
        raise ConfigurationError(
            "Either 'clusters' or 'k' must be provided to compute the score."
        )

    if np.isinf(beta):
        return (beta, float("inf"), 0.0)

    if beta == 0:
        return (beta, float("inf"), 0.0)

    # Optional stabilization: floor denominator to avoid late-k explosion
    # when beta approaches 0 on large trees.
    beta_floor_frac = float(getattr(pc, "score_beta_floor_frac", 0.0) or 0.0)
    beta_floor_abs = float(getattr(pc, "score_beta_floor_abs", 0.0) or 0.0)
    beta_floor = max(beta_floor_abs, beta_floor_frac * float(pc.beta_1))
    beta_denom = max(float(beta), beta_floor)

    beta_ratios = (pc.beta_1 - beta) / beta_denom
    norm_ratios = (
        (num_terminals - num_clusters) / float(num_clusters)
        if num_clusters
        else float("inf")
    )

    if not np.isfinite(beta_ratios) or not np.isfinite(norm_ratios):
        score = float("inf")
    else:
        score = beta_ratios * norm_ratios

    return (beta, beta_ratios, score)


def _vectorised_dp_row_scores(pc) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Vectorised equivalent of looping `_single_cluster_score(pc, k=k)`
    over k = 1..max_k. Returns (betas, beta_ratios, scores) of length
    max_k. Same edge-case semantics as the per-k version:

      - k beyond DP row length → (inf, inf, inf)
      - beta non-finite or zero → (beta, inf, 0)
      - otherwise               → (beta, br, br * norm)
    """
    active_tree = pc._tree_wo_outgroup if pc.outgroup else pc.tree
    root_id = pc.node_to_id[active_tree.root]

    use_pen = getattr(pc, "use_penalized_beta_for_scoring", False)
    dp_row = pc.dp_table[root_id] if use_pen else pc.raw_dp_table[root_id]
    if dp_row is None:
        raise MissingDPTableError("Root DP row missing.")

    dp_row = np.asarray(dp_row, dtype=float)
    pc.beta_1 = float(dp_row[0])
    num_terminals = pc.num_terminals
    max_k = int(pc.max_k)

    beta_floor_frac = float(getattr(pc, "score_beta_floor_frac", 0.0) or 0.0)
    beta_floor_abs = float(getattr(pc, "score_beta_floor_abs", 0.0) or 0.0)
    beta_floor = max(beta_floor_abs, beta_floor_frac * pc.beta_1)

    n_in = min(max_k, len(dp_row))
    betas = np.full(max_k, np.inf, dtype=float)
    betas[:n_in] = dp_row[:n_in]

    ks = np.arange(1, max_k + 1, dtype=float)
    norm_ratios = (num_terminals - ks) / ks

    edge = ~np.isfinite(betas) | (betas == 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        beta_denom = np.maximum(betas, beta_floor)
        beta_ratios_full = (pc.beta_1 - betas) / beta_denom
        scores_full = beta_ratios_full * norm_ratios

    beta_ratios = np.where(edge, np.inf, beta_ratios_full)
    scores = np.where(
        edge | ~np.isfinite(scores_full),
        0.0,
        scores_full,
    )

    # k beyond dp_row length: original returned (inf, inf, inf).
    # Edge handling above set scores=0 for those slots (since betas=inf).
    # Restore inf so the downstream invalid_mask filters them out.
    if n_in < max_k:
        scores[n_in:] = np.inf

    return betas, beta_ratios, scores


def calculate_scores(pc, plot: bool = False) -> None:
    if pc.k is not None:
        from .dp import cluster_map

        cmap = cluster_map(pc, pc.k)
        beta, br, sc = _single_cluster_score(pc, clusters=cmap)
        beta_values = np.array([beta], dtype=float)
        den_list = np.array([br], dtype=float)
        scores = np.array([sc], dtype=float)
    else:
        if not pc.max_k or pc.max_k <= 0:
            raise ConfigurationError(
                "max_k must be set and positive to compute DP-based scores."
            )
        beta_values, den_list, scores = _vectorised_dp_row_scores(pc)

    scores[scores < 0] = 0
    beta_values[beta_values < 0] = 0
    beta_values = np.nan_to_num(beta_values, nan=0.0, posinf=0.0, neginf=0.0)

    n = len(beta_values)
    elbow_scores = np.zeros(n, dtype=float)
    if n >= 3:
        # Use drop magnitudes (not signed diffs) and suppress unstable tail ratios
        # where beta improvements become numerically tiny.
        drops = np.maximum(0.0, -np.diff(beta_values))
        prev_drop = drops[:-1]
        next_drop = drops[1:]

        drop_eps = 1e-10
        stable = (prev_drop > drop_eps) & (next_drop > drop_eps)
        ratios = np.zeros_like(prev_drop)
        ratios[stable] = prev_drop[stable] / next_drop[stable]

        # Prevent tiny-denominator explosions from dominating peak ranking.
        ratios = np.clip(ratios, 0.0, 50.0)
        elbow_scores[1 : n - 1] = ratios

    invalid_mask = (
        np.isnan(scores)
        | np.isinf(scores)
        | np.isnan(elbow_scores)
        | np.isinf(elbow_scores)
    )
    valid_mask = ~invalid_mask

    scores_valid = scores[valid_mask]
    beta_valid = beta_values[valid_mask]
    den_valid = den_list[valid_mask]
    elbow_valid = elbow_scores[valid_mask]

    if len(scores_valid) == 0:
        pc.scores = np.array([], dtype=float)
        pc.beta_values = np.array([], dtype=float)
        pc.norm_ratios = np.array([], dtype=float)
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

    pc.scores = combined_scores
    pc.beta_values = beta_valid
    pc.norm_ratios = den_valid


def _rank_peaks(
    peak_data: list,
    ranking_mode: str,
    lambda_weight: float,
) -> list:
    """Rank peaks by prominence and score; returns list of dicts sorted best-first."""
    if ranking_mode not in ("raw", "adjusted"):
        raise ConfigurationError("ranking_mode must be either 'raw' or 'adjusted'.")

    all_prom = [x[1] for x in peak_data]
    all_sc = [x[2] for x in peak_data]
    prom_min, prom_max = min(all_prom), max(all_prom)
    score_min, score_max = min(all_sc), max(all_sc)

    ranked_data = []
    for pk, prom, sc in peak_data:
        if ranking_mode == "raw":
            prom_norm = prom
            score_norm = sc
            base_metric = prom
            combined_metric = base_metric
        else:
            prom_norm = (
                (prom - prom_min) / (prom_max - prom_min)
                if prom_max > prom_min
                else 1.0
            )
            score_norm = (
                (sc - score_min) / (score_max - score_min)
                if score_max > score_min
                else 1.0
            )
            base_metric = lambda_weight * prom_norm + (1 - lambda_weight) * score_norm
            combined_metric = base_metric

        ranked_data.append(
            {
                "k": pk,
                "prominence": prom,
                "score": sc,
                "prom_norm": prom_norm if ranking_mode == "adjusted" else None,
                "score_norm": score_norm if ranking_mode == "adjusted" else None,
                "base_metric": base_metric,
                "combined_metric": combined_metric,
            }
        )

    ranked_data.sort(key=lambda x: x["combined_metric"], reverse=True)
    return ranked_data


def _plot_raw(
    peaks_to_plot: list,
    *,
    pc,
    plot: bool,
    scores: np.ndarray,
    k_end: int,
    resolution_on: bool,
    num_bins: int,
) -> None:
    """Build and assign score plots onto pc.plot_of_scores / pc.plot_of_raw_scores."""
    if not plot:
        return

    scores_cfg = getattr(getattr(pc, "plot_config", None), "scores", None)
    prefer_unsmoothed_primary = bool(
        getattr(scores_cfg, "prefer_unsmoothed_primary", True)
    )
    show_secondary_score_plot = bool(
        getattr(scores_cfg, "show_secondary_score_plot", False)
    )

    if resolution_on:
        primary_arr = scores[1:k_end].copy()
        pc.plot_of_scores = _plot_scores(
            pc,
            scores_subset=primary_arr,
            peaks=peaks_to_plot,
            k_start=2,
            k_end=k_end,
            resolution_on=True,
            num_bins=num_bins,
        )
        if show_secondary_score_plot:
            secondary_arr = scores[2:k_end].copy()
            pc.plot_of_raw_scores = _plot_scores(
                pc,
                scores_subset=secondary_arr,
                peaks=peaks_to_plot,
                k_start=3,
                k_end=k_end,
                resolution_on=False,
                num_bins=num_bins,
            )
        else:
            pc.plot_of_raw_scores = None
        return

    if prefer_unsmoothed_primary:
        primary_arr = scores[2:k_end].copy()
        pc.plot_of_scores = _plot_scores(
            pc,
            scores_subset=primary_arr,
            peaks=peaks_to_plot,
            k_start=3,
            k_end=k_end,
            resolution_on=False,
            num_bins=num_bins,
        )
        if show_secondary_score_plot:
            secondary_arr = scores[1:k_end].copy()
            pc.plot_of_raw_scores = _plot_scores(
                pc,
                scores_subset=secondary_arr,
                peaks=peaks_to_plot,
                k_start=2,
                k_end=k_end,
                resolution_on=resolution_on,
                num_bins=num_bins,
            )
        else:
            pc.plot_of_raw_scores = None
    else:
        primary_arr = scores[1:k_end].copy()
        pc.plot_of_scores = _plot_scores(
            pc,
            scores_subset=primary_arr,
            peaks=peaks_to_plot,
            k_start=2,
            k_end=k_end,
            resolution_on=resolution_on,
            num_bins=num_bins,
        )
        if show_secondary_score_plot:
            secondary_arr = scores[2:k_end].copy()
            pc.plot_of_raw_scores = _plot_scores(
                pc,
                scores_subset=secondary_arr,
                peaks=peaks_to_plot,
                k_start=3,
                k_end=k_end,
                resolution_on=False,
                num_bins=num_bins,
            )
        else:
            pc.plot_of_raw_scores = None


def find_score_peaks(
    pc,
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
    from scipy.signal import find_peaks

    cfg = peak_config or PeakConfig()

    # Unpack config
    min_k = cfg.min_k
    min_prominence = cfg.min_prominence
    use_log_peak_input = bool(getattr(cfg, "use_log_peak_input", False))
    log_peak_offset = float(getattr(cfg, "log_peak_offset", 1e-12))
    use_relative_prominence = bool(getattr(cfg, "use_relative_prominence", False))
    min_relative_prominence = getattr(cfg, "min_relative_prominence", None)
    prominence_k_power = getattr(cfg, "prominence_k_power", 0.0)
    ranking_mode = cfg.ranking_mode
    lambda_weight = cfg.lambda_weight
    boundary_window_size = cfg.boundary_window_size
    boundary_ratio_threshold = cfg.boundary_ratio_threshold
    resolution_fallback_mode = cfg.resolution_fallback_mode

    if scores is None:
        scores = pc.scores

    if scores is None or len(scores) == 0:
        pc.peaks_by_rank = []
        pc.resolution_info = None
        pc.peaks_by_resolution = None

        if plot:
            pc.plot_of_scores = _plot_scores(
                pc,
                scores_subset=np.array([], dtype=float),
                peaks=[],
                k_start=1,
                resolution_on=resolution_on,
                num_bins=num_bins,
            )

        logger.info("PhytClust did not find any clusters.")
        return []

    scores = np.asarray(scores, dtype=float)
    scores = np.nan_to_num(scores, nan=0.0, posinf=0.0, neginf=0.0)

    # Detect the first k that splits a zero-length branch
    zero_split_k = _find_zero_length_split_k(pc, len(scores))
    pc.zero_length_split_k = zero_split_k
    if zero_split_k is not None:
        logger.warning(
            "At k=%d, the partition first splits terminals connected by a "
            "zero-length branch. Peaks at k>=%d may not be meaningful.",
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
    # Defensive normalization for callers that pass k values instead of indices,
    # or when the feasible score vector is very short under strict constraints.
    k_start = int(max(0, min(k_start, len(scores) - 2)))
    k_end = int(max(k_start + 1, min(k_end, len(scores))))

    # Convenience: bind shared keyword args for _plot_raw calls in this function.
    def _emit_plot(peaks_to_plot):
        _plot_raw(
            peaks_to_plot,
            pc=pc,
            plot=plot,
            scores=scores,
            k_end=k_end,
            resolution_on=resolution_on,
            num_bins=num_bins,
        )

    # Special case: only meaningful score is at k=2
    if len(scores) > 1 and nonzero_idx.size == 1 and nonzero_idx[0] == 1:
        pc.peaks_by_rank = [2]
        if not resolution_on:
            pc.resolution_info = None
            pc.peaks_by_resolution = None
        else:
            pc.resolution_info = {
                "special_case": [(2, float(scores[1]), float(scores[1]), 1.0)]
            }
            pc.peaks_by_resolution = {"special_case": [2]}
        _emit_plot([2])
        return pc.peaks_by_rank

    # ----------------------------------------------------------------
    # 1) Boundary candidate: k=2 via right-window median comparison
    # ----------------------------------------------------------------
    peak_data = []  # list of (k, prominence, score)

    if len(scores) > 2 and min_k <= 2:
        score_k2 = scores[1]  # index 1 = k=2
        w = min(boundary_window_size, len(scores) - 2)  # don't exceed array
        right_window = scores[2 : 2 + w]  # k=3 onward

        if len(right_window) > 0 and score_k2 > right_window[0]:
            window_baseline = float(np.median(right_window))
            window_ratio = score_k2 / (window_baseline + eps)

            if window_ratio > boundary_ratio_threshold:
                # Use the excess over baseline as a prominence-like measure
                boundary_prom = score_k2 - window_baseline
                if use_relative_prominence:
                    boundary_prom = (score_k2 + eps) / (window_baseline + eps)

                if min_relative_prominence is None or boundary_prom >= float(
                    min_relative_prominence
                ):
                    peak_data.append((2, boundary_prom, score_k2))

    # ----------------------------------------------------------------
    # 2) Interior peaks: find_peaks on k>=3 subarray (index 2 onward)
    # ----------------------------------------------------------------
    interior_start = 2  # array index for k=3
    interior_scores = scores[interior_start:k_end].astype(float)
    interior_scores = np.nan_to_num(interior_scores, nan=0.0, posinf=0.0, neginf=0.0)

    # Trim trailing zeros
    if len(interior_scores) > 1 and interior_scores[-1] <= 0:
        peak_input = interior_scores[:-1]
    else:
        peak_input = interior_scores

    if use_log_peak_input:
        # log-domain peak search suppresses multiplicative tail inflation.
        # Scores are non-negative; clamp defensively for numerical safety.
        peak_input = np.log(np.maximum(peak_input, 0.0) + max(log_peak_offset, 1e-15))

    if len(peak_input) >= 2:
        score_range = float(np.max(peak_input) - np.min(peak_input))
        if min_prominence is None:
            auto_prom = 0.01 * score_range
        else:
            auto_prom = min_prominence

        # Prepend a sentinel strictly below peak_input so k=3 (original
        # index 0) can be detected by scipy.signal.find_peaks, which
        # otherwise excludes the first sample as a boundary.
        sentinel = float(np.min(peak_input)) - 1.0
        padded_input = np.concatenate([[sentinel], peak_input])
        peaks_idx_padded, props = find_peaks(padded_input, prominence=auto_prom)
        peaks_idx = peaks_idx_padded - 1
        prominences = props["prominences"]

        for i, pidx in enumerate(peaks_idx):
            pk = pidx + 3  # pidx=0 corresponds to k=3
            if pk >= min_k:
                prom_raw = float(prominences[i])
                score_at_k = float(scores[pk - 1])
                baseline = max(score_at_k - prom_raw, 0.0)
                prom_used = prom_raw

                if use_relative_prominence:
                    # Fold-change over local baseline: log-like behavior while
                    # keeping detection in linear-score space.
                    prom_used = (score_at_k + eps) / (baseline + eps)

                if use_relative_prominence:
                    if min_relative_prominence is None or prom_used >= float(
                        min_relative_prominence
                    ):
                        peak_data.append((pk, prom_used, score_at_k))
                elif min_prominence is None or prominence_k_power <= 0:
                    peak_data.append((pk, prom_used, score_at_k))
                else:
                    min_prom_k = float(min_prominence) * (
                        float(pk) ** float(prominence_k_power)
                    )
                    if prom_used >= min_prom_k:
                        peak_data.append((pk, prom_used, score_at_k))

    # Deduplicate (keep highest prominence per k)
    dedup = {}
    for pk, prom, sc in peak_data:
        if pk not in dedup or prom > dedup[pk][0]:
            dedup[pk] = (prom, sc)
    peak_data = [(pk, prom, sc) for pk, (prom, sc) in dedup.items()]

    if len(peak_data) == 0:
        # No peaks found at all — fallback to k=2 if it beats k=3
        if len(scores) > 2 and scores[1] > scores[2]:
            pc.peaks_by_rank = [2]
            if not resolution_on:
                pc.resolution_info = None
                pc.peaks_by_resolution = None
            else:
                pc.resolution_info = {
                    "fallback_k2": [(2, float(scores[1]), float(scores[1]), 1.0)]
                }
                pc.peaks_by_resolution = {"fallback_k2": [2]}
            _emit_plot([2])
            return pc.peaks_by_rank

        logger.info("PhytClust did not find any clusters.")
        pc.peaks_by_rank = []
        pc.resolution_info = None
        pc.peaks_by_resolution = None
        _emit_plot([])
        return pc.peaks_by_rank

    if not resolution_on:
        pc.resolution_info = None
        pc.peaks_by_resolution = None

        ranked_data = _rank_peaks(peak_data, ranking_mode, lambda_weight)

        chosen = ranked_data[:global_peaks]
        final_peaks = [int(x["k"]) for x in chosen]
        pc.peaks_by_rank = final_peaks
        pc.peak_ranking_details = ranked_data

    else:
        bin_ranges = _define_bins(pc, num_bins=num_bins, k_lo=min_k, k_hi=k_end)
        pc.bin_ranges_current = bin_ranges
        pc.resolution_info = {}
        pc.peaks_by_resolution = {}

        ranked_data = _rank_peaks(peak_data, ranking_mode, lambda_weight)
        pc.peak_ranking_details = ranked_data

        all_picked_peaks = []
        for i, (start_k, end_k) in enumerate(bin_ranges, start=1):
            bin_label = f"Bin {i}: {start_k}-{end_k}"
            candidates = [item for item in ranked_data if start_k <= item["k"] <= end_k]
            chosen = candidates[:peaks_per_bin]

            # Bin diagnostics: why empty bins can happen.
            bin_k_lo = max(int(start_k), int(min_k), 2)
            bin_k_hi = min(int(end_k), len(scores))
            bin_best_k = None
            bin_best_score = None
            if bin_k_lo <= bin_k_hi:
                bin_ks = np.arange(bin_k_lo, bin_k_hi + 1, dtype=int)
                bin_vals = scores[bin_ks - 1]
                if bin_vals.size > 0:
                    best_idx = int(np.argmax(bin_vals))
                    bin_best_k = int(bin_ks[best_idx])
                    bin_best_score = float(bin_vals[best_idx])

            used_fallback = False
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
                    used_fallback = True

            chosen_kvals = [x["k"] for x in chosen]
            pc.resolution_info[bin_label] = chosen
            pc.peaks_by_resolution[bin_label] = chosen_kvals
            all_picked_peaks.extend(chosen_kvals)

        final_peaks = sorted({int(pk) for pk in all_picked_peaks})
        pc.peaks_by_rank = final_peaks

    _emit_plot(pc.peaks_by_rank)
    return pc.peaks_by_rank
