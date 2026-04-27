import logging
import warnings
from dataclasses import dataclass, field
from math import ceil
from typing import Any, Optional

import numpy as np
from pathlib import Path
from io import StringIO

from Bio import Phylo

from ..exceptions import (
    InvalidKError,
    ConfigurationError,
    MissingDPTableError,
    ValidationError,
)
from Bio.Phylo.BaseTree import Tree

from ..algo.dp import (
    validate_args,
    prepare_tree,
    compute_dp_table,
    backtrack,
)
from ..algo.scoring import calculate_scores, find_score_peaks
from ..config import OutlierConfig, PeakConfig, RuntimeConfig
from ..metrics.indices import cluster_alpha

logger = logging.getLogger("phytclust")

IntMap = dict[Any, int]


def _coerce_to_tree(obj: Any) -> Tree:
    """
    Accepts:
      - Bio.Phylo.BaseTree.Tree  -> returned as-is
      - pathlib.Path             -> read as Newick
      - str:
          * if it looks like a file path and exists -> read as Newick file
          * otherwise -> treat as a Newick string

    Raises ValidationError if the object cannot be interpreted as a tree.
    """
    if isinstance(obj, Tree):
        return obj

    if isinstance(obj, Path):
        return Phylo.read(str(obj), "newick")

    if isinstance(obj, str):
        try:
            candidate = Path(obj)
            if candidate.exists():
                return Phylo.read(str(candidate), "newick")
        except OSError:
            # String too long to be a valid path (Errno 36) — treat as Newick
            pass

        handle = StringIO(obj)
        return Phylo.read(handle, "newick")

    raise ValidationError(
        f"Unsupported tree input type: {type(obj)!r}. "
        "Expected a Bio.Phylo Tree, a Newick string, or a path to a Newick file."
    )


@dataclass
class PhytClust:
    """
    Dynamic-programming phylogenetic clustering.

    Parameters
    ----------
    tree : Bio.Phylo.BaseTree.Tree
        The input tree (you can pass a Newick string/file upstream).
    outgroup : str | None, default=None
        Taxon to exclude from all clusters (treated as outgroup).
    min_cluster_size : int, default=1
        Hard constraint: final clusters smaller than this are disallowed.
    k : int | None, default=None
        Fixed number of clusters (only used when you call `run(k=...)`
        or `get_clusters(k)` explicitly).
    max_k : int | None, default=None
        Upper bound on k for scoring / peak search. If None, derived
        from `max_k_limit * num_terminals`.
    max_k_limit : float, default=0.9
        When `max_k` is not set, `max_k` = ceil(max_k_limit * num_terminals).
    num_bins : int, default=3
        Number of log-resolution bins for best_by_resolution.

    Support / branch-length tuning
    ------------------------------
    use_branch_support : bool, default=False
        If True, branch supports are incorporated into effective branch
        lengths and internal split penalties.
    min_support : float, default=0.05
        Minimal support used when normalizing (avoid division by 0).
    support_weight : float, default=1.0
        Weight of the support-derived penalty in `_eff_length`.

    Outlier handling
    ----------------
    outlier : OutlierConfig
        Controls outlier detection thresholds, DP penalty, and tie-breaking.
        See :class:`~phytclust.config.OutlierConfig`.

    Other flags
    -----------
    optimize_polytomies : bool, default=True
        If True, use native DP over multifurcations.
    polytomy_mode : {"hard", "soft"}, default="hard"
        Hard mode forbids cross-child partial merges; soft mode allows them.
    soft_polytomy_max_degree : int, default=18
        Degree guardrail for soft mode's exponential subset DP.
    preserve_dp_tables : bool, default=False
        Keep intermediate child DP rows instead of freeing them after merge.
        Useful for debugging/inspection; increases memory usage.
    compute_all_clusters : bool, default=False
        If True in best_global, compute and cache all k clusterings up to max_k.
    """

    tree: Any
    outgroup: Optional[str] = None
    root_taxon: Optional[str] = None
    min_cluster_size: int = 1
    k: Optional[int] = None
    max_k: Optional[int] = None
    max_k_limit: float = 0.9
    num_bins: int = 3

    # tunables
    use_branch_support: bool = False
    min_support: float = 0.05
    support_weight: float = 1.0

    # outlier handling
    outlier: OutlierConfig = field(default_factory=OutlierConfig)
    use_penalized_beta_for_scoring: bool = False

    # zero-length edge handling
    no_split_zero_length: bool = False
    optimize_polytomies: bool = True
    polytomy_mode: str = "hard"
    soft_polytomy_max_degree: int = 18
    preserve_dp_tables: bool = False

    compute_all_clusters: bool = False
    runtime_config: RuntimeConfig = field(default_factory=RuntimeConfig)
    peak_config: PeakConfig = field(default_factory=PeakConfig)

    def __post_init__(self) -> None:
        self.tree = _coerce_to_tree(self.tree)

        self.name_leaves_per_node = {}
        self.num_leaves_per_node = {}
        self.backptr = {}
        self.dp_table = None
        self.postorder_nodes = None
        self.node_to_id = None
        self.num_terminals = 0
        self._tree_wo_outgroup = None

        self.scores = None
        self.peaks_by_rank = None
        self.alpha_by_k: dict[int, dict[str, Any]] = {}

        self._dp_ready = False
        # Full cache key for the DP table: tree hash + every parameter that
        # changes the contents of dp_table/raw_dp_table/backptr.
        self._dp_cache_sig: Optional[tuple] = None
        # Cache key for pc.scores: depends on DP + max_k + scoring flags.
        self._scores_cache_sig: Optional[tuple] = None
        self.clusters: dict[int, IntMap] = {}
        self._last_result: Optional[dict[str, Any]] = None

        prepare_tree(self)

    def _hash_tree(self) -> int:
        """Tree fingerprint, used to detect modifications."""
        try:
            return hash(self.tree.format("newick"))
        except Exception:
            return hash(repr(self.tree))

    def _dp_signature(self) -> tuple:
        """
        Full cache key for the DP table.

        Anything that can change the contents of dp_table / raw_dp_table /
        backptr must be included here, otherwise ``_ensure_dp`` will reuse
        stale results after a parameter change on the same instance.
        """
        return (
            self._hash_tree(),
            self.outgroup,
            self.min_cluster_size,
            self.no_split_zero_length,
            getattr(self, "zero_length_eps", 1e-12),
            self.use_branch_support,
            self.min_support,
            self.support_weight,
            self.optimize_polytomies,
            self.polytomy_mode,
            self.soft_polytomy_max_degree,
            self.preserve_dp_tables,
            self.outlier.size_threshold,
            self.outlier.prefer_fewer,
        )

    def _scores_signature(self) -> tuple:
        """Cache key for ``pc.scores`` (depends on DP sig + scoring knobs)."""
        return (
            self._dp_cache_sig,
            self.max_k,
            bool(getattr(self, "use_penalized_beta_for_scoring", False)),
        )

    def _ensure_dp(self) -> None:
        current = self._dp_signature()

        if self._dp_ready and (self._dp_cache_sig == current):
            logger.debug("DP exists, not recalculating")
            return

        # DP is being rebuilt; every downstream cache must be invalidated.
        self.clusters = {}
        self.scores = None
        self.peaks_by_rank = None
        self.alpha_by_k = {}
        self._scores_cache_sig = None

        validate_args(self)
        compute_dp_table(self)

        self._dp_ready = True
        self._dp_cache_sig = current

        if self.max_k is None or self.max_k < 1:
            self.max_k = max(2, ceil(self.num_terminals * self.max_k_limit))

    def _effective_max_k(self, max_k: Optional[int] = None) -> int:
        """Resolve max_k without mutating self.

        Precedence is:
        1. explicit method argument ``max_k``
        2. instance attribute ``self.max_k`` if already set
        3. derived cap from ``max_k_limit``
        """
        if max_k is not None:
            return min(self.num_terminals, max_k)
        if self.max_k is not None:
            return min(self.num_terminals, self.max_k)
        return max(2, ceil(self.num_terminals * self.max_k_limit))

    def _build_run_result(
        self,
        *,
        mode: str,
        clusters: list[IntMap],
        selected_ks: list[int],
        exact_k: Optional[int] = None,
        include_scores: bool = True,
    ) -> dict[str, Any]:
        """Build a stable run() result payload.

        Canonical keys are ``k_values`` (list[int]) and ``selected_k`` (int|None).
        Legacy keys (``ks``, ``peaks``, ``k``) are retained for compatibility.
        """
        selected_k = (
            int(exact_k)
            if exact_k is not None
            else (int(selected_ks[0]) if selected_ks else None)
        )
        alpha_info = self._compute_and_log_alphas(selected_ks, clusters)
        result: dict[str, Any] = {
            "mode": mode,
            "k_values": list(selected_ks),
            "selected_k": selected_k,
            "ks": list(selected_ks),
            "clusters": clusters,
            "scores": (
                None
                if (not include_scores or self.scores is None)
                else self.scores.copy()
            ),
            "peaks": list(selected_ks),
            "alphas": [info["alpha"] for info in alpha_info],
            "alpha_details": alpha_info,
        }
        if exact_k is not None:
            result["k"] = int(exact_k)
        return result

    def _compute_and_log_alphas(
        self,
        selected_ks: list[int],
        clusters: list[IntMap],
    ) -> list[dict[str, Any]]:
        """Compute alpha per selected k, cache on ``self.alpha_by_k``, and print it.

        Alpha = mean extra-cluster branch length / mean intra-cluster branch length.
        """
        active_tree = (
            self._tree_wo_outgroup
            if (self.outgroup and self._tree_wo_outgroup is not None)
            else self.tree
        )
        details: list[dict[str, Any]] = []
        for k_val, cmap in zip(selected_ks, clusters):
            info = cluster_alpha(active_tree, cmap)
            info_with_k = {"k": int(k_val), **info}
            self.alpha_by_k[int(k_val)] = info_with_k
            details.append(info_with_k)
            logger.info(
                "alpha(k=%d) = %.6g  (avg_extra=%.6g, avg_intra=%.6g, "
                "n_extra=%d, n_intra=%d)",
                k_val,
                info["alpha"],
                info["avg_extra_branch_length"],
                info["avg_intra_branch_length"],
                info["n_extra_nodes"],
                info["n_intra_nodes"],
            )
            print(
                f"[phytclust] alpha(k={k_val}) = {info['alpha']:.6g}  "
                f"(avg_extra={info['avg_extra_branch_length']:.6g}, "
                f"avg_intra={info['avg_intra_branch_length']:.6g}, "
                f"n_extra={info['n_extra_nodes']}, "
                f"n_intra={info['n_intra_nodes']})"
            )
        return details

    @property
    def plot_config(self):
        """Convenience accessor for plot-specific runtime config."""
        return self.runtime_config.plot

    # ------------------------------------------------------------------ #
    #  Single code path for retrieving / computing a k-partition          #
    # ------------------------------------------------------------------ #

    def get_clusters(self, k: int, *, verbose: bool = False) -> IntMap:
        """Return the exact k-cluster partition (cached after first call)."""
        if k is None:
            raise InvalidKError("Please provide k")
        if k < 1:
            raise InvalidKError("k must be >= 1")
        self._ensure_dp()

        if k in self.clusters:
            return self.clusters[k]

        cmap = backtrack(self, k, verbose=verbose)
        self.clusters[k] = cmap
        return cmap

    # ------------------------------------------------------------------ #
    #  Peak-search modes                                                  #
    # ------------------------------------------------------------------ #

    def _no_peaks_fallback(self) -> list[IntMap]:
        """When scores are empty or too short, fall back gracefully."""
        logger.info("No meaningful peaks found.")
        self.k = None
        self.peaks_by_rank = []
        return []

    def _run_peak_mode(
        self,
        *,
        resolution_on: bool,
        top_n: int = 1,
        num_bins: int = 3,
        max_k: Optional[int] = None,
        plot_scores: bool = True,
        compute_all_clusters: bool = False,
        peak_config: Optional[PeakConfig] = None,
    ) -> list[IntMap]:
        """Internal shared implementation for global and resolution peak modes."""
        self._ensure_dp()

        eff_max_k = self._effective_max_k(max_k)
        self.max_k = eff_max_k

        if eff_max_k < 4:
            raise InvalidKError("max_k must be at least 4.")

        # Scores are a pure function of (DP table, max_k, scoring flags).
        # Recompute only if the signature has changed since the last call;
        # this is what lets back-to-back `run(top_n=...)` calls be cheap.
        scores_sig = self._scores_signature()
        if self.scores is None or self._scores_cache_sig != scores_sig:
            calculate_scores(self, plot=plot_scores)
            self._scores_cache_sig = scores_sig
        else:
            logger.debug("Scores cached, not recalculating")
        if self.scores is None or len(self.scores) == 0:
            return self._no_peaks_fallback()

        score_k_count = min(eff_max_k, len(self.scores))
        if score_k_count < 4:
            return self._no_peaks_fallback()

        active_peak_config = peak_config or self.peak_config

        if resolution_on:
            # Keep existing behavior for small trees: fallback to global peaks.
            if score_k_count < 50:
                top = max(1, min(3, score_k_count - 1))
                return self._run_peak_mode(
                    resolution_on=False,
                    top_n=top,
                    max_k=eff_max_k,
                    plot_scores=plot_scores,
                    compute_all_clusters=False,
                    peak_config=active_peak_config,
                )

            score_len = score_k_count - 1
            find_score_peaks(
                self,
                resolution_on=True,
                num_bins=num_bins,
                peaks_per_bin=1,
                k_start=2,
                k_end=score_len,
                plot=plot_scores,
                peak_config=active_peak_config,
            )
        else:
            score_len = min(eff_max_k, len(self.scores))
            if score_len <= 2:
                return self._no_peaks_fallback()

            find_score_peaks(
                self,
                global_peaks=top_n,
                resolution_on=False,
                k_start=1,
                k_end=score_len,
                plot=plot_scores,
                peak_config=active_peak_config,
            )

        # NOTE: do NOT clear ``self.clusters`` here. ``get_clusters(k)`` is a
        # per-k backtrack cache, and the DP table it was built from hasn't
        # changed (``_ensure_dp`` already cleared it on a true DP miss). Keeping
        # the cache means back-to-back runs with different ``top_n`` only run
        # backtrack for k values that are actually new this round.
        if (not resolution_on) and compute_all_clusters:
            for k_val in range(1, eff_max_k + 1):
                try:
                    self.get_clusters(k_val)
                except MissingDPTableError:
                    continue
        else:
            for k_val in self.peaks_by_rank or []:
                self.get_clusters(k_val)

        self.k = None
        return [
            self.clusters[kv]
            for kv in (self.peaks_by_rank or [])
            if kv in self.clusters
        ]

    def best_global(
        self,
        *,
        top_n: int = 1,
        max_k: Optional[int] = None,
        plot_scores: bool = True,
        compute_all_clusters: bool = False,
        peak_config: Optional[PeakConfig] = None,
    ) -> list[IntMap]:
        """
        Cluster-validity index-based global peak search.

        Returns a list of cluster maps in peak-rank order.
        """
        if top_n < 1:
            raise InvalidKError("`top_n` must be >= 1.")
        return self._run_peak_mode(
            resolution_on=False,
            top_n=top_n,
            max_k=max_k,
            plot_scores=plot_scores,
            compute_all_clusters=compute_all_clusters,
            peak_config=peak_config,
        )

    def best_by_resolution(
        self,
        *,
        num_bins: int = 3,
        max_k: Optional[int] = None,
        plot_scores: bool = True,
        peak_config: Optional[PeakConfig] = None,
    ) -> list[IntMap]:
        return self._run_peak_mode(
            resolution_on=True,
            num_bins=num_bins,
            max_k=max_k,
            plot_scores=plot_scores,
            peak_config=peak_config,
        )

    # ------------------------------------------------------------------ #
    #  Unified entry point                                                #
    # ------------------------------------------------------------------ #

    def run(
        self,
        *,
        k: Optional[int] = None,
        top_n: int = 1,
        by_resolution: bool = False,
        num_bins: Optional[int] = None,
        max_k: Optional[int] = None,
        max_k_limit: Optional[float] = None,
        plot_scores: bool = True,
        peak_config: Optional[PeakConfig] = None,
    ) -> dict[str, Any]:
        """
        Unified high-level entry point.

        Modes
        -----
        1. Exact k:
            pc.run(k=5)

        2. Global peaks:
            pc.run(top_n=3)

        3. Multi-resolution peaks (one per log-bin):
            pc.run(by_resolution=True, num_bins=3)

        All modes accept ``peak_config`` for tuning peak detection::

            from phytclust.config import PeakConfig
            pc.run(top_n=3, peak_config=PeakConfig(lambda_weight=0.5))

        Returns
        -------
        dict with keys:
            mode : str — "k", "global", or "resolution"
            k_values : list[int] — canonical selected k values
            selected_k : int | None — first selected k for convenience
            ks : list[int] — selected k values
            clusters : list[dict] — cluster maps in rank order
            scores : ndarray | None — score vector
            peaks : list[int] — same as ks (for convenience)
        """
        # Apply max_k_limit override for this run only
        saved_limit = self.max_k_limit
        if max_k_limit is not None:
            self.max_k_limit = max_k_limit

        try:
            result = self._run_inner(
                k=k,
                top_n=top_n,
                by_resolution=by_resolution,
                num_bins=num_bins,
                max_k=max_k,
                plot_scores=plot_scores,
                peak_config=peak_config,
            )
        finally:
            self.max_k_limit = saved_limit

        self._last_result = result
        return result

    def _run_inner(
        self,
        *,
        k: Optional[int],
        top_n: int,
        by_resolution: bool,
        num_bins: Optional[int],
        max_k: Optional[int],
        plot_scores: bool,
        peak_config: Optional[PeakConfig],
    ) -> dict[str, Any]:
        # If the caller explicitly asks for resolution mode, ignore any
        # lingering `self.k` from a previous run — otherwise a sequence like
        # `pc.run(k=5); pc.run(by_resolution=True, num_bins=3)` raises
        # "Cannot combine k with by_resolution=True" because `self.k` is
        # silently substituted. Explicit `by_resolution=True` is a clear
        # mode signal that should override the fallback.
        k_val = k if k is not None else (None if by_resolution else self.k)
        if k_val is not None:
            if k_val < 1:
                raise InvalidKError("k must be >= 1.")
            if by_resolution:
                raise ConfigurationError(
                    "Cannot combine `k` with `by_resolution=True`."
                )
            if top_n != 1:
                raise ConfigurationError("`top_n` is meaningless when `k` is given.")

            self._ensure_dp()
            cmap = self.get_clusters(k_val)

            self.k = int(k_val)
            self.peaks_by_rank = [int(k_val)]
            return self._build_run_result(
                mode="k",
                clusters=[cmap],
                selected_ks=[int(k_val)],
                exact_k=int(k_val),
                include_scores=False,
            )

        if top_n < 1:
            raise InvalidKError("`top_n` must be >= 1.")

        if by_resolution:
            clusters = self.best_by_resolution(
                num_bins=num_bins or self.num_bins,
                max_k=max_k,
                plot_scores=plot_scores,
                peak_config=peak_config or self.peak_config,
            )
            return self._build_run_result(
                mode="resolution",
                clusters=clusters,
                selected_ks=list(self.peaks_by_rank or []),
            )

        # global peak mode
        clusters = self.best_global(
            top_n=top_n,
            max_k=max_k,
            plot_scores=plot_scores,
            compute_all_clusters=self.compute_all_clusters,
            peak_config=peak_config or self.peak_config,
        )
        return self._build_run_result(
            mode="global",
            clusters=clusters,
            selected_ks=list(self.peaks_by_rank or []),
        )

    # ------------------------------------------------------------------ #
    #  Convenience wrappers                                               #
    # ------------------------------------------------------------------ #

    def plot(self, results_dir: Optional[str] = None, **kwargs) -> None:
        """Plot clustering results (requires matplotlib)."""
        from ..viz.cluster import plot_clusters

        plot_clusters(self, results_dir=results_dir, **kwargs)

    def save(
        self,
        results_dir: str,
        top_n: int = 1,
        filename: str = "phytclust_results.tsv",
        outlier: bool = True,
        k: Optional[int] = None,
        n: Optional[int] = None,
        output_all: bool = False,
    ) -> Optional[str]:
        """Save clustering results to file (requires pandas).

        ``k`` is the standard selector name. ``n`` is retained as a backwards-
        compatible alias and is only used if ``k`` is not provided.
        """
        if n is not None and k is None:
            warnings.warn(
                "Parameter 'n' is deprecated; use 'k' instead.",
                DeprecationWarning,
                stacklevel=2,
            )
        from ..io.save import save_clusters

        return save_clusters(
            self,
            results_dir=results_dir,
            top_n=top_n,
            filename=filename,
            outlier=outlier,
            n=(k if k is not None else n),
            output_all=output_all,
        )
