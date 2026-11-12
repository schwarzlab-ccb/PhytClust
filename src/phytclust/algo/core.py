import logging
import warnings
from dataclasses import dataclass, field
from math import ceil
from typing import Any, Optional

from pathlib import Path
from io import StringIO

from Bio import Phylo
from Bio.Phylo.NewickIO import NewickError

from ..exceptions import (
    InvalidKError,
    ConfigurationError,
    InvalidClusteringError,
    MissingDPTableError,
    ValidationError,
)
from Bio.Phylo.BaseTree import Tree

from ..algo.dp import (
    validate_clustering_parameters,
    prepare_tree,
    compute_dp_table,
    backtrack,
)
from ..algo.scoring import (
    calculate_scores,
    find_score_peaks,
    zero_distance_groups,
    zero_distance_splits,
)
from ..config import OutlierConfig, PeakConfig, RuntimeConfig
from ..metrics.indices import AlphaIndex, cluster_alpha
from ..utils.traversal import tree_fingerprint
from ..validation import positive_int

logger = logging.getLogger("phytclust")

IntMap = dict[Any, int]


def _coerce_to_tree(obj: Any) -> Tree:
    """Read a tree from a Tree object, Newick string, or file path.

    A string naming an existing file is read as a file; any other string must
    contain a parenthesised Newick tree, so a mistyped path is an error rather
    than a one-leaf tree. Raises ValidationError for other types and for read
    or parse errors.
    """
    if isinstance(obj, Tree):
        return obj

    def read_newick(source):
        try:
            return Phylo.read(source, "newick")
        except (OSError, ValueError, NewickError) as exc:
            raise ValidationError(f"Cannot read tree as Newick: {exc}") from exc

    if isinstance(obj, Path):
        return read_newick(obj)

    if isinstance(obj, str):
        if not obj.strip():
            raise ValidationError("Tree input is an empty string.")
        try:
            candidate = Path(obj)
            if candidate.is_file():
                return read_newick(candidate)
        except (OSError, ValueError):
            pass

        if "(" not in obj:
            raise ValidationError(
                f"No such file: {obj!r}. Pass an existing Newick file path or a "
                "Newick string such as '((A:1,B:1):1,C:2);'."
            )
        handle = StringIO(obj)
        return read_newick(handle)

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
    tree : Bio.Phylo.BaseTree.Tree, str, or pathlib.Path
        Input tree, Newick string, or path to a Newick file. A ``Tree`` object
        is prepared in place, not copied: unnamed and duplicate nodes are
        renamed, negative branch lengths are set to 0, missing lengths are
        filled in, single-child nodes are merged, and ``root_taxon`` re-roots
        it. Returned cluster maps are keyed by the leaves of this tree (with an
        outgroup, of an internal copy). Pass ``copy.deepcopy(tree)`` to keep
        the original unchanged.
    outgroup : str | None, default=None
        Name of the node whose subtree is excluded from clustering.
    root_taxon : str | None, default=None
        Node name to root on, or "midpoint" for midpoint rooting.
    min_cluster_size : int, default=1
        Minimum number of leaves allowed in a final cluster.
    k : int | None, default=None
        Default number of clusters for ``run()`` when no mode argument is
        given. ``run(k=...)`` overrides it for that call; ``top_n > 1`` or
        ``by_resolution=True`` run a peak search instead. ``run()`` never
        changes this attribute; the mode of the last run is recorded in
        ``exact_k`` (the fixed k, or None after a peak search).
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
        Weight of the branch-support penalty in effective branch lengths.

    Outlier handling
    ----------------
    outlier : OutlierConfig
        Controls outlier detection thresholds, DP penalty, and tie-breaking.
        See :class:`~phytclust.config.OutlierConfig`.

    Other flags
    -----------
    polytomy_mode : {"hard", "soft"}, default="soft"
        Soft mode treats a multifurcation as unresolved zero-length branching,
        so any group of two or more of its children may form a cluster. Hard
        mode forbids cross-child partial merges. Soft mode is bounded by
        `soft_polytomy_max_degree`; use hard mode for larger multifurcations.
    soft_polytomy_max_degree : int, default=12
        Maximum number of children allowed in soft mode. Computation grows
        exponentially with the number of children.
    preserve_dp_tables : bool, default=False
        Keep intermediate child DP rows instead of freeing them after merge.
        Useful for debugging/inspection; increases memory usage.
    compute_all_clusters : bool, default=False
        If True in best_global, compute and cache all k clusterings up to max_k.
    save_tied_optima : bool, default=False
        Record partitions that tie for the optimum at each backtracked k,
        in `tied_optima`. Ties are common under `polytomy_mode="soft"`, where
        regrouping children across a zero-length implied branch can leave the
        cost unchanged. Forces `preserve_dp_tables`, since the enumeration
        reads the per-child tables the DP would otherwise free.
    max_tied_optima : int, default=100
        Cap on partitions recorded per k. The optimal set can be exponentially
        large, so enumeration stops at the cap and flags the result truncated.
    """

    tree: Any
    outgroup: Optional[str] = None
    root_taxon: Optional[str] = None
    min_cluster_size: int = 1
    k: Optional[int] = None
    max_k: Optional[int] = None
    max_k_limit: float = 0.9
    num_bins: int = 3

    use_branch_support: bool = False
    min_support: float = 0.05
    support_weight: float = 1.0

    outlier: OutlierConfig = field(default_factory=OutlierConfig)
    use_penalized_beta_for_scoring: bool = False

    no_split_zero_length: bool = False
    zero_length_eps: float = 1e-12
    polytomy_mode: str = "soft"
    soft_polytomy_max_degree: int = 12
    preserve_dp_tables: bool = False

    save_tied_optima: bool = False
    max_tied_optima: int = 100

    dp_float32: bool = False

    compute_all_clusters: bool = False
    runtime_config: RuntimeConfig = field(default_factory=RuntimeConfig)
    peak_config: PeakConfig = field(default_factory=PeakConfig)

    def __post_init__(self) -> None:
        if self.runtime_config is None:
            self.runtime_config = RuntimeConfig()
        if self.peak_config is None:
            self.peak_config = PeakConfig()
        if self.outlier is None:
            self.outlier = OutlierConfig()
        self.tree = _coerce_to_tree(self.tree)

        self.name_leaves_per_node = {}
        self.num_leaves_per_node = {}
        self.backptr = {}
        self.dp_table = None
        self.postorder_nodes = None
        self.node_to_id = None
        self.dp_children = None
        self.num_terminals = 0
        self._tree_wo_outgroup = None

        self.scores = None
        self.peaks_by_rank = None
        self.peaks_by_resolution = None
        self.resolution_info = None
        self.bin_ranges_current = None
        self.peak_ranking_details: list = []
        self.tied_optima: dict[int, dict[str, Any]] = {}
        self.score_plot = None
        self.secondary_score_plot = None
        self._support_scale: Optional[float] = None
        self.exact_k: Optional[int] = self.k
        self.resolution_fell_back = False
        self.zero_length_splits: dict[int, int] = {}
        self.alpha_by_k: dict[int, dict[str, Any]] = {}
        self._alpha_index: Optional[AlphaIndex] = None

        self._dp_ready = False
        self._dp_cache_sig: Optional[tuple] = None
        self._dp_cap: Optional[int] = None
        self._score_arrays_cache: Optional[tuple] = None
        self._score_arrays_cap: Optional[int] = None
        self._score_arrays_sig: Optional[tuple] = None
        self._scores_cache_sig: Optional[tuple] = None
        self.cluster_maps: dict[int, IntMap] = {}

        user_max_k = self.max_k
        prepare_tree(self)
        self._auto_max_k: Optional[int] = self.max_k if user_max_k is None else None
        self._prepared_tree_sig = (self._hash_tree(), self.outgroup)

    def __repr__(self) -> str:
        parts = [f"terminals={self.num_terminals}"]
        if self.outgroup:
            parts.append(f"outgroup={self.outgroup!r}")
        max_k = self.max_k
        if max_k is not None:
            parts.append(f"max_k={max_k}")
        peaks = self.peaks_by_rank
        if peaks:
            parts.append(f"peaks={list(peaks)}")
        elif self.exact_k is not None:
            parts.append(f"k={self.exact_k}")
        if not self._dp_ready:
            parts.append("dp=uncomputed")
        return f"PhytClust({', '.join(parts)})"

    def _hash_tree(self) -> int:
        """Return a fingerprint of topology, names, lengths, and confidence.

        Ignore child order. Backtracking uses the order stored in ``dp_children``.
        """
        return tree_fingerprint(self.tree.root)

    def _dp_signature(self) -> tuple:
        """
        Return the cache key for the tree and parameters that affect clustering.
        """
        return (
            self._hash_tree(),
            self.outgroup,
            self.min_cluster_size,
            self.no_split_zero_length,
            self.zero_length_eps,
            self.use_branch_support,
            self.min_support,
            self.support_weight,
            self.polytomy_mode,
            self.soft_polytomy_max_degree,
            self.preserve_dp_tables,
            self.save_tied_optima,
            self.dp_float32,
            self.outlier.size_threshold,
            self.outlier.prefer_fewer,
            self.outlier.penalty_enabled,
            self.outlier.ratio_weight,
            self.outlier.ratio_mode,
        )

    def _scores_signature(self) -> tuple:
        """Return the cache key for scores."""
        return (
            self._dp_cache_sig,
            self.max_k,
            self.use_penalized_beta_for_scoring,
        )

    def _resolve_dp_cap(self, required_cap: Optional[int]) -> int:
        """Choose the maximum k to compute, limited to the number of leaves.

        Include a caller's ``max_k`` even if the current request needs fewer clusters.
        """
        if required_cap is not None:
            cap = required_cap
            if self._max_k_is_user_set():
                cap = max(cap, self.max_k)
        elif self.max_k is not None:
            cap = self.max_k
        else:
            cap = max(2, ceil(self.num_terminals * self.max_k_limit))
        return max(1, min(self.num_terminals, int(cap)))

    def _ensure_dp(self, required_cap: Optional[int] = None) -> None:
        """Reuse valid DP tables or rebuild them for the current tree and settings."""
        tree_sig = (self._hash_tree(), self.outgroup)
        if tree_sig != self._prepared_tree_sig:
            user_max_k = self.max_k if self._max_k_is_user_set() else None
            root_taxon = self.root_taxon
            self.max_k = user_max_k
            self.root_taxon = None
            try:
                prepare_tree(self)
            finally:
                self.root_taxon = root_taxon
            if user_max_k is None:
                self._auto_max_k = self.max_k
            self._prepared_tree_sig = (self._hash_tree(), self.outgroup)
            self._dp_ready = False
        current = self._dp_signature()
        needed_cap = self._resolve_dp_cap(required_cap)

        sig_matches = self._dp_ready and (self._dp_cache_sig == current)
        if sig_matches and self._dp_cap is not None and needed_cap <= self._dp_cap:
            logger.debug("Reusing cached DP tables.")
            return

        new_cap = needed_cap
        if sig_matches and self._dp_cap is not None:
            grown = 2 * self._dp_cap
            if self.max_k is not None:
                grown = min(grown, self.max_k)
            new_cap = max(needed_cap, min(grown, self.num_terminals))

        # Cluster maps are keyed by leaves of the rebuilt tree and are discarded.
        # When only the cap grows, the DP is exact and smaller-k costs are
        # unchanged, so the previous scores and peaks remain valid.
        self.cluster_maps = {}
        self.alpha_by_k = {}
        self._alpha_index = None
        if not sig_matches:
            self.scores = None
            self.peaks_by_rank = None
            self._scores_cache_sig = None
            self._score_arrays_cache = None
            self._score_arrays_cap = None
            self._score_arrays_sig = None

        self._dp_ready = False
        validate_clustering_parameters(self)
        self._dp_cap = new_cap
        compute_dp_table(self)

        self._dp_ready = True
        self._dp_cache_sig = current

        if self.max_k is None or self.max_k < 1:
            self.max_k = max(2, ceil(self.num_terminals * self.max_k_limit))
            self._auto_max_k = self.max_k

    def _max_k_is_user_set(self) -> bool:
        """True if ``self.max_k`` was set by the caller, not derived by a run."""
        return self.max_k is not None and self.max_k != self._auto_max_k

    def _effective_max_k(self, max_k: Optional[int] = None) -> int:
        """Resolve max_k without mutating self.

        Precedence is:
        1. explicit method argument ``max_k`` (applies to that call only)
        2. ``self.max_k`` if the caller set it
        3. derived cap from the current ``max_k_limit``
        """
        if max_k is not None:
            return min(self.num_terminals, max_k)
        if self._max_k_is_user_set():
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
        """Build the result dictionary returned by ``run()``.

        Canonical keys are ``k_values`` (list[int]) and ``selected_k`` (int|None).
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

    def alpha_info(self, k: int, cluster_map: Optional[IntMap] = None) -> dict[str, Any]:
        """Return alpha and branch statistics for a partition.

        Results for get_clusters(k) are cached until the DP is rebuilt; a
        supplied map is always computed.
        """
        k = positive_int(k, "k", error=InvalidKError)
        self._ensure_dp(required_cap=k)
        cache_result = cluster_map is None or cluster_map is self.cluster_maps.get(int(k))
        if cluster_map is None:
            cached = self.alpha_by_k.get(int(k))
            if cached is not None:
                return cached
            cluster_map = self.get_clusters(k)
        active_tree = (
            self._tree_wo_outgroup
            if (self.outgroup and self._tree_wo_outgroup is not None)
            else self.tree
        )
        if self._alpha_index is None:
            self._alpha_index = AlphaIndex(active_tree)
        info = {
            "k": int(k),
            **cluster_alpha(active_tree, cluster_map, index=self._alpha_index),
        }
        if cache_result:
            self.alpha_by_k[int(k)] = info
        return info

    def _compute_and_log_alphas(
        self,
        selected_ks: list[int],
        clusters: list[IntMap],
    ) -> list[dict[str, Any]]:
        """Calculate and log alpha for each selected partition.

        Alpha = mean extra-cluster branch length / mean intra-cluster branch length.
        """
        details: list[dict[str, Any]] = []
        for k, cluster_map in zip(selected_ks, clusters):
            info = self.alpha_info(k, cluster_map)
            details.append(info)
            logger.info(
                "alpha(k=%d) = %.6g  (avg_extra=%.6g, avg_intra=%.6g, "
                "n_extra=%d, n_intra=%d)",
                k,
                info["alpha"],
                info["avg_extra_branch_length"],
                info["avg_intra_branch_length"],
                info["n_extra_nodes"],
                info["n_intra_nodes"],
            )
        return details

    @property
    def plot_config(self):
        """Return the plotting settings."""
        return self.runtime_config.plot

    def get_clusters(self, k: int, *, verbose: bool = False) -> IntMap:
        """Return the exact k-cluster partition (cached after first call)."""
        k = positive_int(k, "k", error=InvalidKError)
        self._ensure_dp(required_cap=k)
        return self._clusters(k, verbose=verbose)

    def _clusters(self, k: int, *, verbose: bool = False) -> IntMap:
        """Return a cached partition or backtrack through the DP tables.

        The caller must first call ``_ensure_dp`` for the largest required k.
        """
        if k in self.cluster_maps:
            return self.cluster_maps[k]

        cluster_map = backtrack(self, k, verbose=verbose)
        self.cluster_maps[k] = cluster_map
        return cluster_map

    def _check_zero_distance_splits(self, cluster_counts: list[int]) -> None:
        """Warn when a selected partition separates leaves at zero distance.

        Records ``{k: number of split groups}`` in ``zero_length_splits``.
        """
        self.zero_length_splits = {}
        if not cluster_counts:
            return
        groups = zero_distance_groups(self)
        if not groups:
            return
        details = []
        for count in cluster_counts:
            cluster_map = self.cluster_maps.get(count)
            if cluster_map is None:
                continue
            split = zero_distance_splits(groups, cluster_map)
            if split:
                self.zero_length_splits[count] = len(split)
                names = sorted(str(leaf.name) for leaf in split[0])[:2]
                details.append(
                    f"k={count}: {len(split)} group(s), e.g. {names[0]!r} and {names[1]!r}"
                )
        if details:
            logger.warning(
                "Selected partitions separate leaves at zero distance from each "
                "other (for example identical sequences): %s. Interpret these "
                "clusters with caution.",
                "; ".join(details),
            )

    def _no_peaks_fallback(self) -> list[IntMap]:
        """Reset selection state and return no clusters."""
        logger.warning("No score peaks found; no partition was selected.")
        self.exact_k = None
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
        self.exact_k = None
        if resolution_on:
            self.resolution_fell_back = False
        from_user = max_k is None and self._max_k_is_user_set()
        eff_max_k = self._effective_max_k(max_k)
        self.max_k = eff_max_k
        if not from_user:
            self._auto_max_k = eff_max_k

        if eff_max_k < 4:
            if self.num_terminals < 4:
                raise InvalidKError(
                    f"Automatic selection of k needs at least 4 leaves; this tree "
                    f"has {self.num_terminals}. Use run(k=...) for a fixed k."
                )
            if max_k is not None or from_user:
                raise InvalidKError(f"max_k must be at least 4 (got {eff_max_k}).")
            raise InvalidKError(
                f"max_k_limit={self.max_k_limit} gives max_k={eff_max_k} for "
                f"{self.num_terminals} leaves; automatic selection needs at least 4. "
                "Increase max_k_limit or pass max_k."
            )

        self._ensure_dp(required_cap=eff_max_k)

        scores_sig = self._scores_signature()
        if self.scores is None or self._scores_cache_sig != scores_sig:
            calculate_scores(self)
            self._scores_cache_sig = scores_sig
        else:
            logger.debug("Reusing cached scores.")
        if self.scores is None or len(self.scores) == 0:
            return self._no_peaks_fallback()

        score_k_count = min(eff_max_k, len(self.scores))
        if score_k_count < 4:
            return self._no_peaks_fallback()

        active_peak_config = peak_config or self.peak_config

        if resolution_on:
            min_scores = active_peak_config.resolution_min_scores
            if score_k_count < min_scores:
                top = max(1, min(max(top_n, 3), score_k_count - 1))
                logger.warning(
                    "Resolution mode needs at least %d scored cluster counts; "
                    "this run has %d (max_k=%d). Running a global search for the "
                    "top %d peaks instead. Lower PeakConfig.resolution_min_scores "
                    "to use resolution mode on this tree.",
                    min_scores,
                    score_k_count,
                    eff_max_k,
                    top,
                )
                self.resolution_fell_back = True
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
                peaks_per_bin=max(1, int(top_n)),
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

        if (not resolution_on) and compute_all_clusters:
            self._ensure_dp(required_cap=eff_max_k)
            for k in range(1, eff_max_k + 1):
                try:
                    self._clusters(k)
                except (MissingDPTableError, InvalidClusteringError):
                    continue
        else:
            for k in self.peaks_by_rank or []:
                self.get_clusters(k)

        self._check_zero_distance_splits(list(self.peaks_by_rank or []))
        return [
            self.cluster_maps[k]
            for k in (self.peaks_by_rank or [])
            if k in self.cluster_maps
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
        Select the highest-ranked peaks in the clustering scores.

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
        top_n: int = 1,
        max_k: Optional[int] = None,
        plot_scores: bool = True,
        peak_config: Optional[PeakConfig] = None,
    ) -> list[IntMap]:
        """Select up to ``top_n`` score peaks per logarithmic resolution bin.

        With fewer than ``PeakConfig.resolution_min_scores`` scores (default 50),
        warn and use a global search for at least three peaks.
        """
        if top_n < 1:
            raise InvalidKError("top_n must be at least 1.")
        if num_bins < 1:
            raise ConfigurationError("num_bins must be at least 1.")
        return self._run_peak_mode(
            resolution_on=True,
            num_bins=num_bins,
            top_n=top_n,
            max_k=max_k,
            plot_scores=plot_scores,
            peak_config=peak_config,
        )

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
        Run clustering for a fixed k or select k values from score peaks.

        Modes
        -----
        1. Exact k:
            pc.run(k=5)

        2. Global peaks:
            pc.run(top_n=3)

        3. Multi-resolution peaks (one per log-bin):
            pc.run(by_resolution=True, num_bins=3)

        Peak modes accept ``peak_config`` to configure peak detection::

            from phytclust.config import PeakConfig
            pc.run(top_n=3, peak_config=PeakConfig(prominence_weight=0.5))

        Returns
        -------
        dict with keys:
            mode : str — "k", "global", or "resolution"
            k_values : list[int] — canonical selected k values
            selected_k : int | None — first selected k for convenience
            ks : list[int] — selected k values
            clusters : list[dict] — global rank order, or increasing k in resolution mode
            scores : ndarray | None — score vector
        """
        saved_limit = self.max_k_limit
        if max_k_limit is not None:
            self.max_k_limit = max_k_limit

        try:
            result = self._select_mode(
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

        return result

    def _select_mode(
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
        use_default_k = k is None and not by_resolution and top_n == 1
        k = k if k is not None else (self.k if use_default_k else None)
        if k is not None:
            k = positive_int(k, "k", error=InvalidKError)
            if by_resolution:
                raise ConfigurationError(
                    "Cannot combine `k` with `by_resolution=True`."
                )
            if top_n != 1:
                raise ConfigurationError("top_n must be 1 when k is specified.")

            self._ensure_dp(required_cap=int(k))
            cluster_map = self.get_clusters(k)
            self._check_zero_distance_splits([int(k)])

            self.exact_k = int(k)
            self.peaks_by_rank = [int(k)]
            return self._build_run_result(
                mode="k",
                clusters=[cluster_map],
                selected_ks=[int(k)],
                exact_k=int(k),
                include_scores=False,
            )

        if top_n < 1:
            raise InvalidKError("`top_n` must be >= 1.")

        if by_resolution:
            clusters = self.best_by_resolution(
                num_bins=self.num_bins if num_bins is None else num_bins,
                top_n=top_n,
                max_k=max_k,
                plot_scores=plot_scores,
                peak_config=peak_config or self.peak_config,
            )
            return self._build_run_result(
                mode="global" if self.resolution_fell_back else "resolution",
                clusters=clusters,
                selected_ks=list(self.peaks_by_rank or []),
            )

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

    def plot(self, results_dir: Optional[str] = None, **kwargs) -> None:
        """Plot clustering results (requires matplotlib)."""
        from ..viz.cluster import plot_clusters

        plot_clusters(self, results_dir=results_dir, **kwargs)

    def save(
        self,
        results_dir: str,
        top_n: int = 1,
        filename: Optional[str] = None,
        outlier: Optional[bool] = None,
        k: Optional[int] = None,
        n: Optional[int] = None,
        output_all: bool = False,
    ) -> Optional[str]:
        """Save clustering results as a tab-separated file.

        ``filename`` and ``outlier`` default to ``runtime_config.save``
        (``tsv_name`` and ``outlier``). ``n`` is an alias for ``k``.
        """
        if filename is None:
            filename = self.runtime_config.save.tsv_name
        if outlier is None:
            outlier = self.runtime_config.save.outlier
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
