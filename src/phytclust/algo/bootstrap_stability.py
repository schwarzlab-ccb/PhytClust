"""Measure how consistently taxon pairs group across bootstrap trees."""

import copy
from numbers import Integral
from typing import Any, Optional

import numpy as np
from scipy import sparse

from .core import PhytClust
from .dp.table import _copy_tree
from ..exceptions import ConfigurationError, DataError, InvalidKError
from ..utils.traversal import terminals


def _positive_integer(value, name, error=ConfigurationError):
    if (
        isinstance(value, (bool, np.bool_))
        or not isinstance(value, Integral)
        or value < 1
    ):
        raise error(f"{name} must be a positive integer.")


def _common_taxon_names(trees):
    """Return sorted leaf names present in every tree; require unique named leaves."""
    if not trees:
        raise DataError("No bootstrap trees supplied.")
    taxon_sets = []
    for replicate_index, tree in enumerate(trees):
        names = [leaf.name for leaf in terminals(tree.root)]
        if any(not isinstance(name, str) or not name for name in names):
            raise DataError(f"Bootstrap tree {replicate_index} has unnamed leaves.")
        if len(names) != len(set(names)):
            raise DataError(
                f"Bootstrap tree {replicate_index} has duplicate leaf names."
            )
        taxon_sets.append(set(names))
    common_taxa = set.intersection(*taxon_sets)
    if not common_taxa:
        raise DataError("No common taxa across bootstrap trees.")
    return sorted(common_taxa)


def _cluster_labels_for_taxa(cluster_map, taxa):
    """Return cluster IDs in taxon order, using -1 for unassigned taxa."""
    name_to_cluster = {
        leaf.name: cluster_id for leaf, cluster_id in cluster_map.items()
    }
    return np.array([name_to_cluster.get(name, -1) for name in taxa], dtype=int)


def _pairwise_coassociation(labels):
    """Return same-cluster frequencies for a replicate-by-taxon label matrix.

    Negative labels are unassigned. Exclude those replicates for each pair.
    Use zero for unobserved pairs and one on the diagonal.
    """
    labels = np.asarray(labels)
    if labels.ndim != 2 or not np.issubdtype(labels.dtype, np.integer):
        raise DataError("Cluster labels must be a two-dimensional integer array.")
    replicate_count, taxon_count = labels.shape
    valid = labels >= 0
    valid_values = valid.astype(float)
    valid_pair_counts = valid_values.T @ valid_values
    if valid.any():
        replicate_indices, taxon_indices = np.nonzero(valid)
        cluster_keys = np.column_stack((replicate_indices, labels[valid]))
        _, column_indices = np.unique(cluster_keys, axis=0, return_inverse=True)
        indicator = sparse.csr_matrix(
            (np.ones(len(column_indices)), (taxon_indices, column_indices)),
            shape=(taxon_count, int(column_indices.max()) + 1),
        )
        same_cluster_counts = (indicator @ indicator.T).toarray()
    else:
        same_cluster_counts = np.zeros((taxon_count, taxon_count))
    coassociation = np.zeros((taxon_count, taxon_count))
    np.divide(
        same_cluster_counts,
        valid_pair_counts,
        out=coassociation,
        where=valid_pair_counts > 0,
    )
    np.fill_diagonal(coassociation, 1.0)
    return coassociation


def _cluster_bootstrap_replicate(arguments):
    """Cluster supplied k values on an independent tree, reusing its DP tables."""
    tree, k_values, outgroup, min_cluster_size, clustering_options, taxa = arguments
    clustering = PhytClust(
        tree=_copy_tree(tree),
        outgroup=outgroup,
        min_cluster_size=min_cluster_size,
        **copy.deepcopy(clustering_options),
    )
    clustering._ensure_dp(required_cap=max(k_values))
    return {
        k: _cluster_labels_for_taxa(clustering.get_clusters(k), taxa) for k in k_values
    }


def _bootstrap_labels(
    trees, k_values, outgroup, min_cluster_size, clustering_options, n_jobs
):
    """Collect replicate labels for all candidate counts using one worker pool."""
    trees = list(trees)
    taxa = _common_taxon_names(trees)
    for k in k_values:
        _positive_integer(k, "k", InvalidKError)
    _positive_integer(min_cluster_size, "min_cluster_size")
    if n_jobs is not None:
        if (
            isinstance(n_jobs, (bool, np.bool_))
            or not isinstance(n_jobs, Integral)
            or (n_jobs != -1 and n_jobs < 1)
        ):
            raise ConfigurationError("n_jobs must be a positive integer, -1, or None.")
    options = {} if clustering_options is None else dict(clustering_options)
    if "tree" in options or "outgroup" in options or "min_cluster_size" in options:
        raise ConfigurationError(
            "Set tree, outgroup, and min_cluster_size through their named arguments."
        )
    payloads = [
        (tree, k_values, outgroup, min_cluster_size, options, taxa) for tree in trees
    ]
    labels_by_k = {k: np.full((len(trees), len(taxa)), -1, dtype=int) for k in k_values}
    if n_jobs == 1 or len(trees) == 1:
        results = map(_cluster_bootstrap_replicate, payloads)
        for replicate_index, result in enumerate(results):
            for k in k_values:
                labels_by_k[k][replicate_index] = result[k]
    else:
        from concurrent.futures import ProcessPoolExecutor

        max_workers = None if n_jobs in (-1, None) else int(n_jobs)
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            for replicate_index, result in enumerate(
                executor.map(_cluster_bootstrap_replicate, payloads)
            ):
                for k in k_values:
                    labels_by_k[k][replicate_index] = result[k]
    return taxa, labels_by_k


def compute_coassoc_for_k(
    trees,
    k: int,
    *,
    outgroup: Optional[str] = None,
    min_cluster_size: int = 1,
    pc_kwargs: Optional[dict[str, Any]] = None,
    n_jobs: Optional[int] = 1,
):
    """Return common taxa, replicate labels, and pairwise co-association for k.

    Preserve input trees. Labels have shape (replicates, taxa); -1 means
    unassigned. Co-association has shape (taxa, taxa).
    Use n_jobs=1 for serial execution, or -1/None for available worker processes.
    """
    taxa, labels_by_k = _bootstrap_labels(
        trees, [k], outgroup, min_cluster_size, pc_kwargs, n_jobs
    )
    labels = labels_by_k[k]
    return taxa, labels, _pairwise_coassociation(labels)


def stability_for_k(coassoc: np.ndarray, counts: Optional[np.ndarray] = None) -> float:
    """Return mean pair consistency: 2 * abs(co-association - 0.5).

    Both consistent grouping and consistent separation score one; a pair
    together in half the replicates scores zero. Exclude unobserved pairs
    when counts are supplied. Without counts, every off-diagonal pair is used.
    Return zero when no pairs remain. This measure alone does not identify
    a useful number of clusters.
    """
    coassociation = np.asarray(coassoc, dtype=float)
    if coassociation.ndim != 2 or coassociation.shape[0] != coassociation.shape[1]:
        raise DataError("Co-association must be a square matrix.")
    if not np.all(np.isfinite(coassociation)) or np.any(
        (coassociation < 0) | (coassociation > 1)
    ):
        raise DataError(
            "Co-association values must be finite and between zero and one."
        )
    if not np.allclose(coassociation, coassociation.T):
        raise DataError("Co-association must be symmetric.")
    pair_indices = np.triu_indices(len(coassociation), k=1)
    pair_frequencies = coassociation[pair_indices]
    if counts is not None:
        counts = np.asarray(counts, dtype=float)
        if (
            counts.shape != coassociation.shape
            or not np.all(np.isfinite(counts))
            or np.any(counts < 0)
            or not np.allclose(counts, counts.T)
        ):
            raise DataError(
                "Valid-pair counts must match the matrix and be finite, symmetric, and non-negative."
            )
        pair_frequencies = pair_frequencies[counts[pair_indices] > 0]
    return (
        float(np.mean(2 * np.abs(pair_frequencies - 0.5)))
        if pair_frequencies.size
        else 0.0
    )


def choose_k_by_stability(
    trees,
    k_values: list[int],
    *,
    outgroup: Optional[str] = None,
    min_cluster_size: int = 1,
    pc_kwargs: Optional[dict[str, Any]] = None,
    n_jobs: Optional[int] = 1,
):
    """Compare pair consistency for supplied candidate cluster counts.

    Supply candidates from clustering-score peaks or another selection rule.
    Consistency alone can favor trivial partitions. Return best_k, scores,
    coassoc matrices, and taxon order. Equal scores select the smaller k.
    Reuse each replicate's DP across candidates.
    """
    k_values = list(k_values)
    if not k_values:
        raise DataError("k_values is empty; nothing to evaluate.")
    taxa, labels_by_k = _bootstrap_labels(
        trees, k_values, outgroup, min_cluster_size, pc_kwargs, n_jobs
    )
    scores, coassoc_by_k = {}, {}
    for k, labels in labels_by_k.items():
        coassociation = _pairwise_coassociation(labels)
        valid = (labels >= 0).astype(float)
        scores[k] = stability_for_k(coassociation, valid.T @ valid)
        coassoc_by_k[k] = coassociation
    return {
        "best_k": max(scores, key=lambda k: (scores[k], -k)),
        "scores": scores,
        "coassoc": coassoc_by_k,
        "taxa": taxa,
    }
