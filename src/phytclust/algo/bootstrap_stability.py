"""Bootstrap co-association stability analysis."""

from io import StringIO
from typing import Any, Optional

import numpy as np
from Bio import Phylo
from scipy import sparse

from .core import PhytClust
from ..exceptions import DataError


def _taxon_order(trees):
    """Sorted list of taxa present in *all* trees (intersection)."""
    sets = [
        {term.name for term in t.get_terminals()}
        for t in trees
    ]
    common = set.intersection(*sets)
    if not common:
        raise DataError("No common taxa across bootstrap trees.")
    return sorted(common)


def _labels_from_cmap(cmap, taxa, missing_label: int = -1):
    """Convert {leaf_obj -> cluster_id} into a label vector aligned with `taxa`."""
    name_to_cluster = {leaf.name: cid for leaf, cid in cmap.items()}
    labels = np.full(len(taxa), missing_label, dtype=int)
    for i, name in enumerate(taxa):
        if name in name_to_cluster:
            labels[i] = name_to_cluster[name]
    return labels


def _coassoc_from_labels(labels: np.ndarray) -> np.ndarray:
    """Compute co-association matrix from a (B, N) label array.

    For each pair (i, j), the co-association is the fraction of bootstrap
    replicates where both taxa were assigned to the same cluster (ignoring
    replicates where either taxon has label < 0).

    Vectorised: instead of looping over replicates and forming an N x N outer
    product each time (O(B * N^2) Python-level), we build a single sparse
    one-hot indicator ``M`` of shape (N, B * C) whose column ``b * C + c`` marks
    the taxa placed in cluster ``c`` of replicate ``b``. Then ``M @ M.T`` is
    exactly the same-cluster co-occurrence count for every pair at once, and the
    valid-pair normaliser is ``valid.T @ valid``.
    """
    B, N = labels.shape
    valid = labels >= 0
    valid_f = valid.astype(np.float64)

    # counts[i, j] = number of replicates where both taxa are valid.
    counts = valid_f.T @ valid_f

    if valid.any():
        rep_idx, taxon_idx = np.nonzero(valid)  # aligned (replicate, taxon) pairs
        lab = labels[rep_idx, taxon_idx]
        n_clusters = int(lab.max()) + 1  # global upper bound on cluster ids
        cols = rep_idx * n_clusters + lab
        indicator = sparse.csr_matrix(
            (np.ones(rep_idx.size), (taxon_idx, cols)),
            shape=(N, B * n_clusters),
        )
        numer = np.asarray((indicator @ indicator.T).todense(), dtype=np.float64)
    else:
        numer = np.zeros((N, N), dtype=np.float64)

    coassoc = np.zeros((N, N), dtype=np.float64)
    mask = counts > 0
    coassoc[mask] = numer[mask] / counts[mask]
    np.fill_diagonal(coassoc, 1.0)
    return coassoc


def _tree_to_newick(tree) -> str:
    """Serialise a Bio.Phylo tree to a Newick string (cheap, picklable)."""
    buf = StringIO()
    Phylo.write(tree, buf, "newick")
    return buf.getvalue()


def _replicate_labels(args) -> np.ndarray:
    """Cluster one bootstrap replicate and return its taxa-aligned label vector.

    Defined at module scope so it is picklable by ``ProcessPoolExecutor``.
    ``tree`` may be a Bio.Phylo tree (serial path) or a Newick string (parallel
    path); ``PhytClust`` accepts either.
    """
    tree, k, outgroup, min_cluster_size, pc_kwargs, taxa = args
    pc = PhytClust(
        tree=tree,
        outgroup=outgroup,
        min_cluster_size=min_cluster_size,
        **pc_kwargs,
    )
    res = pc.run(k=k, plot_scores=False)
    cmap = res["clusters"][0]
    return _labels_from_cmap(cmap, taxa)


def compute_coassoc_for_k(
    trees,
    k: int,
    *,
    outgroup: Optional[str] = None,
    min_cluster_size: int = 1,
    pc_kwargs: Optional[dict[str, Any]] = None,
    n_jobs: int = 1,
):
    """For a given k, run PhytClust on each bootstrap tree and compute
    co-association.

    Parameters
    ----------
    n_jobs : int, default=1
        Number of worker processes. Each replicate is fully independent, so
        this is embarrassingly parallel. ``1`` runs serially (unchanged
        behaviour); ``-1`` (or ``None``) uses all available cores; any other
        positive integer caps the pool at that many workers.

    Returns
    -------
    taxa : list[str]
        Taxon names in order.
    labels : ndarray of shape (B, N)
        Cluster IDs per replicate (or -1 for missing).
    coassoc : ndarray of shape (N, N)
        Co-association frequencies.
    """
    if pc_kwargs is None:
        pc_kwargs = {}

    taxa = _taxon_order(trees)
    B = len(trees)
    N = len(taxa)
    labels = np.full((B, N), -1, dtype=int)

    if n_jobs == 1 or B <= 1:
        for b, tree in enumerate(trees):
            labels[b, :] = _replicate_labels(
                (tree, k, outgroup, min_cluster_size, pc_kwargs, taxa)
            )
    else:
        from concurrent.futures import ProcessPoolExecutor

        max_workers = None if n_jobs in (-1, None) else int(n_jobs)
        # Newick strings pickle more cheaply than live tree objects.
        # and free of any Bio.Phylo cross-reference pickling quirks.
        payloads = [
            (_tree_to_newick(tree), k, outgroup, min_cluster_size, pc_kwargs, taxa)
            for tree in trees
        ]
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # executor.map preserves input order, so labels[b] stays aligned.
            for b, lab in enumerate(executor.map(_replicate_labels, payloads)):
                labels[b, :] = lab

    coassoc = _coassoc_from_labels(labels)
    return taxa, labels, coassoc


def stability_for_k(
    coassoc: np.ndarray, counts: Optional[np.ndarray] = None
) -> float:
    """Average co-association over off-diagonal pairs.

    If ``counts`` (the per-pair number of replicates in which both taxa were
    valid) is given, pairs that were never simultaneously valid (``counts == 0``,
    e.g. taxa stripped from every replicate) are excluded from the mean instead
    of being counted as structural zeros — otherwise they bias stability low.
    """
    N = coassoc.shape[0]
    triu_idx = np.triu_indices(N, k=1)
    vals = coassoc[triu_idx]
    if counts is not None:
        keep = counts[triu_idx] > 0
        vals = vals[keep]
    if vals.size == 0:
        return 0.0
    return float(np.mean(vals))


def choose_k_by_stability(
    trees,
    k_values: list[int],
    *,
    outgroup: Optional[str] = None,
    min_cluster_size: int = 1,
    pc_kwargs: Optional[dict[str, Any]] = None,
    n_jobs: int = 1,
):
    """For each k in k_values, compute co-association and stability.

    Returns
    -------
    dict with keys:
        best_k : int
        scores : dict[int, float]
        coassoc : dict[int, ndarray]
        taxa : list[str]
    """
    if not k_values:
        raise DataError("k_values is empty; nothing to evaluate.")

    scores = {}
    coassoc_by_k = {}
    taxa_ref = None

    for k in k_values:
        taxa, _labels, coassoc = compute_coassoc_for_k(
            trees,
            k,
            outgroup=outgroup,
            min_cluster_size=min_cluster_size,
            pc_kwargs=pc_kwargs,
            n_jobs=n_jobs,
        )
        if taxa_ref is None:
            taxa_ref = taxa
        elif taxa != taxa_ref:
            raise DataError(
                "Taxon order mismatch across k; this should not happen."
            )

        valid = (_labels >= 0).astype(np.float64)
        counts = valid.T @ valid
        scores[k] = stability_for_k(coassoc, counts)
        coassoc_by_k[k] = coassoc

    # Highest stability wins; ties break toward the smallest (most parsimonious) k.
    best_k = max(scores, key=lambda kv: (scores[kv], -kv))
    return {
        "best_k": best_k,
        "scores": scores,
        "coassoc": coassoc_by_k,
        "taxa": taxa_ref,
    }
