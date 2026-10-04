import numpy as np

from ..exceptions import ConfigurationError, InvalidKError


_EPS_CACHE: dict = {}


def dtype_eps(dtype) -> float:
    """``np.finfo(dtype).eps``, cached.

    finfo builds a machine-limits object on every call; the DP asks for this
    once per k per node, and the answer only depends on the dtype.
    """
    eps = _EPS_CACHE.get(dtype)
    if eps is None:
        eps = float(np.finfo(dtype).eps)
        _EPS_CACHE[dtype] = eps
    return eps


def tie_atol(value: float, dtype) -> float:
    """Absolute tolerance for treating two DP costs as tied.

    A fixed 1e-12 is far below float32 ULP (~1e-2..1e-4 for realistic
    dispersion sums), so under ``dp_float32`` it essentially never fires and the
    outlier tie-break silently degrades. Scale by the dtype epsilon and the
    magnitude, keeping a 1e-12 floor to preserve the previous float64 behaviour.
    """
    return max(1e-12, 8.0 * dtype_eps(dtype) * max(1.0, abs(float(value))))


def validate_args(pc) -> None:
    if pc.k is not None and pc.k < 1:
        raise InvalidKError("k must be ≥ 1 if provided.")
    if not 0 < pc.max_k_limit <= 1:
        raise ConfigurationError("max_k_limit must be between 0 and 1")

    if pc.outlier.ratio_weight < 0:
        raise ConfigurationError("outlier.ratio_weight must be ≥ 0")

    if pc.outlier.ratio_mode not in {"exp", "inverse", "power"}:
        raise ConfigurationError(
            "outlier.ratio_mode must be one of: 'exp', 'inverse', 'power'"
        )



    polytomy_mode = getattr(pc, "polytomy_mode", "hard")
    if polytomy_mode not in {"hard", "soft"}:
        raise ConfigurationError("polytomy_mode must be 'hard' or 'soft'.")

    # Soft mode is exponential in the node degree.
    if polytomy_mode == "soft":
        max_deg = getattr(pc, "soft_polytomy_max_degree", 12)
        if max_deg < 2:
            raise ConfigurationError("soft_polytomy_max_degree must be ≥ 2.")


def penalty_active(pc) -> bool:
    """True iff the outlier cluster-formation penalty should shape the DP."""
    o = pc.outlier
    return bool(
        getattr(o, "penalty_enabled", False)
        and o.size_threshold is not None
        and o.ratio_weight > 0
    )


def cluster_formation_penalty(n_leaves: int, pc) -> float:
    """Additive penalty for forming one cluster of ``n_leaves`` leaves.

    Zero for clusters at or above ``size_threshold``; otherwise
    ``ratio_weight * shape(n_leaves)`` per ``ratio_mode`` (see OutlierConfig).
    Being per-cluster, this term is additive over clusters and therefore
    decomposes exactly over the tree DP.
    """
    o = pc.outlier
    thresh = o.size_threshold
    if thresh is None or n_leaves >= thresh:
        return 0.0
    s = max(1, int(n_leaves))
    mode = o.ratio_mode
    if mode == "inverse":
        shape = 1.0 / s
    elif mode == "exp":
        shape = float(np.expm1(thresh - s))
    else:  # "power": linear in the size deficit
        shape = float(thresh - s)
    return float(o.ratio_weight) * shape


def node_support_factor(node, pc) -> float:
    if not getattr(pc, "use_branch_support", False):
        return 1.0
    raw_support = node.confidence if node.confidence is not None else 100.0
    return max(raw_support / 100.0, pc.min_support)


def state_better(
    cand_raw: float,
    cand_ns: int,
    best_raw: float,
    best_ns: int,
    *,
    use_outlier: bool,
    prefer_fewer: bool,
    atol: float = 1e-12,
) -> bool:
    """Return True iff candidate state is better than current best."""
    if not np.isfinite(cand_raw):
        return False
    if not np.isfinite(best_raw):
        return True

    if use_outlier and prefer_fewer:
        if cand_ns != best_ns:
            return cand_ns < best_ns
        return cand_raw < best_raw - atol

    if cand_raw < best_raw - atol:
        return True

    if use_outlier and abs(cand_raw - best_raw) <= atol:
        return cand_ns < best_ns

    return False


def child_as_atomic_cluster_cost(child, pc) -> float:
    """
    Cost of taking the entire child subtree as ONE cluster and lifting its root
    from `child` up to its parent node.
    """
    child_id = pc.node_to_id[child]
    child_raw0 = pc.raw_dp_table[child_id][0]
    if not np.isfinite(child_raw0):
        return np.inf
    return child_raw0 + pc.num_leaves_per_node[child] * eff_length(child, pc)


def eff_length(node, pc):
    """
    If `pc.use_branch_support` is True, we augment the branch length
    with a penalty derived from the support value:

        eff_len = branch_length + support_weight * (-log(p))

    where p = max(confidence/100, pc.min_support).
    """
    base = node.branch_length or 0.0
    if not pc.use_branch_support:
        return base

    raw = node.confidence if node.confidence is not None else 100.0
    p = max(raw / 100.0, pc.min_support)
    penalty = -np.log(p)
    return base + pc.support_weight * penalty


def subtree_all_zero(node, pc, eps: float = 1e-12):
    """True if every edge strictly below `node` has near-zero effective length."""
    stack = list(node.clades)
    while stack:
        child = stack.pop()
        if eff_length(child, pc) > eps:
            return False
        if not child.is_terminal():
            stack.extend(child.clades)
    return True
