import pathlib

import pytest
from Bio import Phylo

from phytclust.algo import core as core_mod
from phytclust.algo import dp as dp_mod
from phytclust.algo import scoring as scoring_mod
from phytclust.algo.core import PhytClust

TREE_PATH = pathlib.Path(__file__).parent.parent / "examples" / "sample_tree.nwk"


def _load_tree():
    return Phylo.read(TREE_PATH, "newick")


class _Counter:
    """Tiny call counter that also records the last args/kwargs seen."""

    def __init__(self, real):
        self.real = real
        self.calls = 0

    def __call__(self, *args, **kwargs):
        self.calls += 1
        return self.real(*args, **kwargs)


@pytest.fixture
def counted(monkeypatch):
    """
    Replace ``compute_dp_table``, ``calculate_scores``, and ``backtrack`` with
    counting wrappers everywhere they're referenced from ``core.py``.

    Both the original module (``algo.dp`` / ``algo.scoring``) and the names
    re-imported into ``algo.core`` must be patched, because ``core.py`` does
    ``from ..algo.dp import compute_dp_table, backtrack`` at import time.
    """
    dp_counter = _Counter(dp_mod.compute_dp_table)
    bt_counter = _Counter(dp_mod.backtrack)
    sc_counter = _Counter(scoring_mod.calculate_scores)

    monkeypatch.setattr(dp_mod, "compute_dp_table", dp_counter)
    monkeypatch.setattr(scoring_mod, "calculate_scores", sc_counter)
    monkeypatch.setattr(dp_mod, "backtrack", bt_counter)

    monkeypatch.setattr(core_mod, "compute_dp_table", dp_counter)
    monkeypatch.setattr(core_mod, "calculate_scores", sc_counter)
    monkeypatch.setattr(core_mod, "backtrack", bt_counter)

    return {"dp": dp_counter, "scores": sc_counter, "backtrack": bt_counter}


# --------------------------------------------------------------------------- #
#  DP cache: reused when nothing relevant changes                             #
# --------------------------------------------------------------------------- #


def test_dp_reused_across_run_calls(counted):
    pc = PhytClust(tree=_load_tree())

    pc.run(top_n=1, plot_scores=False)
    dp_after_first = counted["dp"].calls

    pc.run(top_n=1, plot_scores=False)

    assert counted["dp"].calls == dp_after_first, (
        "DP table was recomputed even though tree and params did not change."
    )


def test_dp_reused_when_top_n_changes(counted):
    pc = PhytClust(tree=_load_tree())

    pc.run(top_n=1, plot_scores=False)
    dp_after_first = counted["dp"].calls

    pc.run(top_n=3, plot_scores=False)

    assert counted["dp"].calls == dp_after_first, (
        "Changing top_n should not trigger DP recomputation."
    )


# --------------------------------------------------------------------------- #
#  Scores cache: reused across top_n changes, invalidated on DP change        #
# --------------------------------------------------------------------------- #


def test_scores_reused_when_only_top_n_changes(counted):
    pc = PhytClust(tree=_load_tree())

    pc.run(top_n=1, plot_scores=False)
    scores_after_first = counted["scores"].calls

    pc.run(top_n=3, plot_scores=False)

    assert counted["scores"].calls == scores_after_first, (
        "Scores should be cached across calls that differ only in top_n."
    )


def test_scores_recomputed_when_max_k_changes(counted):
    pc = PhytClust(tree=_load_tree())

    pc.run(top_n=1, plot_scores=False)
    scores_after_first = counted["scores"].calls

    # max_k is part of the scores signature: changing it must invalidate.
    pc.run(top_n=1, max_k=4, plot_scores=False)

    assert counted["scores"].calls > scores_after_first, (
        "Scores should be recomputed when max_k changes."
    )


# --------------------------------------------------------------------------- #
#  Backtrack cache: persists across peak-mode runs                            #
# --------------------------------------------------------------------------- #


def test_backtrack_cache_persists_across_runs(counted):
    """Repeating the exact same run must not re-run backtrack at all."""
    pc = PhytClust(tree=_load_tree())

    pc.run(top_n=2, plot_scores=False)
    bt_after_first = counted["backtrack"].calls

    pc.run(top_n=2, plot_scores=False)

    assert counted["backtrack"].calls == bt_after_first, (
        "Repeated identical run() should not re-invoke backtrack; the "
        "per-k cluster cache should short-circuit it."
    )


def test_backtrack_only_for_new_peaks_when_top_n_grows(counted):
    """
    Going from top_n=1 to top_n=3 should backtrack only for the *new* peaks;
    the already-computed k stays cached.
    """
    pc = PhytClust(tree=_load_tree())

    pc.run(top_n=1, plot_scores=False)
    first_peaks = set(pc.peaks_by_rank or [])
    bt_after_first = counted["backtrack"].calls

    pc.run(top_n=3, plot_scores=False)
    second_peaks = set(pc.peaks_by_rank or [])

    new_peaks = second_peaks - first_peaks
    added_backtracks = counted["backtrack"].calls - bt_after_first

    assert added_backtracks == len(new_peaks), (
        f"Expected backtrack to run only for new peaks ({new_peaks}), "
        f"but it ran {added_backtracks} additional times."
    )


def test_backtrack_cache_hit_when_top_n_shrinks(counted):
    """Going from top_n=3 down to top_n=1 should not invoke backtrack at all."""
    pc = PhytClust(tree=_load_tree())

    pc.run(top_n=3, plot_scores=False)
    bt_after_first = counted["backtrack"].calls

    pc.run(top_n=1, plot_scores=False)

    assert counted["backtrack"].calls == bt_after_first, (
        "Shrinking top_n should be a pure cache hit for backtrack."
    )


# --------------------------------------------------------------------------- #
#  DP cache invalidation on DP-affecting parameter changes                    #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "attr,new_value",
    [
        ("min_cluster_size", 2),
        ("no_split_zero_length", True),
        ("use_branch_support", True),
        ("optimize_polytomies", False),
    ],
)
def test_dp_invalidated_on_param_change(counted, attr, new_value):
    """
    Mutating a DP-affecting parameter between runs must trigger a DP rebuild,
    not silently reuse the stale table.
    """
    pc = PhytClust(tree=_load_tree())

    pc.run(top_n=1, plot_scores=False)
    dp_after_first = counted["dp"].calls
    scores_after_first = counted["scores"].calls

    setattr(pc, attr, new_value)
    pc.run(top_n=1, plot_scores=False)

    assert counted["dp"].calls > dp_after_first, (
        f"Changing {attr} did not invalidate the DP cache."
    )
    assert counted["scores"].calls > scores_after_first, (
        f"Changing {attr} did not invalidate the scores cache."
    )


def test_dp_invalidated_on_outlier_threshold_change(counted):
    """Nested config fields on ``self.outlier`` must also invalidate."""
    pc = PhytClust(tree=_load_tree())

    pc.run(top_n=1, plot_scores=False)
    dp_after_first = counted["dp"].calls

    pc.outlier.size_threshold = 3
    pc.run(top_n=1, plot_scores=False)

    assert counted["dp"].calls > dp_after_first, (
        "Mutating pc.outlier.size_threshold in place did not invalidate "
        "the DP cache."
    )
