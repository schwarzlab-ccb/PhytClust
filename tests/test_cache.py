import pathlib

import pytest
from Bio import Phylo

from phytclust.algo import core as core_mod
from phytclust.algo import dp as dp_mod
from phytclust.algo import scoring as scoring_mod
from phytclust.algo.core import PhytClust

TREE_PATH = pathlib.Path(__file__).parent.parent / "examples" / "sample_tree.nwk"


def test_tree_edit_refreshes_leaf_maps_and_cluster_membership():
    from Bio.Phylo.BaseTree import Clade
    pc = PhytClust("(a:1,b:2);")
    pc.get_clusters(1)
    added = Clade(name="c", branch_length=3)
    pc.tree.root.clades.append(added)
    result = pc.get_clusters(1)
    assert pc.num_terminals == 3
    assert set(result) == set(pc.tree.get_terminals())
    assert any(leaf is added for leaf in pc.name_leaves_per_node[pc.tree.root])


def test_tree_edit_refreshes_pruned_outgroup_copy():
    pc = PhytClust("((a:1,b:2):3,out:1);", outgroup="out")
    pc.get_clusters(1)
    pc.tree.find_any(name="a").branch_length = 7
    pc.get_clusters(1)
    assert pc._tree_wo_outgroup.find_any(name="a").branch_length == 7


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


# --------------------------------------------------------------------------- #
#  Swapping subtrees is the same tree                                         #
# --------------------------------------------------------------------------- #

POLYTOMY_TREE = (
    "(((A:1,B:2):1,(C:1,D:3):1,(J:2,K:1):4):1,"
    "((E:1,F:5):2,(G:4,H:1):1):7,I:2);"
)


def _signature(cmap):
    groups = {}
    for terminal, cid in cmap.items():
        groups.setdefault(cid, set()).add(terminal.name)
    return frozenset(frozenset(g) for g in groups.values())


def _rotate(tree, seed):
    import random

    rng = random.Random(seed)
    for clade in tree.find_clades():
        rng.shuffle(clade.clades)
    return tree


def _partitions(pc, ks):
    return {k: _signature(pc.get_clusters(k)) for k in ks}


@pytest.mark.parametrize("newick", [None, POLYTOMY_TREE])
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_dp_reused_after_subtrees_swapped(counted, newick, seed):
    """Rotating children in place keeps the DP and backtracks correctly."""
    from io import StringIO

    def fresh_tree():
        return Phylo.read(StringIO(newick), "newick") if newick else _load_tree()

    ks = range(1, 8)
    pc = PhytClust(tree=fresh_tree(), max_k=8)
    pc._ensure_dp(required_cap=8)
    dp_after_first = counted["dp"].calls

    _rotate(pc.tree, seed)
    pc.clusters = {}
    rotated = _partitions(pc, ks)

    assert counted["dp"].calls == dp_after_first, (
        "Swapping subtrees recomputed the DP although the tree is the same."
    )

    reference = PhytClust(tree=fresh_tree(), max_k=8)
    reference._ensure_dp(required_cap=8)
    assert rotated == _partitions(reference, ks)


def test_dp_invalidated_on_branch_length_change(counted):
    """The order-free fingerprint still sees real edits."""
    pc = PhytClust(tree=_load_tree())
    pc.run(top_n=1, plot_scores=False)
    dp_after_first = counted["dp"].calls

    leaf = pc.tree.get_terminals()[0]
    leaf.branch_length = (leaf.branch_length or 0.0) + 1.0
    pc.run(top_n=1, plot_scores=False)

    assert counted["dp"].calls > dp_after_first


def test_gui_reuses_instance_for_rotated_newick():
    pytest.importorskip("fastapi")
    from io import StringIO

    from phytclust.gui import api

    tree = Phylo.read(StringIO(POLYTOMY_TREE), "newick")
    buf = StringIO()
    Phylo.write(_rotate(tree, 3), buf, "newick")
    rotated_newick = buf.getvalue().strip()
    assert rotated_newick != POLYTOMY_TREE

    api.LAST_PC = api.LAST_CONSTRUCTION_KEY = None
    first = api._run_phytclust(api.PhytclustRequest(newick=POLYTOMY_TREE))
    pc = api.LAST_PC
    second = api._run_phytclust(api.PhytclustRequest(newick=rotated_newick))

    assert api.LAST_PC is pc
    assert second["notes"] == []
    assert [t.name for t in pc.tree.get_terminals()] == [
        t.name for t in Phylo.read(StringIO(rotated_newick), "newick").get_terminals()
    ]
    assert first["clusters"] and second["clusters"]


def test_dp_cap_doubles_when_k_grows_one_at_a_time(counted):
    """Asking for k = 2, 3, ... reruns the DP O(log k) times, not once per k."""
    pc = PhytClust(tree=_load_tree(), max_k=None)
    for k in range(2, pc.num_terminals + 1):
        pc.get_clusters(k)
    assert counted["dp"].calls <= pc.num_terminals.bit_length() + 1


def test_user_max_k_builds_full_range_once(counted):
    """A max_k the caller set is built on the first call, whatever k asked."""
    pc = PhytClust(tree=_load_tree(), max_k=6)
    for k in range(2, 7):
        pc.get_clusters(k)
    assert counted["dp"].calls == 1
    assert pc._dp_cap == 6
