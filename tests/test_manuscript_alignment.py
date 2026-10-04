"""Behaviour the manuscript's Methods specify, checked directly.

Each test names the part of the Methods it pins down.
"""

import itertools
import pathlib
import random
import subprocess
import sys
from io import StringIO

import numpy as np
import pytest
from Bio import Phylo

from phytclust import PeakConfig, PhytClust
from phytclust.algo import core as core_mod
from phytclust.metrics.indices import cluster_alpha

SAMPLE_TREE = pathlib.Path(__file__).parent.parent / "examples" / "sample_tree.nwk"
CLI_TREE = pathlib.Path(__file__).parent / "test_tree.nwk"


def _tree(newick):
    return Phylo.read(StringIO(newick), "newick")


def _signature(cmap):
    groups = {}
    for terminal, cid in cmap.items():
        groups.setdefault(cid, set()).add(terminal.name)
    return frozenset(frozenset(g) for g in groups.values())


def _random_newick(seed, n_leaves):
    rng = random.Random(seed)
    names = [f"t{i}" for i in range(n_leaves)]

    def build(n):
        length = rng.uniform(0.05, 3.0)
        if n == 1:
            return f"{names.pop()}:{length:.6f}"
        left = rng.randint(1, n - 1)
        return f"({build(left)},{build(n - left)}):{length:.6f}"

    return build(n_leaves) + ";"


def _scaled(newick, factor):
    tree = _tree(newick)
    for clade in tree.find_clades():
        if clade.branch_length is not None:
            clade.branch_length *= factor
    return tree


def test_k2_is_a_candidate_by_default():
    """Scoring selects k* from k in {2, ..., N-1}."""
    assert PeakConfig().exclude_k2 is False
    newick = "(((A:1,B:1):1,(C:1,D:1):1):40,((E:1,F:1):1,(G:1,H:1):1):40);"

    result = PhytClust(tree=_tree(newick)).run(top_n=1, plot_scores=False)
    assert result["k_values"] == [2]

    excluded = PhytClust(tree=_tree(newick)).run(
        top_n=1, plot_scores=False, peak_config=PeakConfig(exclude_k2=True)
    )
    assert 2 not in excluded["k_values"]


@pytest.mark.parametrize("factor", [2.0**-30, 1e-9, 1e-12, 1e6])
@pytest.mark.parametrize("seed", [1, 2, 3])
def test_invariant_to_global_rescaling(seed, factor):
    """Multiplying every branch length by a positive constant changes nothing."""
    newick = _random_newick(seed, 60)
    base = PhytClust(tree=_tree(newick))
    scaled = PhytClust(tree=_scaled(newick, factor))

    base_result = base.run(top_n=3, plot_scores=False)
    scaled_result = scaled.run(top_n=3, plot_scores=False)
    assert base_result["k_values"] == scaled_result["k_values"]
    np.testing.assert_allclose(base.scores, scaled.scores, rtol=1e-9)

    for k in (2, 5, 10, 25, 50):
        assert _signature(base.get_clusters(k)) == _signature(scaled.get_clusters(k))


@pytest.mark.parametrize("factor", [2.0**-30, 1e-9, 1e-12, 1e6])
def test_rescaling_keeps_tie_breaks_on_integer_tree(factor):
    base = PhytClust(tree=_tree(SAMPLE_TREE.read_text()))
    scaled = PhytClust(tree=_scaled(SAMPLE_TREE.read_text(), factor))
    for k in range(1, base.num_terminals + 1):
        assert _signature(base.get_clusters(k)) == _signature(scaled.get_clusters(k))


def test_alpha_ignores_the_root_branch():
    """Eq. 18 averages edges; the root has no incoming edge."""
    tree = _tree("((A:1,B:1):1,(C:1,D:1):1):5;")
    pc = PhytClust(tree=tree)
    info = cluster_alpha(pc.tree, pc.get_clusters(2))
    assert info["n_extra_nodes"] == 2
    assert info["n_intra_nodes"] == 4
    assert info["alpha"] == pytest.approx(1.0)


def test_alpha_after_outgroup_pruning_ignores_new_root_length():
    tree = _tree("(((A:1,B:1):1,(C:1,D:1):1):3,O:1);")
    pc = PhytClust(tree=tree, outgroup="O")
    result = pc.run(k=2, plot_scores=False)
    assert result["alpha_details"][0]["alpha"] == pytest.approx(1.0)


def _hard_allocations(pc, k):
    """Brute force: every allocation (m_1..m_q) at the root, with its cost."""
    rows = [pc.raw_dp_table[pc.node_to_id[c]] for c in pc.tree.root.clades]
    out = []
    for alloc in itertools.product(*[range(1, len(r)) for r in rows]):
        if sum(alloc) != k:
            continue
        cost = sum(float(r[m - 1]) for r, m in zip(rows, alloc))
        if np.isfinite(cost):
            out.append((cost, sum(m * m for m in alloc)))
    return out


@pytest.mark.parametrize("seed", range(12))
def test_hard_polytomy_ties_minimise_sum_of_squares(seed):
    """Among cost-tied allocations at a node, the smallest sum of m_j^2 wins."""
    rng = random.Random(seed)
    names = [f"t{i}" for i in range(60)]

    def sub(n):
        length = rng.choice([1, 1, 2])
        if n == 1:
            return f"{names.pop()}:{length}"
        left = rng.randint(1, n - 1)
        return f"({sub(left)},{sub(n - left)}):{length}"

    sizes = [rng.choice([1, 2, 3, 5, 9]) for _ in range(rng.randint(3, 5))]
    newick = "(" + ",".join(sub(s) for s in sizes) + ");"
    pc = PhytClust(
        tree=_tree(newick), polytomy_mode="hard", preserve_dp_tables=True
    )
    n = pc.num_terminals
    pc.get_clusters(n)
    leafsets = [set(pc.name_leaves_per_node[c]) for c in pc.tree.root.clades]

    for k in range(len(sizes), n + 1):
        cmap = pc.get_clusters(k)
        chosen = [len({cmap[t] for t in leaves}) for leaves in leafsets]
        candidates = _hard_allocations(pc, k)
        best = min(c for c, _ in candidates)
        tied = [sq for c, sq in candidates if abs(c - best) <= 1e-9 * max(1.0, best)]
        assert sum(m * m for m in chosen) == min(tied)


def _brute_force_costs(tree):
    """Optimal cost for every k by enumerating monophyletic partitions."""
    def depth_sum(node):
        return sum(
            sum(step.branch_length or 0.0 for step in node.get_path(t))
            for t in node.get_terminals()
        )

    def partitions(node):
        options = [[depth_sum(node)]]
        if node.clades:
            child_sets = [partitions(c) for c in node.clades]
            for combo in itertools.product(*child_sets):
                options.append([x for part in combo for x in part])
        return options

    best = {}
    for part in partitions(tree.root):
        k = len(part)
        best[k] = min(best.get(k, np.inf), sum(part))
    return best


@pytest.mark.parametrize(
    "newick",
    [
        SAMPLE_TREE.read_text(),
        "((A:1,B:2,C:1):1,(D:1,(E:1,F:3):1):2,G:4);",
    ],
)
def test_hard_dp_cost_matches_brute_force(newick):
    expected = _brute_force_costs(_tree(newick))
    pc = PhytClust(tree=_tree(newick), polytomy_mode="hard")
    pc.get_clusters(pc.num_terminals)
    root_row = pc.raw_dp_table[pc.node_to_id[pc.tree.root]]
    for k in range(1, pc.num_terminals + 1):
        want = expected.get(k, np.inf)
        assert root_row[k - 1] == pytest.approx(want)


def test_max_k_limit_is_honoured_and_not_sticky():
    pc = PhytClust(tree=_tree(SAMPLE_TREE.read_text()))
    pc.run(top_n=1, max_k_limit=0.5, plot_scores=False)
    assert pc.max_k == 4
    pc.run(top_n=1, max_k=4, plot_scores=False)
    pc.run(top_n=1, plot_scores=False)
    assert pc.max_k == 8


def test_explicit_max_k_on_the_instance_is_kept():
    pc = PhytClust(tree=_tree(SAMPLE_TREE.read_text()), max_k=5)
    pc.run(top_n=1, plot_scores=False)
    pc.run(top_n=1, plot_scores=False)
    assert pc.max_k == 5


def test_tied_optima_backtracks_once(monkeypatch):
    calls = {"n": 0}
    real = core_mod.backtrack

    def counted(*args, **kwargs):
        calls["n"] += 1
        return real(*args, **kwargs)

    monkeypatch.setattr(core_mod, "backtrack", counted)
    pc = PhytClust(
        tree=_tree("((A:1,B:1):1,(C:1,D:1):1,(E:1,F:1):1,(G:1,H:1):1);"),
        max_k=8,
        save_tied_optima=True,
    )
    pc.get_clusters(2)
    assert calls["n"] == 1
    assert pc.tied_optima[2]["n_solutions"] > 1


def test_cli_accepts_trivial_k(tmp_path):
    result = subprocess.run(
        [sys.executable, "-m", "phytclust.cli", str(CLI_TREE), "-k", "1",
         "--out-dir", str(tmp_path)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr


def test_resolution_bins_start_at_one():
    """Clade-level bins span 1..K, so k=3 and the next coarse level separate."""
    pc = PhytClust(tree=_tree(_random_newick(7, 300)))
    pc.run(by_resolution=True, num_bins=5, plot_scores=False)
    bins = pc.bin_ranges_current
    assert bins[0][0] == 1
    assert bins[0][1] <= 4
    for (_, hi), (lo, _) in zip(bins, bins[1:]):
        assert lo == hi + 1
