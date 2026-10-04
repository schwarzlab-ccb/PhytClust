import itertools
from io import StringIO

import pytest
from Bio import Phylo

from phytclust import PhytClust
from phytclust.algo.dp import enumerate_tied_optima
from phytclust.exceptions import ConfigurationError

# Degree-4 polytomy of four identical cherries: every way of pairing the
# children costs the same, so soft mode has genuine ties at k=2 and k=3.
BALANCED_POLYTOMY = "((A:1,B:1):1,(C:1,D:1):1,(E:1,F:1):1,(G:1,H:1):1);"


def _tree(newick):
    return Phylo.read(StringIO(newick), "newick")


def _signature(cluster_map):
    groups = {}
    for terminal, cluster_id in cluster_map.items():
        groups.setdefault(cluster_id, set()).add(terminal.name)
    return frozenset(frozenset(names) for names in groups.values())


def _pc(newick=BALANCED_POLYTOMY, **kwargs):
    kwargs.setdefault("max_k", 8)
    kwargs.setdefault("save_tied_optima", True)
    return PhytClust(tree=_tree(newick), **kwargs)


def test_disabled_by_default():
    pc = PhytClust(tree=_tree(BALANCED_POLYTOMY), max_k=8, polytomy_mode="soft")
    pc.get_clusters(2)
    assert not getattr(pc, "tied_optima", {})


def test_soft_ties_are_recorded():
    pc = _pc(polytomy_mode="soft")
    pc.get_clusters(2)

    result = pc.tied_optima[2]
    assert result["k"] == 2
    assert not result["truncated"]
    # {AB|CD}{EF|GH}, {AB|EF}{CD|GH}, {AB|GH}{CD|EF}, and the same three read
    # the other way round collapse to three distinct pairings... plus the
    # split that leaves one cherry alone is more expensive, so exactly the
    # pairings tie.
    assert result["n_solutions"] == len(result["solutions"])
    assert result["n_solutions"] > 1

    signatures = {_signature(s) for s in result["solutions"]}
    assert len(signatures) == result["n_solutions"], "solutions must be distinct"
    for signature in signatures:
        assert len(signature) == 2
        assert frozenset().union(*signature) == frozenset("ABCDEFGH")


def test_hard_mode_has_no_subset_merges():
    """Hard mode cannot group a proper subset of children, so k=2 is infeasible."""
    pc = _pc(polytomy_mode="hard")
    with pytest.raises(Exception):
        pc.get_clusters(2)


def test_first_solution_matches_backtrack():
    pc = _pc(polytomy_mode="soft")
    clusters = pc.get_clusters(3)
    recorded = pc.tied_optima[3]["solutions"]
    assert _signature(clusters) == _signature(recorded[0])


@pytest.mark.parametrize("cap", [1, 2])
def test_capped_results_keep_selected_partition_first(cap):
    pc = _pc(polytomy_mode="soft", max_tied_optima=cap)
    chosen = pc.get_clusters(3)
    result = pc.tied_optima[3]
    assert _signature(result["solutions"][0]) == _signature(chosen)
    assert result["n_solutions"] <= cap
    assert result["truncated"]


def test_binary_cost_tie_selects_even_cluster_allocation():
    pc = _pc("(((A:0,B:0):0,(C:0,D:0):0):0,((E:0,F:0):0,(G:0,H:0):0):0);")
    # Keep zero lengths: preparation substitutes unit lengths for an all-zero tree.
    for node in pc.tree.find_clades():
        node.branch_length = 0
    # One positive root length leaves all edges below it at zero after normalization.
    pc.tree.root.clades[0].branch_length = 1
    pc.tree.root.clades[1].branch_length = 1
    chosen = pc.get_clusters(4)
    left = pc.tree.root.clades[0].get_terminals()
    right = pc.tree.root.clades[1].get_terminals()
    assert len({chosen[leaf] for leaf in left}) == 2
    assert len({chosen[leaf] for leaf in right}) == 2


def test_binary_tree_unique_optimum():
    pc = _pc("(((A:1,B:1):5,C:3):2,(D:1,E:4):1);", max_k=4)
    for k in (1, 2, 3):
        pc.get_clusters(k)
        assert pc.tied_optima[k]["n_solutions"] >= 1


def test_cap_truncates_and_flags():
    pc = _pc(polytomy_mode="soft", max_tied_optima=2)
    pc.get_clusters(3)
    result = pc.tied_optima[3]
    assert result["truncated"]
    assert result["n_solutions"] <= 2


def test_prefer_fewer_is_rejected():
    pc = _pc(polytomy_mode="soft")
    pc.outlier.size_threshold = 2
    pc.outlier.prefer_fewer = True
    with pytest.raises(ConfigurationError):
        enumerate_tied_optima(pc, 2)


def test_matches_brute_force_on_polytomy():
    """Independent enumeration of every antichain partition, soft candidates
    including subsets of >=2 polytomy children joined at a zero-length node."""
    newick = BALANCED_POLYTOMY

    def clade_cost(node):
        return sum(
            sum(step.branch_length or 0.0 for step in node.get_path(t))
            for t in node.get_terminals()
        )

    def leafset(node):
        return frozenset(t.name for t in node.get_terminals())

    tree = _tree(newick)
    candidates = [(leafset(n), clade_cost(n)) for n in tree.find_clades()]
    for node in tree.find_clades():
        children = list(node.clades)
        if len(children) <= 2:
            continue
        for size in range(2, len(children)):
            for combo in itertools.combinations(children, size):
                candidates.append(
                    (
                        frozenset().union(*[leafset(c) for c in combo]),
                        sum(
                            clade_cost(c)
                            + len(leafset(c)) * (c.branch_length or 0.0)
                            for c in combo
                        ),
                    )
                )

    all_leaves = leafset(tree.root)
    for k in (2, 3, 5):
        best, best_parts = float("inf"), []
        for combo in itertools.combinations(candidates, k):
            sets = [c[0] for c in combo]
            if sum(len(s) for s in sets) != len(all_leaves):
                continue
            if frozenset().union(*sets) != all_leaves:
                continue
            cost = sum(c[1] for c in combo)
            if cost < best - 1e-9:
                best, best_parts = cost, [frozenset(sets)]
            elif abs(cost - best) <= 1e-9:
                best_parts.append(frozenset(sets))

        pc = _pc(polytomy_mode="soft", max_tied_optima=10_000)
        pc.get_clusters(k)
        result = pc.tied_optima[k]

        assert result["score"] == pytest.approx(best)
        assert {_signature(s) for s in result["solutions"]} == set(best_parts)
