"""Check balanced child allocations among tied polytomy partitions."""

from itertools import product
import pytest
from Bio.Phylo.BaseTree import Clade, Tree
from phytclust import PhytClust


def child_subtree(prefix):
    leaves = [Clade(name=f"{prefix}{index}", branch_length=0) for index in range(4)]
    return Clade(
        branch_length=1,
        clades=[
            Clade(branch_length=0, clades=leaves[:2]),
            Clade(branch_length=0, clades=leaves[2:]),
        ],
    )


@pytest.mark.parametrize("mode", ["hard", "soft"])
@pytest.mark.parametrize("cluster_count", [6, 7, 8])
@pytest.mark.parametrize("reverse_children", [False, True])
@pytest.mark.parametrize("use_float32", [False, True])
def test_polytomy_chooses_even_child_cluster_counts(
    mode, cluster_count, reverse_children, use_float32
):
    children = [child_subtree(prefix) for prefix in ["A", "B", "C"]]
    if reverse_children:
        children.reverse()
    pc = PhytClust(
        Tree(root=Clade(clades=children)),
        polytomy_mode=mode,
        max_k=9,
        dp_float32=use_float32,
        save_tied_optima=True,
        max_tied_optima=1,
    )
    chosen = pc.get_clusters(cluster_count)
    allocations = []
    child_cluster_ids = []
    for child in children:
        cluster_ids = {chosen[leaf] for leaf in child.get_terminals()}
        allocations.append(len(cluster_ids))
        child_cluster_ids.append(cluster_ids)
    # These zero-cost optima keep clusters inside each child subtree.
    assert sum(allocations) == cluster_count
    assert not any(
        a & b
        for index, a in enumerate(child_cluster_ids)
        for b in child_cluster_ids[index + 1 :]
    )
    feasible_allocations = [
        counts
        for counts in product(range(1, 5), repeat=3)
        if sum(counts) == cluster_count
    ]
    assert sum(count**2 for count in allocations) == min(
        sum(count**2 for count in counts) for counts in feasible_allocations
    )
    assert max(allocations) - min(allocations) <= 1
    result = pc.tied_optima[cluster_count]
    assert result["n_solutions"] == 1
    assert result["solutions"][0] == chosen
