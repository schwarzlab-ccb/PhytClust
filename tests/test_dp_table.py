"""Regression checks for tree preparation and partition reconstruction."""

import pytest
from Bio.Phylo.BaseTree import Clade, Tree
from phytclust import PhytClust
from phytclust.algo.dp.table import _copy_tree, _assign_subtrees_to_cluster, backtrack
from phytclust.exceptions import InvalidClusteringError, InvalidKError


def test_outgroup_copy_preserves_precision_and_metadata():
    length = 0.123456789012345
    leaf = Clade(name="a", branch_length=length, confidence=88.123456789)
    leaf.extra = {"values": [1, 2]}
    tree = Tree(
        root=Clade(
            clades=[
                Clade(clades=[leaf, Clade(name="b", branch_length=2)]),
                Clade(name="out", branch_length=1),
            ]
        )
    )
    pc = PhytClust(tree, outgroup="out")
    copied = pc._tree_wo_outgroup.find_any(name="a")
    assert copied.branch_length == length
    assert copied.confidence == leaf.confidence
    assert copied.extra == leaf.extra
    copied.extra["values"].append(3)
    assert leaf.extra == {"values": [1, 2]}
    assert tree.find_any(name="out") is not None


def test_copy_handles_deep_topology_and_metadata_references():
    leaf = Clade(name="leaf")
    root = leaf
    for _ in range(1500):
        root = Clade(clades=[root])
    root.reference = leaf
    copied = _copy_tree(Tree(root=root))
    node = copied.root
    for _ in range(1500):
        node = node.clades[0]
    assert node.name == "leaf"
    assert node is not leaf
    assert copied.root.reference is node


def test_fixed_k_preserves_explicit_maximum():
    pc = PhytClust("(a:1,b:1,c:1,d:1);", k=2, max_k=4)
    assert pc.max_k == 4
    assert len(set(pc.get_clusters(2).values())) == 2


def test_duplicate_leaf_assignment_is_rejected():
    leaf = Clade(name="a")
    with pytest.raises(InvalidClusteringError, match="more than once"):
        _assign_subtrees_to_cluster([leaf, leaf], {}, 0)


def test_missing_leaf_assignment_is_rejected():
    pc = PhytClust("(a:1,b:1);")
    pc.get_clusters(1)
    pc.name_leaves_per_node = {pc.tree.root: [pc.tree.root.clades[0]]}
    with pytest.raises(InvalidClusteringError, match="every active leaf"):
        backtrack(pc, 1)


@pytest.mark.parametrize("k", [True, 1.5])
def test_backtracking_rejects_non_integer_counts(k):
    pc = PhytClust("(a:1,b:1);")
    pc.get_clusters(1)
    with pytest.raises(InvalidKError, match="positive integer"):
        backtrack(pc, k)
