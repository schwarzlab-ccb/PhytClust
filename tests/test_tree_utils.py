"""Tests for tree utility functions."""

from Bio.Phylo.BaseTree import Clade, Tree
import numpy as np
import pytest

from phytclust.exceptions import ConfigurationError, InvalidTreeError
from phytclust.utils.tree import (
    get_parent, get_pairwise_distances, rename_internal_nodes,
    find_all_min_indices, count_branches_in_clusters,
)


def test_get_parent_handles_root_and_descendants():
    leaf = Clade(name="leaf")
    internal = Clade(clades=[leaf])
    root = Clade(clades=[internal])
    tree = Tree(root=root)

    assert get_parent(tree, internal) is root
    assert get_parent(tree, leaf) is internal
    assert get_parent(tree, root) is None
    assert get_parent(tree, Clade(name="outside")) is None


@pytest.mark.parametrize("mode,names", [
    ("terminals", ["A", "B"]),
    ("nonterminals", ["X.", "Y"]),
    ("all", ["A", "B", "X.", "Y"]),
])
def test_distances_honor_mode_with_subtree(mode, names):
    subtree = Clade(name="X.", clades=[
        Clade(name="A", branch_length=2),
        Clade(name="Y", branch_length=3, clades=[Clade(name="B", branch_length=4)]),
    ])
    tree = Tree(root=Clade(clades=[subtree, Clade(name="outside")]))
    frame = get_pairwise_distances(tree, mode, True, "X.")
    assert list(frame.index) == names
    assert list(frame.columns) == names
    for a in names:
        for b in names:
            assert frame.loc[a, b] == tree.distance(a, b)


def test_subtree_names_are_literal_and_unique_and_mode_is_validated():
    tree = Tree(root=Clade(clades=[Clade(name="Xa"), Clade(name="same"), Clade(name="same")]))
    with pytest.raises(InvalidTreeError, match="No node"):
        get_pairwise_distances(tree, mrca="X.")
    with pytest.raises(InvalidTreeError, match="ambiguous"):
        get_pairwise_distances(tree, mrca="same")
    with pytest.raises(ConfigurationError):
        get_pairwise_distances(tree, mode="bad", mrca="Xa")


def test_generated_dataframe_labels_do_not_collide():
    tree = Tree(root=Clade(clades=[Clade(), Clade(name="node_0")]))
    frame = get_pairwise_distances(tree, as_dataframe=True)
    assert list(frame.index) == ["node_0_", "node_0"]


def test_internal_names_do_not_collide_with_leaf_names():
    leaf = Clade(name="internal_1")
    tree = Tree(root=Clade(clades=[Clade(clades=[leaf]), Clade(name="internal_3")]))
    rename_internal_nodes(tree)
    names = [node.name for node in tree.find_clades()]
    assert len(names) == len(set(names))
    assert leaf.name == "internal_1"


def test_minimum_handles_ties_empty_arrays_and_nan():
    assert find_all_min_indices([2, 1, 1]) == ([1, 2], 1)
    assert find_all_min_indices(np.array([])) == ([], float("inf"))
    assert find_all_min_indices([float("inf"), float("inf")]) == ([0, 1], float("inf"))
    with pytest.raises(ValueError, match="NaN"):
        find_all_min_indices([1, float("nan")])


def test_binary_branch_estimate_excludes_singletons():
    assert count_branches_in_clusters({0: [], 1: [Clade()], 2: [Clade(), Clade(), Clade()]}) == 4
