"""Representative species must come from their assigned cluster."""

from io import StringIO

from Bio import Phylo

from phytclust.selection.representatives import select_representative_species


def _tree(newick):
    return Phylo.read(StringIO(newick), "newick")


def test_representative_excludes_other_leaves_below_mrca():
    tree = _tree("((a:10,b:1):1,c:10);")
    cluster = {"a": 0, "c": 0}

    assert select_representative_species(tree, cluster, strategy="central") == ["c"]
    assert select_representative_species(tree, cluster, strategy="divergent") == ["a"]
    assert select_representative_species(tree, cluster, strategy="medoid") == ["a"]


def test_medoid_uses_only_cluster_members_for_distance_sum():
    tree = _tree("((a:10,b:1):1,(c:2,d:3):1);")

    assert select_representative_species(
        tree, {"a": 0, "c": 0, "d": 0}, strategy="medoid"
    ) == ["c"]


import pytest
from phytclust.exceptions import ConfigurationError, InvalidTreeError
from phytclust.selection.representatives import (
    compute_species_distance,
    rank_terminal_nodes,
    maximize_pd,
)


@pytest.mark.parametrize("options", [{"mode": "bad"}, {"distance_ref": "bad"}])
def test_singleton_validates_options(options):
    with pytest.raises(ConfigurationError):
        select_representative_species(_tree("(a:1,b:2);"), {"a": 0}, **options)


def test_missing_singleton_rejected():
    with pytest.raises(ConfigurationError, match="absent"):
        select_representative_species(_tree("(a:1,b:2);"), {"missing": 0})


@pytest.mark.parametrize("count", [0, -1, 1.5, True])
def test_species_count_requires_positive_integer(count):
    with pytest.raises(ConfigurationError, match="num_species"):
        rank_terminal_nodes(_tree("(a:1,b:2);"), num_species=count)


def test_names_match_literally_and_duplicates_rejected():
    tree = _tree("('a.':1,ab:2);")
    assert compute_species_distance(tree, ["a.", "ab"], distance_ref="mrca") == {
        "a.": 1,
        "ab": 2,
    }
    with pytest.raises(ConfigurationError, match="absent"):
        compute_species_distance(tree, ["a.*"])
    with pytest.raises(InvalidTreeError, match="duplicated"):
        rank_terminal_nodes(_tree("(a:1,a:2);"))


def test_distance_sums_match_pairwise_reference():
    tree = _tree("((a:1,b:2):3,c:4);")
    expected = {
        leaf.name: sum(tree.distance(leaf, other) for other in tree.get_terminals())
        for leaf in tree.get_terminals()
    }
    assert compute_species_distance(tree, list(expected)) == expected


def test_deep_tree_uses_iterative_traversal():
    from Bio.Phylo.BaseTree import Clade, Tree

    root = Clade(name="a", branch_length=1)
    for _ in range(1500):
        root = Clade(clades=[root], branch_length=1)
    tree = Tree(root=Clade(clades=[root, Clade(name="b", branch_length=2)]))
    assert rank_terminal_nodes(tree) == [("a", 1503), ("b", 1503)]
    assert select_representative_species(tree, {"a": 0, "b": 0}) == ["b"]


def test_summary_describes_distance_scores():
    summary = maximize_pd(_tree("(a:1,b:2);"), 1)[0][0]
    assert "distance scores" in summary
    assert "Maximizing PD" not in summary
