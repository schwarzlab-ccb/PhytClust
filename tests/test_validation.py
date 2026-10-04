"""Tests for validation module."""

import pathlib
from collections.abc import Mapping
import pytest
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

from phytclust.validation import (
    validate_and_set_outgroup,
    prune_outgroup,
    ensure_branch_lengths,
    merge_single_child_clades,
)
from phytclust.exceptions import InvalidTreeError

TREE_PATH = pathlib.Path(__file__).parent.parent / "examples" / "sample_tree.nwk"


@pytest.mark.parametrize("parent_confidence", [60, None])
def test_unary_collapse_retains_parent_confidence_and_warns(caplog, parent_confidence):
    leaves = [Clade(name="A"), Clade(name="B")]
    child = Clade(name="child", confidence=95, branch_length=2, clades=leaves)
    parent = Clade(name="parent", confidence=parent_confidence, branch_length=1,
                   clades=[child])
    tree = Tree(root=parent)

    merge_single_child_clades(tree)

    assert tree.root is parent
    assert parent.confidence == parent_confidence
    assert parent.branch_length == 3
    assert parent.clades == leaves
    assert len(caplog.records) == 1
    assert caplog.records[0].levelname == "WARNING"
    assert "retaining parent confidence" in caplog.text
    assert "discarding child confidence 95" in caplog.text


def _load_tree() -> Tree:
    """Load test tree."""
    return Phylo.read(TREE_PATH, "newick")


@pytest.mark.parametrize("nested", [False, True])
@pytest.mark.parametrize("parent_confidence,child_confidence", [(60, 95), (None, 95), (60, 60), (60, None)])
def test_pruning_retains_parent_confidence(caplog, nested, parent_confidence, child_confidence):
    child = Clade(name="remaining", confidence=child_confidence, branch_length=2,
                  clades=[Clade(name="A", branch_length=1), Clade(name="B", branch_length=1)])
    parent = Clade(name="parent", confidence=parent_confidence, branch_length=3,
                   clades=[Clade(name="outgroup", branch_length=1), child])
    root = Clade(clades=[parent, Clade(name="C", branch_length=4)]) if nested else parent
    tree = Tree(root=root)
    pairs = [("A", "B")] + ([("A", "C"), ("B", "C")] if nested else [])
    distances = [tree.distance(*pair) for pair in pairs]

    leaves, counts = prune_outgroup(tree, "outgroup")

    assert child.confidence == parent_confidence
    assert tree.root is (root if nested else child)
    assert [tree.distance(*pair) for pair in pairs] == pytest.approx(distances)
    assert counts[tree.root] == (3 if nested else 2)
    assert "outgroup" not in {leaf.name for leaf in leaves[tree.root]}
    should_warn = child_confidence is not None and child_confidence != parent_confidence
    assert len(caplog.records) == int(should_warn)
    if should_warn:
        assert caplog.records[0].levelname == "WARNING"
        assert "retaining parent confidence" in caplog.text


class TestValidateAndSetOutgroup:
    """Test outgroup validation and setting."""

    def test_with_valid_outgroup(self):
        """Test with a valid outgroup name."""
        tree = _load_tree()
        terminals = [t.name for t in tree.get_terminals()]
        if terminals:
            outgroup = terminals[0]
            result_tree, result_outgroup = validate_and_set_outgroup(tree, outgroup)

            assert isinstance(result_tree, Tree)
            assert result_outgroup == outgroup
            assert result_tree.root is not None

    def test_with_none_outgroup(self):
        """Test with None outgroup (tree should remain unchanged)."""
        tree = _load_tree()
        result_tree, result_outgroup = validate_and_set_outgroup(tree, None)

        assert result_tree is not None
        assert result_outgroup is None

    def test_with_invalid_outgroup_raises_error(self):
        """Test that invalid outgroup name raises InvalidTreeError."""
        tree = _load_tree()
        with pytest.raises(InvalidTreeError):
            validate_and_set_outgroup(tree, "nonexistent_taxon_xyz")


class TestPruneOutgroup:
    """Test outgroup pruning."""

    def test_prune_without_outgroup_returns_all_nodes(self):
        """Test pruning without outgroup returns all nodes."""
        tree = _load_tree()
        node_terminals, terminal_count = prune_outgroup(tree, None)

        assert isinstance(node_terminals, Mapping)
        assert isinstance(terminal_count, dict)
        assert len(node_terminals) > 0
        assert len(terminal_count) == len(node_terminals)

        for node, count in terminal_count.items():
            assert count >= 1
            assert len(node_terminals[node]) == count

    def test_prune_with_outgroup_reduces_nodes(self):
        """Test that pruning with outgroup reduces node count."""
        tree = _load_tree()
        terminals = [t.name for t in tree.get_terminals()]

        if terminals:
            node_terminals_no_prune, _ = prune_outgroup(tree, None)

            node_terminals_with_prune, _ = prune_outgroup(tree, terminals[0])
            assert len(node_terminals_with_prune) < len(node_terminals_no_prune)
            assert terminals[0] not in {n.name for n in tree.get_terminals()}



class TestEnsureBranchLengths:
    """Test branch length validation."""

    def test_ensure_branch_lengths_adds_missing_lengths(self):
        """Test that ensure_branch_lengths adds missing branch lengths."""
        tree = _load_tree()

        # Set some branch lengths to None to test
        for clade in tree.find_clades():
            if clade != tree.root:
                clade.branch_length = None

        ensure_branch_lengths(tree)

        # After ensuring, all non-root clades should have branch lengths
        for clade in tree.find_clades():
            if clade != tree.root:
                assert clade.branch_length is not None
                assert clade.branch_length >= 0

    def test_ensure_branch_lengths_preserves_existing(self):
        """Test that existing branch lengths are preserved."""
        tree = _load_tree()

        # Store original branch lengths
        original_lengths = {}
        for clade in tree.find_clades():
            if clade != tree.root:
                clade.branch_length = 1.5
                original_lengths[id(clade)] = 1.5

        ensure_branch_lengths(tree)

        # Check that set lengths are preserved
        for clade in tree.find_clades():
            if id(clade) in original_lengths:
                assert clade.branch_length == original_lengths[id(clade)]


def _newick(text):
    from io import StringIO

    return Phylo.read(StringIO(text), 'newick')


@pytest.mark.parametrize('length', [float('nan'), float('inf'), -float('inf'), 'bad'])
@pytest.mark.parametrize('at_root', [False, True])
def test_invalid_branch_lengths_are_rejected(length, at_root):
    tree = _newick('(a:1,b:1);')
    node = tree.root if at_root else tree.root.clades[0]
    node.branch_length = length
    with pytest.raises(InvalidTreeError, match='finite number'):
        ensure_branch_lengths(tree)


def test_negative_lengths_are_clamped_before_unary_collapse():
    tree, _ = validate_and_set_outgroup(_newick('((a:2):-1,b:1);'), None)
    assert tree.find_any(name='a').branch_length == 2
    assert tree.distance('a', 'b') == 3


@pytest.mark.parametrize('newick', ['(a,b,c);', '(a:1,b,c:2);', '(a:0,b:0,c:0);'])
def test_midpoint_rooting_normalizes_missing_lengths(newick):
    tree, _ = validate_and_set_outgroup(_newick(newick), None, 'midpoint')
    assert {n.name for n in tree.get_terminals()} == {'a', 'b', 'c'}
    depths = tree.depths()
    assert max(depths[n] for n in tree.get_terminals()) > 0


def test_midpoint_preserves_pairwise_distances():
    tree = _newick('((a:1,b:3):2,c:4);')
    pairs = [('a', 'b'), ('a', 'c'), ('b', 'c')]
    before = [tree.distance(*pair) for pair in pairs]
    validate_and_set_outgroup(tree, None, 'midpoint')
    assert [tree.distance(*pair) for pair in pairs] == pytest.approx(before)


@pytest.mark.parametrize('option', ['outgroup', 'root_taxon'])
def test_requested_names_are_literal_and_unique(option):
    tree = _newick('(a:1,ab:1,c:1);')
    with pytest.raises(InvalidTreeError, match='not found'):
        validate_and_set_outgroup(tree, **{option: 'a.'}, **({'outgroup': None} if option == 'root_taxon' else {}))
    tree = _newick('(a:1,a:1,c:1);')
    with pytest.raises(InvalidTreeError, match='ambiguous'):
        validate_and_set_outgroup(tree, **{option: 'a'}, **({'outgroup': None} if option == 'root_taxon' else {}))
    tree = _newick("('a.':1,ab:1,c:1);")
    validate_and_set_outgroup(tree, **{option: 'a.'}, **({'outgroup': None} if option == 'root_taxon' else {}))


@pytest.mark.parametrize('newick', [
    '(((a:1,b:1)X:1,c:2):3,d:4);',
    '((a:1,b:1)X:1,c:2,d:4);',
    '((a:1,b:1)X:1,(c:2,d:4):3);',
    '(((a:1)X:1,c:2):3,d:4);',
])
def test_internal_outgroup_pruning_preserves_remaining_distances(newick):
    tree, _ = validate_and_set_outgroup(_newick(newick), 'X')
    distance = tree.distance('c', 'd')
    leaves, counts = prune_outgroup(tree, 'X')
    assert {n.name for n in tree.get_terminals()} == {'c', 'd'}
    assert tree.distance('c', 'd') == pytest.approx(distance)
    assert counts[tree.root] == 2
    assert {n.name for n in leaves[tree.root]} == {'c', 'd'}


def test_single_leaf_midpoint_has_clear_error():
    with pytest.raises(InvalidTreeError, match='at least two leaves'):
        validate_and_set_outgroup(_newick('a;'), None, 'midpoint')


def test_duplicate_names_do_not_collide_with_existing_suffixes():
    tree, _ = validate_and_set_outgroup(_newick('(a:1,a:1,a_1:1);'), None)
    names = [n.name for n in tree.find_clades()]
    assert len(names) == len(set(names))
    assert {n.name for n in tree.get_terminals()} == {'a', 'a_1', 'a_2'}


def test_zero_length_branches_are_preserved_without_warning(caplog):
    tree = _newick('((a:0,b:0):1,c:0);')
    original = [node.branch_length for node in tree.find_clades()]
    with caplog.at_level('WARNING'):
        ensure_branch_lengths(tree)
    assert [node.branch_length for node in tree.find_clades()] == original
    assert not caplog.records


def test_missing_lengths_are_reported_separately_from_zero_lengths(caplog):
    tree = _newick('((a:0,b):1,c:0);')
    with caplog.at_level('WARNING'):
        ensure_branch_lengths(tree)
    assert '1 of 4 branches have missing lengths' in caplog.text
    assert 'Provide' not in caplog.text
    assert tree.find_any(name='a').branch_length == 0
    assert tree.find_any(name='b').branch_length is None
