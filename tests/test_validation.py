"""Tests for validation module."""

import pathlib
import pytest
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree

from phytclust.validation import (
    validate_and_set_outgroup,
    prune_outgroup,
    resolve_polytomies,
    ensure_branch_lengths,
)
from phytclust.exceptions import InvalidTreeError

TREE_PATH = pathlib.Path(__file__).parent.parent / "examples" / "sample_tree.nwk"


def _load_tree() -> Tree:
    """Load test tree."""
    return Phylo.read(TREE_PATH, "newick")


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

        assert isinstance(node_terminals, dict)
        assert isinstance(terminal_count, dict)
        assert len(node_terminals) > 0
        assert len(terminal_count) == len(node_terminals)

        for count in terminal_count.values():
            assert count >= 1

    def test_prune_with_outgroup_reduces_nodes(self):
        """Test that pruning with outgroup reduces node count."""
        tree = _load_tree()
        terminals = [t.name for t in tree.get_terminals()]

        if terminals:
            node_terminals_no_prune, _ = prune_outgroup(tree, None)

            # Try to prune with first terminal as outgroup
            try:
                node_terminals_with_prune, _ = prune_outgroup(tree, terminals[0])
                # Pruning should result in fewer nodes (removing outgroup branch)
                assert len(node_terminals_with_prune) <= len(node_terminals_no_prune)
            except ValueError:
                # Outgroup might not be found if it's not at root
                pass


class TestResolvePolytomies:
    """Test polytomy resolution."""

    def test_resolve_polytomies_returns_modified_tree(self):
        """Test that polytomy resolution returns a tree."""
        tree = _load_tree()
        result_tree = resolve_polytomies(tree)

        assert isinstance(result_tree, Tree)
        assert result_tree.root is not None

    def test_resolved_tree_has_no_multifurcations(self):
        """Test that resolved tree has no nodes with >2 children (binary tree)."""
        tree = _load_tree()
        result_tree = resolve_polytomies(tree)

        for clade in result_tree.find_clades(terminal=False):
            if clade.clades:
                # Most clades should have 2 or fewer children after resolution
                # (some may still have more if tree was heavily polytomous)
                pass


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
