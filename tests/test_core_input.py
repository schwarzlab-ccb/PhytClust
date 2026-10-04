"""Regression tests for tree input handling."""

import pytest
from Bio.Phylo.BaseTree import Tree

from phytclust.algo.core import _coerce_to_tree, PhytClust
from phytclust.exceptions import ValidationError


def test_tree_input_preserves_existing_tree():
    tree = Tree()
    assert _coerce_to_tree(tree) is tree


def test_newick_strings_and_paths(tmp_path):
    path = tmp_path / "tree.nwk"
    path.write_text("(a:1,b:2);")
    for source in [path, str(path), path.read_text()]:
        assert len(_coerce_to_tree(source).get_terminals()) == 2


def test_invalid_newick_and_missing_files_raise_validation_error(tmp_path):
    for source in ["(a,b", tmp_path / "missing.nwk", 123]:
        with pytest.raises(ValidationError):
            _coerce_to_tree(source)


def test_fingerprint_failure_is_not_hidden(monkeypatch):
    pc = PhytClust("(a:1,b:1);")
    def fail(root):
        raise ValueError("invalid metadata")
    monkeypatch.setattr("phytclust.algo.core.tree_fingerprint", fail)
    with pytest.raises(ValueError, match="invalid metadata"):
        pc.get_clusters(1)


def test_alpha_refreshes_after_tree_edit():
    pc = PhytClust("((a:1,b:2):3,c:4);")
    before = pc.alpha_info(2)
    pc.tree.find_any(name="a").branch_length = 10
    after = pc.alpha_info(2)
    assert after is not before


def test_custom_alpha_map_does_not_return_cached_partition():
    pc = PhytClust("((a:1,b:2):3,c:4);")
    cached = pc.alpha_info(2)
    custom = {leaf: 0 for leaf in pc.tree.get_terminals()}
    result = pc.alpha_info(2, custom)
    assert result is not cached
    assert pc.alpha_by_k[2] is cached


def test_resolution_rejects_zero_bins():
    from phytclust.exceptions import ConfigurationError
    pc = PhytClust("(a:1,b:1,c:1,d:1);")
    with pytest.raises(ConfigurationError, match="num_bins"):
        pc.run(by_resolution=True, num_bins=0, plot_scores=False)
