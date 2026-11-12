import pathlib
from collections import Counter

import numpy as np
import pytest
from Bio.Phylo.BaseTree import Clade, Tree

from phytclust import PhytClust
from phytclust.algo.core import _coerce_to_tree
from phytclust.algo.dp.table import _copy_tree
from phytclust.config import PeakConfig
from phytclust.exceptions import (
    ConfigurationError,
    InvalidClusteringError,
    InvalidKError,
    InvalidTreeError,
    ValidationError,
)

EXAMPLES = pathlib.Path(__file__).resolve().parent.parent / "examples"
POLYTOMY = EXAMPLES / "sample_polytomy.newick"
SAMPLE = EXAMPLES / "sample_tree.nwk"


def _sizes(cluster_map):
    return sorted(Counter(cluster_map.values()).values())


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
    pc = PhytClust("(a:1,b:1,c:1,d:1);")
    with pytest.raises(ConfigurationError, match="num_bins"):
        pc.run(by_resolution=True, num_bins=0, plot_scores=False)


@pytest.mark.parametrize("cached", [False, True])
@pytest.mark.parametrize(
    "method,value",
    [
        (method, value)
        for method in ["get_clusters", "run", "alpha_info"]
        for value in [2.0, 2.7, True, np.bool_(True), "2", None, 0, -1, np.nan]
        if not (method == "run" and value is None)
    ],
)
def test_invalid_cluster_counts_do_not_depend_on_cached_partitions(
    cached, method, value
):
    clustering = PhytClust("((a:1,b:2):3,(c:4,d:5):6);")
    if cached:
        clustering.get_clusters(2)
        clustering.get_clusters(1)
    kwargs = {"k": value}
    if method == "run":
        kwargs["plot_scores"] = False
    with pytest.raises(InvalidKError, match="positive integer"):
        getattr(clustering, method)(**kwargs)


@pytest.mark.parametrize("count", [2, np.int64(2)])
def test_integer_cluster_counts_work_before_and_after_caching(count):
    clustering = PhytClust("((a:1,b:2):3,(c:4,d:5):6);")
    first = clustering.get_clusters(count)
    second = clustering.get_clusters(count)
    assert len(set(first.values())) == 2
    assert first == second
    assert clustering.run(k=count, plot_scores=False)["k_values"] == [2]
    assert clustering.alpha_info(count)["k"] == 2


def test_small_children_can_form_one_valid_parent_cluster():
    pc = PhytClust("((a:1,b:1):1,c:1);", min_cluster_size=2)

    assert _sizes(pc.get_clusters(1)) == [3]
    with pytest.raises(InvalidClusteringError):
        pc.get_clusters(2)


def test_binary_minimum_size_allows_only_feasible_partitions():
    pc = PhytClust("((a:1,b:2):1,(c:2,d:1):1);", min_cluster_size=2)

    assert _sizes(pc.get_clusters(1)) == [4]
    assert _sizes(pc.get_clusters(2)) == [2, 2]
    with pytest.raises(InvalidClusteringError):
        pc.get_clusters(3)

    root = pc.node_to_id[pc.tree.root]
    assert pc.unpenalised_dp_table[root][0] == pytest.approx(10.0)
    assert pc.unpenalised_dp_table[root][1] == pytest.approx(6.0)


def test_soft_polytomy_merges_small_children_into_valid_clusters():
    pc = PhytClust("(a:1,b:1,c:1,d:1);", min_cluster_size=2)

    assert _sizes(pc.get_clusters(1)) == [4]
    assert _sizes(pc.get_clusters(2)) == [2, 2]
    with pytest.raises(InvalidClusteringError):
        pc.get_clusters(3)


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


def test_internal_clade_outgroup_counts_all_removed_leaves():
    pc = PhytClust("((a:1,b:1)X:1,(c:1,d:1):1);", outgroup="X")

    assert pc.num_terminals == 2
    assert pc.max_k == 2
    assert {leaf.name for leaf in pc.get_clusters(2)} == {"c", "d"}


def test_entire_tree_cannot_be_pruned_as_outgroup():
    with pytest.raises(InvalidTreeError, match="entire tree"):
        PhytClust("((a:1,b:1):1,c:1)ROOT;", outgroup="ROOT")


@pytest.mark.parametrize("text", ["nonexistent_file.nwk", "not a tree", "  "])
def test_mistyped_path_is_not_read_as_one_leaf_tree(text):
    with pytest.raises(ValidationError):
        PhytClust(text)


def test_none_configs_get_defaults():
    pc = PhytClust(SAMPLE, runtime_config=None, peak_config=None, outlier=None)
    assert pc.plot_config is not None
    assert isinstance(pc.peak_config, PeakConfig)


def test_small_tree_error_names_the_leaf_count():
    with pytest.raises(InvalidKError, match="this tree has 3"):
        PhytClust("((A:1,B:1):1,C:1);").run(plot_scores=False)


def test_small_max_k_limit_error_names_the_setting():
    with pytest.raises(InvalidKError, match="max_k_limit"):
        PhytClust(SAMPLE, max_k_limit=0.3).run(plot_scores=False)


def test_outgroup_only_child_leaves_no_phantom_leaf():
    tree = "(((A:1,B:1):1,(O:1)U:1):1,(E:1,F:1):2);"
    pc = PhytClust(tree, outgroup="O")
    names = sorted(leaf.name for leaf in pc.get_clusters(2))
    assert names == ["A", "B", "E", "F"]


def test_growing_dp_keeps_last_peaks(tmp_path):
    pc = PhytClust(POLYTOMY)
    peaks = pc.run(top_n=3, plot_scores=False)["k_values"]
    pc.get_clusters(pc.num_terminals)
    assert pc.peaks_by_rank == peaks
    assert pc.save(str(tmp_path), top_n=3) is not None
