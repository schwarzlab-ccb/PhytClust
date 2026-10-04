"""Check input filtering and row alignment in the Clonetrac notebook pipeline."""

import importlib.util
from io import StringIO
from pathlib import Path

from Bio import Phylo
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest


spec = importlib.util.spec_from_file_location(
    "clonetrac_analysis", Path(__file__).parents[1] / "clonetrac" / "run_analysis.py"
)
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


def profiles():
    return pd.DataFrame(
        [
            {
                "sample_id": cell,
                "chrom": chromosome,
                "start": 1,
                "end": 10,
                "major": value,
            }
            for cell, values in [("a", [1, 2]), ("b", [3, 4]), ("c", [5, 6])]
            for chromosome, value in zip(["chr2", "chr1"], values)
        ]
    )


def test_discovery_uses_only_original_top_level_pairs(tmp_path):
    (tmp_path / "trees").mkdir()
    (tmp_path / "trees" / "NB02_final_tree.new").touch()
    (tmp_path / "NB01_no_normals.new").touch()
    (tmp_path / "NB01_final_tree.new").touch()
    (tmp_path / "NB01_final_cn_profiles_without_ecdna.tsv").touch()
    samples = analysis.discover_samples(tmp_path)
    assert [sample[0] for sample in samples] == ["NB01"]
    assert samples[0][2].name.endswith("without_ecdna.tsv")


def test_normal_cells_are_removed_from_copies_of_both_inputs(tmp_path):
    tree_path = tmp_path / "input.new"
    original_tree = "((CTR_NB01_A1_A1.final.bam:1,a:2):3,diploid:4);"
    tree_path.write_text(original_tree)
    profile_path = tmp_path / "input.tsv"
    frame = profiles().replace({"a": "CTR_NB01_A1_A1.final.bam"})
    frame.to_csv(profile_path, sep="\t", index=False)
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    tree, filtered, original_count, removed = analysis.prepare_inputs(
        tree_path, profile_path, {"NB01_A1_A1"}, output_dir
    )
    assert original_count == 3
    assert removed == ["CTR_NB01_A1_A1.final.bam"]
    assert [leaf.name for leaf in tree.get_terminals()] == ["a", "diploid"]
    assert not filtered.sample_id.str.contains("NB01").any()
    assert tree_path.read_text() == original_tree
    assert len(pd.read_csv(profile_path, sep="\t")) == len(frame)


def test_profiles_follow_tree_order_and_genomic_order():
    values, bins, edges = analysis.profile_matrix(profiles(), ["c", "a", "b"])
    np.testing.assert_array_equal(values, [[6, 5], [2, 1], [4, 3]])
    assert [segment[0] for segment in bins] == ["chr1", "chr2"]
    np.testing.assert_array_equal(edges, [0, 10, 20])


@pytest.mark.parametrize("problem", ["missing_cell", "missing_segment", "duplicate"])
def test_incomplete_or_duplicate_profiles_fail_clearly(problem):
    frame = profiles()
    if problem == "missing_cell":
        frame = frame[frame.sample_id != "c"]
    elif problem == "missing_segment":
        frame = frame.drop(index=0)
    else:
        frame = pd.concat([frame, frame.iloc[:1]])
    with pytest.raises(ValueError, match="missing|finite|Duplicate"):
        analysis.profile_matrix(frame, ["a", "b", "c"])


def test_tree_cluster_strip_and_profile_rows_have_identical_coordinates():
    tree = Phylo.read(StringIO("((c:1,a:2):1,b:3);"), "newick")
    leaves = tree.get_terminals()
    clusters = dict(zip(leaves, [1, 0, 1]))
    figure, leaf_order = analysis.plot_tree_and_profiles(
        tree, clusters, profiles(), "NB_test", 2
    )
    try:
        assert leaf_order == ["c", "a", "b"]
        tree_axis, cluster_axis, profile_axis, _ = figure.axes
        assert (
            tree_axis.get_ylim()
            == cluster_axis.get_ylim()
            == profile_axis.get_ylim()
            == (2.5, -0.5)
        )
        np.testing.assert_array_equal(
            profile_axis.collections[0].get_array(), [[6, 5], [2, 1], [4, 3]]
        )
        assert len(figure.legends[0].texts) == 2
    finally:
        plt.close(figure)
