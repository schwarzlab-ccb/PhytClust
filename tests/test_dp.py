import itertools
import numpy as np
import pytest
from Bio.Phylo.BaseTree import Clade

from conftest import SAMPLE_TREE, leaf_groups, random_newick, tree_from_newick
from phytclust import PhytClust
from phytclust.algo.dp import merge
from phytclust.algo.dp.costs import (
    floating_point_precision,
    is_better_partition,
    small_cluster_penalty,
)
from phytclust.algo.dp.table import _assign_subtrees_to_cluster, backtrack
from phytclust.config import OutlierConfig
from phytclust.exceptions import (
    ConfigurationError,
    InvalidClusteringError,
    InvalidKError,
    InvalidTreeError,
)


def reference(
    left,
    right,
    left_total,
    right_total,
    left_counts,
    right_counts,
    cap,
    settings,
    left_balance=None,
):
    outputs = [
        np.full(cap + 1, np.inf),
        np.full(cap + 1, np.inf),
        np.zeros(cap + 1, dtype=int),
        np.full(cap + 1, -1, dtype=int),
    ]
    for k in range(1, cap + 1):
        low = max(0, k - len(right))
        candidates = []
        for i in range(low, min(k, len(left))):
            j = k - 1 - i
            cost = left[i] + right[j]
            total = left_total[i] + right_total[j]
            if not np.isfinite(cost) or not np.isfinite(total):
                continue
            count = (
                left_counts[i] + right_counts[j]
                if settings.outlier_detection_enabled
                else 0
            )
            candidates.append(
                (
                    cost,
                    total,
                    count,
                    ((i + 1) ** 2 if left_balance is None else left_balance[i])
                    + (j + 1) ** 2,
                    i,
                )
            )
        if not candidates:
            continue

        def filter_cost(index):
            nonlocal candidates
            best = min(c[index] for c in candidates)
            candidates = [
                c
                for c in candidates
                if c[index] <= best + settings.relative_cost_tolerance * abs(best)
            ]

        def filter_count():
            nonlocal candidates
            best = min(c[2] for c in candidates)
            candidates = [c for c in candidates if c[2] == best]

        if settings.small_cluster_penalty_enabled:
            filter_cost(1)
            filter_count()
            filter_cost(0)
        elif settings.outlier_detection_enabled and settings.prioritize_fewer_outliers:
            filter_count()
            filter_cost(0)
        else:
            filter_cost(0)
            if settings.outlier_detection_enabled:
                filter_count()
        chosen = min(candidates, key=lambda c: (c[3], c[4]))
        outputs[0][k], outputs[1][k], outputs[2][k], outputs[3][k] = (
            chosen[0],
            chosen[1] if settings.small_cluster_penalty_enabled else chosen[0],
            chosen[2],
            chosen[4] - low,
        )
    return outputs


@pytest.mark.parametrize("lengths", [(2, 12), (12, 2), (9, 13), (13, 9)])
@pytest.mark.parametrize("mode", ["cost", "outliers", "fewer", "penalty"])
@pytest.mark.parametrize("prefix_balance", [False, True])
@pytest.mark.parametrize("batched", [False, True])
def test_merge_matches_reference(monkeypatch, lengths, mode, batched, prefix_balance):
    monkeypatch.setattr(merge, "_BATCH_MIN_STATES", 8 if batched else 10000)
    rng = np.random.default_rng(41)
    left, right = [rng.integers(0, 5, n).astype(float) for n in lengths]
    if lengths[0] == 2:
        left[1] = np.inf
    if lengths[1] == 2:
        right[1] = np.inf
    left_total, right_total = left + 2, right + 1
    left_counts, right_counts = [rng.integers(0, 3, n) for n in lengths]
    settings = merge.MergeSettings(
        mode != "cost", mode == "penalty", mode == "fewer", 0.0, np.float64
    )
    cap = 10
    left_balance = (
        np.arange(len(left), dtype=np.int64) ** 2 + 7 if prefix_balance else None
    )
    balance_output = np.zeros(cap + 1, dtype=np.int64) if prefix_balance else None
    expected = reference(
        left,
        right,
        left_total,
        right_total,
        left_counts,
        right_counts,
        cap,
        settings,
        left_balance=left_balance,
    )
    outputs = [
        np.full(cap + 1, np.inf),
        np.full(cap + 1, np.inf),
        np.zeros(cap + 1, dtype=int),
        np.full(cap + 1, -1, dtype=int),
    ]
    with np.errstate(invalid="raise"):
        merge.merge_child_tables(
            left,
            right,
            left_total,
            right_total,
            left_counts if mode != "cost" else None,
            right_counts if mode != "cost" else None,
            cap,
            settings,
            *outputs,
            left_cluster_count_sum_of_squares=left_balance,
            cluster_count_sum_of_squares_out=balance_output,
        )
    for actual, wanted in zip(outputs, expected):
        np.testing.assert_array_equal(actual, wanted)


@pytest.mark.parametrize("batched", [False, True])
@pytest.mark.parametrize("lengths", [(9, 13), (13, 9)])
@pytest.mark.parametrize("mode", ["cost", "outliers", "fewer", "penalty"])
@pytest.mark.parametrize("cluster_count", [8, 9])
def test_cost_ties_choose_smallest_child_state_difference(
    monkeypatch, batched, lengths, mode, cluster_count
):
    monkeypatch.setattr(merge, "_BATCH_MIN_STATES", 8 if batched else 10000)
    left_costs, right_costs = [np.zeros(length) for length in lengths]
    left_counts, right_counts = [np.zeros(length, dtype=int) for length in lengths]
    settings = merge.MergeSettings(
        mode != "cost", mode == "penalty", mode == "fewer", 0.0, np.float64
    )
    outputs = [
        np.full(cluster_count + 1, np.inf),
        np.full(cluster_count + 1, np.inf),
        np.zeros(cluster_count + 1, dtype=int),
        np.full(cluster_count + 1, -1, dtype=int),
    ]
    merge.merge_child_tables(
        left_costs, right_costs, left_costs, right_costs,
        left_counts if mode != "cost" else None,
        right_counts if mode != "cost" else None,
        cluster_count, settings, *outputs,
    )
    left_state, right_state = merge.decode_split_offset(
        int(outputs[3][cluster_count]), cluster_count, len(right_costs)
    )
    candidates = [
        (i, cluster_count - 1 - i)
        for i in range(len(left_costs))
        if 0 <= cluster_count - 1 - i < len(right_costs)
    ]
    assert (left_state, right_state) == min(
        candidates, key=lambda pair: (abs(pair[0] - pair[1]), pair[0])
    )


def test_precision_is_relative_spacing():
    assert floating_point_precision(np.float64) == np.finfo(np.float64).eps


@pytest.mark.parametrize('weight', [-1, float('nan'), float('inf')])
def test_support_settings_checked_only_when_enabled(weight):
    pc = PhytClust('(a:1,b:2);', support_weight=weight)
    pc.tree.root.confidence = 200
    pc.get_clusters(1)
    pc.use_branch_support = True
    with pytest.raises(ConfigurationError, match='support_weight'):
        pc.get_clusters(1)


@pytest.mark.parametrize('confidence', [-1, 101, float('nan'), float('inf')])
def test_invalid_confidence_rejected_when_support_enabled(confidence):
    pc = PhytClust('(a:1,b:2);', use_branch_support=True)
    pc.tree.root.confidence = confidence
    with pytest.raises(InvalidTreeError, match='Confidence'):
        pc.get_clusters(1)


@pytest.mark.parametrize('name,value', [('min_cluster_size', 1.5), ('soft_polytomy_max_degree', True), ('min_support', 0)])
def test_invalid_parameters(name, value):
    pc = PhytClust('(a:1,b:2);', use_branch_support=True, **{name: value})
    with pytest.raises(ConfigurationError):
        pc.get_clusters(1)


def test_overflowing_penalty_has_clear_error():
    pc = PhytClust('(a:1,b:2);', outlier=OutlierConfig(size_threshold=1000, penalty_enabled=True))
    with pytest.raises(ConfigurationError, match='exponential outlier penalty'):
        small_cluster_penalty(1, pc)
    pc.outlier.penalty_enabled = False
    assert small_cluster_penalty(1, pc) == 0


def test_comparison_priorities_and_relative_ties():
    assert is_better_partition(10, 0, 1, 1, outlier_detection_enabled=True, prioritize_fewer_outliers=True)
    assert not is_better_partition(10, 0, 1, 1, outlier_detection_enabled=True, prioritize_fewer_outliers=False)
    for scale in [1, 1000]:
        assert is_better_partition(1.000001 * scale, 0, 1 * scale, 1,
                                  outlier_detection_enabled=True, prioritize_fewer_outliers=False,
                                  relative_cost_tolerance=1e-5)


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


@pytest.mark.parametrize(
    "newick",
    [
        SAMPLE_TREE.read_text(),
        "((A:1,B:2,C:1):1,(D:1,(E:1,F:3):1):2,G:4);",
    ],
)
def test_hard_dp_cost_matches_brute_force(newick):
    expected = _brute_force_costs(tree_from_newick(newick))
    pc = PhytClust(tree=tree_from_newick(newick), polytomy_mode="hard")
    pc.get_clusters(pc.num_terminals)
    root_row = pc.unpenalised_dp_table[pc.node_to_id[pc.tree.root]]
    for k in range(1, pc.num_terminals + 1):
        want = expected.get(k, np.inf)
        assert root_row[k - 1] == pytest.approx(want)


def _brute_force_costs(tree):
    """Optimal cost for every k by enumerating monophyletic partitions."""
    def depth_sum(node):
        return sum(
            sum(step.branch_length or 0.0 for step in node.get_path(t))
            for t in node.get_terminals()
        )

    def partitions(node):
        options = [[depth_sum(node)]]
        if node.clades:
            child_sets = [partitions(c) for c in node.clades]
            for combo in itertools.product(*child_sets):
                options.append([x for part in combo for x in part])
        return options

    best = {}
    for part in partitions(tree.root):
        k = len(part)
        best[k] = min(best.get(k, np.inf), sum(part))
    return best


@pytest.mark.parametrize("factor", [2.0**-30, 1e-9, 1e-12, 1e6])
@pytest.mark.parametrize("seed", [1, 2, 3])
def test_invariant_to_global_rescaling(seed, factor):
    """Multiplying every branch length by a positive constant changes nothing."""
    newick = random_newick(seed, 60)
    base = PhytClust(tree=tree_from_newick(newick))
    scaled = PhytClust(tree=_scaled(newick, factor))

    base_result = base.run(top_n=3, plot_scores=False)
    scaled_result = scaled.run(top_n=3, plot_scores=False)
    assert base_result["k_values"] == scaled_result["k_values"]
    np.testing.assert_allclose(base.scores, scaled.scores, rtol=1e-9)

    for k in (2, 5, 10, 25, 50):
        assert leaf_groups(base.get_clusters(k)) == leaf_groups(scaled.get_clusters(k))


@pytest.mark.parametrize("factor", [2.0**-30, 1e-9, 1e-12, 1e6])
def test_rescaling_keeps_tie_breaks_on_integer_tree(factor):
    base = PhytClust(tree=tree_from_newick(SAMPLE_TREE.read_text()))
    scaled = PhytClust(tree=_scaled(SAMPLE_TREE.read_text(), factor))
    for k in range(1, base.num_terminals + 1):
        assert leaf_groups(base.get_clusters(k)) == leaf_groups(scaled.get_clusters(k))


def _scaled(newick, factor):
    tree = tree_from_newick(newick)
    for clade in tree.find_clades():
        if clade.branch_length is not None:
            clade.branch_length *= factor
    return tree
