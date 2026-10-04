"""Check soft-polytomy costs against exhaustive leaf partitions."""

from itertools import product
import numpy as np
import pytest
from phytclust import PhytClust
from phytclust.config import OutlierConfig
from phytclust.algo.dp.polytomy import compute_polytomy_dp
from phytclust.exceptions import ConfigurationError


def partitions(items):
    if not items:
        yield []
        return
    first, *rest = items
    for partition in partitions(rest):
        yield [[first], *partition]
        for index in range(len(partition)):
            yield partition[:index] + [[first, *partition[index]]] + partition[
                index + 1 :
            ]


@pytest.mark.parametrize(
    "minimum_size,support_enabled,prefer_fewer",
    list(product([1, 2], [False, True], [False, True])),
)
def test_soft_star_matches_exhaustive_partitions(
    minimum_size, support_enabled, prefer_fewer
):
    pc = PhytClust(
        "(a:1,b:2,c:4,d:8,e:16);",
        min_cluster_size=minimum_size,
        use_branch_support=support_enabled,
        max_k=5,
        outlier=OutlierConfig(size_threshold=3, prefer_fewer=prefer_fewer),
    )
    pc.tree.root.confidence = 50
    pc._ensure_dp(required_cap=5)
    leaves = pc.tree.get_terminals()
    root_id = pc.node_to_id[pc.tree.root]
    root_costs = pc.raw_dp_table[root_id]
    for count in range(1, 6):
        candidates = []
        for partition in partitions(list(range(5))):
            if len(partition) != count or any(
                len(group) < minimum_size for group in partition
            ):
                continue
            cost = sum(
                sum(leaves[i].branch_length for i in group)
                for group in partition
                if len(group) > 1
            )
            if count == 1 and support_enabled:
                cost /= 0.5
            outliers = sum(len(group) < 3 for group in partition)
            candidates.append((cost, outliers))
        if not candidates:
            assert np.isinf(root_costs[count - 1])
            continue
        expected = min(
            candidates, key=lambda item: (item[1], item[0]) if prefer_fewer else item
        )
        assert root_costs[count - 1] == expected[0]
        assert pc._outlier_counts[root_id][count - 1] == expected[1]
        clusters = pc.get_clusters(count)
        assert set(clusters) == set(leaves)
        assert len(set(clusters.values())) == count


def test_dispatch_rejects_unknown_mode():
    pc = PhytClust("(a:1,b:2,c:3);")
    pc.polytomy_mode = "unknown"
    with pytest.raises(ConfigurationError, match="polytomy_mode"):
        compute_polytomy_dp(pc.tree.root, pc, 3, 1, None, np.float64)
