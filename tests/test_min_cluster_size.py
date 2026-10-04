"""Regression checks for minimum cluster sizes in the DP."""

from collections import Counter

import pytest

from phytclust import PhytClust
from phytclust.exceptions import InvalidClusteringError


def _sizes(cmap):
    return sorted(Counter(cmap.values()).values())


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
    assert pc.raw_dp_table[root][0] == pytest.approx(10.0)
    assert pc.raw_dp_table[root][1] == pytest.approx(6.0)


def test_soft_polytomy_merges_small_children_into_valid_clusters():
    pc = PhytClust("(a:1,b:1,c:1,d:1);", min_cluster_size=2)

    assert _sizes(pc.get_clusters(1)) == [4]
    assert _sizes(pc.get_clusters(2)) == [2, 2]
    with pytest.raises(InvalidClusteringError):
        pc.get_clusters(3)
