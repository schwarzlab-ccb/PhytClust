from types import MappingProxyType

import numpy as np
import pytest

from phytclust import PhytClust
from phytclust.algo.dp.costs import small_cluster_penalty
from phytclust.config import (
    OutlierConfig,
)
from phytclust.config.runtime import build_runtime_config
from phytclust.exceptions import ConfigurationError


@pytest.mark.parametrize(
    "settings, field_name",
    [
        ({"size_threshold": 0}, "size_threshold"),
        ({"size_threshold": 2.5}, "size_threshold"),
        ({"size_threshold": True}, "size_threshold"),
        ({"prefer_fewer": "false"}, "prefer_fewer"),
        ({"penalty_enabled": "false"}, "penalty_enabled"),
        (
            {"penalty_enabled": True, "size_threshold": 3, "ratio_weight": np.nan},
            "ratio_weight",
        ),
        (
            {"penalty_enabled": True, "size_threshold": 3, "ratio_weight": np.inf},
            "ratio_weight",
        ),
        (
            {"penalty_enabled": True, "size_threshold": 3, "ratio_weight": -1},
            "ratio_weight",
        ),
        (
            {"penalty_enabled": True, "size_threshold": 3, "ratio_weight": True},
            "ratio_weight",
        ),
        (
            {"penalty_enabled": True, "size_threshold": 3, "ratio_mode": "other"},
            "ratio_mode",
        ),
    ],
)
def test_invalid_active_settings(settings, field_name):
    with pytest.raises(ConfigurationError, match=field_name):
        OutlierConfig(**settings)


@pytest.mark.parametrize(
    "settings",
    [
        {"ratio_weight": np.nan, "ratio_mode": "unused"},
        {"size_threshold": 3, "ratio_weight": -1, "ratio_mode": "unused"},
        {"penalty_enabled": True, "ratio_weight": np.inf, "ratio_mode": "unused"},
        {
            "penalty_enabled": True,
            "size_threshold": 3,
            "ratio_weight": 0,
            "ratio_mode": "unused",
        },
    ],
)
def test_unused_penalty_settings_do_not_block_clustering(settings):
    clustering = PhytClust("(a:1,b:2);", outlier=OutlierConfig(**settings))
    assert len(clustering.get_clusters(1)) == 2
    assert small_cluster_penalty(1, clustering) == 0


def test_settings_mutated_after_dp_are_checked():
    clustering = PhytClust("(a:1,b:2);", outlier=OutlierConfig(size_threshold=3))
    clustering.get_clusters(1)
    clustering.outlier.ratio_weight = np.nan
    clustering.outlier.penalty_enabled = True
    with pytest.raises(ConfigurationError, match="ratio_weight"):
        clustering.get_clusters(1)


@pytest.mark.parametrize(
    "mode, expected", [("exp", np.expm1(2) * 2), ("inverse", 2.0), ("power", 4.0)]
)
def test_penalty_shapes_are_preserved(mode, expected):
    config = OutlierConfig(
        size_threshold=np.int64(3),
        penalty_enabled=True,
        ratio_weight=np.float64(2),
        ratio_mode=mode,
    )
    clustering = PhytClust("(a:1,b:2);", outlier=config)
    assert small_cluster_penalty(1, clustering) == pytest.approx(expected)
    assert small_cluster_penalty(3, clustering) == 0


@pytest.mark.parametrize("mode", ["soft", "hard"])
@pytest.mark.parametrize("newick", ["((a:1,b:2):1,(c:1,d:2):3);", "(a:1,b:2,c:3,d:4);"])
def test_prefer_fewer_defaults_to_singletons(newick, mode):
    default = OutlierConfig(prefer_fewer=True)
    assert default.size_threshold is None
    assert default.counting_threshold == 2
    assert OutlierConfig().counting_threshold is None
    implicit = PhytClust(newick, outlier=default, polytomy_mode=mode)
    explicit = PhytClust(
        newick,
        outlier=OutlierConfig(size_threshold=2, prefer_fewer=True),
        polytomy_mode=mode,
    )
    counts = (1, 4) if mode == "hard" and newick == "(a:1,b:2,c:3,d:4);" else (2, 3)
    for count in counts:
        first = {
            leaf.name: cluster for leaf, cluster in implicit.get_clusters(count).items()
        }
        second = {
            leaf.name: cluster for leaf, cluster in explicit.get_clusters(count).items()
        }
        assert first == second


def test_mutating_prefer_fewer_uses_default_singleton_threshold():
    config = OutlierConfig()
    config.prefer_fewer = True
    config.validate()
    assert config.counting_threshold == 2
    config.size_threshold = 3
    assert config.counting_threshold == 3


def test_nested_overrides_keep_other_defaults_and_input():
    overrides = {
        "plot": {"scores": {"fig_width": 7.5, "log_base": None}},
        "save": {"outlier": False},
    }
    config = build_runtime_config(overrides)
    assert config.plot.scores.fig_width == 7.5
    assert config.plot.scores.log_base is None
    assert config.plot.cluster.width_scale == 2.0
    assert config.save.outlier is False
    assert overrides["plot"]["scores"] == {"fig_width": 7.5, "log_base": None}


def test_configs_have_independent_nested_defaults():
    first, second = build_runtime_config(), build_runtime_config({})
    first.plot.scores.fig_width = 20
    first.save.tsv_name = "other.tsv"
    assert second.plot.scores.fig_width == 9
    assert second.save.tsv_name == "phytclust_results.tsv"


def test_unknown_setting_warns_and_valid_setting_applies():
    with pytest.warns(UserWarning, match=r"plot.scores.fig_wdith.*ignored"):
        config = build_runtime_config(
            {"plot": {"scores": {"fig_wdith": 20, "fig_height": 6}}}
        )
    assert config.plot.scores.fig_width == 9
    assert config.plot.scores.fig_height == 6


def test_dataclass_attributes_are_not_settings():
    with pytest.warns(UserWarning, match="__dict__"):
        config = build_runtime_config({"__dict__": {"save": "broken"}})
    assert config.save.tsv_name == "phytclust_results.tsv"


@pytest.mark.parametrize(
    "overrides, name",
    [
        ({"plot": 3}, "plot"),
        ({"plot": {"scores": None}}, "plot.scores"),
        ({"plot": {"cluster": "wrong"}}, "plot.cluster"),
        ({"save": []}, "save"),
        ({"save": {"tsv_name": {}}}, "save.tsv_name"),
    ],
)
def test_malformed_sections_raise_clear_errors(overrides, name):
    with pytest.raises(ConfigurationError, match=name):
        build_runtime_config(overrides)


@pytest.mark.parametrize("overrides", [[], "", False, 0])
def test_invalid_top_level_overrides(overrides):
    with pytest.raises(ConfigurationError, match="Runtime overrides must be a mapping"):
        build_runtime_config(overrides)


def test_read_only_mappings_are_supported():
    config = build_runtime_config(
        MappingProxyType(
            {"plot": MappingProxyType({"scores": MappingProxyType({"fig_height": 6})})}
        )
    )
    assert config.plot.scores.fig_height == 6
