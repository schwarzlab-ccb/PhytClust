"""Nested runtime overrides and their error messages."""

from types import MappingProxyType

import pytest

from phytclust.config.runtime import build_runtime_config
from phytclust.exceptions import ConfigurationError


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
