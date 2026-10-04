"""Peak settings and relative prominence on logarithmic inputs."""

from types import SimpleNamespace

import numpy as np
import pytest

from phytclust.algo import scoring
from phytclust.config import PeakConfig
from phytclust.exceptions import ConfigurationError


def test_relative_prominence_uses_original_score_baseline(monkeypatch):
    monkeypatch.setattr(scoring, "_first_zero_length_pair_split", lambda *args: None)
    clustering = SimpleNamespace(scores=np.array([0.0, 0.0, 1.0, 10.0, 1.0, 1.0]))
    peaks = scoring.find_score_peaks(
        clustering,
        plot=False,
        peak_config=PeakConfig(
            use_log_peak_input=True,
            use_relative_prominence=True,
            min_relative_prominence=5,
        ),
    )
    assert peaks == [4]
    assert clustering.peak_ranking_details[0]["prominence"] == pytest.approx(10)


@pytest.mark.parametrize(
    "settings, field_name",
    [
        ({"ranking_mode": "other"}, "ranking_mode"),
        ({"resolution_fallback_mode": "other"}, "resolution_fallback_mode"),
        ({"prominence_weight": np.nan}, "prominence_weight"),
        ({"prominence_weight": 1.1}, "prominence_weight"),
        ({"prominence_weight": "0.7"}, "prominence_weight"),
        ({"boundary_window_size": 2.5}, "boundary_window_size"),
        ({"boundary_window_size": True}, "boundary_window_size"),
        ({"boundary_window_size": 0}, "boundary_window_size"),
        ({"boundary_ratio_threshold": np.inf}, "boundary_ratio_threshold"),
        ({"min_k": 0}, "min_k"),
        ({"min_k": 2.5}, "min_k"),
        ({"min_prominence": -1}, "min_prominence"),
        ({"use_log_peak_input": True, "log_peak_offset": 0}, "log_peak_offset"),
        ({"use_log_peak_input": True, "log_peak_offset": np.nan}, "log_peak_offset"),
        (
            {"use_relative_prominence": True, "min_relative_prominence": -1},
            "min_relative_prominence",
        ),
        ({"min_prominence": 1.0, "prominence_k_power": np.inf}, "prominence_k_power"),
        ({"use_relative_prominence": "false"}, "use_relative_prominence"),
    ],
)
def test_invalid_active_settings(settings, field_name):
    with pytest.raises(ConfigurationError, match=field_name):
        PeakConfig(**settings)


@pytest.mark.parametrize(
    "settings",
    [
        {"ranking_mode": "raw", "prominence_weight": np.nan},
        {
            "exclude_k2": True,
            "boundary_window_size": 0,
            "boundary_ratio_threshold": np.nan,
        },
        {"min_k": 3, "boundary_window_size": 0},
        {"use_log_peak_input": False, "log_peak_offset": np.nan},
        {"use_log_peak_input": False, "log_peak_offset": "unused"},
        {"use_relative_prominence": False, "min_relative_prominence": np.nan},
        {"min_prominence": None, "prominence_k_power": np.nan},
        {
            "use_relative_prominence": True,
            "min_prominence": 0.0,
            "prominence_k_power": np.nan,
        },
    ],
)
def test_unused_settings_are_ignored(settings, monkeypatch):
    monkeypatch.setattr(scoring, "_first_zero_length_pair_split", lambda *args: None)
    config = PeakConfig(**settings)
    clustering = SimpleNamespace(scores=np.array([0.0, 0.0, 1.0, 10.0, 1.0, 1.0]))
    assert scoring.find_score_peaks(clustering, peak_config=config, plot=False) == [4]


def test_changes_after_creation_are_validated(monkeypatch):
    config = PeakConfig()
    config.min_k = 0
    with pytest.raises(ConfigurationError, match="min_k"):
        scoring.find_score_peaks(
            SimpleNamespace(scores=[]), peak_config=config, plot=False
        )


def test_numpy_numbers_and_zero_thresholds():
    PeakConfig(
        min_k=np.int64(2),
        boundary_window_size=np.int64(3),
        min_prominence=0.0,
        prominence_weight=np.float64(0.5),
        boundary_ratio_threshold=0.0,
    )


@pytest.mark.parametrize("log_input", [False, True])
def test_relative_ratio_uses_higher_of_the_two_bases(monkeypatch, log_input):
    monkeypatch.setattr(scoring, "_first_zero_length_pair_split", lambda *args: None)
    clustering = SimpleNamespace(scores=np.array([0.0, 0.0, 2.0, 10.0, 4.0, 4.0]))
    scoring.find_score_peaks(
        clustering,
        plot=False,
        peak_config=PeakConfig(
            use_relative_prominence=True, use_log_peak_input=log_input
        ),
    )
    assert clustering.peak_ranking_details[0]["prominence"] == pytest.approx(2.5)


def test_unused_relative_threshold_does_not_filter_boundary_peak(monkeypatch):
    monkeypatch.setattr(scoring, "_first_zero_length_pair_split", lambda *args: None)
    clustering = SimpleNamespace(scores=np.array([0.0, 10.0, 1.0, 5.0, 1.0, 1.0]))
    peaks = scoring.find_score_peaks(
        clustering,
        plot=False,
        global_peaks=2,
        peak_config=PeakConfig(min_relative_prominence=np.nan),
    )
    assert 2 in peaks
