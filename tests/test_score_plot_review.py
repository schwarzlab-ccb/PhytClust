"""Score ranges, axis choices, and resolution labels."""

from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import pytest

from phytclust.config.runtime import ScorePlotConfig
from phytclust.exceptions import ConfigurationError, InvalidClusteringError
from phytclust.viz.scores import plot_scores


@pytest.fixture
def clustering():
    yield SimpleNamespace(
        scores=np.array([0.0, 2.0, 5.0, 3.0, 4.0, 1.0]),
        plot_config=SimpleNamespace(scores=ScorePlotConfig()),
        bin_ranges_current=[(1, 2), (3, 4), (5, 6)],
    )
    plt.close("all")


def test_subset_positions_numpy_peaks_and_title(clustering):
    figure = plot_scores(
        clustering,
        scores_subset=np.array([8.0, 4.0, 6.0]),
        k_start=3,
        peaks=np.array([3, 5, 3]),
        title="My scores",
        log_scale_y=False,
    )
    axis = figure.axes[0]
    np.testing.assert_array_equal(axis.lines[0].get_xdata(), [3, 4, 5])
    assert [text.get_text() for text in axis.texts] == ["3", "5"]
    assert axis.get_title(loc="left") == "My scores"
    np.testing.assert_array_equal(clustering.scores, [0, 2, 5, 3, 4, 1])


def test_log_base_and_zero_warning(clustering):
    with pytest.warns(UserWarning, match="omitted"):
        axis = plot_scores(clustering, log_base=2, peaks=[1, 3]).axes[0]
    assert axis.xaxis._scale.base == axis.yaxis._scale.base == 2
    assert np.isnan(axis.lines[0].get_ydata()[0])
    assert [text.get_text() for text in axis.texts] == ["3"]


def test_small_scores_and_all_zero(clustering):
    axis = plot_scores(
        clustering, scores_subset=np.array([1e-20, 2e-20]), log_scale_y=False
    ).axes[0]
    assert axis.get_ylim()[1] < 1e-19
    axis = plot_scores(clustering, scores_subset=np.array([1e-20, 2e-20])).axes[0]
    assert axis.get_ylim()[0] < 1e-20
    axis = plot_scores(clustering, scores_subset=np.zeros(3)).axes[0]
    assert axis.get_yscale() == "linear"


def test_negative_scores_can_be_shown(clustering):
    clustering.plot_config.scores.clamp_negative_to_zero = False
    axis = plot_scores(clustering, scores_subset=np.array([-2.0, 0.0, 3.0])).axes[0]
    assert axis.get_yscale() == "symlog"
    np.testing.assert_array_equal(axis.lines[0].get_ydata(), [-2, 0, 3])


def test_resolution_bands_only_cover_visible_counts(clustering):
    axis = plot_scores(
        clustering, k_start=4, k_end=5, resolution_on=True, log_scale_y=False
    ).axes[0]
    assert [text.get_text() for text in axis.texts] == ["CL2", "CL3"]
    assert [text.get_text() for text in axis.get_xticklabels()] == ["[4–4]", "[5–5]"]
    assert len(axis.patches) == 2


def test_explicit_empty_subset_is_supported(clustering):
    assert (
        plot_scores(clustering, scores_subset=np.array([])).axes[0].texts[0].get_text()
        == "No scores to plot"
    )
    clustering.scores = None
    with pytest.raises(InvalidClusteringError, match="Scores are empty"):
        plot_scores(clustering)


@pytest.mark.parametrize(
    "options",
    [
        {"k_start": 0},
        {"k_start": True},
        {"k_start": 4, "k_end": 2},
        {"log_base": 1},
        {"log_base": np.inf},
        {"fig_width": 0},
        {"x_axis_mode": "other"},
        {"resolution_on": True, "num_bins": 0},
    ],
)
def test_invalid_options_do_not_create_figures(clustering, options):
    before = plt.get_fignums()
    with pytest.raises(ConfigurationError):
        plot_scores(clustering, **options)
    assert plt.get_fignums() == before


@pytest.mark.parametrize(
    "values", [np.array([np.nan]), np.array([np.inf]), np.ones((2, 2))]
)
def test_invalid_scores(clustering, values):
    with pytest.raises(InvalidClusteringError, match="finite one-dimensional"):
        plot_scores(clustering, scores_subset=values)


def test_custom_sizes_and_empty_title(clustering):
    clustering.plot_config.scores.title_fontsize = 16
    figure = plot_scores(clustering, k_start=2, fig_width=7.5, title="")
    assert figure.get_size_inches()[0] == 7.5
    assert figure.axes[0].get_title(loc="left") == ""


def test_peak_label_boxes_do_not_overlap(clustering):
    figure = plot_scores(
        clustering,
        scores_subset=np.ones(12),
        peaks=list(range(1, 13)),
        x_axis_mode="linear",
    )
    figure.canvas.draw()
    boxes = [
        text.get_window_extent(figure.canvas.get_renderer())
        for text in figure.axes[0].texts
    ]
    assert not any(
        first.overlaps(second)
        for index, first in enumerate(boxes)
        for second in boxes[index + 1 :]
    )


def test_large_peak_labels_near_right_edge_stay_visible(clustering):
    figure = plot_scores(clustering, scores_subset=np.ones(95), peaks=[83, 89, 91])
    figure.canvas.draw()
    axis = figure.axes[0]
    boxes = [
        text.get_window_extent(figure.canvas.get_renderer()) for text in axis.texts
    ]
    assert [text.get_text() for text in axis.texts] == ["83", "89", "91"]
    assert all(
        axis.bbox.contains(box.x0, box.y0) and axis.bbox.contains(box.x1, box.y1)
        for box in boxes
    )
    assert not any(
        first.overlaps(second)
        for index, first in enumerate(boxes)
        for second in boxes[index + 1 :]
    )
