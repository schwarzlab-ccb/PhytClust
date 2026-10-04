import json
import re
import shutil
import subprocess
from pathlib import Path

import matplotlib.colors as mcolors
import numpy as np
import pytest

from phytclust.exceptions import ConfigurationError
from phytclust.viz.palette import BASE_HEX, expand_palette, get_cmap


@pytest.mark.parametrize("count", [-1, True, 1.5, float("nan")])
def test_invalid_color_counts(count):
    with pytest.raises(ConfigurationError, match="Colour count"):
        expand_palette(count)


def test_empty_palette_and_zero_count():
    assert expand_palette(0) == []
    with pytest.raises(ConfigurationError, match="base palette"):
        expand_palette(2, base=[])
    with pytest.raises(ConfigurationError, match="greater than zero"):
        get_cmap(0)


def test_base_order_and_numpy_count_are_preserved():
    assert expand_palette(np.int64(16)) == [
        mcolors.to_rgba(color) for color in BASE_HEX
    ]
    assert get_cmap(n=20).N == 20


def test_custom_alpha_is_preserved_during_expansion():
    colors = expand_palette(3, base=[(0.2, 0.4, 0.6, 0.5)])
    assert [color[3] for color in colors] == pytest.approx([0.5, 0.44, 0.38])
    assert colors[1][:3] == pytest.approx([0.312, 0.484, 0.656])


@pytest.mark.parametrize("color", ["not-a-color", (float("nan"), 0, 0, 1)])
def test_invalid_colors_are_reported(color):
    with pytest.raises(ConfigurationError, match="colour"):
        expand_palette(1, base=[color])


def test_large_palette_stops_before_white_and_reports_repeats():
    with pytest.warns(UserWarning, match="repeat beyond 112"):
        colors = expand_palette(160)
    assert all(any(channel < 1 for channel in color[:3]) for color in colors)
    assert colors[112] == colors[96]


@pytest.mark.skipif(
    shutil.which("node") is None, reason="Node is needed to compare GUI colours"
)
def test_gui_and_python_generate_matching_colors():
    path = Path(__file__).resolve().parents[1] / "src/phytclust/gui/static/js/colors.js"
    script = """const fs = require('fs');
const vm = require('vm');
const path = require('path');
const d3 = require(path.join(path.dirname(process.argv[1]), '../vendor/d3.v7.min.js'));
const source = fs.readFileSync(process.argv[1], 'utf8').replace(/^import .*$/gm, '').replace(/export /g, '');
const colors = vm.runInNewContext(source + '; generateClusterColors(64)', { d3 });
process.stdout.write(JSON.stringify(colors));"""
    output = subprocess.run(
        ["node", "-e", script, str(path)],
        capture_output=True,
        text=True,
        check=True,
        timeout=10,
    )
    colors = json.loads(output.stdout)
    parsed = []
    for color in colors:
        if color.startswith("#"):
            parsed.append(mcolors.to_rgba(color))
        else:
            channels = [float(value) for value in re.findall(r"[\d.]+", color)]
            parsed.append(
                (
                    *[value / 255 for value in channels[:3]],
                    channels[3] if len(channels) > 3 else 1,
                )
            )
    # CSS output rounds RGB channels to the nearest 8-bit value.
    expected = np.asarray(expand_palette(64))
    np.testing.assert_allclose(
        np.asarray(parsed)[:, :3], expected[:, :3], atol=0.5 / 255, rtol=0
    )
    np.testing.assert_allclose(np.asarray(parsed)[:, 3], expected[:, 3], atol=1e-14)
