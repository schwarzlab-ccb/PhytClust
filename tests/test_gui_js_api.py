"""Run the JavaScript GUI JavaScript checks."""

from pathlib import Path
import shutil
import subprocess

import pytest


@pytest.mark.skipif(
    shutil.which("node") is None, reason="Node is needed for GUI JavaScript checks"
)
@pytest.mark.parametrize(
    "filename",
    [
        "api.test.cjs",
        "utils.test.cjs",
        "state.test.cjs",
        "session.test.cjs",
        "colors.test.cjs",
    ],
)
def test_gui_api_requests_and_results(filename):
    test_file = Path(__file__).parent / "js" / filename
    result = subprocess.run(
        ["node", str(test_file)], capture_output=True, text=True, timeout=15
    )
    assert result.returncode == 0, result.stdout + result.stderr
