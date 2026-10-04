"""CLI configuration errors, precedence, and output failures."""

from copy import deepcopy
from pathlib import Path
import logging

import pytest
from PIL import Image

from phytclust.cli import app
from phytclust.config import build_runtime_config

TREE = Path(__file__).with_name("test_tree.nwk")


@pytest.fixture(autouse=True)
def restore_library_logging():
    logger = logging.getLogger("phytclust")
    original_handlers = list(logger.handlers)
    original_level, original_propagation = logger.level, logger.propagate
    yield
    for handler in logger.handlers:
        if handler not in original_handlers:
            handler.close()
    logger.handlers[:] = original_handlers
    logger.setLevel(original_level)
    logger.propagate = original_propagation


def test_runtime_sections_merge_without_losing_settings():
    config = {
        "plot": {"marker_size": 12},
        "scores_plot": {"fig_width": 7},
        "save": {"tsv_name": "custom.tsv"},
        "runtime": {"plot": {"scores": {"fig_height": 4}}, "save": {"outlier": False}},
    }
    original = deepcopy(config)
    runtime = build_runtime_config(app._runtime_overrides_from_config(config))
    assert runtime.plot.cluster.marker_size == 12
    assert runtime.plot.scores.fig_width == 7
    assert runtime.plot.scores.fig_height == 4
    assert runtime.save.tsv_name == "custom.tsv"
    assert runtime.save.outlier is False
    assert config == original


@pytest.mark.parametrize(
    "section", ["plot", "scores_plot", "save", "runtime", "algorithm", "peak"]
)
def test_malformed_sections_report_cli_error(tmp_path, capsys, section):
    config_path = tmp_path / "bad.yaml"
    config_path.write_text(f"{section}: wrong\n")
    assert (
        app.main([str(TREE), "--config", str(config_path), "--no-tsv", "--no-color"])
        == 2
    )
    assert f"{section} must be a mapping" in capsys.readouterr().err


def test_bad_nested_runtime_reports_cli_error(tmp_path, capsys):
    config_path = tmp_path / "bad.yaml"
    config_path.write_text("runtime:\n  plot:\n    scores: 3\n")
    assert (
        app.main([str(TREE), "--config", str(config_path), "--no-tsv", "--no-color"])
        == 2
    )
    assert "plot.scores must be a mapping" in capsys.readouterr().err


def test_cli_peak_override_takes_priority():
    config = {"peak": {"prominence_weight": 0.2, "exclude_k2": False}}
    settings = app._peak_config_from_config(
        config, {"prominence_weight": 0.8, "exclude_k2": True}
    )
    assert settings.prominence_weight == 0.8
    assert settings.exclude_k2 is True


@pytest.mark.parametrize(
    "flags", [["-k", "2", "--resolution"], ["-k", "2", "--top-n", "2"]]
)
def test_conflicting_modes_fail_before_running(flags):
    with pytest.raises(SystemExit) as error:
        app.main([str(TREE), *flags, "--no-tsv"])
    assert error.value.code == 2


def test_output_directory_errors_are_reported(tmp_path, capsys):
    occupied = tmp_path / "file"
    occupied.write_text("occupied")
    assert (
        app.main([str(TREE), "-k", "2", "--out-dir", str(occupied), "--no-color"]) == 2
    )
    assert "Cannot create output directory" in capsys.readouterr().err


def test_plot_failure_returns_failure_and_still_saves_tsv(
    monkeypatch, tmp_path, capsys
):
    def fail_plot(*args, **kwargs):
        raise OSError("example plot failure")

    monkeypatch.setattr(app.PhytClust, "plot", fail_plot)
    assert (
        app.main(
            [
                str(TREE),
                "-k",
                "2",
                "--save-fig",
                "--out-dir",
                str(tmp_path),
                "--no-color",
            ]
        )
        == 1
    )
    assert (tmp_path / "phytclust_results.tsv").is_file()
    assert "Could not create the requested plots" in capsys.readouterr().err


def test_tree_export_honours_dpi_and_quiet(tmp_path, capsys):
    assert (
        app.main(
            [
                str(TREE),
                "-k",
                "2",
                "--save-fig",
                "--dpi",
                "72",
                "--out-dir",
                str(tmp_path),
                "--no-color",
                "-q",
            ]
        )
        == 0
    )
    with Image.open(tmp_path / "tree_k2.png") as image:
        assert image.info["dpi"][0] == pytest.approx(72, abs=0.02)
    assert capsys.readouterr().out == ""
    assert not (tmp_path / "scores.png").exists()


def test_input_directories_are_rejected(tmp_path):
    with pytest.raises(SystemExit) as error:
        app.build_parser().parse_args([str(tmp_path)])
    assert error.value.code == 2


def test_deprecated_save_key_still_works():
    with pytest.warns(DeprecationWarning, match="save.csv_name"):
        config = app._runtime_overrides_from_config({"save": {"csv_name": "old.tsv"}})
    assert config["save"]["tsv_name"] == "old.tsv"


def test_peak_method_name_cannot_be_overridden(monkeypatch):
    messages = []
    monkeypatch.setattr(app.LOG, "warning", lambda *args: messages.append(args))
    settings = app._peak_config_from_config({"peak": {"validate": "broken"}})
    assert callable(settings.validate)
    assert messages == [("config: unknown key 'peak.%s' ignored.", "validate")]


def test_python_warnings_use_cli_logging(tmp_path, capsys):
    config_path = tmp_path / "unknown.yaml"
    config_path.write_text("scores_plot:\n  misspelled_setting: 12\n")
    assert (
        app.main([str(TREE), "--config", str(config_path), "--no-tsv", "--no-color"])
        == 0
    )
    error_output = capsys.readouterr().err
    assert (
        "WARNING: Unknown runtime setting 'plot.scores.misspelled_setting'"
        in error_output
    )
    assert ".py:" not in error_output


@pytest.mark.parametrize("quiet_flag", ["-q", "-qq"])
def test_errors_remain_visible_in_quiet_mode(tmp_path, capsys, quiet_flag):
    config_path = tmp_path / "bad.yaml"
    config_path.write_text("runtime: wrong\n")
    assert (
        app.main(
            [
                str(TREE),
                "--config",
                str(config_path),
                "--no-tsv",
                "--no-color",
                quiet_flag,
            ]
        )
        == 2
    )
    assert "ERROR: runtime must be a mapping" in capsys.readouterr().err


def test_user_warnings_are_suppressed_in_quiet_mode(tmp_path, capsys):
    config_path = tmp_path / "unknown.yaml"
    config_path.write_text("scores_plot:\n  misspelled_setting: 12\n")
    assert (
        app.main(
            [str(TREE), "--config", str(config_path), "--no-tsv", "--no-color", "-q"]
        )
        == 0
    )
    captured = capsys.readouterr()
    assert captured.out == captured.err == ""


def test_log_markup_is_kept_literal(tmp_path, capsys, monkeypatch):
    config_path = tmp_path / "bad[red].yaml"
    config_path.write_bytes(b"\xff")
    monkeypatch.setattr(app.sys.stderr, "isatty", lambda: True)
    assert app.main([str(TREE), "--config", str(config_path), "--no-tsv"]) == 2
    assert "bad[red].yaml" in "".join(
        app.Text.from_ansi(capsys.readouterr().err).plain.split()
    )
