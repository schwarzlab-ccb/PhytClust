#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import pathlib
import sys
import textwrap
import logging
import time
import warnings
from typing import Any
from dataclasses import fields
from copy import deepcopy

from Bio import Phylo
from ..algo.core import PhytClust
from ..config import OutlierConfig, PeakConfig, build_runtime_config
from ..exceptions import ConfigurationError
from pathlib import Path
from importlib import resources

from importlib.metadata import version, PackageNotFoundError

try:
    from rich.console import Console
    from rich.text import Text
except ImportError:
    Console = None
    Text = None

console = Console(stderr=True) if Console else None
output_console = Console() if Console else None

LOG = logging.getLogger("phytclust.cli")


def print_banner():
    try:
        text = (
            resources.files("phytclust")
            .joinpath("assets", "ascii_logo.txt")
            .read_text(encoding="utf-8")
        )
        print(text)
    except (OSError, ModuleNotFoundError):
        LOG.debug("ascii_logo.txt not found; skipping banner.")


def _positive_int(value: str) -> int:
    integer_value = int(value)
    if integer_value < 1:
        raise argparse.ArgumentTypeError("value must be ≥ 1")
    return integer_value


def _min_int(minimum_value: int):
    def _check(value: str) -> int:
        integer_value = int(value)
        if integer_value < minimum_value:
            raise argparse.ArgumentTypeError(f"value must be ≥ {minimum_value}")
        return integer_value

    return _check


def _existing_path_or_stdin(value: str) -> pathlib.Path | str:
    if value == "-":
        return value
    path = pathlib.Path(value)
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"file not found: {value}")
    return path


def _existing_path(value: str) -> pathlib.Path:
    path = pathlib.Path(value)
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"file not found: {value}")
    return path


def _package_version() -> str:
    if version is None:
        return "unknown"
    try:
        return version("phytclust")
    except PackageNotFoundError:
        return "unknown"


def _load_config(path: pathlib.Path | None) -> dict:
    """Read a YAML or JSON configuration file."""
    if not path:
        return {}
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        raise ConfigurationError(f"Could not read config {path}: {exc}") from exc

    import yaml

    try:
        data = yaml.safe_load(text)
    except yaml.YAMLError as exc:
        raise ConfigurationError(f"Could not parse config {path}: {exc}") from exc

    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ConfigurationError(
            f"Config {path} must contain a mapping at the top level, "
            f"got {type(data).__name__}."
        )
    return data


def _warn_deprecated_key(old: str, new: str) -> None:
    """Report a deprecated setting and its replacement."""
    message = f"config key '{old}' is deprecated, use '{new}' instead."
    warnings.warn(message, DeprecationWarning, stacklevel=3)
    LOG.warning("config: %s", message)


def _config_section(config: dict[str, Any], name: str) -> dict[str, Any]:
    """Read a configuration section and check that it is a mapping."""
    section = config.get(name, {})
    if not isinstance(section, dict):
        raise ConfigurationError(f"{name} must be a mapping of settings.")
    return section


def _merge_settings(base: dict[str, Any], overrides: dict[str, Any]) -> dict[str, Any]:
    """Merge nested settings without changing either input."""
    merged = deepcopy(base)
    for name, value in overrides.items():
        if isinstance(value, dict) and isinstance(merged.get(name), dict):
            merged[name] = _merge_settings(merged[name], value)
        else:
            merged[name] = deepcopy(value)
    return merged


def _runtime_overrides_from_config(config: dict[str, Any]) -> dict[str, Any]:
    """Combine plotting and save sections with nested runtime overrides."""
    overrides: dict[str, Any] = {}
    cluster_settings = _config_section(config, "plot")
    score_settings = _config_section(config, "scores_plot")
    save_settings = dict(_config_section(config, "save"))
    if cluster_settings:
        overrides.setdefault("plot", {})["cluster"] = dict(cluster_settings)
    if score_settings:
        overrides.setdefault("plot", {})["scores"] = dict(score_settings)
    if "csv_name" in save_settings:
        _warn_deprecated_key("save.csv_name", "save.tsv_name")
        if save_settings.get("tsv_name") is None:
            save_settings["tsv_name"] = save_settings["csv_name"]
        del save_settings["csv_name"]
    if save_settings:
        overrides["save"] = save_settings
    return _merge_settings(overrides, _config_section(config, "runtime"))


def _peak_config_from_config(
    config: dict[str, Any],
    cli_overrides: dict[str, Any] | None = None,
) -> PeakConfig:
    """Create peak settings; explicit CLI flags override file settings."""
    peak_config = PeakConfig()

    algorithm_settings = _config_section(config, "algorithm")
    peak_block = (
        _config_section(config, "peak")
        if "peak" in config
        else _config_section(algorithm_settings, "peak")
    )
    field_names = {config_field.name for config_field in fields(PeakConfig)}
    if isinstance(peak_block, dict):
        for key, val in peak_block.items():
            if key in field_names:
                setattr(peak_config, key, val)
            else:
                LOG.warning("config: unknown key 'peak.%s' ignored.", key)

    for key, val in (cli_overrides or {}).items():
        if val is not None and key in field_names:
            setattr(peak_config, key, val)

    peak_config.validate()
    return peak_config


def _add_common_run_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--plot",
        dest="plot",
        default=False,
        action=argparse.BooleanOptionalAction,
        help="Open interactive plot windows after running; use --save-fig to write PNGs instead.",
    )
    parser.add_argument(
        "--save-fig",
        action="store_true",
        help="Write coloured tree PNGs and, in peak modes, scores.png to --out-dir.",
    )
    parser.add_argument(
        "--save-tree",
        action="store_true",
        help="Explicitly write coloured tree PNG(s) (in addition to --save-fig).",
    )
    parser.add_argument(
        "--save-all-k",
        action="store_true",
        help="Write tsv rows for *every* k from 1..max_k (can be large!).",
    )
    parser.add_argument(
        "--tsv-name",
        default=None,
        help="Output tsv filename (default: from config or 'phytclust_results.tsv').",
    )
    parser.add_argument(
        "--no-tsv", action="store_true", help="Skip writing the results tsv."
    )
    parser.add_argument(
        "--dpi", type=_positive_int, default=150, help="PNG resolution (DPI)."
    )


def _configure_logging(
    verbosity: int, quiet: int, no_color: bool, log_format: str
) -> None:
    """Set CLI logging level and output format."""
    level = logging.WARNING
    if verbosity > 0:
        level = logging.INFO
    if verbosity > 1:
        level = logging.DEBUG
    if quiet > 0:
        level = logging.ERROR

    root = logging.getLogger("phytclust")
    root.handlers.clear()
    root.setLevel(level)
    root.propagate = False

    handler = None
    if not no_color and sys.stderr.isatty():
        try:
            from rich.logging import RichHandler
            from rich.highlighter import NullHighlighter

            handler = RichHandler(
                level=level,
                console=console,
                markup=False,
                highlighter=NullHighlighter(),
                rich_tracebacks=False,
                show_time=False,
                show_level=True,
                show_path=False,
            )
            formatter = logging.Formatter(log_format)
            handler.setFormatter(formatter)
        except ImportError:
            handler = None

    if handler is None:
        handler = logging.StreamHandler()
        plain_format = (
            "%(levelname)s: %(message)s" if log_format == "%(message)s" else log_format
        )
        formatter = logging.Formatter(plain_format)
        handler.setFormatter(formatter)

    root.addHandler(handler)


class _HelpFormatter(
    argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter
):
    """Preserve example formatting and omit empty defaults."""

    def _get_help_string(self, action):
        text = action.help or ""
        if "%(default)" in text or "default:" in text or "default " in text:
            return text
        if action.default in (None, False, argparse.SUPPRESS):
            return text
        return super()._get_help_string(action)


class RunPhase:
    """Show progress when enabled and log elapsed time."""

    def __init__(self, enabled: bool, label: str):
        self.enabled = bool(enabled and console is not None)
        self.label = label
        self.started_at = None
        self.status = None

    def __enter__(self):
        self.started_at = time.perf_counter()
        if self.enabled:
            self.status = console.status(f"[bold]{self.label}…", spinner="dots")
            self.status.start()
        return self

    def __exit__(self, exc_type, exc, tb):
        elapsed = time.perf_counter() - self.started_at
        if self.enabled and self.status:
            self.status.stop()
        LOG.info("%s completed in %.3fs", self.label, elapsed)


def build_parser() -> argparse.ArgumentParser:
    """Define clustering flags and command examples."""
    epilog = textwrap.dedent(
        """\
        examples:
          # exact k=5
          phytclust data/tree.nwk -k 5 --save-fig

          # global validity-index peaks: top 3, cap k at 200
          phytclust data/tree.nwk --top-n 3 --max-k 200 --save-fig --out-dir out

          # one peak per 4 log-bins, headless save
          phytclust data/tree.nwk --resolution --bins 4 --no-plot --save-all-k --save-fig

          # launch the web GUI (see `phytclust gui --help`)
          phytclust gui --port 8000
        """
    )

    parser = argparse.ArgumentParser(
        prog="phytclust",
        description="Monophyletic clustering of phylogenetic trees (dynamic programming).",
        formatter_class=_HelpFormatter,
        epilog=epilog,
    )

    parser.add_argument(
        "tree", type=_existing_path_or_stdin, help="Newick file or '-' for stdin."
    )

    parser.add_argument(
        "-o",
        "--out-dir",
        type=pathlib.Path,
        default=pathlib.Path("results"),
        help="Directory for PNGs / tsv (created if needed).",
    )
    parser.add_argument("--outgroup", help="Taxon to exclude from all clusters.")
    parser.add_argument(
        "--root-taxon",
        dest="root_taxon",
        help="Taxon to root tree on (or 'midpoint' for midpoint rooting). If not specified, tree is assumed rooted.",
    )
    parser.add_argument(
        "--no-outlier",
        action="store_true",
        help="Do NOT mark singleton clusters as -1 in tsv.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase logging verbosity (-vv for DEBUG).",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="count",
        default=0,
        help="Reduce logging verbosity (-qq for CRITICAL).",
    )
    parser.add_argument(
        "--no-color", action="store_true", help="Disable colorized logging."
    )
    parser.add_argument(
        "--log-format",
        default="%(message)s",
        help="Logging format (plain logging.Formatter style).",
    )
    parser.add_argument(
        "--time",
        action="store_true",
        help="Add the total runtime to the final summary.",
    )
    parser.add_argument(
        "--progress",
        action="store_true",
        help="Show progress spinners even when stderr is not a terminal (requires rich).",
    )
    parser.add_argument(
        "--config",
        type=_existing_path,
        help="Optional JSON/YAML config with plotting/saving options.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"phytclust {_package_version()}",
    )

    parser.add_argument(
        "-k",
        "--k",
        type=_min_int(1),
        help="Exact k-way partition. Cannot be combined with --resolution or --top-n greater than one.",
    )
    parser.add_argument(
        "--top-n",
        type=_positive_int,
        default=1,
        help="Number of clustering validity index peaks to report (per bin with --resolution). Cannot be combined with -k.",
    )
    parser.add_argument(
        "--resolution",
        action="store_true",
        help="Multi-resolution mode: one peak per log-spaced bin.",
    )
    parser.add_argument(
        "--bins",
        type=_positive_int,
        default=3,
        help="Number of log bins when using --resolution.",
    )
    parser.add_argument(
        "--max-k",
        type=_min_int(4),
        help="Upper bound for k in peak search (global or resolution modes).",
    )
    parser.add_argument(
        "--max-k-limit",
        type=float,
        help=(
            "If --max-k is not given, set max_k = ceil(max_k_limit * num_leaves). "
            "Ignored when -k is specified."
        ),
    )
    parser.add_argument(
        "--prominence-weight",
        type=float,
        dest="prominence_weight",
        default=None,
        metavar="W",
        help="How peaks are ranked: 1.0 = purely by prominence (how much a "
        "peak stands out from its neighbours), 0.0 = purely by absolute "
        "score. Only used when peak.ranking_mode is 'adjusted' (the "
        "default); ignored for 'raw'. "
        "(default: 0.7, from config file or PeakConfig if unset).",
    )
    parser.add_argument(
        "--exclude-k2",
        dest="exclude_k2",
        action="store_true",
        default=False,
        help="Drop k=2 from automatic peak selection; top_n / resolution "
        "mode then pick the next-ranked peak instead. An explicit -k 2 "
        "is always honoured.",
    )
    parser.add_argument(
        "--include-k2",
        dest="include_k2_deprecated",
        action="store_true",
        default=False,
        help=argparse.SUPPRESS,
    )

    parser.add_argument(
        "--min-cluster-size",
        type=_min_int(1),
        default=1,
        dest="min_cluster_size",
        help="Minimum number of leaves allowed per cluster (hard constraint).",
    )
    parser.add_argument(
        "--outlier-size-threshold",
        type=_min_int(1),
        dest="outlier_size_threshold",
        help="Count clusters smaller than this as outliers; also use this threshold for TSV labels.",
    )
    parser.add_argument(
        "--prefer-fewer-outliers",
        action="store_true",
        dest="prefer_fewer_outliers",
        default=False,
        help="At fixed k, minimise outlier count first, then break ties by raw cost.",
    )
    parser.add_argument(
        "--polytomy-mode",
        choices=("hard", "soft"),
        default=None,
        help="Polytomy DP mode (default: soft). Soft lets any group of a "
        "multifurcation's children form a cluster; hard forbids "
        "cross-child subset merges and has no degree limit.",
    )
    parser.add_argument(
        "--soft-polytomy-max-degree",
        type=_min_int(2),
        default=None,
        dest="soft_polytomy_max_degree",
        help="Maximum polytomy degree allowed in soft mode (guardrail for exponential complexity).",
    )
    parser.add_argument(
        "--no-split-zero-length",
        action="store_true",
        dest="no_split_zero_length",
        default=False,
        help="Prevent splitting zero-length branches (blocks DP from separating clades with all-zero edges).",
    )

    _add_common_run_flags(parser)

    return parser


def main(argv=None) -> int:
    """Run the CLI with Python warnings sent through its logger."""
    with warnings.catch_warnings():

        def report_warning(message, category, filename, lineno, file=None, line=None):
            if not issubclass(category, DeprecationWarning):
                LOG.warning("%s", message)

        warnings.showwarning = report_warning
        return _run_cli(argv)


def _run_cli(argv=None) -> int:
    """Parse arguments, run clustering, and report results."""
    if argv is None:
        argv = sys.argv[1:]

    if argv and argv[0] == "gui":
        from ..gui.launch import main as gui_main

        return gui_main(argv[1:])

    parser = build_parser()
    args = parser.parse_args(argv)
    global console, output_console
    console = Console(stderr=True, no_color=args.no_color) if Console else None
    output_console = Console(no_color=args.no_color) if Console else None
    if args.k is not None and (args.resolution or args.top_n != 1):
        parser.error(
            "-k cannot be combined with --resolution or --top-n greater than one."
        )
    if args.exclude_k2 and args.include_k2_deprecated:
        parser.error("--exclude-k2 cannot be combined with --include-k2.")
    try:
        _configure_logging(args.verbose, args.quiet, args.no_color, args.log_format)
    except ValueError as exc:
        print(f"Invalid logging format: {exc}", file=sys.stderr)
        return 2

    if args.include_k2_deprecated:
        LOG.warning(
            "--include-k2 is deprecated: k=2 is now a candidate by default. "
            "Use --exclude-k2 to drop it."
        )

    show_phase_ui = bool(console is not None and (sys.stderr.isatty() or args.progress))

    if sys.stdout.isatty() and not args.quiet:
        _emit_header(args)

    if (
        args.plot
        and sys.platform.startswith("linux")
        and not (os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))
    ):
        LOG.warning(
            "No display found; ignoring --plot. Use --save-fig to write figures."
        )
        args.plot = False

    try:
        config = _load_config(args.config)
        runtime_overrides = _runtime_overrides_from_config(config)
        runtime_config = build_runtime_config(runtime_overrides)
        algorithm_config = _config_section(config, "algorithm")
        polytomy_mode = args.polytomy_mode or algorithm_config.get(
            "polytomy_mode", "soft"
        )
        soft_polytomy_max_degree = (
            args.soft_polytomy_max_degree
            if args.soft_polytomy_max_degree is not None
            else algorithm_config.get("soft_polytomy_max_degree", 12)
        )
        cli_peak_overrides = {}
        if args.prominence_weight is not None:
            cli_peak_overrides["prominence_weight"] = args.prominence_weight
        if args.exclude_k2:
            cli_peak_overrides["exclude_k2"] = True
        peak_config = (
            _peak_config_from_config(config, cli_peak_overrides)
            if args.k is None
            else None
        )
    except ConfigurationError as exc:
        LOG.error("%s", exc)
        return 2

    plot_settings = {
        "width_scale": runtime_config.plot.cluster.width_scale,
        "marker_size": runtime_config.plot.cluster.marker_size,
        "show_branch_lengths": runtime_config.plot.cluster.show_branch_lengths,
        "hide_internal_nodes": runtime_config.plot.cluster.hide_internal_nodes,
    }

    if "height_scale" in runtime_overrides.get("plot", {}).get("cluster", {}):
        plot_settings["height_scale"] = runtime_config.plot.cluster.height_scale
    plot_settings["cmap"] = runtime_config.plot.cluster.cmap

    save_default_filename = runtime_config.save.tsv_name
    save_default_outlier = runtime_config.save.outlier

    started_at = time.perf_counter()

    will_write_anything = bool(args.save_fig or args.save_tree or (not args.no_tsv))
    if will_write_anything:
        out_dir = Path(args.out_dir) if args.out_dir is not None else Path("results")
        try:
            out_dir.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            LOG.error("Cannot create output directory %s: %s", out_dir, exc)
            return 2
        args.out_dir = str(out_dir)

    try:
        with RunPhase(show_phase_ui, "load tree"):
            handle: Any = sys.stdin if args.tree == "-" else args.tree
            tree = Phylo.read(handle, "newick")
    except Exception as exc:
        LOG.error("Cannot read tree: %s", exc)
        return 2

    try:
        with RunPhase(show_phase_ui, "initialization"):
            outlier_config = OutlierConfig(
                size_threshold=args.outlier_size_threshold,
                prefer_fewer=args.prefer_fewer_outliers,
            )
            pc = PhytClust(
                tree=tree,
                outgroup=args.outgroup,
                root_taxon=args.root_taxon,
                min_cluster_size=args.min_cluster_size,
                outlier=outlier_config,
                polytomy_mode=polytomy_mode,
                soft_polytomy_max_degree=soft_polytomy_max_degree,
                no_split_zero_length=args.no_split_zero_length,
                runtime_config=runtime_config,
            )
        LOG.info("Tree terminals (after outgroup handling): %d", pc.num_terminals)
    except Exception as exc:
        LOG.error("Failed to initialize PhytClust: %s", exc)
        return 1

    try:
        with RunPhase(show_phase_ui, "clustering"):
            LOG.info(
                "Running PhytClust (k=%s, top_n=%d, resolution=%s, max_k=%s)…",
                args.k or "auto",
                args.top_n,
                args.resolution,
                args.max_k or "auto",
            )

            pc.peak_config = peak_config
            run_kwargs = dict(
                k=args.k,
                top_n=args.top_n,
                by_resolution=args.resolution,
                num_bins=args.bins,
                max_k=args.max_k,
                max_k_limit=args.max_k_limit,
                plot_scores=args.k is None and (args.plot or args.save_fig),
                peak_config=peak_config,
            )

            run_started_at = time.perf_counter()
            result = pc.run(**run_kwargs)
            run_elapsed = time.perf_counter() - run_started_at

            LOG.debug("Clustering took %.3fs", run_elapsed)

            selected_cluster_counts = []
            if isinstance(result, dict):
                if "k_values" in result and result["k_values"] is not None:
                    selected_cluster_counts = [int(x) for x in result["k_values"]]
                elif "ks" in result and result["ks"] is not None:
                    selected_cluster_counts = [int(x) for x in result["ks"]]
                elif "k" in result and result["k"] is not None:
                    selected_cluster_counts = [int(result["k"])]

            if selected_cluster_counts:
                LOG.info(
                    "Selected k: %s",
                    selected_cluster_counts
                    if len(selected_cluster_counts) > 1
                    else selected_cluster_counts[0],
                )

    except (ConfigurationError, ValueError) as e:
        LOG.error("Clustering failed: %s", e)
        return 1
    except Exception as exc:
        LOG.error("Unexpected error during clustering: %s", exc)
        return 1

    plotting_failed = False
    try:
        if args.plot or args.save_fig or args.save_tree:
            with RunPhase(show_phase_ui, "render/plot"):
                pc.plot(
                    results_dir=(
                        args.out_dir if (args.save_fig or args.save_tree) else None
                    ),
                    save=(args.save_fig or args.save_tree),
                    k=args.k,
                    top_n=max(1, len(selected_cluster_counts)),
                    dpi=args.dpi,
                    **plot_settings,
                )
                scores_fig = getattr(pc, "plot_of_scores", None)
                if args.save_fig and scores_fig is not None:
                    scores_fig.savefig(
                        Path(args.out_dir) / "scores.png",
                        dpi=args.dpi,
                        bbox_inches="tight",
                    )
    except Exception as exc:
        LOG.error("Could not create the requested plots: %s", exc)
        plotting_failed = True

    if not args.no_tsv:
        try:
            with RunPhase(show_phase_ui, "write tsv"):
                pc.save(
                    results_dir=args.out_dir,
                    filename=(
                        args.tsv_name
                        if args.tsv_name is not None
                        else save_default_filename
                    ),
                    outlier=(False if args.no_outlier else save_default_outlier),
                    output_all=args.save_all_k,
                    top_n=max(1, len(selected_cluster_counts)),
                )
        except Exception as exc:
            LOG.error("Failed to write tsv: %s", exc)
            return 1

    if will_write_anything and not args.quiet:
        _emit_artifacts(args.out_dir)

    total = time.perf_counter() - started_at
    LOG.info("Total runtime: %.3fs", total)
    if not args.quiet:
        _emit_summary(
            pc,
            selected_cluster_counts,
            args.out_dir if will_write_anything else None,
            total if args.time else None,
        )

    if not plotting_failed:
        LOG.info("Done.")
    return 1 if plotting_failed else 0


def _emit_summary(
    pc: Any, selected_cluster_counts: list[int], out_dir: Any, total: Any
) -> None:
    """Print selected partitions and their output location."""
    rows = []
    threshold = (
        pc.outlier.size_threshold if pc.outlier.size_threshold is not None else 2
    )
    for rank, count in enumerate(selected_cluster_counts, start=1):
        cluster_map = pc.clusters.get(count)
        if not cluster_map:
            continue
        sizes: dict[int, int] = {}
        for cluster_id in cluster_map.values():
            sizes[cluster_id] = sizes.get(cluster_id, 0) + 1
        outliers = sum(size < threshold for size in sizes.values())
        alpha = pc.alpha_info(count, cluster_map)["alpha"]
        rows.append(
            (str(rank), str(count), str(len(sizes)), str(outliers), f"{alpha:.3g}")
        )
    if output_console is not None and sys.stdout.isatty():
        from rich import box
        from rich.table import Table

        table = Table(
            title="Selected partitions",
            title_style="bold",
            box=box.SIMPLE_HEAD,
            header_style="bold",
            padding=(0, 2),
        )
        for name in ("Rank", "k", "Clusters", "Outliers", "Alpha"):
            table.add_column(name, justify="right")
        for row in rows:
            table.add_row(*row)
        if rows:
            output_console.print(table)
        else:
            output_console.print("No partitions selected.")
        if out_dir is not None:
            output_console.print(Text(f"Results: {Path(out_dir).resolve()}"))
        if total is not None:
            output_console.print(Text(f"Elapsed: {total:.2f}s"))
    else:
        if rows:
            print("Rank  k  Clusters  Outliers  Alpha")
            for row in rows:
                print("  ".join(row))
        else:
            print("No partitions selected.")
        if out_dir is not None:
            print(f"Results: {Path(out_dir).resolve()}")
        if total is not None:
            print(f"Elapsed: {total:.2f}s")


def _emit_header(args: argparse.Namespace) -> None:
    """Show the input and selection mode before a terminal run."""
    mode = (
        f"exact k = {args.k}"
        if args.k is not None
        else ("resolution peaks" if args.resolution else "global peaks")
    )
    if console is not None and Text is not None:
        header = Text(f"PhytClust {_package_version()}", style="bold")
        header.append(f"  ·  {mode}")
        console.print(header)
        console.print(Text(f"Tree: {args.tree}"))
    else:
        print(f"PhytClust {_package_version()} | {mode}", file=sys.stderr)
        print(f"Tree: {args.tree}", file=sys.stderr)


def _emit_artifacts(out_dir: str) -> None:
    """List output files in terminals; keep redirected summaries compact."""
    if not sys.stdout.isatty():
        return
    try:
        directory = Path(out_dir)
        files = sorted(path.name for path in directory.iterdir() if path.is_file())
    except OSError as exc:
        LOG.debug("Could not list output files: %s", exc)
        return
    if files:
        message = "Files: " + ", ".join(files)
        if output_console is not None and Text is not None:
            output_console.print(Text(message))
        else:
            print(message)
