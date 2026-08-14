#!/usr/bin/env python3
from __future__ import annotations

import argparse
import pathlib
import sys
import textwrap
import logging
import time
import warnings
from typing import Any

from Bio import Phylo
from ..algo.core import PhytClust
from ..config import OutlierConfig, PeakConfig, build_runtime_config
from ..exceptions import ConfigurationError
from pathlib import Path
from importlib import resources

try:
    from importlib.metadata import version, PackageNotFoundError
except Exception:
    try:
        from importlib_metadata import version, PackageNotFoundError
    except Exception:
        version = None
        PackageNotFoundError = Exception

try:
    from rich.console import Console
    from rich.status import Status
    from rich.text import Text
    from rich.panel import Panel
except Exception:
    Console = None
    Status = None
    Text = None
    Panel = None

console = Console() if Console else None

LOG = logging.getLogger("phytclust.cli")


def print_banner():
    try:
        text = (
            resources.files("phytclust")
            .joinpath("assets", "ascii_logo.txt")
            .read_text(encoding="utf-8")
        )
        print(text)
    except Exception:
        LOG.debug("ascii_logo.txt not found; skipping banner.")


def _positive_int(value: str) -> int:
    ivalue = int(value)
    if ivalue < 1:
        raise argparse.ArgumentTypeError("value must be ≥ 1")
    return ivalue


def _min_int(min_value: int):
    def _check(value: str) -> int:
        ivalue = int(value)
        if ivalue < min_value:
            raise argparse.ArgumentTypeError(f"value must be ≥ {min_value}")
        return ivalue

    return _check


def _existing_path_or_stdin(p: str) -> pathlib.Path | str:
    if p == "-":
        return p
    path = pathlib.Path(p)
    if not path.exists():
        raise argparse.ArgumentTypeError(f"file not found: {p}")
    return path


def _existing_path(p: str) -> pathlib.Path:
    path = pathlib.Path(p)
    if not path.exists():
        raise argparse.ArgumentTypeError(f"file not found: {p}")
    return path


def _package_version() -> str:
    if version is None:
        return "unknown"
    try:
        return version("phytclust")
    except PackageNotFoundError:
        return "unknown"


def _load_config(path: pathlib.Path | None) -> dict:
    """Load a YAML/JSON config file into a dict.

    Raises ConfigurationError (rather than silently returning {}) on a missing,
    unparseable, or non-mapping config, so a typo can't silently fall back to a
    full-defaults run. yaml.safe_load parses JSON too, so no separate branch.
    """
    if not path:
        return {}
    try:
        text = path.read_text()
    except OSError as exc:
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
    """Announce a renamed config key on both channels.

    ``DeprecationWarning`` is what programmatic callers and tests can catch;
    Python hides it from end users by default, so the log line is what actually
    reaches someone running the CLI.
    """
    message = f"config key '{old}' is deprecated, use '{new}' instead."
    warnings.warn(message, DeprecationWarning, stacklevel=3)
    LOG.warning("config: %s", message)


def _runtime_overrides_from_cfg(raw_cfg: dict[str, Any]) -> dict[str, Any]:
    """Map user config file keys to RuntimeConfig overrides."""
    out: dict[str, Any] = {}

    plot_block = dict(raw_cfg.get("plot") or {})
    scores_block = dict(raw_cfg.get("scores_plot") or {})
    save_block = dict(raw_cfg.get("save") or {})

    if plot_block:
        cluster = {
            "width_scale": plot_block.get("width_scale"),
            "height_scale": plot_block.get("height_scale"),
            "marker_size": plot_block.get("marker_size"),
            "show_branch_lengths": plot_block.get("show_branch_lengths"),
            "hide_internal_nodes": plot_block.get("hide_internal_nodes"),
        }
        cluster = {k: v for k, v in cluster.items() if v is not None}
        if cluster:
            out.setdefault("plot", {})["cluster"] = cluster

    if scores_block:
        out.setdefault("plot", {})["scores"] = dict(scores_block)

    if save_block:
        tsv_name = save_block.get("tsv_name")
        if tsv_name is None and save_block.get("csv_name") is not None:
            # Renamed: the output was always tab-separated.
            _warn_deprecated_key("save.csv_name", "save.tsv_name")
            tsv_name = save_block.get("csv_name")
        save_map = {
            "tsv_name": tsv_name,
            "outlier": save_block.get("outlier"),
        }
        save_map = {k: v for k, v in save_map.items() if v is not None}
        if save_map:
            out["save"] = save_map

    # New-style nested runtime overrides can be passed directly.
    runtime_block = raw_cfg.get("runtime")
    if isinstance(runtime_block, dict):
        out.update(runtime_block)

    return out


#: Deprecated config-file keys under `peak`, mapped old -> new.
_DEPRECATED_PEAK_KEYS = {"lambda_weight": "prominence_weight"}


def _peak_config_from_cfg(
    raw_cfg: dict[str, Any],
    cli_overrides: dict[str, Any] | None = None,
) -> PeakConfig:
    """Build a PeakConfig with precedence: explicit CLI flags > config file > defaults.

    Start from dataclass defaults, overlay the config file's ``peak`` block, then
    overlay only the CLI flags the user actually set (``cli_overrides``). This is
    the standard precedence; the previous version let the config file clobber
    explicit command-line flags.
    """
    peak_cfg = PeakConfig()

    peak_block = raw_cfg.get("peak")
    if not isinstance(peak_block, dict):
        peak_block = (raw_cfg.get("algorithm") or {}).get("peak")
    if isinstance(peak_block, dict):
        for key, val in peak_block.items():
            if key in _DEPRECATED_PEAK_KEYS:
                new_key = _DEPRECATED_PEAK_KEYS[key]
                if new_key in peak_block:
                    continue  # new key wins; no need to warn twice
                _warn_deprecated_key(f"peak.{key}", f"peak.{new_key}")
                key = new_key
            if hasattr(peak_cfg, key):
                setattr(peak_cfg, key, val)
            else:
                LOG.warning("config: unknown key 'peak.%s' ignored.", key)

    for key, val in (cli_overrides or {}).items():
        if val is not None and hasattr(peak_cfg, key):
            setattr(peak_cfg, key, val)

    # Overlays bypass __post_init__, so re-check here — this reports a bad
    # config value before the DP runs rather than partway through it.
    peak_cfg.validate()
    return peak_cfg


def _add_common_run_flags(sp: argparse.ArgumentParser) -> None:
    sp.add_argument(
        "--plot",
        dest="plot",
        default=True,
        action=argparse.BooleanOptionalAction,
        help="Display interactive plots after running.",
    )
    sp.add_argument(
        "--save-fig",
        action="store_true",
        help="Write scores.png and coloured tree PNG(s) to --out-dir.",
    )
    sp.add_argument(
        "--save-tree",
        action="store_true",
        help="Explicitly write coloured tree PNG(s) (in addition to --save-fig).",
    )
    sp.add_argument(
        "--save-all-k",
        action="store_true",
        help="Write tsv rows for *every* k from 1..max_k (can be large!).",
    )
    sp.add_argument(
        "--tsv-name",
        default=None,
        help="Output tsv filename (default: from config or 'phytclust_results.tsv').",
    )
    sp.add_argument(
        "--no-tsv", action="store_true", help="Skip writing the results tsv."
    )
    sp.add_argument(
        "--dpi", type=_positive_int, default=150, help="PNG resolution (DPI)."
    )


def _configure_logging(
    verbosity: int, quiet: int, no_color: bool, log_format: str
) -> None:
    """
    Configure a single handler for the 'phytclust' logger hierarchy so that
    all submodules (phytclust.*, using logging.getLogger(__name__)) inherit it.
    """
    level = logging.WARNING
    if verbosity > 0:
        level = logging.INFO
    if verbosity > 1:
        level = logging.DEBUG
    if quiet > 0:
        level = logging.ERROR
    if quiet > 1:
        level = logging.CRITICAL

    root = logging.getLogger("phytclust")
    root.handlers.clear()
    root.setLevel(level)
    root.propagate = False

    handler = None
    if not no_color:
        try:
            from rich.logging import RichHandler

            handler = RichHandler(
                level=level,
                markup=True,
                rich_tracebacks=True,
                show_time=True,
                show_level=True,
                show_path=False,
            )
            formatter = logging.Formatter(log_format)
            handler.setFormatter(formatter)
        except Exception:
            handler = None

    if handler is None:
        handler = logging.StreamHandler()
        formatter = logging.Formatter(log_format)
        handler.setFormatter(formatter)

    root.addHandler(handler)


class phase:
    """
    Spinner + timing in one place. Shows a 'rich' spinner if --progress is set
    and rich is installed. Always logs elapsed time at INFO.
    """

    def __init__(self, enabled: bool, label: str):
        self.enabled = bool(enabled and console is not None)
        self.label = label
        self.t0 = None
        self.status = None

    def __enter__(self):
        self.t0 = time.perf_counter()
        if self.enabled:
            self.status = console.status(f"[bold]{self.label}…", spinner="dots")
            self.status.start()
        return self

    def __exit__(self, exc_type, exc, tb):
        dt = time.perf_counter() - self.t0
        if self.enabled and self.status:
            self.status.stop()
        LOG.info("⏱ %s: %.3fs", self.label, dt)


def build_parser() -> argparse.ArgumentParser:
    epilog = textwrap.dedent(
        """\
        examples:
          # exact k=5
          phytclust data/tree.nwk -k 5 --save-fig

          # global CalBow peaks: top 3, cap k at 200
          phytclust data/tree.nwk --top-n 3 --max-k 200 --save-fig --out-dir out

          # one peak per 4 log-bins, headless save
          phytclust data/tree.nwk --resolution --bins 4 --no-plot --save-all-k --save-fig

          # launch the web GUI (see `phytclust gui --help`)
          phytclust gui --port 8000
        """
    )

    p = argparse.ArgumentParser(
        prog="phytclust",
        description="Monophyletic clustering of phylogenetic trees (dynamic programming).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=epilog,
    )

    p.add_argument(
        "tree", type=_existing_path_or_stdin, help="Newick file or '-' for stdin."
    )

    p.add_argument(
        "-o",
        "--out-dir",
        type=pathlib.Path,
        default=pathlib.Path("results"),
        help="Directory for PNGs / tsv (created if needed).",
    )
    p.add_argument("--outgroup", help="Taxon to exclude from all clusters.")
    p.add_argument(
        "--root-taxon",
        dest="root_taxon",
        help="Taxon to root tree on (or 'midpoint' for midpoint rooting). If not specified, tree is assumed rooted.",
    )
    p.add_argument(
        "--no-outlier",
        action="store_true",
        help="Do NOT mark singleton clusters as -1 in tsv.",
    )
    p.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase logging verbosity (-vv for DEBUG).",
    )
    p.add_argument(
        "-q",
        "--quiet",
        action="count",
        default=0,
        help="Reduce logging verbosity (-qq for CRITICAL).",
    )
    p.add_argument("--no-color", action="store_true", help="Disable colorized logging.")
    p.add_argument(
        "--log-format",
        default="%(message)s",
        help="Logging format (plain logging.Formatter style).",
    )
    p.add_argument(
        "--time",
        action="store_true",
        help="Report elapsed time per phase and total.",
    )
    p.add_argument(
        "--progress",
        action="store_true",
        help="Show progress spinners for major phases (requires rich).",
    )
    p.add_argument(
        "--config",
        type=_existing_path,
        help="Optional JSON/YAML config with plotting/saving options.",
    )
    p.add_argument(
        "--version",
        action="version",
        version=f"phytclust { _package_version() }",
    )

    p.add_argument(
        "-k",
        "--k",
        type=_min_int(2),
        help="Exact k-way partition. If given, overrides --top-n/--resolution.",
    )
    p.add_argument(
        "--top-n",
        type=_positive_int,
        default=1,
        help="Number of global CalBow peaks to return (ignored if -k is set or --resolution is used).",
    )
    p.add_argument(
        "--resolution",
        action="store_true",
        help="Multi-resolution mode: one peak per log-spaced bin.",
    )
    p.add_argument(
        "--bins",
        type=_positive_int,
        default=3,
        help="Number of log bins when using --resolution.",
    )
    p.add_argument(
        "--max-k",
        type=_min_int(4),
        help="Upper bound for k in peak search (global or resolution modes).",
    )
    p.add_argument(
        "--max-k-limit",
        type=float,
        help=(
            "If --max-k is not given, set max_k = ceil(max_k_limit * num_leaves). "
            "Ignored when -k is specified."
        ),
    )
    p.add_argument(
        "--prominence-weight",
        "--lambda-weight",  # deprecated alias, scheduled for removal in 2.0.0
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
    p.add_argument(
        "--include-k2",
        dest="exclude_k2",
        action="store_false",
        default=True,
        help="Allow k=2 as a peak. By default k=2 is excluded; top_n / "
             "resolution mode then pick the next-ranked peak instead.",
    )

    # Outlier and cluster size constraints
    p.add_argument(
        "--min-cluster-size",
        type=_min_int(1),
        default=1,
        dest="min_cluster_size",
        help="Minimum number of leaves allowed per cluster (hard constraint).",
    )
    p.add_argument(
        "--outlier-size-threshold",
        type=_min_int(1),
        dest="outlier_size_threshold",
        help="Clusters smaller than this are marked as outliers in output and penalized in DP.",
    )
    p.add_argument(
        "--prefer-fewer-outliers",
        action="store_true",
        dest="prefer_fewer_outliers",
        default=False,
        help="At fixed k, minimise outlier count first, then break ties by raw cost.",
    )
    p.add_argument(
        "--no-optimize-polytomies",
        action="store_false",
        dest="optimize_polytomies",
        default=True,
        help="Use greedy dummy-node polytomy resolution instead of optimal sequential convolution.",
    )
    p.add_argument(
        "--polytomy-mode",
        choices=("hard", "soft"),
        default=None,
        help="Polytomy DP mode: hard forbids cross-child subset merges; soft allows subset merges.",
    )
    p.add_argument(
        "--soft-polytomy-max-degree",
        type=_min_int(2),
        default=None,
        dest="soft_polytomy_max_degree",
        help="Maximum polytomy degree allowed in soft mode (guardrail for exponential complexity).",
    )
    p.add_argument(
        "--no-split-zero-length",
        action="store_true",
        dest="no_split_zero_length",
        default=False,
        help="Prevent splitting zero-length branches (blocks DP from separating clades with all-zero edges).",
    )

    _add_common_run_flags(p)

    return p


def main(argv=None) -> int:
    if argv is None:
        argv = sys.argv[1:]

    # `phytclust gui [...]` launches the web GUI. Intercepted before the main
    # (tree-clustering) parser, whose required positional `tree` would otherwise
    # reject the subcommand.
    if argv and argv[0] == "gui":
        from ..gui.launch import main as gui_main

        return gui_main(argv[1:])

    args = build_parser().parse_args(argv)
    _configure_logging(args.verbose, args.quiet, args.no_color, args.log_format)

    # argparse accepts the deprecated alias silently, so flag it here.
    if any(a == "--lambda-weight" or a.startswith("--lambda-weight=") for a in argv):
        warnings.warn(
            "--lambda-weight is deprecated, use --prominence-weight instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        LOG.warning(
            "--lambda-weight is deprecated, use --prominence-weight instead."
        )

    show_phase_ui = bool(console is not None and sys.stderr.isatty())

    print_banner()
    _emit_header(args)

    try:
        cfg = _load_config(args.config)
    except ConfigurationError as exc:
        LOG.error("%s", exc)
        sys.exit(2)
    runtime_cfg = build_runtime_config(_runtime_overrides_from_cfg(cfg))
    algorithm_cfg = (
        cfg.get("algorithm") if isinstance(cfg.get("algorithm"), dict) else {}
    )

    polytomy_mode = args.polytomy_mode or algorithm_cfg.get("polytomy_mode", "hard")
    soft_polytomy_max_degree = (
        args.soft_polytomy_max_degree
        if args.soft_polytomy_max_degree is not None
        else int(algorithm_cfg.get("soft_polytomy_max_degree", 12))
    )

    plot_cfg = {
        "width_scale": runtime_cfg.plot.cluster.width_scale,
        "height_scale": runtime_cfg.plot.cluster.height_scale,
        "marker_size": runtime_cfg.plot.cluster.marker_size,
        "show_branch_lengths": runtime_cfg.plot.cluster.show_branch_lengths,
        "hide_internal_nodes": runtime_cfg.plot.cluster.hide_internal_nodes,
    }

    save_default_filename = runtime_cfg.save.tsv_name
    save_default_outlier = runtime_cfg.save.outlier

    t0_total = time.perf_counter()

    will_write_anything = bool(args.save_fig or args.save_tree or (not args.no_tsv))
    if will_write_anything:
        out_dir = Path(args.out_dir) if args.out_dir is not None else Path("results")
        out_dir.mkdir(parents=True, exist_ok=True)
        args.out_dir = str(out_dir)

    try:
        with phase(show_phase_ui, "load tree"):
            handle: Any = sys.stdin if args.tree == "-" else args.tree
            tree = Phylo.read(handle, "newick")
    except Exception as exc:
        LOG.error("Cannot read tree: %s", exc)
        sys.exit(2)

    try:
        with phase(show_phase_ui, "initialization"):
            outlier_cfg = OutlierConfig(
                size_threshold=args.outlier_size_threshold,
                prefer_fewer=args.prefer_fewer_outliers,
            )
            pc = PhytClust(
                tree=tree,
                outgroup=args.outgroup,
                root_taxon=args.root_taxon,
                min_cluster_size=args.min_cluster_size,
                outlier=outlier_cfg,
                optimize_polytomies=args.optimize_polytomies,
                polytomy_mode=polytomy_mode,
                soft_polytomy_max_degree=soft_polytomy_max_degree,
                no_split_zero_length=args.no_split_zero_length,
                runtime_config=runtime_cfg,
            )
        LOG.info("Tree terminals (after outgroup handling): %d", pc.num_terminals)
    except Exception as exc:
        LOG.error("Failed to initialize PhytClust: %s", exc)
        sys.exit(1)

    try:
        with phase(show_phase_ui, "clustering"):
            LOG.info(
                "Running PhytClust (k=%s, top_n=%d, resolution=%s, max_k=%s)…",
                args.k or "auto",
                args.top_n,
                args.resolution,
                args.max_k or "auto",
            )

            # Explicit CLI flags win over config/defaults. prominence_weight
            # defaults to None (unset); exclude_k2 is only reachable as False
            # via --include-k2, so treat False as an explicit override.
            cli_peak_overrides: dict[str, Any] = {}
            if args.prominence_weight is not None:
                cli_peak_overrides["prominence_weight"] = args.prominence_weight
            if not args.exclude_k2:
                cli_peak_overrides["exclude_k2"] = False
            peak_cfg = _peak_config_from_cfg(cfg, cli_peak_overrides)
            pc.peak_config = peak_cfg
            run_kwargs = dict(
                k=args.k,
                top_n=args.top_n,
                by_resolution=args.resolution,
                num_bins=args.bins,
                max_k=args.max_k,
                max_k_limit=args.max_k_limit,
                plot_scores=args.plot,
                peak_config=peak_cfg,
            )

            t0_run = time.perf_counter()
            result = pc.run(**run_kwargs)
            dt_run = time.perf_counter() - t0_run

            LOG.info("⏱ phytclust.run(): %.3fs", dt_run)

            ks = []
            if isinstance(result, dict):
                if "k_values" in result and result["k_values"] is not None:
                    ks = [int(x) for x in result["k_values"]]
                elif "ks" in result and result["ks"] is not None:
                    ks = [int(x) for x in result["ks"]]
                elif "k" in result and result["k"] is not None:
                    ks = [int(result["k"])]

            if ks:
                LOG.info("Selected k: %s", ks if len(ks) > 1 else ks[0])

    except (ConfigurationError, ValueError) as e:
        LOG.error("Clustering failed: %s", e)
        sys.exit(1)
    except Exception as exc:
        LOG.error("Unexpected error during clustering: %s", exc)
        sys.exit(1)

    try:
        if args.plot or args.save_fig or args.save_tree:
            with phase(show_phase_ui, "render/plot"):
                pc.plot(
                    results_dir=(
                        args.out_dir if (args.save_fig or args.save_tree) else None
                    ),
                    save=(args.save_fig or args.save_tree),
                    k=args.k,
                    **plot_cfg,
                )
    except Exception as exc:
        LOG.warning("Plotting encountered a problem: %s", exc)

    if not args.no_tsv:
        try:
            with phase(show_phase_ui, "write tsv"):
                pc.save(
                    results_dir=args.out_dir,
                    filename=(
                        args.tsv_name
                        if args.tsv_name is not None
                        else save_default_filename
                    ),
                    outlier=(False if args.no_outlier else save_default_outlier),
                    output_all=args.save_all_k,
                )
        except Exception as exc:
            LOG.error("Failed to write tsv: %s", exc)
            sys.exit(1)

    if will_write_anything:
        _emit_artifacts(args.out_dir)

    LOG.info("Total runtime: %.3fs", time.perf_counter() - t0_total)

    LOG.info("Done.")
    return 0


def _emit_header(args: argparse.Namespace) -> None:
    if console is None or Panel is None or Text is None:
        LOG.info(
            "phytclust | k=%s | top_n=%s | resolution=%s | out_dir=%s",
            args.k or "auto",
            args.top_n,
            args.resolution,
            args.out_dir,
        )
        return

    mode = (
        "k" if args.k is not None else ("resolution" if args.resolution else "global")
    )
    body = Text()
    body.append("Mode: ", style="bold")
    body.append(f"{mode}\n")
    body.append("Tree: ", style="bold")
    body.append(f"{args.tree}\n")
    body.append("Output: ", style="bold")
    body.append(f"{args.out_dir}\n")
    body.append("Flags: ", style="bold")
    body.append(
        f"plot={args.plot} save_fig={args.save_fig} save_tree={args.save_tree} tsv={not args.no_tsv}"
    )

    console.print(Panel(body, title="phytclust", expand=False))


def _emit_artifacts(out_dir: str) -> None:
    try:
        p = Path(out_dir)
        if not p.exists():
            return
        files = sorted([x.name for x in p.iterdir() if x.is_file()])
        if not files:
            return

        if console is None:
            LOG.info("Wrote %d file(s) to %s", len(files), out_dir)
            return

        from rich.table import Table  # type: ignore

        t = Table(title="Output files", show_lines=False)
        t.add_column("File", style="cyan")
        for f in files:
            t.add_row(f)
        console.print(t)
    except Exception:
        return
