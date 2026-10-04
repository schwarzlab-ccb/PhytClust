# Changelog

All notable changes to _PhytClust_ will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/) and the
project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [1.0.0] – 2026-10-04

First stable release, accompanying the PhytClust manuscript.

Aligns the algorithm with the manuscript Methods.

### Added

- Optional peak ranking by fewer singleton outliers or balanced cluster sizes,
  through `partition_preference` and `partition_weight`. This reorders detected
  candidates without changing score curves or the DP objective.
- CLI flags `--peak-partition-preference` and `--peak-partition-weight`.

### Changed

- Balanced split selection for tied binary and polytomy partitions.
- Larger score-plot labels and clearer default cluster plots.
- `OutlierConfig(prefer_fewer=True)` counts singleton clusters when no size
  threshold is supplied.
- **Soft polytomy is now the default** (`polytomy_mode="soft"`, CLI and GUI
  alike). Nodes above `soft_polytomy_max_degree` (12) raise; use
  `polytomy_mode="hard"` / `--polytomy-mode hard` for them.
- **k = 2 is a candidate for automatic selection by default**
  (`PeakConfig.exclude_k2=False`). The CLI flag is now `--exclude-k2`;
  `--include-k2` is accepted with a deprecation warning.
- Equal-cost splits at a hard polytomy minimise the sum of squared cluster
  counts over all its children, not pairwise along the fold.
- Tie tolerances and the EL1 drop guard are relative, so a global rescaling of
  branch lengths no longer changes which partition or *k* is selected.
- α (`cluster_alpha`) no longer counts the root as an extra-cluster edge. After
  outgroup pruning this could add the old root branch to the backbone mean.
- `-k 1` is accepted on the CLI.
- Resolution-mode clade-level bins start at k = 1 instead of the smallest
  candidate (2). Starting at 2 made the first bin run to about
  2·(K/2)^(1/b), which on the avian tree (b = 5) put k = 3 and k = 6 in the
  same bin, so the second coarse level was never reported.
- The GUI has a polytomy-mode selector in place of the removed "Optimize
  polytomies" toggle, and reports configuration errors (such as a polytomy
  above the soft-mode degree limit) instead of a generic server error.

### Fixed

- `run(max_k_limit=...)` was ignored after the first run, and a `max_k` passed
  to one `run()` was reused by later runs.
- `save_tied_optima=True` re-entered `backtrack` recursively (hundreds of times
  per *k*) before the recursion limit stopped it.
- Documentation and the demo notebook no longer use the removed
  `optimize_polytomies` argument.

### Added

- `phytclust gui` command to launch the web GUI. The GUI itself shipped in
  0.1.2; this adds the launcher subcommand, so `uvicorn` is no longer needed.
- `PhytClust.__repr__` for readable REPL and notebook summaries.
- `strategy` argument (`central` / `divergent` / `medoid`) for
  `select_representative_species`.
- Warnings for duplicate taxa and for partial or negative branch lengths.
- `RANKING_MODES` and `RESOLUTION_FALLBACK_MODES` exported as the accepted value
  sets for the corresponding config fields.

### Changed

- **Renamed** `PeakConfig.lambda_weight` to `prominence_weight`, and the CLI flag
  `--lambda-weight` to `--prominence-weight`. The name now says what the value
  does: `1.0` ranks peaks by prominence, `0.0` by absolute score height. Default
  remains `0.7`.
- **Renamed** `SaveConfig.csv_name` to `tsv_name`, matching the long-standing
  `--tsv-name` flag. The output was always tab-separated.
- `ranking_mode` and `resolution_fallback_mode` now validate their values and
  raise `ConfigurationError` on an unrecognised one, instead of silently falling
  back to the default.
- Unknown keys in a config file's `peak` block now log a warning instead of being
  dropped silently.
- CLI flags now override config-file values.
- Invalid or missing `--config` files now raise an error instead of silently
  using defaults.
- Faster Colless index, representative-distance, and cluster/outgroup setup.
- Cluster MRCAs are now found in a single postorder pass rather than one
  `Bio.Phylo` `common_ancestor` call per cluster, which re-walked the tree each
  time. This dominated total runtime on trees of any size — roughly 7x faster
  end to end at 600 leaves and 14x at 1200. Results are unchanged.
- Unfiltered tree walks use a plain traversal instead of `Bio.Phylo`'s
  `find_clades`, which runs its attribute matcher on every node even when no
  filter is given. Same nodes, same order; about 1.15x faster end to end.
- The DP cache fingerprint hashes topology, names and branch lengths directly
  instead of serialising the whole tree to a Newick string on every check.
- The outlier-aware DP inner loop now pairs children with contiguous slices
  instead of building index arrays and gathering through them once per *k*,
  matching what the no-outlier branch already did. Since `size_threshold`
  is set, this reduces temporary allocations. Selected *k*, cluster
  assignments and DP costs are unchanged.

#### Changes that alter output

These change results for code that is otherwise unmodified.

- `select_representative_species` now defaults to the central leaf rather than the
  most divergent one. Pass `strategy="divergent"` to restore the old behaviour.
- `soft_polytomy_max_degree` default lowered from 18 to 12. Polytomies of degree
  13 to 18 now require hard mode or an explicitly increased degree limit.
- Cluster palette reordered so that adjacent clusters stay visually distinct.
  Plots and the GUI now agree, but figures regenerated with this release differ in
  colour from earlier ones.
- `ClusterPlotConfig` defaults changed from `width_scale=2.5` / `height_scale=0.30`
  to `2.0` / `0.1`, so that all three plotting entry points share one set of
  defaults. Trees written by the CLI are correspondingly less tall; figures drawn
  from Python are unchanged, since `plot_clusters` already used these values.
  Set the fields explicitly to restore the old proportions.
- `plot_cluster` defaults brought into line with the same source of truth:
  `height_scale` `0.4` → `0.1` and `marker_size` `50` → `40`. This affects direct
  calls to `plot_cluster` only; `plot_clusters` has always passed both explicitly.

### Deprecated

- `SaveConfig.csv_name` still works but emits a `DeprecationWarning`. Scheduled
  for removal in 2.0.0.

### Removed

- The old `lambda_weight` config key and `--lambda-weight` CLI alias. Use
  `prominence_weight` and `--prominence-weight`.
- `CoreConfig`, which was unused.

### Fixed

- `ClusterPlotConfig` set via `runtime_config` is now honoured by `pc.plot()` and
  `plot_clusters()`. It was previously read only by the CLI, so cluster-plot
  settings passed from Python were silently ignored. Precedence is explicit
  argument → `ClusterPlotConfig` → built-in default, matching how
  `ScorePlotConfig` already behaved.
- `--no-outlier` is now respected alongside an outlier size threshold.
- Bootstrap stability: empty-`k` crash, tie-breaking, and score dilution.
- `RecursionError` when plotting very deep trees.
- GUI optimal-*k* axis toggle redrawing stale data.
- Output filename now defaults to `.tsv`, matching the tab-separated content.

### Security

- GUI: fixed a backslash-triggered infinite loop (denial of service) and a
  path-traversal file write.
- GUI: `/api/save` confines its `results_dir` to the server's working
  directory. The output *filename* was already guarded, but the directory was
  taken from the request unchecked, so a caller could create and write to any
  path the server user could reach. Restart the GUI in the target directory to
  write elsewhere.

## [0.1.2] – 2026-05-02

Reconstructed from `git log v0.1.2`; this release was tagged and uploaded without
a changelog entry.

### Added

- Web GUI (beta): FastAPI application with tree visualisation, a redesigned
  sidebar, and the Compare tab. Launched at this point via
  `uvicorn phytclust.gui.api:app`; the `phytclust gui` command came later, in
  1.0.0.
- Polytomy handling, including the convolution-based multifurcation DP, a worked
  example, and its test suite.
- Documentation site built with MkDocs.
- Visualisation examples.

### Changed

- CLI flag `--alpha` renamed to `--lambda-weight`. (Renamed again in 1.0.0, to
  `--prominence-weight`.)
- Faster DP and node lookup.
- Better parsing and error messages for `--max-k`.

### Fixed

- Repeated calls to `pc.run()` no longer fail; tree preparation moved into
  `__post_init__`.
- GUI API: NaN and infinite alpha scores are sanitised before JSON serialisation.
- Logo rendering bug.

## [0.1.1] – 2026-03-04

### Changed

- Plotting updates.

## [0.1.0] – 2026-03-04

- Initial release.

[1.0.0]: https://github.com/schwarzlab-ccb/PhytClust/releases/tag/v1.0.0
[0.1.2]: https://github.com/schwarzlab-ccb/PhytClust/releases/tag/v0.1.2
[0.1.1]: https://pypi.org/project/phytclust/0.1.1/
[0.1.0]: https://pypi.org/project/phytclust/0.1.0/
