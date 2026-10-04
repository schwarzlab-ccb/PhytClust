# Reference

CLI, config schema, and Python API for `phytclust`. For what the score is
and how *k* is selected, see [concepts](../concepts.md). Walkthroughs:
[CLI](../tutorials/cli-end-to-end.md) · [Python](../tutorials/python-end-to-end.md).

---

## CLI

```bash
phytclust TREE [options]     # cluster a tree (TREE = Newick file, or - for stdin)
phytclust gui [options]      # launch the web GUI
```

### Common flags

| Flag | Default | Meaning |
|------|---------|---------|
| `-k, --k` |: | Exactly this many clusters (≥ 1) |
| `--top-n` | `1` | Number of peaks taken from the global score curve |
| `--resolution` | off | One peak per log bin |
| `--max-k` |: | Upper bound on *k*; defaults to `ceil(max_k_limit · n_leaves)` |
| `--outgroup` |: | Outgroup to drop before clustering |
| `--root-taxon` |: | Re-root on this taxon, or `midpoint` |
| `-o, --out-dir` | `results/` | Output directory |
| `--save-fig` | off | Write score and tree PNGs |
| `--tsv-name` | `phytclust_results.tsv` | Output table |
| `--config` |: | YAML/JSON config file |
| `--plot` | off | Open interactive plot windows (blocks until they are closed) |

### Other flags

| Flag | Default | Meaning |
|------|---------|---------|
| `--bins` | `3` | Bins for `--resolution` |
| `--max-k-limit` | `0.9` | *k* ceiling as a fraction of `n_leaves` |
| `--prominence-weight` | `0.7` | Peak-rank blend: 1 = rank by prominence, 0 = by score height |
| `--peak-partition-preference` | `none` | Rank detected peaks using `none`, `fewer_outliers`, or `balanced` |
| `--peak-partition-weight` | `0.5` | Blend peak strength with partition quality; `1` uses partition quality alone |
| `--exclude-k2` | off | Drop *k* = 2 from automatic selection (`--include-k2` is accepted but deprecated; *k* = 2 is a candidate by default) |
| `--min-cluster-size` | `1` | Hard minimum cluster size |
| `--outlier-size-threshold` |: | Clusters below this size are marked `-1` |
| `--prefer-fewer-outliers` | off | Make the DP minimise outlier count before cost; counts singletons when no threshold is supplied |
| `--polytomy-mode` | `soft` | `soft` = any group of a multifurcation's children may form a cluster; `hard` = each child goes wholly into one cluster |
| `--soft-polytomy-max-degree` | `12` | Soft mode raises at nodes above this degree; use `hard` for them |
| `--no-split-zero-length` | off | Forbid splitting zero-length edges |
| `--save-tree` `--save-all-k` `--no-tsv` `--dpi` |: | Tree-only PNGs · every *k* · skip TSV · PNG dpi (`150`) |
| `-v` `-q` `--time` `--progress` `--no-color` `--version` |: | Verbosity · timing · spinner · colour · version |

CLI flags override config-file values.

**`phytclust gui`** requires the `[gui]` extra. Flags: `--host` (`127.0.0.1`),
`--port` (`8000`), `--reload`, `--no-browser`.

```bash
phytclust tree.nwk --k 5 --save-fig            # exact k
phytclust tree.nwk --top-n 3 --max-k 120       # global peaks
phytclust tree.nwk --resolution --bins 4       # resolution mode
cat tree.nwk | phytclust - --k 5               # read from stdin
```

---

## Configuration

Two groups: **`PeakConfig`** (peak detection) and **`RuntimeConfig`**
(plots and output). Set either in Python or in a YAML file passed with `--config`.

```yaml
peak:
  prominence_weight: 0.7
  min_prominence: 5.0
runtime:
  plot:
    scores:  { log_scale_y: false, fig_width: 22 }
    cluster: { height_scale: 0.25, marker_size: 55 }
  save:      { tsv_name: phytclust_results.tsv, outlier: true }
```

**`PeakConfig`**

| Field | Default | Meaning |
|-------|---------|---------|
| `prominence_weight` | `0.7` | Blend when `ranking_mode="adjusted"`: `1` ranks by prominence only, `0` by score height only |
| `ranking_mode` | `"adjusted"` | `"adjusted"` min–max normalises prominence and score, then blends them by `prominence_weight`; `"raw"` ranks by absolute prominence alone and ignores it |
| `partition_preference` | `"none"` | Keep existing ranking, prefer fewer singleton cells (`"fewer_outliers"`), or also prefer even cluster sizes (`"balanced"`) |
| `partition_weight` | `0.5` | Weight from 0 to 1; `0` keeps existing ranking, `1` uses partition quality alone. Ignored for `"none"` |
| `min_prominence` | `None` | Minimum peak prominence; `None` derives it from the score range |
| `exclude_k2` | `False` | Drop *k* = 2 from automatic selection ([details](../concepts.md#k-2-and-the-boundary-test)). An explicit `k=2` is always honoured |
| `min_k` | `2` | Ignore peaks below this *k* |
| `use_relative_prominence` | `False` | Rank by fold-change rather than absolute prominence |

Rarely changed: `boundary_window_size` (`5`), `boundary_ratio_threshold` (`1.5`),
`use_log_peak_input` (`False`), `prominence_k_power` (`0.0`).

`resolution_fallback_mode` (`"none"`) applies to `--resolution` mode only. It
decides what happens to a bin in which no peak was detected:

| Value | Effect |
|-------|--------|
| `"none"` *(default)* | Leave the bin empty, so fewer *k* values come back than there are bins |
| `"max_score"` | Fall back to the highest-scoring *k* in that bin, so every bin returns a *k* |

Any other value raises `ConfigurationError`. `ranking_mode` likewise accepts only
`"raw"` or `"adjusted"`. Both sets are importable as
`phytclust.config.RESOLUTION_FALLBACK_MODES` and `RANKING_MODES`.

**`RuntimeConfig`** holds `plot.cluster`, `plot.scores`, and `save`.

`ClusterPlotConfig`

| Field | Default | Meaning |
|-------|---------|---------|
| `cmap` | `"phytclust"` | Colormap (project palette) |
| `width_scale` / `height_scale` | `2.0` / `0.1` | Horizontal / vertical stretch |
| `marker_size` | `40` | Leaf marker size |
| `show_branch_lengths` | `False` | Label branch lengths |
| `hide_internal_nodes` | `True` | Hide internal labels |

`ScorePlotConfig`: `fig_width` / `fig_height` (`9` / `5`), `log_scale_y` (`True`),
`x_axis_mode` (`"log"`), `clamp_negative_to_zero` (`True`), `colorblind_palette`,
and font sizes (title `20`, axis `18`, tick `14`, peak `16`, bin `14`).
The default title font is Liberation Sans.

`SaveConfig`: `tsv_name` (`"phytclust_results.tsv"`), `outlier` (`True` marks
outlier clusters `-1`).

---

## Python API

```python
from phytclust import PhytClust
from phytclust.config import OutlierConfig

pc = PhytClust(
    tree,                          # Bio.Phylo tree, Newick string, or file path
    outgroup=None,                 # taxon to exclude
    root_taxon=None,               # taxon to root on, or "midpoint"
    min_cluster_size=1,
    outlier=OutlierConfig(),       # size threshold / DP penalty / tie-breaking
    polytomy_mode="soft",          # "soft" | "hard"
    soft_polytomy_max_degree=12,
    no_split_zero_length=False,
    dp_float32=False,              # halve DP memory on very large trees
    runtime_config=None,           # RuntimeConfig; defaults applied if None
    peak_config=None,              # PeakConfig; defaults applied if None
)
```

### Large trees

The DP uses float64 by default. `dp_float32=True` halves the storage used by
floating-point DP tables, but reduces precision and may change selected peaks.
Use it when memory is limiting, and compare with float64 on a smaller tree
when possible.

Two other settings affect memory use: `max_k` (or `max_k_limit`) bounds the DP
directly, and `preserve_dp_tables=False` (the default) frees each child's table
as soon as its parent is computed.

Outlier settings go through `OutlierConfig`. With `prefer_fewer=True` and no
threshold, singleton clusters count as outliers. For example,
`OutlierConfig(size_threshold=3, prefer_fewer=True)` counts clusters of one or
two leaves instead.

Peak partition preferences only reorder detected candidates; they leave the DP
objective and score curve unchanged. Their singleton counts do not depend on
`OutlierConfig.size_threshold`.

### `run()` → dict

The mode is inferred from the arguments:

- `k` selects exact-*k*; `by_resolution=True` selects resolution mode; otherwise
  the global search runs and returns `top_n` peaks.
- `k` cannot be combined with `by_resolution` or with `top_n ≠ 1`.

```python
pc.run(k=5)                              # exact k
pc.run(top_n=3, max_k=120)               # global; top-n peaks
pc.run(by_resolution=True, num_bins=4)   # peak from each bin containing one
```

`run()` also accepts `peak_config=`, which overrides the `PeakConfig` given to
the constructor for that call only. Both forms are supported.

Returns `{mode, k_values, selected_k, clusters, scores}`. `clusters` is always a
list of leaf-to-cluster-ID maps, one entry per selected *k*, including in
exact-*k* mode where the list has length one. `scores` is `None` in exact-*k* mode.

### Other methods

```python
pc.get_clusters(k=7)   # {Clade: id}: keyed by Bio.Phylo leaf objects; use clade.name
pc.plot(results_dir="results/", save=True, dpi=150)
pc.save(results_dir="results/", filename="phytclust_results.tsv")
```

Both `run()` and `get_clusters()` return cluster maps keyed by `Bio.Phylo`
leaf objects. Use `leaf.name` to retrieve each label. Returned maps retain
cluster IDs; saving with `outlier=True` marks small clusters as `-1` in the
exported table. With no threshold, export marking applies to singletons.

### Plotting functions

`pc.plot()` is the convenience wrapper. The `phytclust.viz` functions take the
full keyword set, documented on the
[visualisation page](../visualisation.md), which is their reference.

```python
from phytclust.viz.cluster import plot_clusters, plot_multiple_k
from phytclust.viz.plots import plot_cluster

plot_clusters(pc, k=3, save=True, results_dir="figures")
plot_multiple_k(pc, k_values=[3, 5, 8], results_dir="figures", save=True)
```

`plot_clusters` backtracks the requested *k* from the DP table, so it can plot
any *k* up to `max_k`, not only the ones the run selected.

Options left unset resolve from the object's `ClusterPlotConfig`, so the CLI,
`pc.plot()` and `plot_clusters()` agree. Precedence is explicit argument →
`ClusterPlotConfig` → built-in default.

### Imports

```python
from phytclust import PeakConfig, RuntimeConfig, PlotConfig, ScorePlotConfig, ClusterPlotConfig
from phytclust.exceptions import (
    PhytClustError, ValidationError, InvalidTreeError, ConfigurationError,
    InvalidKError, InvalidClusteringError, MissingDPTableError, ComputationError, DataError,
)
from phytclust.selection.representatives import select_representative_species  # strategy="central"|"divergent"|"medoid"
from phytclust.metrics.indices import colless_index_calc, normalized_colless
```
