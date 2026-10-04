# PhytClust Config Tutorial

This page walks you through configuring PhytClust in a practical, step-by-step way.

If you are new to configuration, start with the quick path below and then copy one of the ready-made templates.

---

## What You Will Learn

By the end of this tutorial, you will know how to:

1. Tune peak selection for `run(top_n=...)` and `run(by_resolution=True, ...)`
2. Customize score and cluster plotting defaults
3. Set save defaults for output tables
4. Use the same settings from Python or a YAML config file for CLI runs

---

## 1. Quick Start (Python API)

All config classes are importable from `phytclust`:

```python
from phytclust import (
    PhytClust,
    PeakConfig,
    RuntimeConfig,
    PlotConfig,
    ScorePlotConfig,
    ClusterPlotConfig,
)

# Build your peak-selection behavior
peak_cfg = PeakConfig(
    prominence_weight=0.5,
    min_prominence=5.0,
    resolution_fallback_mode="max_score",
)

# Build your plotting behavior
runtime_cfg = RuntimeConfig(
    plot=PlotConfig(
        scores=ScorePlotConfig(log_scale_y=False, fig_width=24),
        cluster=ClusterPlotConfig(height_scale=0.2, marker_size=60),
    )
)

pc = PhytClust(tree, runtime_config=runtime_cfg)
result = pc.run(top_n=3, peak_config=peak_cfg)
```

Use this pattern when you want reproducible runs from notebooks or scripts.

---

## 2. Tune Peak Selection

`PeakConfig` controls how candidate k values are detected and ranked.

```python
from phytclust import PeakConfig

peak_cfg = PeakConfig(
    prominence_weight=0.7,
    ranking_mode="adjusted",
    min_prominence=None,
    min_k=2,
)
```

### Most Important Parameters

| Field                      | Default      | When to change it                                                              |
| -------------------------- | ------------ | ------------------------------------------------------------------------------ |
| `prominence_weight`            | `0.7`        | Lower toward `0` to rank peaks by score height; raise toward `1` to rank by prominence. |
| `ranking_mode`             | `"adjusted"` | `"adjusted"` normalises then blends the two signals; `"raw"` ranks by absolute prominence alone. |
| `min_prominence`           | `None`       | Increase to suppress minor local peaks.                                        |
| `min_k`                    | `2`          | Increase if you want to ignore very small k values.                            |
| `resolution_fallback_mode` | `"none"`     | One of `"none"` \| `"max_score"`. Set to `"max_score"` so each resolution bin always returns a k. |

### Boundary Controls for k=2

These options specifically gate whether the boundary candidate at `k=2` is accepted:

| Field                      | Default | Meaning                                           |
| -------------------------- | ------- | ------------------------------------------------- |
| `boundary_window_size`     | `5`     | Right-side window size used for comparison.       |
| `boundary_ratio_threshold` | `1.5`   | Minimum ratio vs right-window median to keep k=2. |

---

## 3. Customize Plots and Saved Outputs

Use `RuntimeConfig` to keep plotting and saving settings in one place.

```python
from phytclust import RuntimeConfig, PlotConfig, ScorePlotConfig, ClusterPlotConfig

runtime_cfg = RuntimeConfig(
    plot=PlotConfig(
        scores=ScorePlotConfig(
            fig_width=18,
            fig_height=10,
            log_scale_y=True,
            x_axis_mode="log",
        ),
        cluster=ClusterPlotConfig(
            cmap="tab20",
            width_scale=2.5,
            height_scale=0.30,
            marker_size=40,
        ),
    )
)
```

### Score Plot Tips

- Set `log_scale_y=False` if small peaks are getting visually compressed.
- Set `x_axis_mode="linear"` for easier interpretation on small k ranges.
- Keep `prefer_unsmoothed_primary=True` if you want raw k>=3 behavior as the main exported curve.

### Cluster Plot Tips

- Increase `height_scale` for crowded trees with many leaves.
- Increase `marker_size` if terminal nodes are hard to see in saved figures.
- Enable `show_branch_lengths=True` when branch-length interpretation is important.

### Save Defaults

`SaveConfig` controls defaults used by `pc.save()`:

| Field      | Default                   | Meaning                                        |
| ---------- | ------------------------- | ---------------------------------------------- |
| `tsv_name` | `"phytclust_results.tsv"` | Default output filename (output is tab-separated). |
| `outlier`  | `True`                    | Mark outlier clusters as `-1` in saved output. |

---

## 4. Use the Same Settings in CLI with YAML

You can define equivalent settings in a config file and pass it via `--config`.

```yaml
# phytclust.config.yaml
peak:
  prominence_weight: 0.5
  min_prominence: 5.0
  resolution_fallback_mode: max_score

plot:
  scores:
    log_scale_y: false
    fig_width: 24
  cluster:
    height_scale: 0.2
    cmap: tab10

save:
  tsv_name: results.tsv
```

```bash
phytclust tree.nwk --top-n 3 --config phytclust.config.yaml --out-dir out
```

This is the easiest way to share exact run settings with collaborators.

---

## 5. Ready-to-Use Configuration Recipes

### Recipe A: Conservative Peak Calls

Use when you want fewer, stronger peaks.

```yaml
peak:
  ranking_mode: adjusted
  prominence_weight: 0.8
  min_prominence: 10.0
  min_k: 3
```

### Recipe B: Explore More Candidate Peaks

Use when exploring noisy trees.

```yaml
peak:
  ranking_mode: raw
  prominence_weight: 0.2
  min_prominence: 1.0
  min_k: 2
  resolution_fallback_mode: max_score
```

### Recipe C: Publication-Ready Plot Scaling

Use when exported figures need larger typography.

```yaml
plot:
  scores:
    fig_width: 24
    fig_height: 12
    title_fontsize: 56
    axis_label_fontsize: 38
    tick_labelsize: 32
  cluster:
    marker_size: 60
    height_scale: 0.35
```

---

## 6. Config Object Nesting (Mental Model)

```text
RuntimeConfig
  plot: PlotConfig
    cluster: ClusterPlotConfig
    scores: ScorePlotConfig
  save: SaveConfig
```

Think of `PeakConfig` as run-time scoring logic, and `RuntimeConfig` as plotting/saving defaults.

---

## 7. Full Reference Tables

Use this section when you need exact field names and defaults.

### PeakConfig

| Field                      | Type            | Default      | Description                                                                                        |
| -------------------------- | --------------- | ------------ | -------------------------------------------------------------------------------------------------- |
| `prominence_weight`            | `float`         | `0.7`        | Blend of the two ranking signals, `0`–`1`: `1` = peak prominence only, `0` = absolute score only. Ignored when `ranking_mode="raw"`. |
| `ranking_mode`             | `str`           | `"adjusted"` | One of `"raw"` \| `"adjusted"`. `"adjusted"` min–max normalises prominence and score across peaks, then blends by `prominence_weight`; `"raw"` ranks by absolute prominence alone. |
| `boundary_window_size`     | `int`           | `5`          | Right-window size for evaluating the `k=2` boundary candidate.                                     |
| `boundary_ratio_threshold` | `float`         | `1.5`        | Minimum score ratio vs right-window median for `k=2` to pass.                                      |
| `min_prominence`           | `float \| None` | `None`       | Minimum prominence for `scipy.signal.find_peaks`; `None` auto-sets to 1% of score range.           |
| `min_k`                    | `int`           | `2`          | Smallest cluster count considered for automatic peak selection.                                                                    |
| `resolution_fallback_mode` | `str`           | `"none"`     | One of `"none"` \| `"max_score"`. Resolution mode only: `"none"` leaves a peak-less bin empty; `"max_score"` picks the max-score k in that bin. Any other value raises `ConfigurationError`. |
| `exclude_k2`               | `bool`          | `False`      | Drop `k=2` from automatic peak selection. By default `k=2` is a candidate, judged by the boundary test. Exact `k=2` requests are always honoured. |
| `use_log_peak_input` | `bool` | `False` | Detect peaks in the logarithm of nonnegative scores plus an offset. |
| `log_peak_offset` | `float` | `1e-12` | Positive offset before taking logarithms; used only with log peak input. |
| `use_relative_prominence` | `bool` | `False` | Rank prominence as the original peak score divided by its baseline score. |
| `min_relative_prominence` | `float \| None` | `None` | Minimum peak-to-baseline ratio; used only with relative prominence. |
| `prominence_k_power` | `float` | `0.0` | Multiply an explicit minimum prominence by cluster count raised to this power; used only without relative prominence. |

Settings are checked again before peak detection, so changes made after construction are also validated. Weights and thresholds for disabled options are ignored. Relative prominence uses original scores even when peak detection uses logarithms.

### OutlierConfig

| Field | Type | Default | Description |
| --- | --- | --- | --- |
| `size_threshold` | `int \| None` | `None` | Clusters with fewer leaves than this count are outliers. None counts singletons when prefer_fewer=True; otherwise it disables outlier counting in the DP. |
| `prefer_fewer` | `bool` | `False` | Prioritize outlier count before partition cost. With no threshold, count singleton clusters. |
| `penalty_enabled` | `bool` | `False` | Add a small-cluster penalty when a size threshold is set. Supported for binary trees. |
| `ratio_weight` | `float` | `10.0` | Multiply the small-cluster penalty by this nonnegative weight. Zero disables the penalty. |
| `ratio_mode` | `str` | `"exp"` | Penalty shape: exp(size shortfall) − 1, inverse cluster size, or size shortfall for the legacy "power" mode. |

Penalty weights and modes are checked only when the penalty is enabled and a size threshold is set. With zero weight, the mode is ignored. DP validation repeats these checks after settings change. When `prefer_fewer=False`, cost takes priority and outlier count breaks cost ties.

Saving with outlier labels enabled uses the configured size threshold, or a threshold of two when none is set. This labels singleton clusters as −1 by default, even when DP outlier counting is disabled.

### ScorePlotConfig

| Field                       | Type            | Default    | Description                                  |
| --------------------------- | --------------- | ---------- | -------------------------------------------- |
| `title_fontsize`            | `int`           | `20`       | Title font size.                             |
| `title_fontfamily` | `str` | `"Liberation Sans"` | Score title font; falls back to DejaVu Sans. |
| `axis_label_fontsize`       | `int`           | `18`       | Axis label font size.                        |
| `tick_labelsize`            | `int`           | `14`       | Tick label font size.                        |
| `peak_labelsize`            | `int`           | `16`       | Peak annotation font size.                   |
| `peak_marker`               | `str`           | `"o"`      | Matplotlib marker for peaks.                 |
| `peak_markersize`           | `int`           | `7`       | Peak marker size; also sets label clearance. |
| `bin_labelsize`             | `int`           | `14`       | Bin label font size (resolution mode).       |
| `fig_width`                 | `float`           | `9`       | Figure width in inches.                      |
| `fig_height`                | `float`           | `5`       | Figure height in inches.                     |
| `clamp_negative_to_zero`    | `bool`          | `True`     | Show negative scores as zero.          |
| `log_scale_y`               | `bool`          | `True`     | Use log scale on y-axis.                     |
| `x_axis_mode`               | `str`           | `"log"`    | `"log"` or `"linear"` x-axis for k.          |
| `log_base`                  | `float \| None` | `None`     | Log base for axes; `None` means base ten. |
| `prefer_unsmoothed_primary` | `bool`          | `True`     | Save k>=3 unsmoothed score curve as primary. |
| `show_secondary_score_plot` | `bool`          | `False`    | Generate a secondary companion figure.       |
| `colorblind_palette`        | `list[str] \| None`     | `None` | Resolution band colours; None uses the built-in palette.                |

### ClusterPlotConfig

| Field                 | Type    | Default   | Description                             |
| --------------------- | ------- | --------- | --------------------------------------- |
| `cmap`                | `str`   | `"phytclust"` | Cluster palette name. |
| `width_scale`         | `float` | `2.0`     | Horizontal tree scaling factor.         |
| `height_scale`        | `float` | `0.1`    | Vertical scaling factor per leaf.       |
| `marker_size`         | `int`   | `40`      | Terminal node marker size.              |
| `show_branch_lengths` | `bool`  | `False`   | Annotate edges with branch lengths.     |
| `hide_internal_nodes` | `bool`  | `True`    | Hide internal node markers.             |

### SaveConfig

| Field      | Type   | Default                   | Description                                    |
| ---------- | ------ | ------------------------- | ---------------------------------------------- |
| `tsv_name` | `str`  | `"phytclust_results.tsv"` | Default output filename (output is tab-separated). |
| `outlier`  | `bool` | `True`                    | Annotate outlier clusters with `-1` in output. |

`build_runtime_config` accepts nested mappings matching `RuntimeConfig`. Unknown settings produce a warning with the full setting name and are ignored. The `plot`, `plot.cluster`, `plot.scores`, and `save` sections must be mappings.

### Rank peaks using cluster sizes

`PeakConfig(partition_preference="balanced", partition_weight=0.5)` reorders
detected peaks without changing the DP partitions, score curve, or candidate
peaks. Keep `OutlierConfig(prefer_fewer=False)` to retain cost-first partitions.
The DP's `prefer_fewer` option remains separate: it changes the optimization
objective and can change the score curve.

| Field | Default | Behavior |
| --- | --- | --- |
| `partition_preference` | `"none"` | `none` retains existing ranking; `fewer_outliers` favors fewer singleton cells; `balanced` also favors more even cluster sizes. |
| `partition_weight` | `0.5` | Weight from 0 to 1. Zero retains existing ranking; one ranks by partition quality alone. Ignored for `none`. |

For each candidate, singleton avoidance is `1 - singleton_count / cell_count`.
Size balance is `cell_count**2 / (cluster_count * sum(size**2))`: one for equal
cluster sizes and lower for uneven sizes. Balanced quality multiplies the two
values, so a partition consisting entirely of singletons has zero quality.
The new rank is `(1 - partition_weight) * peak_strength + partition_weight * quality`.
Peak strength is the existing ranking metric divided by its maximum among
candidates; it is one for all candidates when all existing metrics are zero.

`peak_ranking_details` retains the original score, prominence, and combined
metric and includes cluster sizes, singleton counts, size balance, partition
quality, and the new ranking metric when the preference is active.
Global and resolution modes both use the new order. It does not add peaks at
requested cluster counts. Backtracking is only needed for detected candidates
when the option is active and its weight is nonzero.

CLI equivalents are `--peak-partition-preference balanced` and
`--peak-partition-weight 0.5`.
