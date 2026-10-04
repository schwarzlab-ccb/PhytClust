# Python tutorial

The same ground as the [CLI tutorial](cli-end-to-end.md), covered through the
Python API.

---

## 1. Load a tree and construct a PhytClust object

PhytClust works on a [Biopython](https://biopython.org/) `Phylo` tree. The
constructor accepts a path to a Newick file, a Newick string, or an
already-parsed tree object.

```python
from Bio import Phylo
from phytclust import PhytClust
from phytclust.config import OutlierConfig

tree = Phylo.read("examples/sample_tree.nwk", "newick")

pc = PhytClust(
    tree,
    min_cluster_size=2,
    outlier=OutlierConfig(size_threshold=3, prefer_fewer=True),
    polytomy_mode="soft",
)
```

Outlier settings are passed as an `OutlierConfig`, as above. The constructor has
no `outlier_size_threshold=` or `prefer_fewer_outliers=` parameters: those names
exist only as the CLI flags `--outlier-size-threshold` and
`--prefer-fewer-outliers`, which the CLI routes into `OutlierConfig` for you.

## 2. Fixed-*k* clustering

With the number of clusters known in advance:

```python
result = pc.run(k=5)

print(result["mode"])        # "k"
print(result["selected_k"])  # 5
print(result["clusters"][0]) # Bio.Phylo leaf object -> cluster ID, for k = 5
```

`clusters` is always a list of leaf-to-cluster-ID maps, one per selected *k*. In
exact-*k* mode the list holds a single entry. `scores` is `None` in this mode.

## 3. Global peak search

When *k* is not known in advance, the global search scores every *k* up to
`max_k` and returns up to `top_n` peaks in rank order:

```python
result = pc.run(top_n=3, max_k=120)

print(result["mode"])           # "global"
print(result["k_values"])       # e.g. [7, 23, 45], ranked
print(len(result["clusters"]))  # up to 3, one map per selected k
```

`result["scores"]` holds the score vector for plotting or further analysis.
`result["k_values"]` is ordered by the configured ranking metric.

## 4. Multi-resolution mode

Select a peak from each logarithmic bin that contains one:

```python
result = pc.run(by_resolution=True, num_bins=4, max_k=120)

print(result["mode"])      # "resolution"
print(result["k_values"])  # up to four values, e.g. [3, 12, 38, 95]
```

## 5. Tuning peak selection

Peak detection is controlled by `PeakConfig`:

```python
from phytclust import PeakConfig

peak_cfg = PeakConfig(
    prominence_weight=0.7,    # 1 = rank by prominence, 0 = by score height
    ranking_mode="adjusted",  # "raw" | "adjusted"
    min_prominence=5.0,       # ignore small bumps in the score curve
    min_k=2,                  # ignore peaks below this k
    resolution_fallback_mode="max_score",   # "none" | "max_score": if a bin
                                            # holds no peak, take its max-score k
)

pc = PhytClust(tree, peak_config=peak_cfg)
result = pc.run(top_n=3, max_k=120)
print(result["k_values"])
```

`resolution_fallback_mode` applies only in resolution mode; it has no effect on
the global search shown here.

To rank detected peaks by balanced cluster sizes without changing the DP
objective or scores:

```python
peak_cfg = PeakConfig(partition_preference="balanced", partition_weight=1.0)
result = pc.run(top_n=3, peak_config=peak_cfg)
```

`partition_weight=1.0` uses partition quality alone. Lower weights blend it
with peak strength. Balanced quality rewards even cluster sizes and penalises
singleton cells, independently of the outlier size threshold.

## 6. Plot and output settings

Plot defaults live in `RuntimeConfig`, set once on the `PhytClust` object:

```python
from phytclust import (
    RuntimeConfig, PlotConfig,
    ScorePlotConfig, ClusterPlotConfig,
)

runtime_cfg = RuntimeConfig(
    plot=PlotConfig(
        scores=ScorePlotConfig(
            log_scale_y=False,
            fig_width=22,
            fig_height=11,
        ),
        cluster=ClusterPlotConfig(
            height_scale=0.25,
            marker_size=60,
        ),
    )
)

pc = PhytClust(tree, runtime_config=runtime_cfg)
result = pc.run(top_n=3, max_k=100)

pc.plot(results_dir="results/python", save=True)
pc.save(results_dir="results/python", filename="phytclust_results.tsv")
```

To retrieve a single partition without re-running, use `pc.get_clusters(k=7)`.
Both this method and `run()` return maps keyed by `Bio.Phylo` leaf objects.
Use `leaf.name` to retrieve each label.

## 7. Comparing several *k*

`plot_multiple_k` renders one tree image per *k*, whether the values came from
peak detection or from your own list:

```python
from phytclust.viz.cluster import plot_multiple_k

plot_multiple_k(
    pc,
    k_values=[3, 5, 8],
    results_dir="results/python_multi",
    save=True,
)
```

Comparing the images side by side shows how clusters subdivide as granularity
increases.

Full keyword arguments for `plot_multiple_k`, `plot_clusters` and `plot_cluster`
are documented on the [visualisation page](../visualisation.md), which is their
reference.

## 8. Batch processing

Running the same analysis over a directory of trees and collecting the results:

```python
from pathlib import Path
from Bio import Phylo
from phytclust import PhytClust, PeakConfig
from phytclust.config import OutlierConfig
import pandas as pd

peak_cfg = PeakConfig(prominence_weight=0.7, min_prominence=5.0)

rows = []
for tree_file in Path("data/trees/").glob("*.nwk"):
    tree = Phylo.read(str(tree_file), "newick")
    pc = PhytClust(
        tree,
        min_cluster_size=2,
        outlier=OutlierConfig(size_threshold=3),
        peak_config=peak_cfg,
    )
    result = pc.run(top_n=1, max_k=100)

    rows.append({
        "tree": tree_file.stem,
        "best_k": result["selected_k"],
        "n_leaves": len(list(tree.get_terminals())),
    })

    # optional: write per-tree outputs
    pc.save(results_dir=f"results/{tree_file.stem}")

summary = pd.DataFrame(rows)
print(summary)
```

## 9. Working with results

Each entry of `clusters` corresponds to the *k* at the same position in
`k_values`:

```python
import pandas as pd

result = pc.run(top_n=3, max_k=100)

for k, clusters in zip(result["k_values"], result["clusters"]):
    print(f"k={k}: {len(set(clusters.values()))} clusters")

    df = pd.DataFrame([
        {"leaf": leaf.name, "cluster": cid}
        for leaf, cid in clusters.items()
    ])

    threshold = pc.outlier.size_threshold
    if threshold is None:
        threshold = 2  # singleton marking used by the default export
    sizes = df.groupby("cluster")["leaf"].transform("size")
    outlier_leaves = df[sizes < threshold]
    print(f"  {len(outlier_leaves)} leaves in clusters below size {threshold}")
```

Returned partitions retain their cluster IDs. `pc.save(outlier=True)` marks
small clusters as `-1` in the exported table without changing those maps.
