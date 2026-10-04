# CLI tutorial

This walks one tree through the main modes of the `phytclust` command line tool:
exact *k*, global peak search, resolution mode, size constraints, polytomy
handling, and config files. Later sections assume the earlier ones.

Examples use the trees bundled under `examples/`. Substitute your own Newick file
anywhere `sample_tree.nwk` appears.

---

## 1. Exact-*k* clustering

Pass `--k` when the number of groups is fixed in advance — by prior biology, or
by a downstream tool that expects a set number of clusters.

```bash
phytclust examples/sample_tree.nwk --k 5 --save-fig --out-dir results/cli_baseline
```

The output directory then contains:

- **`phytclust_results.tsv`** — one row per leaf with its cluster assignment
- **`tree_k5.png`** — the tree, coloured by cluster
- **`scores.png`** — the score curve across all *k*

## 2. Global peak search

When *k* is not known in advance, the global search scans the score curve over
the full *k* range and returns the highest-ranked peaks.

```bash
phytclust examples/sample_tree.nwk \
  --top-n 3 \
  --max-k 120 \
  --save-fig \
  --out-dir results/cli_global
```

Three files carry the result:

- **`peaks_by_rank.txt`** — the selected *k* values ordered by prominence
- **`scores.png`** — the curve the peaks were taken from. Widely separated peaks
  indicate distinct structural scales in the tree; adjacent peaks usually mean
  one scale with an uncertain boundary
- **`phytclust_results.tsv`** — assignments for every selected *k*, side by side

Without `--max-k`, the scan runs up to 90% of the leaf count. On large trees,
setting it explicitly to a plausible range is faster.

## 3. Resolution mode: one *k* per scale

Resolution mode splits the *k* range into logarithmic bins and takes the best *k*
in each, rather than returning neighbouring peaks from one region of the curve.

```bash
phytclust examples/sample_tree.nwk \
  --resolution \
  --bins 4 \
  --save-fig \
  --out-dir results/cli_resolution
```

This returns four solutions spread from broad to fine-grained.

## 4. Constraining cluster size

### Hard minimum

To exclude clusters smaller than three leaves:

```bash
phytclust examples/sample_tree.nwk \
  --top-n 3 \
  --min-cluster-size 3 \
  --out-dir results/cli_min_size
```

The DP enforces this during optimisation rather than filtering afterwards. Any
*k* that cannot be reached without violating the constraint does not appear in
the results at all.

### Outlier marking

The alternative is to let small clusters exist and mark them, optionally biasing
the ranking toward solutions that produce fewer of them:

```bash
phytclust examples/sample_tree.nwk \
  --top-n 3 \
  --outlier-size-threshold 3 \
  --prefer-fewer-outliers \
  --out-dir results/cli_outlier
```

Outlier clusters are written as `-1` in the output TSV.

## 5. Polytomies and zero-length edges

### Polytomies

Internal nodes with more than two children are handled natively by a DP that
reasons over all children of the node at once:

```bash
phytclust examples/sample_polytomy.newick --k 5
```

Two polytomy modes are available:

- **`soft`** (default) — the multifurcation is treated as unresolved
  zero-length branching, so any group of two or more children may form a cluster
- **`hard`** — each child of a polytomy goes entirely into one cluster, unless
  all of them form one cluster together

```bash
phytclust tree.nwk --k 5 --polytomy-mode hard
```

Soft mode is exponential in node degree. `--soft-polytomy-max-degree`
(default 12) is a guardrail: a polytomy above that degree stops the run with an
error. Use `--polytomy-mode hard` for such trees, or raise the limit
deliberately.

### Zero-length edges

Internal edges of length zero usually come from collapsed uncertain nodes. By
default the DP is free to split at them. To keep zero-length clades together:

```bash
phytclust tree.nwk --k 5 --no-split-zero-length
```

This matters on trees with many collapsed or near-identical leaves, where a split
at a zero-length branch is arbitrary.

## 6. Config files

Settled settings go in a YAML file so runs are reproducible:

```yaml
# my_config.yaml
peak:
  prominence_weight: 0.6
  ranking_mode: adjusted          # "raw" | "adjusted"
  min_prominence: 5.0
  resolution_fallback_mode: max_score   # "none" | "max_score"

runtime:
  plot:
    scores:
      log_scale_y: false
      fig_width: 20
    cluster:
      height_scale: 0.25
      marker_size: 55
  save:
    tsv_name: phytclust_results.tsv
```

Pass it with `--config`:

```bash
phytclust examples/sample_tree.nwk \
  --top-n 3 \
  --max-k 120 \
  --config my_config.yaml \
  --save-fig \
  --out-dir results/cli_config
```

CLI flags override config-file values, so the file can serve as a baseline that
individual runs adjust.

The curve referred to in sections 2 and 3, and plotted in `scores.png`, is the
product of a Calinski-Harabasz-style validity score and an elbow score, evaluated
at every *k*. Both terms, and the reasons the criterion departs from textbook CH,
are defined in
[the score curve and choosing *k*](../concepts.md#the-score-curve-and-choosing-k).
