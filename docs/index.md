# PhytClust

**Threshold-free clustering within phylogenetic trees**

PhytClust splits the leaves of a rooted phylogenetic tree into monophyletic
clusters. On binary trees, each cluster is a complete clade. Soft polytomy
mode also allows groups of children at an unresolved node to form a cluster.

The algorithm is dynamic programming directly on the tree topology. It finds the
partition of leaves into *k* groups that minimises within-cluster dispersion,
measured as summed leaf-to-MRCA distance. With the default objective, it finds
the lowest-cost valid partition for each *k*. Outlier preferences and penalties
can change that objective.

---

## Who it is for

Use PhytClust to group leaves in a rooted phylogenetic tree for analysis or
visualisation. Examples include:

- Grouping clones from tumour phylogenies into subclonal populations
- Partitioning gene or species trees into coherent groups for comparative analysis
- Summarising metagenomic diversity at several resolutions in one pass

## Three modes

Which mode applies depends on how much is known about the tree in advance.

- **Exact *k***: return a partition with the requested number of clusters.
- **Global peak search**: return up to `top_n` detected peaks, ordered by the
  configured ranking method.
- **Multi-resolution**: split the *k* range into logarithmic bins and select
  a peak from each bin that contains one.

### Output

The CLI saves assignment tables by default. Add `--save-fig` to save plots.

| File | Contents |
|------|----------|
| `phytclust_results.tsv` | Leaf-to-cluster assignments, for one or several *k* |
| `scores.png` | Score curve with annotated peaks, in automatic-selection modes |
| `tree_k{K}.png` | Coloured tree for each selected *k* |
| `peaks_by_rank.txt` | Selected *k* values in rank order |

---

## Where to start

In order:

1. **[Getting started](getting-started.md)**
2. **[How it works](concepts.md)**
3. **[CLI tutorial](tutorials/cli-end-to-end.md)**
4. **[Python tutorial](tutorials/python-end-to-end.md)**

The [reference](reference/index.md) covers the [CLI](reference/index.md#cli),
[configuration](reference/index.md#configuration), and
[Python API](reference/index.md#python-api) on one page.
