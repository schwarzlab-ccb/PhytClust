# PhytClust

**Threshold-free clustering within phylogenetic trees**

PhytClust splits the leaves of a rooted phylogenetic tree into monophyletic
clusters, every cluster is a complete clade.

The algorithm is dynamic programming directly on the tree topology. It finds the
partition of leaves into *k* groups that minimises within-cluster dispersion,
measured as summed leaf-to-MRCA distance. The DP explores the full space of valid
partitions for each *k*, so the result is exact rather than heuristic.

---

## Who it is for

PhytClust is for anybody that needs taxa/leaves grouped into interpretable units for
downstream analysis, visualisation, or summarising large trees. It take as an input any tree-like structure. It can be used in various fields. For example:

- Grouping clones from tumour phylogenies into subclonal populations
- Partitioning gene or species trees into coherent groups for comparative analysis
- Summarising metagenomic diversity at several resolutions in one pass

## Three modes

Which mode applies depends on how much is known about the tree in advance.

- **Exact *k*** — the number of clusters is fixed beforehand.
- **Global peak search** — every *k* in range is scored and the values where the
  tree has natural breakpoints are returned, ranked by prominence.
- **Multi-resolution** — the *k* range is split into logarithmic bins and one
  representative *k* is returned per bin, giving a coarse-to-fine panel from a
  single run.

### Output

| File | Contents |
|------|----------|
| `phytclust_results.tsv` | Leaf-to-cluster assignments, for one or several *k* |
| `scores.png` | Score curve with annotated peaks |
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
