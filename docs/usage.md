# Usage

## Install

Requires Python 3.10 or later.

```bash
pip install phytclust
```

Replace `tree.nwk` with the path to your rooted Newick tree. Python examples reuse
`pc` from the first section.

## Choose the number of clusters

=== "CLI"

    ```bash
    phytclust tree.nwk --k 5 --save-fig --out-dir results/k5
    ```

=== "Python"

    ```python
    from phytclust import PhytClust

    pc = PhytClust("tree.nwk")
    result = pc.run(k=5)
    clusters = result["clusters"][0]
    ```

The CLI saves `phytclust_results.tsv` and `tree_k5.png` under `results/k5/`.
Singleton clusters are marked `-1` in the table by default. In Python,
`clusters` maps leaf objects to cluster IDs; use `leaf.name` for their labels.

## Find candidate partitions

If *k* is unknown, return up to three ranked peaks from the score curve:

=== "CLI"

    ```bash
    phytclust tree.nwk --top-n 3 --max-k 120 --save-fig
    ```

=== "Python"

    ```python
    result = pc.run(top_n=3, max_k=120)
    print(result["k_values"])
    ```

Set `max_k` to a suitable bound for your tree. Without it, the search extends
to 90% of the leaf count, rounded up.

To compare coarse and fine partitions, use resolution mode. It selects one
peak per logarithmic bin:

=== "CLI"

    ```bash
    phytclust tree.nwk --resolution --bins 4 --max-k 120 --save-fig
    ```

=== "Python"

    ```python
    result = pc.run(by_resolution=True, num_bins=4, max_k=120)
    ```

Bins without peaks remain empty. With fewer than 50 scored *k* values,
resolution mode falls back to a global search.

Automatic selection with `--save-fig` also saves the score curve as
`scores.png`.

## Save results from Python

```python
pc.plot(results_dir="results/python", save=True, top_n=3)
pc.save(results_dir="results/python", top_n=3)
```

Both methods use only the first selected partition unless `top_n` is given.

## Common options

| To… | CLI option | Python constructor option |
|-----|------------|---------------------------|
| Root at the midpoint | `--root-taxon midpoint` | `root_taxon="midpoint"` |
| Require at least three leaves per cluster | `--min-cluster-size 3` | `min_cluster_size=3` |
| Mark clusters smaller than three as outliers | `--outlier-size-threshold 3` | `outlier=OutlierConfig(size_threshold=3)` |
| Preserve the clades at polytomies strictly | `--polytomy-mode hard` | `polytomy_mode="hard"` |
| Keep all-zero subtrees together | `--no-split-zero-length` | `no_split_zero_length=True` |

Import `OutlierConfig` from `phytclust` when using it in Python. A minimum
cluster size can make some *k* values infeasible. Outlier marking only changes
exported IDs; it keeps the small clusters in the partition.

## Web interface

```bash
pip install "phytclust[gui]"
phytclust gui
```

Open <http://127.0.0.1:8000> and upload your Newick file.

See the [reference](reference.md) for configuration and plotting options, or
[How PhytClust works](concepts.md) for the objective and score.
