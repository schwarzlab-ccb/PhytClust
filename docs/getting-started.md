# Getting started

## Requirements

- Python 3.10 or later
- A rooted tree in **Newick** format

## Install

```bash
pip install phytclust
```

From source:

```bash
git clone https://github.com/schwarzlab-ccb/PhytClust.git
cd PhytClust
pip install -e ".[dev]"
```

Confirm the install with `phytclust --version`.

## First run

Cluster the bundled sample tree into five groups:

```bash
phytclust examples/sample_tree.nwk --k 5 --save-fig --out-dir results/quickstart
```

`results/quickstart/` then contains:

- `phytclust_results.tsv`: which leaf is in which cluster
- `tree_k5.png`: the tree coloured by cluster, one colour per clade

Automatic selection also saves `scores.png` when `--save-fig` is enabled.

## When *k* is unknown

`--top-n 3` returns up to three ranked peaks. `--resolution` selects a peak
from each bin that contains one. The [CLI tutorial](tutorials/cli-end-to-end.md) covers both.

## Rooting

PhytClust requires a rooted tree. An unrooted tree can be rooted at run time:

```bash
phytclust tree.nwk --k 5 --root-taxon "species_A"   # on a named taxon
phytclust tree.nwk --k 5 --root-taxon midpoint      # midpoint rooting
```

To drop a taxon before clustering rather than root on it, use `--outgroup`.

## Web GUI (experimental)

```bash
pip install "phytclust[gui]"
phytclust gui
```

This opens <http://127.0.0.1:8000>, where a Newick string can be pasted and
explored interactively. Use the CLI or the Python API for reproducible analysis.

## Next

- [How it works](concepts.md)
- [CLI](tutorials/cli-end-to-end.md) / [Python](tutorials/python-end-to-end.md)
- [Reference](reference/index.md)
