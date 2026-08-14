# PhytClust

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19731229.svg)](https://doi.org/10.5281/zenodo.19731229)
[![PyPI](https://img.shields.io/pypi/v/phytclust.svg)](https://pypi.org/project/phytclust/)
[![License: GPL-3.0](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

Clustering for rooted phylogenetic trees.

PhytClust groups the tips of a rooted tree into monophyletic clusters using
dynamic programming that minimises the summed distance from each tip in a cluster
to that cluster's most recent common ancestor (MRCA). For a fixed number of
clusters the result is the exact global optimum, not a heuristic.

It runs as a command-line tool, a Python API, and an experimental web interface.

## Requirements

Python 3.10 or later, and a rooted tree in Newick format.

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

## Usage

```bash
# cluster a tree into exactly 5 clades
phytclust examples/sample_tree.nwk --k 5 --save-fig

# let PhytClust choose the 3 best values of k
phytclust examples/sample_tree.nwk --top-n 3 --save-fig

# one k per resolution scale
phytclust examples/sample_tree.nwk --resolution --bins 4 --save-fig
```

From Python:

```python
from phytclust import PhytClust

pc = PhytClust("examples/sample_tree.nwk")
result = pc.run(k=5)
```

## Web interface (experimental)

![PhytClust web GUI: sample tree and interface](docs/img/gui_screenshot.png)

```bash
pip install "phytclust[gui]"
phytclust gui
```

This opens <http://127.0.0.1:8000>. Run `phytclust gui --help` for options
(`--host`, `--port`, `--reload`, `--no-browser`), or launch
`uvicorn phytclust.gui.api:app` directly.

Use the CLI or the Python API for reproducible analysis.

## Documentation

Tutorials, reference, and background:
[schwarzlab-ccb.github.io/PhytClust](https://schwarzlab-ccb.github.io/PhytClust/).

Release history and upgrade notes: [CHANGELOG.md](CHANGELOG.md).

To build the docs locally:

```bash
pip install -e ".[docs]"
mkdocs serve --dev-addr 127.0.0.1:8001
```

The default MkDocs port is 8000, the same one the GUI uses, so `--dev-addr`
avoids a collision if both are running.

## Citation

If you use PhytClust, please cite the paper:

> K. Ganesan, E. Billard, T.L. Kaufmann, C.B. Strange, M.C. Cwikla,
> A.M. Altenhoff, C. Dessimoz, R.F. Schwarz. *PhytClust: efficient and optimal
> monophyletic partitioning of rooted phylogenetic trees*. bioRxiv (2026).
> https://doi.org/10.64898/2025.12.11.693738

and, to cite a specific version of the software, the Zenodo record:

> https://doi.org/10.5281/zenodo.19731229

<!-- TODO(kat): swap the preprint for the journal reference on acceptance, and
     check this block matches CITATION.cff — which currently says 0.1.1 and needs
     updating to 1.0.0 with the current year before tagging. -->

## License

GPL-3.0 — see [LICENSE](LICENSE).
