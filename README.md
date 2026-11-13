<p align="center">
  <img src="https://raw.githubusercontent.com/schwarzlab-ccb/PhytClust/master/docs/img/phytclust_logo.png" alt="PhytClust" width="360">
</p>

[![PyPI](https://img.shields.io/pypi/v/phytclust.svg)](https://pypi.org/project/phytclust/)
[![CI](https://github.com/schwarzlab-ccb/PhytClust/actions/workflows/ci.yml/badge.svg)](https://github.com/schwarzlab-ccb/PhytClust/actions/workflows/ci.yml)
[![License: GPL-3.0](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://github.com/schwarzlab-ccb/PhytClust/blob/master/LICENSE)
[![Manuscript data](https://zenodo.org/badge/DOI/10.5281/zenodo.19731229.svg)](https://doi.org/10.5281/zenodo.19731229)

PhytClust partitions rooted phylogenetic trees into monophyletic clusters.
Specify the number of clusters or let PhytClust suggest candidate partitions.
Available as a command-line tool, Python API, and experimental web interface.

## Install

Requires Python 3.10 or later.

```bash
pip install phytclust
```

## Usage

Use a rooted tree in Newick format:

```bash
# Five clusters
phytclust tree.nwk --k 5 --save-fig

# Up to three automatically selected partitions
phytclust tree.nwk --top-n 3 --save-fig
```

Assignments and plots are saved under `results/`.

```python
from phytclust import PhytClust

pc = PhytClust("tree.nwk")
result = pc.run(k=5)
```

## Web interface (experimental)

```bash
pip install "phytclust[gui]"
phytclust gui
```

Open <http://127.0.0.1:8000> and upload a Newick file.

## Documentation

[Usage, reference, and algorithm](https://schwarzlab-ccb.github.io/PhytClust/).

## Citation

If you use PhytClust, please cite the paper:

> K. Ganesan, E. Billard, T.L. Kaufmann, C.B. Strange, M.C. Cwikla,
> A.M. Altenhoff, C. Dessimoz, R.F. Schwarz. *PhytClust: efficient and optimal
> monophyletic partitioning of rooted phylogenetic trees*. bioRxiv (2026).
> https://doi.org/10.64898/2025.12.11.693738

The manuscript's data and analysis scripts are archived on
[Zenodo](https://doi.org/10.5281/zenodo.19731229). Machine-readable citation
metadata is available in
[CITATION.cff](https://github.com/schwarzlab-ccb/PhytClust/blob/master/CITATION.cff).

PhytClust is developed at the
[Institute for Computational Cancer Biology](https://iccb-cologne.org/home).

<a href="https://iccb-cologne.org/home">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/schwarzlab-ccb/PhytClust/master/docs/img/iccb_logo_dark.png">
    <img src="https://raw.githubusercontent.com/schwarzlab-ccb/PhytClust/master/docs/img/iccb_logo.png" alt="ICCB" width="170">
  </picture>
</a>

## License

See [LICENSE](https://github.com/schwarzlab-ccb/PhytClust/blob/master/LICENSE).
