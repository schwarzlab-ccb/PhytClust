# PhytClust

**Tree-aware clustering for rooted phylogenies.**

PhytClust partitions the leaves of a phylogenetic tree into clusters that are guaranteed to be monophyletic — every cluster is a proper clade. Under the hood it uses dynamic programming on the tree structure itself, so the optimization respects the topology rather than flattening it into a distance matrix.

It ships as a single `pip install` with three interfaces: a **CLI** for quick runs, a **Python API** for scripting and pipelines, and an experimental **web GUI** for interactive exploration.

## Install

```bash
pip install phytclust
```

Or from source:

```bash
git clone https://github.com/ICCB-Cologne/PhytClust.git
cd PhytClust
pip install -e ".[dev]"
```

## Quick taste

```bash
# cluster a tree into exactly 5 clades
phytclust examples/sample_tree.nwk --k 5 --save-fig

# let PhytClust pick the 3 best k values
phytclust examples/sample_tree.nwk --top-n 3 --save-fig

# get one representative k per resolution scale
phytclust examples/sample_tree.nwk --resolution --bins 4 --save-fig
```

## Web GUI [beta]

PhytClust comes with an experimental browser-based interface for interactive exploration.

**Install with GUI extras:**

```bash
pip install "phytclust[gui]"
```

**Launch the server:**

```bash
uvicorn phytclust.gui.api:app --host 127.0.0.1 --port 8000
```

Then open **http://127.0.0.1:8000** in a browser. Paste any Newick string, choose a clustering mode, and explore the results interactively.

For a public deployment (multi-user, capped tree size, no server-side writes):

```bash
PHYTCLUST_PUBLIC_MODE=1 PHYTCLUST_MAX_TIPS=5000 \
  uvicorn phytclust.gui.api:app --host 0.0.0.0 --port 8000
```

## Documentation

Full docs (tutorials, reference, concepts) are hosted at:

**[iccb-cologne.github.io/PhytClust](https://iccb-cologne.github.io/PhytClust/)**

If you want to browse or build the docs locally:

```bash
pip install -e ".[docs]"
mkdocs serve          # opens at http://127.0.0.1:8000
```

## Citation

If PhytClust is useful in your research, please cite:

> K. Ganesan, E. Billard, T.L. Kaufmann, C.B. Strange, M.C. Cwikla, A. Altenhoff, C. Dessimoz, R.F. Schwarz.
> *PhytClust* (2025). [github.com/ICCB-Cologne/PhytClust](https://github.com/ICCB-Cologne/PhytClust)

## License

GPL-3.0 — see [LICENSE](LICENSE) for details.
