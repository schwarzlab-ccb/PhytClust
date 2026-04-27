# PhytClust

**Tree-aware clustering for rooted phylogenies.**

PhytClust takes a rooted phylogenetic tree and splits its leaves into clusters that are guaranteed to be monophyletic — every cluster is a complete clade. Monophyly isn't a post-hoc filter or a heuristic; it's a constraint the optimisation respects from the start.

The core algorithm is dynamic programming directly on the tree topology. It finds the partition of leaves into *k* groups that minimises within-cluster dispersion (the summed leaf-to-MRCA distances), while never breaking monophyly. Because the DP explores the full valid solution space for each *k*, the result is exact — not an approximation, and not sensitive to initial conditions or random seeds.

---

## Who is this for?

PhytClust was built for researchers who work with phylogenies and need to group taxa into meaningful clusters — for downstream analysis, visualisation, or summarising large trees. Concrete examples:

- Grouping clones from tumour phylogenies into subclonal populations.
- Partitioning gene or species trees into coherent groups for comparative analysis.
- Summarising metagenomic diversity at multiple resolutions in one pass.

If you have a rooted Newick tree and you want monophyletic clusters out of it, PhytClust is the right tool.

## What it does

You can run PhytClust in three ways depending on how much you already know about the tree:

- **Exact *k*** — "Give me exactly 5 clusters." When biology or a downstream pipeline already dictates the number.
- **Global peak search** — "Find the best *k* values." PhytClust scores every *k* in range and picks out the values where the tree has natural breakpoints. The output is the top *N* *k* values ranked by how prominent each one is.
- **Multi-resolution** — "Show me the tree at several scales." Splits the *k* range into logarithmic bins and returns one representative *k* per bin, giving you a coarse-to-fine panel from a single run.

Real trees are messy, so PhytClust also handles the standard headaches: minimum cluster size as a hard constraint, optional outlier marking for tiny noise clusters, native support for polytomies (no need to resolve them beforehand), an opt-in guard against splitting at zero-length branches, and outgroup / midpoint rooting on the fly.

There are three ways to talk to it: a **CLI** for quick runs and shell pipelines, a **Python API** for notebooks and batch processing, and an experimental **web GUI** for interactive exploration in the browser.

### What you get out

A normal run produces:

| File | What it is |
|------|-----------|
| `phytclust_results.tsv` | Leaf-to-cluster assignments (one or multiple *k* values) |
| `scores.png` | Score curve with annotated peaks |
| `tree_k{K}.png` | Coloured tree for each selected *k* |
| `peaks_by_rank.txt` | Selected *k* values in rank order |

Plots use a colourblind-safe palette and are configurable down to font sizes, marker sizes, and axis scales.

---

## Where to start

If you're new here, go through these pages in order:

1. **[Getting Started](getting-started.md)** — install, verify, run your first clustering
2. **[How It Works](concepts.md)** — what the algorithm actually does (no math prerequisites)
3. **[CLI Tutorial](tutorials/cli-end-to-end.md)** — practical walkthrough from first run to tuned output
4. **[Python Tutorial](tutorials/python-end-to-end.md)** — same journey, but in scripts and notebooks

When you need to look something up:

- [CLI Flags Reference](reference/cli.md)
- [Configuration Reference](reference/configuration.md)
- [Python API Reference](reference/api.md)
