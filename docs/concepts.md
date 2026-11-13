# How PhytClust works

## The objective

PhytClust partitions a rooted phylogenetic tree into *k* clusters. Each
cluster is a clade containing all descendants of its root, with an allowance
for unresolved branching in soft polytomy mode.

For each cluster, PhytClust sums the distances from its leaves to their most
recent common ancestor (MRCA). The default objective minimises this total
across all clusters, favouring groups whose leaves are close to their common
ancestor.

A dynamic program combines subtree solutions from the leaves to the root,
then recovers the minimum-cost partition for the requested *k*. For binary
trees with *n* leaves, the time complexity is O(*n* · *k*²).

## Choosing *k*

When *k* is unknown, PhytClust uses peaks in a score curve to suggest candidate
partitions. Let β(*k*) be the minimum partition cost and *n* the leaf count:

```text
V(k) = [(β(1) − β(k)) / β(k)] × [(n − k) / k]
E(k) = [β(k−1) − β(k)] / [β(k) − β(k+1)]
score(k) = E(k) × V(k)
```

V measures explained dispersion relative to what remains, weighted by *k*.
E identifies elbows i.e. a large reduction in cost followed by a small one.
It is set to zero when either neighbouring reduction is below
`1e-12 × β(1)`.

Peaks are ranked by a blend of prominence and score height. Global selection
returns the highest-ranked peaks; resolution mode selects one per
logarithmic bin to cover different scales. Use `--exclude-k2` if you want to
exclude the root split from automatic selection.

## Polytomies

A polytomy is a node with more than two children.

- **Soft mode** (default) treats the branching order as unresolved. Subsets
  of child subtrees may form clusters through implied zero-length branches.
- **Hard mode** preserves the represented clades strictly. A cluster must
  contain the whole subtree at that node or stay within one child subtree.

Soft mode becomes expensive as the number of children grows and is limited
to 12 children per node by default. Use `--polytomy-mode hard` for larger
multifurcations.

## Small clusters and zero-length branches

`min_cluster_size` enforces a minimum during optimisation. Outlier marking
instead labels small clusters `-1` in the exported table; it marks singletons
by default. `prefer_fewer_outliers` changes the optimisation to minimise
outlier count before cost.

`--no-split-zero-length` keeps subtrees whose branches are all zero in a
single cluster. Individual zero-length edges elsewhere may still be cluster
boundaries. Branches at or below `1e-12` count as zero by default.

See the [reference](reference.md) for settings and defaults.
