# How PhytClust works

---

## The problem

Given a rooted phylogenetic tree, split its leaves into groups such that every
group is a **clade** — a complete subtree containing all descendants of its
root. This property is monophyly, and it is what makes a cluster interpretable.
A clade is everything descended from one ancestor.

The question is which of the possible partitions of a tree into *k* monophyletic
groups is best.

## Where cuts can happen

Each internal edge is a potential cut point. Cutting an edge separates the
subtree below it from the rest, and a clustering is a set of cuts that divides
the leaves into *k* groups.

Not every combination of cuts yields monophyletic clusters; the cuts have to be
mutually consistent. The dynamic program works bottom-up and only ever
enumerates valid partitions, so monophyly is guaranteed by construction rather
than enforced by filtering afterwards.

## The objective: within-cluster dispersion

A partition is scored by summing, over every cluster, the distance from each
member leaf to that cluster's MRCA — the deepest node containing every leaf of
the cluster. Shorter mean distance to the MRCA means a tighter group. Long
branches at the cuts separate well-differentiated lineages; short distances
inside a cluster keep similar taxa together.

### Why distance-to-MRCA

Three cost definitions are plausible for a cluster containing leaves A, B, C:

1. **Sum of every internal branch in the cluster's subtree.** Ignores which
   leaves are present, so two clusters with the same skeleton and different
   membership score identically.
2. **Sum of pairwise distances between members.** Matches the intuition, but
   costs O(*n*²) per cluster and does not decompose when walking up the tree.
3. **Sum of leaf-to-MRCA distances.** What PhytClust uses.

The third factors recursively. Extending a cluster up by one branch — merging it
with its sibling at the parent — gives a new cost of the old cost plus
`(branch length) × (number of leaves in the cluster)`, because every leaf below
the new branch picks up that branch's length exactly once.

That recurrence is what the DP implements, and it is why the algorithm runs in
O(*n* · *k*²). Pairwise distance forces a more expensive update at every node.

Minimising the sum of leaf-to-MRCA distances is equivalent to minimising average
dispersion around the cluster's centre on the tree. The tree-aware counterpart
of within-cluster sum of squares in flat clustering.

## Dynamic programming on the tree

1. **At the leaves**, each leaf is a cluster of size one.
2. **At each internal node**, the DP decides how to allocate that node's children
   across clusters, recording the best achievable cost of partitioning the
   subtree rooted there into *j* groups, for every valid *j*.
3. **At the root**, the optimal partition for the target *k* is read off.

Each node is visited once and all valid splits at that node are considered, so
the result is exact rather than heuristic. For a binary tree with *n* leaves the
time complexity is O(*n* · *k*²), which is tractable for trees with thousands of
leaves.

## The score curve and choosing *k*

When *k* is not known in advance, PhytClust evaluates every *k* from 2 up to a
maximum and derives a curve over *k*. Peaks in that curve mark the *k* values at
which the tree's own structure most strongly suggests a cluster boundary,
typically where a long branch separates two subtrees.

Three quantities are involved, and only the last of them is "the score curve"
that the rest of the documentation refers to.

**1. The DP cost, β(*k*).** The minimum total within-cluster dispersion
achievable with *k* clusters — the objective defined above. It decreases
monotonically with *k*, so it has no peaks and is never itself peak-detected.

**2. The validity score, V(*k*).** How much of the dispersion has been explained
relative to what remains, weighted by resolution:

```text
V(k) = [ (β(1) − β(k)) / β(k) ] × [ (n − k) / k ]
```

**3. The elbow score, E(*k*).** How sharply the returns flatten immediately after
*k*, as the ratio of the improvement gained reaching *k* to the improvement
available just beyond it:

```text
E(k) = (β(k−1) − β(k)) / (β(k) − β(k+1))
```

E is large exactly at an elbow, where a large gain is followed by a small one. It
is clamped to at most 50 so that a near-zero denominator cannot dominate.

**The score curve is their product**, `score(k) = E(k) · V(k)`. This is what
`pc.scores` holds, what `scores.png` plots, and what `scipy.signal.find_peaks`
is run on. Requiring both factors is what makes a peak meaningful: V alone
favours large *k* on many trees, and E alone fires on any local flattening
however small the cluster structure behind it.

Peak detection runs on the *k* ≥ 3 portion of that curve (see
[skipping *k* = 2](#skipping-k-2)), and peaks are ranked by **prominence**, the
height of a peak relative to its neighbours.

### Relationship to the Calinski-Harabasz index

The validity score V(*k*) above is a **tree-adapted analogue of the
Calinski-Harabasz index, not the CH index itself**. The correspondence is close
enough to be worth stating precisely, because the differences are deliberate.

Textbook CH, for flat data in Euclidean space, is

```text
CH(k) = [ B(k) / (k − 1) ] / [ W(k) / (n − k) ]
      = [ B(k) / W(k) ] × [ (n − k) / (k − 1) ]
```

with `B` the between-group scatter and `W` the within-group scatter, both sums of
squared distances to centroids. PhytClust replaces each term with its tree
counterpart:

| CH term | PhytClust counterpart |
|---------|----------------------|
| Within-group scatter `W(k)` | `β(k)`, the summed leaf-to-MRCA distance under the optimal *k*-clustering |
| Between-group scatter `B(k)` | `β(1) − β(k)`, the dispersion explained by splitting into *k* clusters rather than one |
| Squared distance to a centroid | Path length along the tree to the cluster's MRCA |
| Degrees-of-freedom weight `(n − k)/(k − 1)` | `(n − k)/k` |

The first three substitutions are what make the criterion computable on a tree:
there is no centroid in tree space, so the MRCA takes its place, and `β(1)`
supplies the total-dispersion baseline that `B + W` provides in the flat case.

The last substitution is the one that departs from CH deliberately. The
denominator is ***k*, not *k* − 1**. In this setting the term is not acting as a
degrees-of-freedom correction on a significance test; it is a resolution weight
used to rank peaks against each other, and using *k* keeps it finite at *k* = 1,
where CH is undefined.

This is not a cosmetic change. The two forms differ by a factor of
`(k − 1)/k`, which is `0.5` at *k* = 2, `0.9` at *k* = 10 and approaches `1` as
*k* grows — so relative to CH, the PhytClust weight penalises small *k* and
converges with CH at large *k*. It can therefore reorder peaks, and it is
intended to: the smallest *k* values are the least informative on a tree (see
[skipping *k* = 2](#skipping-k-2)). The practical consequence is that the
reported numbers are not CH values and should not be compared against published
CH thresholds.

### How peaks are ranked

Two properties describe a peak, and they do not always agree:

- **Prominence** — how far the peak rises above the surrounding curve. A modest
  bump in a flat stretch is prominent; a taller bump on an already-high shoulder
  is not.
- **Absolute score** — the peak's height on the *y* axis, ignoring context.

`prominence_weight` sets the balance between them. At `1.0` peaks are ranked by
prominence alone, at `0.0` by absolute score alone, and the default `0.7` leans
on prominence while letting a high absolute score break near-ties. Because the
two quantities live on different scales, each is min–max normalised across the
detected peaks before being blended.

That normalise-then-blend behaviour is what `ranking_mode="adjusted"` (the
default) means. `ranking_mode="raw"` instead ranks by absolute prominence with
no normalisation and ignores `prominence_weight`.

### Skipping *k* = 2

*k* = 2 is the trivial first split. Every binary tree has a clean root-level cut,
and that cut almost always produces a large drop in cost, so left in the running
it dominates the top of the ranking on most trees. Halving the tree at the root
is rarely the reason anyone runs a clustering algorithm.

*k* = 2 is therefore dropped from automatic peak selection by default. In
`top_n` mode the next-ranked peak takes its slot; in resolution mode, a bin that
would have selected *k* = 2 advances to its next candidate.

`exclude_k2=False`, or `--include-k2` on the CLI, puts it back in the running.
An explicit `run(k=2)` always returns *k* = 2 — the flag affects automatic
selection only.

## Polytomies

A polytomy is an internal node with more than two children. Real phylogenies
contain them for two distinct reasons:

- **Soft polytomies** — the data could not resolve the branching order and an
  inference tool collapsed an uncertain bifurcation. The underlying tree is
  binary; the resolution is simply not recoverable from the data.
- **Hard polytomies** — several lineages diverged in rapid succession (rapid
  radiation, gene duplication bursts) and there is no internal binary structure
  to recover.

The clustering decision at a polytomy is genuinely underdetermined. Consider a
node with five children:

```text
          ,-- A
          |-- B
   ,------+-- C
---+      |-- D
   |      `-- E
   `-- F
```

To keep that node's subtree as two clusters, the split could be `{A,B}` and
`{C,D,E}`, or `{A}` and `{B,C,D,E}`, or others. The tree does not say which; the
polytomy is by definition unresolved about that ordering.

Two strategies are available.

**Hard mode** (`--polytomy-mode hard`, the default) treats the polytomy as
all-or-nothing. At the polytomy node the DP has exactly two kinds of option:

- merge **all** of the children into a single cluster rooted at the node, or
- merge **none** of them, so that every cluster lies entirely inside one child
  subtree.

There is no third option. A *proper subset* of the children — some but not all —
can never be grouped into one cluster. This is what makes it a faithful model of
a hard multifurcation: since the tree asserts that these lineages diverged
simultaneously, it provides no grounds for preferring `{A,B}` over `{B,C}`, and
hard mode declines to invent one.

**Soft mode** (`--polytomy-mode soft`) removes exactly that restriction. Any
subset of the children may be grouped, by treating them as though they hung from
a hidden zero-length node inside the polytomy. This is the resolution the data
could not supply, reintroduced as a search over groupings.

So for the five-child node above, with the target of splitting its subtree into
two clusters:

| Partition | Hard | Soft |
|-----------|------|------|
| `{A,B,C,D,E}` as one cluster | ✅ | ✅ |
| Each child its own cluster | ✅ | ✅ |
| `{A,B}` + `{C,D,E}` | ❌ | ✅ |
| `{A,C,E}` + `{B,D}` | ❌ | ✅ |

The third row is the distinction in a single line: **hard mode cannot produce
`{A,B}` + `{C,D,E}` at all.** Asked for that shape it must either pull all five
into one cluster or break all five apart, and the cost of those two options is
what it actually compares.

That expressiveness is expensive. Soft mode's search space is the set of
partitions of the children, which grows as the Bell number of the child count:
12 children give about 4.2 million partitions, 16 give about 10.5 billion, and 18
give roughly 680 billion. PhytClust therefore reverts to hard mode automatically
once the node degree exceeds `soft_polytomy_max_degree` (default 12).

Note that `F` in the diagram is a *sibling* of the polytomy, not one of its
children, so no polytomy setting can group `A` with `F` without also taking the
rest of the polytomy — that grouping is governed by the ordinary DP, not by
`--polytomy-mode`.

Hard mode is the right default. Soft mode is worth trying on polytomies of
moderate degree (5-15) where the resulting clusters look forced.

## Outlier handling

Small clusters are unavoidable in real phylogenies: isolated taxa on long
branches, contaminant sequences, biological noise. Two mechanisms handle them.

- **Hard constraint** (`min_cluster_size`). The DP enforces a minimum cluster
  size during optimisation. A *k* that cannot be reached without violating it is
  not returned at all.
- **Soft marking** (`outlier_size_threshold` with `prefer_fewer_outliers`).
  Clusters below the threshold are marked with ID `-1` in the output but still
  exist. `prefer_fewer_outliers` biases toward solutions that concentrate noise
  into fewer groups rather than scattering it.

`prefer_fewer_outliers` acts **inside the dynamic program**, not on peak
ranking. It changes how two candidate DP states are compared: normally cost is
minimised first and outlier count only breaks ties, whereas with the flag set the
outlier count is minimised first and cost becomes the tie-breaker. The effect is
therefore on which partition is chosen for a given *k*, not on which *k* is
selected. It requires `outlier_size_threshold` to be set, and raises without it.

## Zero-length edges

Internal branches of length zero arise for two common reasons:

- **Identical sequences.** With no observed substitutions between two leaves, the
  inference produces a zero-length branch somewhere in their lineage. The two
  leaves are indistinguishable in the data; the edge is real but uninformative
  for clustering.
- **Collapsed weakly-supported nodes.** Some tools set low-support edges to zero
  length so that downstream consumers treat them as polytomies. The branch exists
  topologically but does not represent evolutionary distance.

Splitting at a zero-length edge is free under the default objective (the branch
contributes nothing to the cost) so the DP will separate two identical
sequences into different clusters if doing so shaves anything off the cost
elsewhere. That is a numerical artefact, not a biological signal: if the data
cannot distinguish A from A', the clustering should not either.

`--no-split-zero-length` forbids a cluster boundary on a zero-length edge, so
identical-sequence neighbours stay together.

Two practical notes:

- The flag applies a numerical tolerance (`zero_length_eps`, default `1e-12`), so
  trees carrying floating-point noise behave as intended.
- If a polytomy sits at the bottom of an all-zero subtree, the flag keeps that
  whole subtree as a single cluster. That is intentional: the subtree carries no
  information about how to split it.

## Summary

| Concept | In PhytClust |
|---------|--------------|
| **Monophyly** | Every cluster is a complete clade, guaranteed by the DP rather than filtered post hoc |
| **Objective** | Minimise total leaf-to-MRCA distance within clusters; decomposes recursively, unlike pairwise distance |
| **DP** | Exact bottom-up algorithm; each node visited once; O(*n* · *k*²) |
| **Score curve** | Evaluated over *k*; peaks mark candidate cluster boundaries |
| **Peak ranking** | Prominence-based, with optional outlier adjustment |
| **Polytomies** | Hard mode keeps each child in one cluster; soft mode allows arbitrary subsets, at Bell-number cost |
| **Outliers** | Hard constraint on size, or soft marking with a threshold |
| **Zero-length edges** | Optional flag forbidding cuts where the branch carries no information |
