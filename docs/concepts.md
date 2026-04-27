# How It Works

This page explains what PhytClust does under the hood. No equations required — just intuition about trees and clustering.

---

## The problem

You have a rooted phylogenetic tree and you want to split its leaves into groups. Sounds simple, but there's a catch: you want each group to be a **clade** — a complete subtree where every descendant is included. This property is called **monophyly**.

Why does monophyly matter? Because a clade has a clear evolutionary interpretation: it's "everything descended from this ancestor." If your clusters aren't monophyletic, they're mixing lineages in ways that are hard to justify biologically.

The question, then, is: *of all possible ways to partition a tree into k monophyletic groups, which one is best?*

## Where cuts can happen

Think of the tree as a branching structure where each internal edge is a potential "cut point." If you cut an edge, you separate the subtree below it from the rest. A clustering is just a set of cuts that divides all leaves into *k* groups.

Not every combination of cuts produces monophyletic clusters — you need the cuts to be *consistent*. PhytClust's dynamic programming guarantees this by working bottom-up through the tree, so it only ever considers valid partitions.

## The objective: minimize within-cluster dispersion

Given a valid partition, PhytClust scores it by adding up, for every cluster, **how far each member is from its cluster's MRCA**. The MRCA (most recent common ancestor) is the deepest node that contains every leaf of the cluster — think of it as the cluster's centre of mass on the tree. Shorter average distance to the MRCA means a tighter, more cohesive group.

Intuitively: long branches *between* clusters (where the cuts land) are good — they separate well-differentiated lineages. Short distances *within* a cluster are good — they keep similar taxa together.

### Why distance-to-MRCA, not something else?

It's worth being precise about what gets summed, because there are several reasonable-sounding alternatives and the choice has consequences.

Take a cluster with three leaves: A, B, C. Three plausible cost definitions:

1. **Sum of every internal branch in the cluster's subtree.** Simple, but it ignores which leaves are present. Two clusters with the same skeleton but different membership would score identically.
2. **Sum of pairwise distances between members.** Faithful to the "spread out" intuition, but it scales as O(*n*²) per cluster and doesn't decompose neatly when you walk up the tree.
3. **Sum of leaf-to-MRCA distances** (what PhytClust uses). Each leaf contributes the length of its path up to the MRCA.

The third option has a property the others don't: it factors recursively. When you extend a cluster up by one branch — say you decide to merge it with its sibling at the parent — the new cost is just the old cost plus `(branch length) × (number of leaves in the cluster)`. Every leaf below the new branch picks up that branch's length once.

That's exactly the recurrence the DP uses, and it's why the algorithm runs in O(*n* · *k*²) instead of something cubic. Pairwise distance would force a more expensive update at every node; subtree-sum would lose track of membership.

There's also a clean interpretation: minimising the sum of leaf-to-MRCA distances is equivalent to minimising **average dispersion around the cluster's centre on the tree**. It's the closest tree-aware analogue to "within-cluster sum of squares" in flat clustering.

## Dynamic programming on the tree

The algorithm works bottom-up:

1. **Start at the leaves.** Each leaf is trivially a cluster of size 1.
2. **At each internal node**, decide how to allocate its *children* across clusters. The DP table records the best achievable cost for "partition the subtree rooted here into *j* groups" for every valid *j*.
3. **At the root**, read off the optimal partition for your target *k*.

Because the DP visits each node once and considers all valid splits at each node, the result is **exact** — not an approximation or heuristic. For a binary tree with *n* leaves, the time complexity is roughly O(*n* · *k*²), which is fast enough for trees with thousands of leaves.

## The score curve: finding good *k* values

If you don't know the right *k*, PhytClust computes the optimal cost for every *k* from 2 up to some maximum. Plotting cost against *k* gives a **score curve**.

The score curve typically drops sharply at *k* values where the tree has natural breakpoints — places where a long branch separates two subtrees. These sharp drops appear as **peaks** in the rate of improvement.

PhytClust uses scipy's peak detection to find these peaks, then ranks them by **prominence** — how much a peak stands out from its neighbors. The most prominent peaks correspond to the *k* values where the tree structure most strongly suggests a cluster boundary.

### Raw vs. adjusted ranking

By default, PhytClust uses "adjusted" ranking, which accounts for how the number of outlier clusters changes with *k*. This prevents the ranker from favoring *k* values that only look good because they produce many tiny noise clusters. You can control the balance between raw and adjusted ranking with the `lambda_weight` parameter.

## Polytomies

A polytomy is an internal node with more than two children. They show up everywhere in real phylogenies, for two very different reasons:

- **Soft polytomies** — the data couldn't resolve the branching order, so an inference tool collapsed an uncertain bifurcation. The "true" tree is binary, you just can't see it.
- **Hard polytomies** — multiple lineages genuinely diverged in rapid succession (rapid radiation, gene duplication bursts), and there is no internal binary structure to recover.

The clustering decision at a polytomy is not obvious. Picture a node with five children:

```
          ┌── A
          ├── B
   ┌──────┤
───┤      ├── C
   │      ├── D
   │      └── E
   └── F
```

If you want to keep that node's subtree as two clusters, where do you draw the line? `{A,B}` and `{C,D,E}`? `{A}` and `{B,C,D,E}`? Any of those? The tree itself doesn't tell you; the polytomy is by definition *unresolved* about that ordering.

PhytClust gives you two strategies:

**Hard mode** (default). Each child of the polytomy goes entirely into one cluster — no child is split between two clusters. The DP just has more children to allocate at the polytomy node, which is still tractable. With five children you have at most a handful of partitions to enumerate per `k`.

**Soft mode** (`--optimize-polytomies`). The DP is allowed to merge any subset of the polytomy's children into the same cluster, including subsets that wouldn't be a clade in any binary resolution. This is strictly more expressive — it can find clusterings that hard mode cannot represent — but the search space is the set of partitions of the children, which grows like the Bell number of the child count. A node with 12 children has roughly 4 million ways to partition its children; a node with 18 has about 10 billion. PhytClust falls back to hard mode automatically once the degree crosses `soft_polytomy_max_degree` (default 18).

When does soft mode actually matter? Take the example above. Suppose the branch leading to `{A,B,C,D,E}` is short, but within that group, `A` and `F` are biologically similar. Hard mode forces you to choose: either `A` is in the same cluster as `B,C,D,E` or in the same cluster as `F` — never both. Soft mode lets the DP put `A` with `F` and leave `{B,C,D,E}` together, by treating the polytomy as if it were resolved into `((F, A), (B,C,D,E))` for the purposes of clustering, even though no such bifurcation appears in the input.

In practice: hard mode is the right default. Try soft mode if you have polytomies of moderate degree (5–15) and the resulting clusters look forced.

## Outlier handling

Small clusters are a fact of life in real phylogenies — isolated taxa on long branches, contaminant sequences, or just biological noise. PhytClust gives you two ways to deal with them:

- **Hard constraint** (`min_cluster_size`): the DP enforces a minimum cluster size during optimization. If a *k* can't be achieved without violating this, it's simply not returned.

- **Soft handling** (`outlier_size_threshold` + `prefer_fewer_outliers`): clusters below the threshold are *marked* as outliers (cluster ID = -1) in the output, but they still exist. With `prefer_fewer_outliers`, the DP slightly favors solutions that concentrate noise into fewer groups rather than scattering it.

## Zero-length edges

Some internal branches in a phylogeny have length zero. There are a couple of common reasons:

- **Identical sequences.** If two leaves have no observed substitutions between them, the inference produces a branch of length zero somewhere in their lineage. Biologically the two leaves are indistinguishable; the zero-length edge is real but uninformative for clustering.
- **Collapsed weakly-supported nodes.** Some tools collapse low-support edges to zero length so downstream consumers treat them as polytomies. The branch is there topologically but has no claim to representing actual evolutionary distance.

Either way, splitting at a zero-length edge is *free* under the default objective — cutting there costs nothing, because the branch itself contributes zero. So the DP will happily put two identical sequences into separate clusters if it shaves a fraction off the cost elsewhere.

That's almost always wrong. If the data can't tell `A` and `A'` apart, neither should the clustering. Putting them in different groups is a numerical artefact, not a biological signal.

`--no-split-zero-length` forbids the DP from placing a cluster boundary on a zero-length edge. Identical-sequence neighbours stay together. Think of it as the cluster equivalent of "don't break a tie by coin flip" — if the tree is genuinely undecided about whether two taxa are distinct, don't pretend you know.

Two practical notes:

- The flag uses a small numerical tolerance (`zero_length_eps`, default `1e-12`) so trees with floating-point noise still behave correctly.
- If a polytomy sits at the bottom of an all-zero subtree, no-split-zero-length effectively keeps that whole subtree as one cluster. That's intentional: the subtree contains no information about how to split it.

## Summary

| Concept | What it means in PhytClust |
|---------|---------------------------|
| **Monophyly** | Every cluster is a complete clade — guaranteed, not post-hoc filtered |
| **Objective** | Minimise total leaf-to-MRCA distance within clusters (decomposes recursively, unlike pairwise distance) |
| **DP** | Exact bottom-up algorithm on the tree; visits each node once |
| **Score curve** | Cost vs. *k*; peaks indicate natural cluster boundaries |
| **Peak ranking** | Prominence-based, with optional outlier adjustment |
| **Polytomies** | Hard mode keeps each child in one cluster; soft mode allows arbitrary subsets at the cost of exponential blow-up |
| **Outliers** | Hard constraint (min size) or soft marking (threshold + preference) |
| **Zero-length edges** | Optional flag to forbid cuts where the branch carries no information |
