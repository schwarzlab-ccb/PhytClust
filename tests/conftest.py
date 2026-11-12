import pathlib
from io import StringIO

import matplotlib

matplotlib.use("Agg")

from Bio import Phylo  # noqa: E402

EXAMPLES = pathlib.Path(__file__).resolve().parent.parent / "examples"
SAMPLE_TREE = EXAMPLES / "sample_tree.nwk"
POLYTOMY_TREE = EXAMPLES / "sample_polytomy.newick"
BOOTSTRAP_TREE = EXAMPLES / "sample_bootstrap.newick"
CANCER_TREE = EXAMPLES / "OV_25_cancer.newick"


def load_tree(path=SAMPLE_TREE):
    return Phylo.read(path, "newick")


def tree_from_newick(newick):
    return Phylo.read(StringIO(newick), "newick")


def leaf_groups(cluster_map):
    """Leaf-name sets, one per cluster, order-independent."""
    groups = {}
    for leaf, cluster_id in cluster_map.items():
        groups.setdefault(cluster_id, set()).add(leaf.name)
    return frozenset(frozenset(group) for group in groups.values())


def random_newick(seed, n_leaves):
    """A random rooted binary tree with positive branch lengths."""
    import random

    rng = random.Random(seed)
    names = [f"t{i}" for i in range(n_leaves)]

    def build(n):
        length = rng.uniform(0.05, 3.0)
        if n == 1:
            return f"{names.pop()}:{length:.6f}"
        left = rng.randint(1, n - 1)
        return f"({build(left)},{build(n - left)}):{length:.6f}"

    return build(n_leaves) + ";"
