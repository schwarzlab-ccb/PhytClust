"""The tree dynamic program: build the cost table, then read a partition out.

- `table`    walks the tree in postorder and fills one cost-vs-k row per node,
             then backtracks a chosen k into leaf -> cluster assignments.
- `merge`    combines two child rows into a parent row. Bifurcating nodes do
             this once; hard multifurcations fold their children with it.
- `polytomy` the multifurcation recurrences, hard and soft.
- `costs`    per-cluster cost terms, penalties, and tie tolerances.
- `tied_partitions` lists cost-tied alternatives, with the selected partition first.
"""

from .table import (
    backtrack,
    cluster_map,
    compute_dp_table,
    prepare_tree,
    validate_clustering_parameters,
)
from .tied_partitions import enumerate_tied_optima

__all__ = [
    "backtrack",
    "cluster_map",
    "compute_dp_table",
    "enumerate_tied_optima",
    "prepare_tree",
    "validate_clustering_parameters",
]
