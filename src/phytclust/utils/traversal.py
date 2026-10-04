"""Fast tree traversal.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, Iterator

import numpy as np


def iter_clades(root: Any, order: str = "preorder") -> Iterator[Any]:
    """Yield every clade at or below ``root``.

    ``order`` is "preorder", "postorder" or "level", matching the orders
    Bio.Phylo's ``find_clades`` produces.
    """
    if order == "preorder":
        stack = [root]
        while stack:
            node = stack.pop()
            yield node
            stack.extend(reversed(node.clades))

    elif order == "postorder":
        stack = [(root, False)]
        while stack:
            node, expanded = stack.pop()
            if expanded:
                yield node
            else:
                stack.append((node, True))
                stack.extend((child, False) for child in reversed(node.clades))

    elif order == "level":
        queue = [root]
        i = 0
        while i < len(queue):
            node = queue[i]
            i += 1
            yield node
            queue.extend(node.clades)

    else:
        raise ValueError(
            f"order must be 'preorder', 'postorder' or 'level', got {order!r}"
        )


def terminals(root: Any) -> Iterator[Any]:
    """Yield the leaves at or below ``root``, in preorder."""
    for node in iter_clades(root, "preorder"):
        if not node.clades:
            yield node


def nonterminals(root: Any) -> Iterator[Any]:
    """Yield the internal nodes at or below ``root``, in preorder."""
    for node in iter_clades(root, "preorder"):
        if node.clades:
            yield node


_HASH_MASK = (1 << 64) - 1


def _fingerprint_pass(
    root: Any, internal_names: bool, keep_all: bool
) -> dict[int, int]:
    """Bottom-up subtree hashes keyed by ``id(clade)``.

    Children are combined by summing their hashes, which is order-free and
    needs no sort. Reversing a stack preorder visits children before parents
    without the postorder generator's per-node bookkeeping.
    """
    order = []
    stack = [root]
    while stack:
        node = stack.pop()
        order.append(node)
        stack.extend(node.clades)
    hashes: dict[int, int] = {}
    take = hashes.get if keep_all else hashes.pop
    for node in reversed(order):
        total = 0
        for child in node.clades:
            total += take(id(child))
        name = node.name if internal_names or not node.clades else None
        hashes[id(node)] = hash(
            (
                name,
                node.branch_length or 0.0,
                node.confidence,
                total & _HASH_MASK,
                len(node.clades),
            )
        )
    return hashes


def clade_fingerprints(
    root: Any, internal_names: bool = True
) -> dict[int, int]:
    """Return ``id(clade) -> subtree hash``, ignoring child order.

    Hashes include topology, branch lengths, confidence values, and names.
    With ``internal_names=False``, only leaf names are included. Missing
    branch lengths count as 0.0.

    Different subtrees may have the same hash because collisions are possible.
    Values are intended for comparisons within one process, not persistent
    identifiers across runs.
    """
    return _fingerprint_pass(root, internal_names, keep_all=True)


def tree_fingerprint(root: Any) -> int:
    """Hash of the tree under ``root`` that ignores the order of children."""
    return _fingerprint_pass(root, True, keep_all=False)[id(root)]


def match_child_order(target: Any, source: Any) -> bool:
    """Reorder ``target``'s children in place to follow ``source``'s order.

    Match subtrees using fingerprints that ignore internal node names.
    Return False without modifying ``target`` if matching fails. Return
    True after applying the matched child order.

    Matching relies on hashes, so hash collisions can cause incorrect matches.
    """
    target_fp = clade_fingerprints(target, internal_names=False)
    source_fp = clade_fingerprints(source, internal_names=False)
    if target_fp[id(target)] != source_fp[id(source)]:
        return False
    orders: list[tuple[Any, list[Any]]] = []
    stack = [(target, source)]
    while stack:
        t_node, s_node = stack.pop()
        pool: dict[int, list[Any]] = {}
        for child in t_node.clades:
            pool.setdefault(target_fp[id(child)], []).append(child)
        ordered = []
        for s_child in s_node.clades:
            bucket = pool.get(source_fp[id(s_child)])
            if not bucket:
                return False
            t_child = bucket.pop()
            ordered.append(t_child)
            stack.append((t_child, s_child))
        orders.append((t_node, ordered))
    for node, ordered in orders:
        node.clades[:] = ordered
    return True


class LeafSpans(Mapping):
    """Map nodes to descendant leaves using views into one shared array.

    Each node stores a ``(start, end)`` span in a depth-first leaf ordering.
    Values are NumPy views containing the original leaf objects.

    Spans describe the tree when this mapping was built and do not update
    after tree edits. Rebuild after changes to topology or child order.
    """

    def __init__(self, leaf_order: np.ndarray, spans: dict[Any, tuple[int, int]]):
        self._leaves = leaf_order
        self._spans = spans

    def __getitem__(self, node: Any) -> np.ndarray:
        start, end = self._spans[node]
        return self._leaves[start:end]

    def __iter__(self) -> Iterator[Any]:
        return iter(self._spans)

    def __len__(self) -> int:
        return len(self._spans)

    def __contains__(self, node: object) -> bool:
        return node in self._spans


def leaf_maps(root: Any) -> tuple[LeafSpans, dict[Any, int]]:
    """Build ``(leaves under each node, leaf count per node)`` in one pass."""
    leaves: list[Any] = []
    spans: dict[Any, tuple[int, int]] = {}
    counts: dict[Any, int] = {}
    for node in iter_clades(root, "postorder"):
        if not node.clades:
            spans[node] = (len(leaves), len(leaves) + 1)
            leaves.append(node)
            counts[node] = 1
        else:
            start = spans[node.clades[0]][0]
            end = spans[node.clades[-1]][1]
            spans[node] = (start, end)
            counts[node] = end - start
    leaf_order = np.empty(len(leaves), dtype=object)
    leaf_order[:] = leaves
    return LeafSpans(leaf_order, spans), counts
