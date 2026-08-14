"""Fast tree traversal.

``Tree.find_clades()`` runs Bio.Phylo's attribute matcher on every node even
when no filter is given: it pops the node's children, runs a search against the
node, then restores them. On a few thousand nodes that overhead is around a
fifth of a PhytClust run, and none of it does any work we need.

``iter_clades`` returns the same nodes in the same order for an unfiltered
walk, without the matcher. Use ``find_clades`` when actually filtering (by
name, terminal, and so on) -- that is what it is for.
"""

from __future__ import annotations

from typing import Any, Iterator


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
