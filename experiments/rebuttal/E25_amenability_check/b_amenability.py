"""Amenability read literally from Perković, Kalisch and Maathuis (UAI 2017).

Definition 3.1: a path p = <X = V0, ..., Vk = Y> is b-possibly causal in the maximal PDAG G iff
no edge Vi <- Vj with 0 <= i < j <= k is in G (edges off p count too).
Definition 4.3 (b-amenability): every proper b-possibly causal path from X to Y starts with a
directed edge out of X.

Written for E25 as a second implementation that shares no code with the submitted
``is_amenable`` (which decides the same property semantically, by one Meek closure per
undirected neighbour of X). It enumerates simple paths depth first, so it is exponential in the
worst case: a budget on path extensions returns None (undecided) instead of an answer.
"""

from __future__ import annotations

from bkrobust.demo.graph import MPDAG


def b_amenable(g: MPDAG, x: str, y: str, max_steps: int = 200_000) -> bool | None:
    """True if b-amenable, False if a proper b-possibly causal path starts undirected, None if
    the step budget ran out first. Single nodes x != y; the path is proper because it is simple."""
    directed = g.directed_edges  # set of (tail, head)
    steps = 0

    def forward(a: str) -> list[str]:
        # neighbours reachable along an edge that is a -> w or a - w
        out = [w for w in g.children(a)]
        out += [w for w in g.neighbors(a)]
        return sorted(set(out))

    starts = sorted(g.neighbors(x))  # undirected edges x - v only: those are the witnesses
    for v in starts:
        # A path x - v ... ; check the pair (x, v): v -> x must not be in G (it is undirected).
        # iterative DFS with explicit path lists (paths are short in practice)
        stack: list[tuple[list[str], int]] = [([x, v], 0)]
        while stack:
            path, _ = stack.pop()
            steps += 1
            if steps > max_steps:
                return None
            last = path[-1]
            if last == y:
                return False
            nodes = set(path)
            for w in forward(last):
                if w in nodes:
                    continue
                # b-possibly causal: no edge w -> Vi for any earlier Vi on the path
                if any((w, vi) in directed for vi in path):
                    continue
                stack.append((path + [w], 0))
    return True
