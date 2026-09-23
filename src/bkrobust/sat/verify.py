"""Turn a SAT answer into a machine-checked certificate.

A solver returning SAT is a claim, not a proof: the claim is only as good as the
encoding, and the encoding is the thing under test. Every witness this module
touches is therefore replayed through the ordinary graph code -- the same
:mod:`bkrobust.demo.meek` and :mod:`bkrobust.search.space_fixed` that the
brute-force search uses -- and re-checked from scratch.

For an E3 walk that means three independent facts, none of which consults the
encoding:

1. every state is a genuine knowledge state of the CPDAG;
2. every consecutive pair differs by **exactly one orientation**, and is hence a
   covering pair unconditionally (nothing lies strictly between two sets that
   differ by a single element -- no appeal to Lemma R or anti-exchange);
3. the final state really does fail.

Together these certify ``d_BFS(G0, ·) <= len(walk) - 1`` with a failure at the
end, without trusting the encoding at all.
"""

from __future__ import annotations

from dataclasses import dataclass

from bkrobust.core.oracle import is_valid
from bkrobust.demo.graph import MPDAG, Edge, canon
from bkrobust.search.space_fixed import is_knowledge_state


@dataclass(frozen=True)
class WalkCertificate:
    """The verdict of replaying an E3 walk, with the first failure located.

    Attributes:
        ok: True iff every check passed.
        length: Number of steps, i.e. the certified distance bound.
        monotone: True iff the walk only ever retracts. Conjecture 2 predicts
            that a shortest walk can always be taken monotone; a *verified*
            non-monotone walk that is shorter than the retraction-only radius is
            what a refutation looks like.
        problems: Human-readable descriptions of every check that failed.
    """

    ok: bool
    length: int
    monotone: bool
    problems: tuple[str, ...]


def to_mpdag(cpdag: MPDAG, oriented: frozenset[Edge]) -> MPDAG:
    """Rebuild an MPDAG from a solver's set of directed edges.

    Args:
        cpdag: The CPDAG, which supplies the skeleton.
        oriented: The ``(tail, head)`` pairs the solver set to true.

    Returns:
        The MPDAG with those edges directed and the rest of the skeleton
        undirected.
    """
    directed = set(oriented)
    und = {
        canon(a, b)
        for a, b in cpdag.skeleton()
        if (a, b) not in directed and (b, a) not in directed
    }
    return MPDAG(sorted(cpdag.nodes), directed=sorted(directed), undirected=sorted(und))


def verify_walk(
    cpdag: MPDAG,
    g0: MPDAG,
    walk: tuple[frozenset[Edge], ...],
    x: str,
    y: str,
    z: frozenset[str],
) -> WalkCertificate:
    """Replay an E3 walk through the ordinary graph code and check every step.

    Args:
        cpdag: The CPDAG.
        g0: The analyst's graph; the walk must start here.
        walk: One orientation set per state, as returned by
            :func:`bkrobust.sat.e3.radius_e3`.
        x: Treatment.
        y: Outcome.
        z: The adjustment set.

    Returns:
        A :class:`WalkCertificate`.
    """
    problems: list[str] = []
    if not walk:
        return WalkCertificate(False, 0, True, ("empty walk",))

    graphs = [to_mpdag(cpdag, s) for s in walk]
    if graphs[0] != g0:
        problems.append("walk does not start at G0")
    for i, g in enumerate(graphs):
        if not is_knowledge_state(cpdag, g):
            problems.append(f"state {i} is not a knowledge state of the CPDAG")

    monotone = True
    for i in range(len(walk) - 1):
        lo, hi = walk[i], walk[i + 1]
        sym = (lo - hi) | (hi - lo)
        if len(sym) != 1:
            problems.append(f"step {i}->{i + 1} changes {len(sym)} orientations, not 1")
        if not hi < lo:
            monotone = False
            if not lo < hi:
                problems.append(f"step {i}->{i + 1} is not comparable")

    if is_valid(z, graphs[-1], x, y):
        problems.append("final state does not actually fail")

    return WalkCertificate(not problems, len(walk) - 1, monotone, tuple(problems))
