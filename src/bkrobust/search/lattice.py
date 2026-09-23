"""Order-theoretic structure of the perturbation space: gradedness, joins, Lemma L.

The route to Conjecture 2 goes through the model-inclusion order rather than
through the failure predicate. Recall the orientation of the order: MORE
knowledge means FEWER represented DAGs means LOWER in the order, so the CPDAG is
the maximum and the DAGs are the minimal elements. "Moving up" is retraction.

Three questions, in the order they must be settled:

1. **Is the order graded?** No -- and this is a recorded dead end, not a result
   to rediscover. Upper semimodularity would imply the Jordan-Dedekind chain
   condition, and that already fails on the three-node chain: orienting
   ``a -> b`` propagates by Meek R1 all the way to a DAG in one covering step,
   while orienting ``b -> a`` propagates nothing and needs two. So maximal
   chains between the same endpoints have different lengths, and any proof of
   Conjecture 2 must use something weaker than semimodularity.
   :func:`chain_length_spread` measures this.

2. **Do joins exist?** The candidate construction is to keep the orientations
   the two graphs agree on and Meek-close: ``G v H = Meek(C, K_G intersect
   K_H)``. This is always an *upper* bound, since every DAG of ``[G]`` satisfies
   every orientation in the intersection; whether it is the *least* one is
   checked by :func:`is_least_upper_bound`.

3. **Does Lemma L hold?** ``d(G0, G v G0) <= d(G0, G)`` -- the join map is
   1-Lipschitz from the base point. This is much sharper than testing
   Conjecture 2 directly, because it mentions neither failure nor adjustment
   sets. With Theorem 1 (failure is upward-closed) it implies Conjecture 2.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from bkrobust.core.spacelib import Space, distances_from
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.search.space_fixed import knowledge_of


def join_candidate(cpdag: MPDAG, g: MPDAG, h: MPDAG) -> MPDAG | None:
    """``Meek(cpdag, K_g intersect K_h)`` -- the candidate least upper bound.

    Keeping only the orientations both graphs assert gives a graph at least as
    permissive as either, so it is an upper bound. Leastness is a separate
    question; see :func:`is_least_upper_bound`.
    """
    shared = sorted(set(knowledge_of(cpdag, g)) & set(knowledge_of(cpdag, h)))
    return apply_orientations(cpdag, shared)


def is_upper_bound(space: Space, w: MPDAG, g: MPDAG, h: MPDAG) -> bool:
    """Whether ``w`` lies above both ``g`` and ``h`` in model inclusion."""
    return space.reps[g] <= space.reps[w] and space.reps[h] <= space.reps[w]


def is_least_upper_bound(space: Space, w: MPDAG, g: MPDAG, h: MPDAG) -> bool:
    """Whether ``w`` is an upper bound of ``g``, ``h`` below every other one."""
    if not is_upper_bound(space, w, g, h):
        return False
    for u in space.elements:
        if u is w:
            continue
        if is_upper_bound(space, u, g, h) and not (space.reps[w] <= space.reps[u]):
            return False
    return True


@dataclass
class JoinReport:
    """Whether the space is a join-semilattice, and how the candidate performs.

    Attributes:
        n_pairs: Ordered-distinct pairs examined.
        n_join_exists: Pairs with a least upper bound in the space.
        n_candidate_is_join: Pairs where the intersect-and-close construction is
            that least upper bound.
        n_candidate_not_upper: Pairs where the candidate is not even an upper
            bound (would indicate the construction is wrong, not merely loose).
        n_candidate_outside_space: Pairs where the candidate is not an element.
        failures: Small sample of offending pairs, for inspection.
    """

    n_pairs: int = 0
    n_join_exists: int = 0
    n_candidate_is_join: int = 0
    n_candidate_not_upper: int = 0
    n_candidate_outside_space: int = 0
    failures: list[dict[str, Any]] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        """JSON-friendly view."""
        return {
            "n_pairs": self.n_pairs,
            "n_join_exists": self.n_join_exists,
            "n_candidate_is_join": self.n_candidate_is_join,
            "n_candidate_not_upper": self.n_candidate_not_upper,
            "n_candidate_outside_space": self.n_candidate_outside_space,
            "frac_candidate_is_join": (
                self.n_candidate_is_join / self.n_pairs if self.n_pairs else None
            ),
            "failures": self.failures[:5],
        }


def check_joins(space: Space) -> JoinReport:
    """Test join existence and the candidate construction over every pair."""
    rep = JoinReport()
    members = {g.edge_string(): g for g in space.elements}
    for i, g in enumerate(space.elements):
        for h in space.elements[i + 1 :]:
            rep.n_pairs += 1
            uppers = [u for u in space.elements if is_upper_bound(space, u, g, h)]
            least = [u for u in uppers if all(space.reps[u] <= space.reps[v] for v in uppers)]
            if least:
                rep.n_join_exists += 1
            cand = join_candidate(space.cpdag, g, h)
            if cand is None or cand.edge_string() not in members:
                rep.n_candidate_outside_space += 1
                if len(rep.failures) < 5:
                    rep.failures.append(
                        {
                            "kind": "candidate_outside_space",
                            "g": g.edge_string(),
                            "h": h.edge_string(),
                            "candidate": None if cand is None else cand.edge_string(),
                        }
                    )
                continue
            cand = members[cand.edge_string()]
            if not is_upper_bound(space, cand, g, h):
                rep.n_candidate_not_upper += 1
                if len(rep.failures) < 5:
                    rep.failures.append(
                        {
                            "kind": "candidate_not_upper",
                            "g": g.edge_string(),
                            "h": h.edge_string(),
                            "candidate": cand.edge_string(),
                        }
                    )
            elif is_least_upper_bound(space, cand, g, h):
                rep.n_candidate_is_join += 1
            elif len(rep.failures) < 5:
                rep.failures.append(
                    {
                        "kind": "candidate_not_least",
                        "g": g.edge_string(),
                        "h": h.edge_string(),
                        "candidate": cand.edge_string(),
                    }
                )
    return rep


@dataclass
class LemmaLReport:
    """Slack in Lemma L: ``d(G0, G) - d(G0, G v G0)``.

    Attributes:
        n_tested: ``(G0, G)`` pairs examined.
        n_violations: Pairs where the join is STRICTLY farther than ``G``.
        slack_hist: Distribution of the slack. Zero slack means the bound is
            tight -- the conjecture holds but with no margin, which is a very
            different state from holding comfortably.
        n_join_missing: Pairs where the join was unavailable, so untestable.
        violations: Sample of offending pairs.
    """

    n_tested: int = 0
    n_violations: int = 0
    slack_hist: dict[int, int] = field(default_factory=dict)
    n_join_missing: int = 0
    violations: list[dict[str, Any]] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        """JSON-friendly view."""
        total = sum(self.slack_hist.values())
        return {
            "n_tested": self.n_tested,
            "n_violations": self.n_violations,
            "n_join_missing": self.n_join_missing,
            "slack_hist": {str(k): v for k, v in sorted(self.slack_hist.items())},
            "frac_zero_slack": (self.slack_hist.get(0, 0) / total) if total else None,
            "violations": self.violations[:5],
        }


def check_lemma_l(space: Space, report: LemmaLReport | None = None) -> LemmaLReport:
    """Test ``d(G0, G v G0) <= d(G0, G)`` over every base point and every element.

    Args:
        space: The enumerated space (use the CORRECTED builder).
        report: Optional report to accumulate into, for sweeping many spaces.

    Returns:
        The accumulated :class:`LemmaLReport`, including the slack distribution.
    """
    rep = report if report is not None else LemmaLReport()
    members = {g.edge_string(): g for g in space.elements}
    for g0 in space.elements:
        dists = distances_from(space, g0)
        for g in space.elements:
            if g is g0:
                continue
            d_g = dists.get(g)
            if d_g is None:
                continue
            cand = join_candidate(space.cpdag, g, g0)
            if cand is None or cand.edge_string() not in members:
                rep.n_join_missing += 1
                continue
            j = members[cand.edge_string()]
            d_j = dists.get(j)
            if d_j is None:
                rep.n_join_missing += 1
                continue
            rep.n_tested += 1
            slack = d_g - d_j
            rep.slack_hist[slack] = rep.slack_hist.get(slack, 0) + 1
            if slack < 0:
                rep.n_violations += 1
                if len(rep.violations) < 5:
                    rep.violations.append(
                        {
                            "cpdag": space.cpdag.edge_string(),
                            "g0": g0.edge_string(),
                            "g": g.edge_string(),
                            "join": j.edge_string(),
                            "d_g": d_g,
                            "d_join": d_j,
                        }
                    )
    return rep


def chain_length_spread(space: Space) -> dict[str, Any]:
    """Non-gradedness: the spread of maximal-chain lengths between the same endpoints.

    Computes, for the top element (the CPDAG) and every minimal element, the
    shortest and longest saturated descending chain. A spread greater than zero
    refutes the Jordan-Dedekind chain condition and hence upper semimodularity.

    Returns:
        ``max_spread`` over endpoint pairs, the number of pairs with a positive
        spread, and a witness.
    """
    lower: dict[MPDAG, list[MPDAG]] = {g: [] for g in space.elements}
    for lo, hi in space.covers:
        lower[hi].append(lo)

    short: dict[MPDAG, int] = {}
    long: dict[MPDAG, int] = {}

    order = sorted(space.elements, key=lambda g: len(space.reps[g]))
    for g in order:
        downs = lower[g]
        if not downs:
            short[g] = 0
            long[g] = 0
        else:
            short[g] = 1 + min(short[d] for d in downs)
            long[g] = 1 + max(long[d] for d in downs)

    top = space.cpdag
    spread = long.get(top, 0) - short.get(top, 0)
    positives = sum(1 for g in space.elements if long[g] - short[g] > 0)
    return {
        "cpdag": top.edge_string(),
        "n_elements": len(space.elements),
        "shortest_maximal_chain_from_top": short.get(top, 0),
        "longest_maximal_chain_from_top": long.get(top, 0),
        "spread": spread,
        "n_elements_with_positive_spread": positives,
        "graded": spread == 0 and positives == 0,
    }
