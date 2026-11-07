"""Adjustment identification distance (Henckel, Wurtzen, Weichwald 2024), oset variant.

Our implementation reads validity and amenability from the submitted GAC module
(bkrobust.gac.mpdag_level), so the AID and the radius are decided by the same
predicate. gadjid, the authors' published implementation, is called beside it
when installed; the two are compared, never merged.

Mistake rules (also in PREREG.md): the guess G answers, for an ordered pair (t, v),
  "zero"  if v is not a possible descendant of t in G,
  "noid"  if G is not amenable for (t, v) or its optimal set differs across [G],
  a set   O(t, v, G) otherwise.
Against a reference T: "zero" is wrong if v is a possible descendant of t in T;
"noid" is wrong if T is amenable for (t, v); a set is wrong if it is not
GAC-valid for (t, v) in T.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from bkrobust.demo.evaluate import optimal_adjustment_set_dag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import enumerate_dag_extensions
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag, why_invalid_gac
from bkrobust.mpdag_criterion.paths import possible_descendants


def guess_answers(g: MPDAG) -> dict[tuple[str, str], Any]:
    """What the guess graph says for every ordered pair: "zero", "noid" or a frozenset."""
    exts = enumerate_dag_extensions(g)
    if not exts:
        raise ValueError("guess graph has no DAG extension")
    nodes = list(g.nodes)
    out: dict[tuple[str, str], Any] = {}
    for t in nodes:
        possde = possible_descendants(g, t) - {t}
        for v in nodes:
            if v == t:
                continue
            if v not in possde:
                out[(t, v)] = "zero"
                continue
            # The reason "not_amenable" is decided before z is looked at, so an
            # empty z reads amenability alone.
            if why_invalid_gac(g, t, v, ()) == "not_amenable":
                out[(t, v)] = "noid"
                continue
            sets = {frozenset(optimal_adjustment_set_dag(d, t, v)) for d in exts}
            out[(t, v)] = sets.pop() if len(sets) == 1 else "noid_oset_undefined"
    return out


def mistakes(answers: dict[tuple[str, str], Any], ref: MPDAG) -> dict[tuple[str, str], int]:
    """1 where the guess's answer is wrong in the reference graph, else 0."""
    possde_cache: dict[str, set[str]] = {}
    out: dict[tuple[str, str], int] = {}
    for (t, v), ans in answers.items():
        if t not in possde_cache:
            possde_cache[t] = possible_descendants(ref, t) - {t}
        if ans == "zero":
            out[(t, v)] = int(v in possde_cache[t])
        elif isinstance(ans, str):  # noid, including an undefined optimal set
            out[(t, v)] = int(why_invalid_gac(ref, t, v, ()) != "not_amenable")
        else:
            out[(t, v)] = int(not is_gac_valid_mpdag(ref, t, v, ans))
    return out


def aid_scopes(answers, ref: MPDAG, x: str, y: str) -> dict[str, Any]:
    """Graph, row (treatment x) and pair (x, y) scopes, counts and normalised values."""
    m = mistakes(answers, ref)
    n = len(ref.nodes)
    graph = sum(m.values())
    row = sum(val for (t, _v), val in m.items() if t == x)
    return {
        "graph_count": graph,
        "graph_norm": graph / (n * (n - 1)),
        "row_count": row,
        "row_norm": row / (n - 1),
        "pair": m[(x, y)],
    }


def answer_summary(answers) -> dict[str, int]:
    kinds = {"zero": 0, "noid": 0, "noid_oset_undefined": 0, "set": 0}
    for ans in answers.values():
        kinds[ans if isinstance(ans, str) else "set"] += 1
    return kinds


# --------------------------------------------------------------------------- gadjid

def _matrix(g: MPDAG, nodes: list[str], symmetric_undirected: bool) -> np.ndarray:
    """Adjacency, row -> column; an undirected edge is 2 (upper triangle, or both)."""
    idx = {v: i for i, v in enumerate(nodes)}
    m = np.zeros((len(nodes), len(nodes)), dtype=np.int8)
    for a, b in g.directed_edges:
        m[idx[a], idx[b]] = 1
    for a, b in g.undirected_edges:
        i, j = sorted((idx[a], idx[b]))
        m[i, j] = 2
        if symmetric_undirected:
            m[j, i] = 2
    return m


def gadjid_probe() -> dict[str, Any]:
    """Import gadjid and find which undirected-edge encoding it accepts.

    The encoding of a CPDAG's undirected edges could not be read from the
    package on the authoring machine, so it is decided here on a known case:
    a three-node undirected chain compared with itself must give zero mistakes.
    """
    try:
        import gadjid  # noqa: F401
    except Exception as exc:  # noqa: BLE001
        return {"available": False, "error": f"import: {exc!r}"}
    import gadjid

    chain = MPDAG(["a", "b", "c"], undirected=[("a", "b"), ("b", "c")])
    errors = []
    for sym in (False, True):
        try:
            m = _matrix(chain, ["a", "b", "c"], sym)
            res = gadjid.oset_aid(m, m, edge_direction="from row to column")
            if int(res[1]) == 0:
                return {"available": True, "symmetric_undirected": sym,
                        "version": getattr(gadjid, "__version__", "unknown")}
            errors.append(f"sym={sym}: self-distance {res}")
        except Exception as exc:  # noqa: BLE001
            errors.append(f"sym={sym}: {exc!r}")
    return {"available": False, "error": "; ".join(errors)}


def gadjid_aids(ref: MPDAG, guess: MPDAG, symmetric_undirected: bool) -> dict[str, Any]:
    """gadjid's three AIDs from guess to ref; an error string where it refuses the input."""
    import gadjid

    nodes = list(ref.nodes)
    a = _matrix(ref, nodes, symmetric_undirected)
    b = _matrix(guess, nodes, symmetric_undirected)
    out: dict[str, Any] = {}
    for name in ("oset_aid", "parent_aid", "ancestor_aid"):
        try:
            norm, count = getattr(gadjid, name)(a, b, edge_direction="from row to column")
            out[f"{name}_norm"] = float(norm)
            out[f"{name}_count"] = int(count)
        except Exception as exc:  # noqa: BLE001
            out[f"{name}_error"] = repr(exc)[:300]
    return out
