r"""Relevant-restricted mean-bias profile mu_rel(d), against the reviewer's
dilution critique of Section 7's uniform mean.

``experiments/mean_table_fullspace.py`` computes

    mu(d) = mean { B(G) : G in G_Chat, d(G0, G) = d }

averaging uniformly over the *whole* enumerated space at each BFS shell. The
reviewer's point: when the CPDAG's undirected part splits into several chain
components and some of them cannot affect the query (X, Y, Z, or the possible
ancestors of X/Y), every orientation of an irrelevant component is a "free"
extra state at whatever distance it sits at, and it dilutes the shell average
even though it contributes nothing to the query's own worst-case bias. A
disconnected pair of nodes at distance 1 from an otherwise-oriented G0, for
instance, divides mu(1) by (q+1) for q such irrelevant pairs.

This script does not re-derive B from scratch. It reuses the exact same
per-instance objects mean_table_fullspace.py builds (same panel, same
--condition knowledge, same seeded SEM via final_table.derive_seed, same
CPDAG parse, same build_space_correct enumeration, same bias_at(...,
method="enumerate")) so that B(G) values here are bit-identical to the
committed results/mean_table_fullspace/ run for every element they share.
Nothing about the SEM or the space construction is changed; only the
*averaging convention* is.

Relevant-component rule (conservative, stated once here and used throughout):
an undirected chain component of Chat is RELEVANT to a query (x, y, z) iff it

    (a) contains x, or
    (b) contains a node of z union {y}, or
    (c) contains a node that is a possible ancestor (in Chat) of x or of y,

where "possible ancestor of v" means: reachable from v backwards along a
possibly-directed path, i.e. there is a sequence v = u_0, u_1, ..., u_k = w
with, for each step, either u_i -> u_{i+1} is NOT present in the wrong
direction -- concretely w is a possible ancestor of v iff there is a path
w = p_0 - p_1 - ... - p_k = v in Chat such that every edge p_i, p_{i+1} is
either undirected or directed p_i -> p_{i+1} (never p_{i+1} -> p_i). This is
the standard "possDe/possAn" notion from the PC/FCI literature, restricted to
a CPDAG.

All other components are IRRELEVANT. Every element of the enumerated space is
grouped by its restriction to the edges of the relevant components (its
"relevant signature"); Task 1 checks, directly on the data, whether B is
constant within every group (the invariance claim: B(G) depends only on the
component containing X, because tau_D = beta(pa_D(X); Sigma) depends only on
X's parents, which are only ever edges incident to X, all of which lie in X's
own undirected component of Chat, and theta_Z is fixed per instance).

mu_rel(d) is then the mean of B over the SLICE of the space that agrees with
G0 on every irrelevant edge (equivalently: the product-space factor for the
relevant components alone), at the same BFS distance d used by the full-space
script (Cor. ties BFS hop count to |K_G Delta K_H|; restricting to a slice
whose irrelevant edges never move makes the restricted distance and the full
distance coincide on that slice by construction).

Usage::

    PYTHONPATH=src .venv/bin/python experiments/mean_table_relevant_v2.py [options]
"""

from __future__ import annotations

import argparse
import json
import signal
import statistics
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.core.spacelib import distances_from  # noqa: E402
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag, random_sem  # noqa: E402
from bkrobust.demo.graph import MPDAG, undirected_components  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.epsilon.bias import bias_at, make_context  # noqa: E402
from bkrobust.search.space_fixed import build_space_correct  # noqa: E402

import final_table as ft  # noqa: E402

DEFAULT_EPSILONS = "0.02,0.10,0.25,0.40,0.60,0.90"
#: Source of the committed full-space run; only its (network, x, y, status,
#: r_val, beta_top) columns are read, to reproduce the exact 39-query subset
#: (status == "ok", r_val finite, beta_top > 0) without re-deriving the filter.
FULLSPACE_INSTANCES = ROOT / "results/mean_table_fullspace/instances.jsonl"


class Timeout(Exception):
    """Raised when one network's space exceeds its build budget."""


def _alarm(_signum: int, _frame: Any) -> None:
    raise Timeout()


# --- relevance -------------------------------------------------------------


def possible_ancestors(cpdag: MPDAG, v: str) -> set[str]:
    """Nodes with a possibly-directed path into ``v`` in ``cpdag`` (exclusive of ``v``).

    ``w`` is a possible ancestor of ``v`` iff there is a path
    ``w = p_0 - p_1 - ... - p_k = v`` where every step is either an undirected
    edge or a directed edge ``p_i -> p_{i+1}`` (never the reverse). Computed by
    reverse BFS from ``v`` over the "allowed forward" relation.

    Args:
        cpdag: The CPDAG.
        v: The node to find possible ancestors of.

    Returns:
        The set of possible ancestors, excluding ``v`` itself.
    """
    # allowed(a, b): a step a -> b is legal on a possibly-directed path, i.e.
    # cpdag has a -> b directed, or a - b undirected.
    def forward_neighbours(a: str) -> set[str]:
        out = set(cpdag.children(a)) | set(cpdag.neighbors(a))
        return out

    # Reverse search: w is a possible ancestor of v iff v is reachable from w
    # via forward_neighbours. Equivalently, search backwards from v using the
    # inverse relation: b is a predecessor of a if a in forward_neighbours(b).
    seen: set[str] = set()
    stack = [n for n in cpdag.nodes if n != v and v in forward_neighbours(n)]
    # The above one-hop seed is expensive to repeat; instead do a proper
    # reverse BFS using an explicit inverse adjacency.
    inverse: dict[str, set[str]] = {n: set() for n in cpdag.nodes}
    for a in cpdag.nodes:
        for b in forward_neighbours(a):
            inverse[b].add(a)
    seen = set()
    stack = list(inverse[v])
    while stack:
        cur = stack.pop()
        if cur in seen or cur == v:
            continue
        seen.add(cur)
        stack.extend(inverse[cur])
    seen.discard(v)
    return seen


def relevant_components(
    cpdag: MPDAG, x: str, y: str, z: frozenset[str], rule: str = "conservative"
) -> tuple[list[frozenset[str]], list[frozenset[str]]]:
    """Split ``undirected_components(cpdag)`` into (relevant, irrelevant).

    A component is relevant iff it contains ``x``; a node of ``z | {y}``; or a
    possible ancestor (in ``cpdag``) of ``x`` or of ``y``.

    Args:
        cpdag: The CPDAG.
        x: Treatment.
        y: Outcome.
        z: The adjustment set.

    Returns:
        ``(relevant, irrelevant)`` lists of components (node frozensets).
    """
    comps = undirected_components(cpdag)
    if rule == "x_component":
        # B(G) depends only on Pa_D(X) over D in [G]; CPDAG chain components
        # are oriented independently, so only X's own component can move it.
        relevant = [c for c in comps if x in c]
        irrelevant = [c for c in comps if x not in c]
        return relevant, irrelevant
    anc_x = possible_ancestors(cpdag, x)
    anc_y = possible_ancestors(cpdag, y)
    must_hit = {x} | set(z) | {y} | anc_x | anc_y
    relevant = [c for c in comps if c & must_hit]
    irrelevant = [c for c in comps if not (c & must_hit)]
    return relevant, irrelevant


def component_edges(comps: list[frozenset[str]]) -> set[tuple[str, str]]:
    """Undirected-edge endpoints (unordered, as stored) covered by a set of components."""
    nodes: set[str] = set()
    for c in comps:
        nodes |= c
    return nodes


# --- profiles ----------------------------------------------------------


def edges_within(cpdag: MPDAG, node_set: set[str]) -> set[tuple[str, str]]:
    """The undirected edges of ``cpdag`` with both endpoints in ``node_set``."""
    return {(a, b) for a, b in cpdag.undirected_edges if a in node_set and b in node_set}


def agrees_outside(g: MPDAG, g0: MPDAG, irrelevant_edges: set[tuple[str, str]]) -> bool:
    """Whether ``g`` orients every edge in ``irrelevant_edges`` exactly as ``g0`` does."""
    for a, b in irrelevant_edges:
        if g0.is_undirected_edge(a, b):
            if not g.is_undirected_edge(a, b):
                return False
        elif g0.is_directed_edge(a, b):
            if not g.is_directed_edge(a, b):
                return False
        else:  # g0 has b -> a
            if not g.is_directed_edge(b, a):
                return False
    return True


def relevant_signature(
    g: MPDAG, relevant_edges: set[tuple[str, str]]
) -> tuple[tuple[str, str, str], ...]:
    """A hashable restriction of ``g`` to the relevant edges, for grouping."""
    out = []
    for a, b in sorted(relevant_edges):
        if g.is_undirected_edge(a, b):
            out.append((a, b, "u"))
        elif g.is_directed_edge(a, b):
            out.append((a, b, "fwd"))
        else:
            out.append((a, b, "bwd"))
    return tuple(out)


def profile_from(
    b: dict[MPDAG, float],
    reachable: list[MPDAG],
    dists: dict[MPDAG, int],
) -> list[dict[str, Any]]:
    """Shell-by-shell mu/beta profile over an arbitrary subset of elements."""
    out: list[dict[str, Any]] = []
    if not reachable:
        return out
    running_max = 0.0
    for d in range(max(dists[g] for g in reachable) + 1):
        sphere = [b[g] for g in reachable if dists[g] == d]
        if not sphere:
            continue
        ball = [b[g] for g in reachable if dists[g] <= d]
        running_max = max(running_max, max(sphere))
        nz = [v for v in sphere if v > 0]
        out.append({
            "d": d,
            "n_sphere": len(sphere),
            "n_ball": len(ball),
            "mu": statistics.fmean(sphere),
            "mu_ball": statistics.fmean(ball),
            "sphere_max": max(sphere),
            "beta": running_max,
            "frac_nonzero": len(nz) / len(sphere),
            "median": statistics.median(sphere),
        })
    return out


def r_onset(profile: list[dict[str, Any]], eps: float, key: str = "mu") -> int:
    for rec in profile:
        if rec[key] > eps:
            return rec["d"]
    return UNREACHED


def persists(profile: list[dict[str, Any]], eps: float, key: str = "mu") -> bool | None:
    hit = r_onset(profile, eps, key)
    if hit == UNREACHED:
        return None
    return all(rec[key] > eps for rec in profile if rec["d"] >= hit)


def load_target_queries() -> set[tuple[str, str, str]]:
    """The committed 39-query subset: status ok, r_val finite, beta_top > 0."""
    out = set()
    for line in FULLSPACE_INSTANCES.open():
        r = json.loads(line)
        if r.get("status") != "ok":
            continue
        if r.get("r_val") in (None, "inf"):
            continue
        if not (r.get("beta_top", 0) > 0):
            continue
        out.add((r["network"], r["x"], r["y"]))
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--condition", default="D_LLM")
    p.add_argument("--panel", default="results/final_table_eps_gt1/instances.jsonl")
    p.add_argument("--epsilons", default=DEFAULT_EPSILONS)
    p.add_argument("--seed", type=int, default=20260917)
    p.add_argument("--budget", type=int, default=900)
    p.add_argument("--max-undirected", type=int, default=10)
    p.add_argument("--out", default="results/mean_table_relevant_v2")
    p.add_argument("--networks", default="", help="comma list; default all target networks")
    p.add_argument("--rule", default="conservative", choices=["conservative", "x_component"],
                   help="which chain components count as relevant to the query")
    args = p.parse_args()

    signal.signal(signal.SIGALRM, _alarm)
    epsilons = [float(v) for v in args.epsilons.split(",")]
    panel = [json.loads(line) for line in (ROOT / args.panel).open()]
    targets = load_target_queries()
    net_filter = {v for v in args.networks.split(",") if v}
    informative = [
        r for r in panel
        if r["status"] == "ok" and (r["network"], r["x"], r["y"]) in targets
        and (not net_filter or r["network"] in net_filter)
    ]
    print(f"[relevant] {len(informative)}/{len(targets)} target queries found in panel", flush=True)

    knowledge = ft.load_knowledge(args.condition)
    want = {r["network"] for r in informative}
    parsed = ft.load_networks(names=want, skip=set(), max_nodes=0)

    out_dir = ROOT / args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, Any]] = []
    t_all = time.perf_counter()

    by_net: dict[str, list[dict[str, Any]]] = {}
    for r in informative:
        by_net.setdefault(r["network"], []).append(r)

    for net in sorted(by_net):
        rows = by_net[net]
        if net not in parsed or net not in knowledge:
            for r in rows:
                records.append({**{k: r[k] for k in ("network", "x", "y", "r_val")},
                                "status": "unavailable"})
            continue
        cpdag, dag = parsed[net]["cpdag"], parsed[net]["dag"]
        k_und = len(cpdag.undirected_edges)
        if k_und > args.max_undirected:
            for r in rows:
                records.append({**{q: r[q] for q in ("network", "x", "y", "r_val")},
                                "status": "space_too_large", "k_undirected": k_und})
            continue
        signal.alarm(args.budget)
        t0 = time.perf_counter()
        try:
            space = build_space_correct(cpdag)
            signal.alarm(0)
        except Timeout:
            signal.alarm(0)
            for r in rows:
                records.append({**{q: r[q] for q in ("network", "x", "y", "r_val")},
                                "status": "space_timeout", "k_undirected": k_und})
            continue
        build_s = time.perf_counter() - t0
        print(f"[relevant] {net}: |space|={len(space):,} k={k_und} built in {build_s:.1f}s "
              f"({len(rows)} instances)", flush=True)

        for r in rows:
            x, y = r["x"], r["y"]
            rec: dict[str, Any] = {
                "network": net, "x": x, "y": y, "r_val": r["r_val"],
                "n_knowledge": r["n_knowledge"], "k_undirected": k_und,
                "space_size": len(space), "space_build_seconds": round(build_s, 2),
            }
            g0 = apply_orientations(cpdag, list(knowledge[net]))
            if g0 is None or g0 not in space.neighbours:
                rec["status"] = "g0_not_in_space"
                records.append(rec); continue
            z = optimal_adjustment_set_mpdag(g0, x, y)
            if z is None:
                rec["status"] = "degenerate"
                records.append(rec); continue
            z = frozenset(z)
            ctx = make_context(
                random_sem(dag, np.random.default_rng(ft.derive_seed(args.seed, net, x, y))),
                cpdag, x, y, z)
            scale = abs(ctx.theta_z)
            if scale < 1e-12:
                rec["status"] = "zero_estimate"
                records.append(rec); continue

            t1 = time.perf_counter()
            dists = distances_from(space, g0)
            reachable = [g for g in space.elements if dists.get(g) is not None]

            # Component split.
            relevant, irrelevant = relevant_components(cpdag, x, y, z, rule=args.rule)
            relevant_nodes = component_edges(relevant)
            irrelevant_nodes = component_edges(irrelevant)
            relevant_edges = edges_within(cpdag, relevant_nodes)
            irrelevant_edges = edges_within(cpdag, irrelevant_nodes)
            assert relevant_edges | irrelevant_edges == set(cpdag.undirected_edges), \
                f"{net} {x}->{y}: component split does not partition the undirected edges"

            # B(G) for every reachable element -- identical evaluator/method to
            # mean_table_fullspace.py, so bit-identical values on overlap.
            b = {g: bias_at(ctx, g, method="enumerate").worst / scale for g in reachable}

            # Task 1: is B constant within each relevant-signature group?
            groups: dict[Any, list[float]] = {}
            for g in reachable:
                sig = relevant_signature(g, relevant_edges)
                groups.setdefault(sig, []).append(b[g])
            spreads = [max(vals) - min(vals) for vals in groups.values() if len(vals) > 1]
            invariance_holds = all(s < 1e-9 for s in spreads)
            max_spread = max(spreads) if spreads else 0.0

            # Full-space profile (identical to mean_table_fullspace.py; used
            # as the "before" column for comparison in this table).
            prof_full = profile_from(b, reachable, dists)

            # Relevant-restricted profile: the slice that agrees with G0 on
            # every irrelevant edge.
            slice_elems = [g for g in reachable if agrees_outside(g, g0, irrelevant_edges)]
            prof_rel = profile_from(b, slice_elems, dists)

            mus_full = [pp["mu"] for pp in prof_full]
            mus_rel = [pp["mu"] for pp in prof_rel]

            rec.update({
                "status": "ok",
                "theta_z": ctx.theta_z,
                "z": sorted(z),
                "n_undirected_components": len(undirected_components(cpdag)),
                "n_relevant_components": len(relevant),
                "n_irrelevant_components": len(irrelevant),
                "irrelevant_component_sizes": sorted(len(c) for c in irrelevant),
                "n_irrelevant_undirected_edges": len(irrelevant_edges),
                "n_relevant_undirected_edges": len(relevant_edges),
                "invariance_holds": invariance_holds,
                "max_signature_spread": max_spread,
                "n_signature_groups": len(groups),
                "max_distance_full": max(pp["d"] for pp in prof_full),
                "max_distance_rel": max(pp["d"] for pp in prof_rel) if prof_rel else None,
                "n_reachable_full": len(reachable),
                "n_reachable_rel": len(slice_elems),
                "mu_peak_full": max(mus_full) if mus_full else 0.0,
                "mu_peak_rel": max(mus_rel) if mus_rel else 0.0,
                "r_mu_full": {str(e): r_onset(prof_full, e) for e in epsilons},
                "r_mu_rel": {str(e): r_onset(prof_rel, e) for e in epsilons},
                "r_mu_full_persists": {str(e): persists(prof_full, e) for e in epsilons},
                "r_mu_rel_persists": {str(e): persists(prof_rel, e) for e in epsilons},
                "profile_full": prof_full,
                "profile_rel": prof_rel,
                "seconds": round(time.perf_counter() - t1, 2),
            })
            records.append(rec)
            print(f"[relevant]   {net} {x}->{y} r_val={r['r_val']} "
                  f"irrel_comps={len(irrelevant)} inv={invariance_holds} "
                  f"r_mu_full={[rec['r_mu_full'][str(e)] for e in epsilons]} "
                  f"r_mu_rel={[rec['r_mu_rel'][str(e)] for e in epsilons]}", flush=True)

    with (out_dir / "instances.jsonl").open("w") as f:
        for rec in records:
            f.write(json.dumps(rec) + "\n")

    ok = [r for r in records if r.get("status") == "ok"]
    n_changed = sum(
        1 for r in ok
        if any(r["r_mu_full"][str(e)] != r["r_mu_rel"][str(e)] for e in epsilons)
    )
    n_bad_invariance = sum(1 for r in ok if not r["invariance_holds"])
    summary = {
        "args": vars(args), "epsilons": epsilons,
        "n_target_queries": len(targets),
        "n_ok": len(ok),
        "n_invariance_violations": n_bad_invariance,
        "n_queries_with_changed_r_mu": n_changed,
        "n_queries_with_irrelevant_components": sum(
            1 for r in ok if r["n_irrelevant_components"] > 0),
        "seconds": round(time.perf_counter() - t_all, 1),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\n[relevant] {summary['n_ok']}/{summary['n_target_queries']} target queries "
          f"processed; invariance violations: {n_bad_invariance}; "
          f"queries with >=1 irrelevant component: {summary['n_queries_with_irrelevant_components']}; "
          f"r_mu changed for {n_changed}; {summary['seconds']}s")


if __name__ == "__main__":
    main()
