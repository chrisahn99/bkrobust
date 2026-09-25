r"""The X-component mean-bias profile, computed by enumerating ONLY X's own
undirected chain component, for the 18 instances that a full-space build could
not reach.

``experiments/mean_table_relevant_v2.py --rule x_component`` establishes
(and this file's own validation step re-derives, see ``--validate``) that

    B(G) depends only on the orientations inside X's own undirected chain
    component of Chat; the other chain components are oriented independently
    by Meek's rules, so a state that agrees with G0 outside X's component has
    the same B as G0 restricted to that component.

That script computes this by building the WHOLE ``3^k`` space (``k`` =
Chat's total undirected-edge count) and then filtering down to the slice that
agrees with G0 outside X's component. For 18 queries on 8 networks (
Kampen_2014, andes, child, ecoli70, magic-irri, magic-niab, paths, win95pts)
``k`` is too large (10-25) for that whole-space build to finish, even though
X's *own* component -- the only part that moves B -- is typically far
smaller.

This script enumerates the slice DIRECTLY instead: fix every edge outside
X's component to G0's own orientation (including "still undirected", when
G0 leaves it that way -- the slice requires *exact* agreement, not just
consistency), range freely over the ``3^{k_x}`` assignments of X's own
``k_x`` undirected edges, Meek-close each against Chat, and deduplicate. This
needs no covering-pairs computation (no ``O(|space|^2)`` step) and no
BFS: by Cor. cor:distance(a) of the paper, ``d(G, H) = |K_G Delta K_H|`` is a
closed form, so shell membership is read straight off the two states'
closure-orientation sets rather than off a neighbour graph.

**Validation** (``--validate``, on by default before the 18 run): the same
slice enumeration, applied to the 39 queries
``results/mean_table_relevant_v2/x_component/instances.jsonl`` already
covers via the full-space filter, must reproduce every field of that file's
"_rel" side exactly (profile, r_mu, r_mu_persists, mu_peak, max_distance,
n_reachable) -- bit for bit, since both are the same slice computed two
different ways.

**The bias evaluator.** ``bias_at(ctx, g, method="enumerate")`` enumerates
ALL DAG extensions of the state it is given -- exponential in every
undirected edge that state still carries, not just the ones at X. A slice
element inherits G0's own leftover undirected edges outside X's component
(whichever ones G0 itself did not resolve), so evaluating B naively on the
slice element itself can be exponential in irrelevant structure even though
B does not depend on it. This script instead evaluates B on a graph where
every OUTSIDE edge is resolved to the ground-truth DAG's own direction (a
valid, Meek-consistent choice, since the DAG is by construction an extension
of Chat) before calling ``bias_at(..., method="enumerate")`` -- provably the
same B, by the invariance fact above, but exponential only in ``k_x``. This
substitution is itself part of what ``--validate`` checks: on the 39 queries
the substituted evaluation must equal both the committed mu values and a
direct (unsubstituted) re-evaluation.

Usage::

    PYTHONPATH=src .venv/bin/python experiments/mean_table_xcomp_extend_v2.py [options]
"""

from __future__ import annotations

import argparse
import itertools
import json
import signal
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag, random_sem  # noqa: E402
from bkrobust.demo.graph import MPDAG, canon, undirected_components  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.epsilon.bias import bias_at, make_context  # noqa: E402
from bkrobust.search.space_fixed import knowledge_of  # noqa: E402

import final_table as ft  # noqa: E402
import mean_table_relevant_v2 as rel  # noqa: E402

DEFAULT_EPSILONS = "0.02,0.10,0.25,0.40,0.60,0.90"
FULLSPACE_INSTANCES = ROOT / "results/mean_table_fullspace/instances.jsonl"
XCOMP_INSTANCES = ROOT / "results/mean_table_relevant_v2/x_component/instances.jsonl"

Edge = tuple[str, str]


class Timeout(Exception):
    """Raised when one instance's slice enumeration exceeds its budget."""


def _alarm(_signum: int, _frame: Any) -> None:
    raise Timeout()


# --- the 18 skipped instances ----------------------------------------------


def load_skipped_queries() -> list[dict[str, Any]]:
    """The instances mean_table_fullspace.py excluded as space_too_large/timeout."""
    out = []
    for line in FULLSPACE_INSTANCES.open():
        r = json.loads(line)
        if r.get("status") in ("space_too_large", "space_timeout"):
            out.append(r)
    return out


# --- direct slice enumeration -----------------------------------------------


def x_component_edges(cpdag: MPDAG, x: str) -> tuple[set[str], list[Edge]]:
    """X's own undirected chain component: its node set and undirected edges.

    Returns:
        ``(nodes, edges)``; both empty when ``x`` has no undirected edge in
        ``cpdag`` (its component is degenerate -- the slice is then just
        ``{G0}``).
    """
    for comp in undirected_components(cpdag):
        if x in comp:
            edges = sorted((a, b) for a, b in cpdag.undirected_edges if a in comp and b in comp)
            return set(comp), edges
    return set(), []


def outside_fixed_orientations(cpdag: MPDAG, g0: MPDAG, x_nodes: set[str]) -> list[Edge]:
    """G0's own orientations, restricted to edges with neither endpoint in X's component.

    Edges G0 leaves undirected are simply omitted (so, since they stay
    undirected in ``cpdag`` too, the slice enumeration below leaves them
    undirected as well -- exact agreement with G0, as the slice requires).
    """
    return sorted(
        (a, b) for a, b in knowledge_of(cpdag, g0) if a not in x_nodes or b not in x_nodes
    )


def resolve_outside_via_dag(
    cpdag: MPDAG, dag: MPDAG, x_nodes: set[str], outside_fixed: list[Edge]
) -> list[Edge]:
    """``outside_fixed`` (G0's own orientations), plus the DAG's direction for
    every OUTSIDE Chat edge G0 itself leaves undirected.

    G0's own choices must be kept verbatim: they need not agree with the
    ground-truth DAG (G0 is the analyst's, possibly imperfect, knowledge), so
    substituting the DAG's direction for an edge G0 *did* orient could
    directly contradict it and make the resulting graph inconsistent. Only
    the edges G0 leaves open are free to pick any consistent resolution --
    B(G) does not depend on outside orientation at all, so this is a real
    choice with no effect on the answer, made solely to remove these edges
    from ``bias_at(..., method="enumerate")``'s DAG-extension enumeration.
    """
    fixed_edges = {canon(a, b) for a, b in outside_fixed}
    ors: list[Edge] = list(outside_fixed)
    for a, b in cpdag.undirected_edges:
        if a in x_nodes and b in x_nodes:
            continue
        if canon(a, b) in fixed_edges:
            continue
        if dag.is_directed_edge(a, b):
            ors.append((a, b))
        elif dag.is_directed_edge(b, a):
            ors.append((b, a))
        else:
            raise ValueError(f"DAG does not resolve Chat edge {(a, b)}")
    return ors


def enumerate_slice(
    cpdag: MPDAG, x_edges: list[Edge], outside_fixed: list[Edge]
) -> dict[str, MPDAG]:
    """Every state agreeing with G0 outside X's component, keyed by edge string.

    Ranges over the ``3^{k_x}`` assignments of X's own undirected edges
    (0: leave undirected, 1: forward, 2: backward -- the same convention
    ``enumerate_space_correct`` uses), Meek-closes each against ``cpdag``
    together with the fixed outside orientations, and deduplicates.
    """
    seen: dict[str, MPDAG] = {}
    for assign in itertools.product((0, 1, 2), repeat=len(x_edges)):
        ors = list(outside_fixed)
        for (a, b), v in zip(x_edges, assign, strict=True):
            if v == 1:
                ors.append((a, b))
            elif v == 2:
                ors.append((b, a))
        g = apply_orientations(cpdag, ors)
        if g is not None:
            seen.setdefault(g.edge_string(), g)
    return seen


def slice_distance(cpdag: MPDAG, g: MPDAG, g0: MPDAG) -> int:
    """``d(G, H) = |K_G Delta K_H|`` (Cor. cor:distance(a)): the closed form."""
    kg = frozenset(knowledge_of(cpdag, g))
    kg0 = frozenset(knowledge_of(cpdag, g0))
    return len(kg ^ kg0)


def bias_worst(ctx: Any, cpdag: MPDAG, outside_resolution: list[Edge], g: MPDAG, scale: float) -> float:
    """B(g)/scale, evaluated on a graph with every outside edge DAG-resolved."""
    g_eval = apply_orientations(cpdag, outside_resolution + list(knowledge_of(cpdag, g)))
    assert g_eval is not None, "resolved graph must stay consistent"
    return bias_at(ctx, g_eval, method="enumerate").worst / scale


# --- profile construction (mirrors mean_table_relevant_v2.profile_from) ----


def build_record(
    net: str,
    x: str,
    y: str,
    r_val: Any,
    n_knowledge: Any,
    cpdag: MPDAG,
    dag: MPDAG,
    g0: MPDAG,
    ctx: Any,
    scale: float,
    epsilons: list[float],
    evaluate_direct: bool = False,
) -> dict[str, Any]:
    """One instance's slice profile, computed directly (no full-space build)."""
    rec: dict[str, Any] = {
        "network": net, "x": x, "y": y, "r_val": r_val, "n_knowledge": n_knowledge,
        "k_undirected": len(cpdag.undirected_edges),
    }
    x_nodes, x_edges = x_component_edges(cpdag, x)
    k_x = len(x_edges)
    rec["k_x"] = k_x
    rec["x_component_size"] = len(x_nodes)

    outside_fixed = outside_fixed_orientations(cpdag, g0, x_nodes)
    if k_x == 0:
        # X has no undirected edge in Chat: the slice is the single point G0.
        slice_elems = {g0.edge_string(): g0}
    else:
        slice_elems = enumerate_slice(cpdag, x_edges, outside_fixed)
    rec["slice_size"] = len(slice_elems)

    outside_resolution = resolve_outside_via_dag(cpdag, dag, x_nodes, outside_fixed)

    dists = {g: slice_distance(cpdag, g, g0) for g in slice_elems.values()}
    b = {g: bias_worst(ctx, cpdag, outside_resolution, g, scale) for g in slice_elems.values()}
    if evaluate_direct:
        # Differential check: evaluate without the DAG-resolution substitution
        # (only ever used on the small, already-enumerable 39 queries).
        b_direct = {g: bias_at(ctx, g, method="enumerate").worst / scale for g in slice_elems.values()}
        rec["_direct_matches"] = all(
            abs(b[g] - b_direct[g]) < 1e-9 for g in slice_elems.values()
        )

    reachable = list(slice_elems.values())
    prof = rel.profile_from(b, reachable, dists)
    mus = [pp["mu"] for pp in prof]

    k_abs = len(knowledge_of(cpdag, g0))
    # B(Chat): the top of the order, i.e. K = empty everywhere, including on
    # X's own component (left fully undirected, exactly as in cpdag). Any
    # consistent resolution of the OUTSIDE edges gives the same value by the
    # same invariance fact, so the DAG's own direction (needs no agreement
    # with G0 here -- this is not a slice element) is used freely.
    b_top_ctx_g = apply_orientations(cpdag, resolve_outside_via_dag(cpdag, dag, x_nodes, []))
    b_top = bias_at(ctx, b_top_ctx_g, method="enumerate").worst / scale

    rec.update({
        "status": "ok",
        "theta_z": ctx.theta_z,
        "z": sorted(ctx.z),
        "n_undirected_components": len(undirected_components(cpdag)),
        "max_distance_rel": max(pp["d"] for pp in prof) if prof else 0,
        "n_reachable_rel": len(reachable),
        "mu_peak_rel": max(mus) if mus else 0.0,
        "mu_peak_rel_at": prof[mus.index(max(mus))]["d"] if mus else None,
        "mu_monotone_rel": all(mus[i] >= mus[i - 1] - 1e-12 for i in range(1, len(mus))),
        "r_mu_rel": {str(e): rel.r_onset(prof, e) for e in epsilons},
        "r_mu_rel_persists": {str(e): rel.persists(prof, e) for e in epsilons},
        "profile_rel": prof,
        "k_abs": k_abs,
        "b_top": b_top,
        "has_finite_r_val": r_val not in (None, -1, "inf"),
        "has_nonzero_bias": b_top > 1e-9,
    })
    return rec


def compare_records(mine: dict[str, Any], theirs: dict[str, Any], epsilons: list[float]) -> list[str]:
    """Field-by-field mismatch report between a freshly computed row and the committed one."""
    problems = []
    if mine["max_distance_rel"] != theirs["max_distance_rel"]:
        problems.append(f"max_distance_rel {mine['max_distance_rel']} != {theirs['max_distance_rel']}")
    if mine["n_reachable_rel"] != theirs["n_reachable_rel"]:
        problems.append(f"n_reachable_rel {mine['n_reachable_rel']} != {theirs['n_reachable_rel']}")
    if abs(mine["mu_peak_rel"] - theirs["mu_peak_rel"]) > 1e-9:
        problems.append(f"mu_peak_rel {mine['mu_peak_rel']} != {theirs['mu_peak_rel']}")
    for e in epsilons:
        k = str(e)
        if mine["r_mu_rel"][k] != theirs["r_mu_rel"][k]:
            problems.append(f"r_mu_rel[{k}] {mine['r_mu_rel'][k]} != {theirs['r_mu_rel'][k]}")
        if mine["r_mu_rel_persists"][k] != theirs["r_mu_rel_persists"][k]:
            problems.append(
                f"r_mu_rel_persists[{k}] {mine['r_mu_rel_persists'][k]} != {theirs['r_mu_rel_persists'][k]}"
            )
    mine_prof = {p["d"]: p for p in mine["profile_rel"]}
    their_prof = {p["d"]: p for p in theirs["profile_rel"]}
    if set(mine_prof) != set(their_prof):
        problems.append(f"profile_rel distances differ: {sorted(mine_prof)} vs {sorted(their_prof)}")
    else:
        for d, mp in mine_prof.items():
            tp = their_prof[d]
            for key in ("n_sphere", "n_ball"):
                if mp[key] != tp[key]:
                    problems.append(f"d={d} {key} {mp[key]} != {tp[key]}")
            for key in ("mu", "mu_ball", "sphere_max", "beta", "frac_nonzero", "median"):
                if abs(mp[key] - tp[key]) > 1e-9:
                    problems.append(f"d={d} {key} {mp[key]} != {tp[key]}")
    return problems


def run_validation(args: argparse.Namespace, epsilons: list[float]) -> bool:
    """Reproduce the committed 39-query x_component table via direct slice enumeration."""
    committed = {}
    for line in XCOMP_INSTANCES.open():
        r = json.loads(line)
        if r.get("status") == "ok":
            committed[(r["network"], r["x"], r["y"])] = r
    print(f"[validate] {len(committed)} committed ok rows to reproduce", flush=True)

    knowledge = ft.load_knowledge(args.condition)
    want = {net for net, _, _ in committed}
    parsed = ft.load_networks(names=want, skip=set(), max_nodes=0)

    all_ok = True
    n_checked = 0
    for (net, x, y), theirs in sorted(committed.items()):
        if net not in parsed or net not in knowledge:
            print(f"[validate] SKIP {net} {x}->{y}: network unavailable", flush=True)
            continue
        cpdag, dag = parsed[net]["cpdag"], parsed[net]["dag"]
        g0 = apply_orientations(cpdag, list(knowledge[net]))
        if g0 is None:
            print(f"[validate] SKIP {net} {x}->{y}: g0 build failed", flush=True)
            continue
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None:
            print(f"[validate] SKIP {net} {x}->{y}: degenerate", flush=True)
            continue
        z = frozenset(z)
        ctx = make_context(
            random_sem(dag, np.random.default_rng(ft.derive_seed(args.seed, net, x, y))),
            cpdag, x, y, z)
        scale = abs(ctx.theta_z)
        if scale < 1e-12:
            print(f"[validate] SKIP {net} {x}->{y}: zero estimate", flush=True)
            continue

        mine = build_record(net, x, y, theirs["r_val"], theirs.get("n_knowledge"),
                             cpdag, dag, g0, ctx, scale, epsilons, evaluate_direct=True)
        n_checked += 1
        problems = compare_records(mine, theirs, epsilons)
        if mine.get("_direct_matches") is False:
            problems.append("DAG-resolved bias eval != direct bias_at on the real slice element")
        if problems:
            all_ok = False
            print(f"[validate] MISMATCH {net} {x}->{y}:", flush=True)
            for p in problems:
                print(f"    {p}", flush=True)
        else:
            print(f"[validate] ok {net} {x}->{y} (k_x={mine['k_x']}, slice={mine['slice_size']})",
                  flush=True)

    print(f"\n[validate] {n_checked}/{len(committed)} checked; "
          f"{'ALL MATCH' if all_ok else 'MISMATCHES FOUND'}", flush=True)
    return all_ok


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--condition", default="D_LLM")
    p.add_argument("--epsilons", default=DEFAULT_EPSILONS)
    p.add_argument("--seed", type=int, default=20260917)
    p.add_argument("--budget", type=int, default=900, help="seconds per instance")
    p.add_argument("--out", default="results/mean_table_relevant_v2/extended")
    p.add_argument("--skip-validate", action="store_true")
    args = p.parse_args()

    signal.signal(signal.SIGALRM, _alarm)
    epsilons = [float(v) for v in args.epsilons.split(",")]
    out_dir = ROOT / args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    if not args.skip_validate:
        ok = run_validation(args, epsilons)
        (out_dir / "validation.json").write_text(json.dumps({"all_match": ok}, indent=2) + "\n")
        if not ok:
            print("[main] validation FAILED; proceeding to the 18-instance run anyway, "
                  "but results should not be trusted until this is fixed", flush=True)

    skipped = load_skipped_queries()
    print(f"\n[xcomp18] {len(skipped)} previously-skipped instances to process", flush=True)
    knowledge = ft.load_knowledge(args.condition)
    want = {r["network"] for r in skipped}
    parsed = ft.load_networks(names=want, skip=set(), max_nodes=0)

    records: list[dict[str, Any]] = []
    t_all = time.perf_counter()
    for r in sorted(skipped, key=lambda r: (r["network"], r["x"], r["y"])):
        net, x, y = r["network"], r["x"], r["y"]
        if net not in parsed or net not in knowledge:
            records.append({"network": net, "x": x, "y": y, "status": "unavailable"})
            continue
        cpdag, dag = parsed[net]["cpdag"], parsed[net]["dag"]
        g0 = apply_orientations(cpdag, list(knowledge[net]))
        if g0 is None:
            records.append({"network": net, "x": x, "y": y, "status": "g0_build_failed"})
            continue
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None:
            records.append({"network": net, "x": x, "y": y, "status": "degenerate"})
            continue
        z = frozenset(z)
        ctx = make_context(
            random_sem(dag, np.random.default_rng(ft.derive_seed(args.seed, net, x, y))),
            cpdag, x, y, z)
        scale = abs(ctx.theta_z)
        if scale < 1e-12:
            records.append({"network": net, "x": x, "y": y, "status": "zero_estimate"})
            continue

        t0 = time.perf_counter()
        signal.alarm(args.budget)
        try:
            rec = build_record(net, x, y, r["r_val"], None, cpdag, dag, g0, ctx, scale, epsilons)
            signal.alarm(0)
            rec["seconds"] = round(time.perf_counter() - t0, 2)
            records.append(rec)
            print(f"[xcomp18]   {net} {x}->{y} r_val={r['r_val']} k_x={rec['k_x']} "
                  f"slice={rec['slice_size']} D={rec['max_distance_rel']} "
                  f"r_mu_rel={[rec['r_mu_rel'][str(e)] for e in epsilons]} "
                  f"b_top={rec['b_top']:.4g} ({rec['seconds']}s)", flush=True)
        except Timeout:
            signal.alarm(0)
            k_x = len(x_component_edges(cpdag, x)[1])
            records.append({
                "network": net, "x": x, "y": y, "status": "timeout",
                "budget": args.budget, "k_x": k_x,
                "seconds": round(time.perf_counter() - t0, 2),
            })
            print(f"[xcomp18]   TIMEOUT {net} {x}->{y} after {args.budget}s (k_x={k_x})",
                  flush=True)

    with (out_dir / "instances.jsonl").open("w") as f:
        for rec in records:
            f.write(json.dumps(rec) + "\n")

    ok_recs = [r for r in records if r.get("status") == "ok"]
    n_finite_nonzero = sum(1 for r in ok_recs if r["has_finite_r_val"] and r["has_nonzero_bias"])
    n_no_failure = sum(1 for r in ok_recs if not r["has_finite_r_val"] or not r["has_nonzero_bias"])
    timeouts = [r for r in records if r.get("status") == "timeout"]
    summary = {
        "args": vars(args), "epsilons": epsilons,
        "n_skipped_total": len(skipped),
        "n_ok": len(ok_recs),
        "n_finite_r_and_nonzero_bias": n_finite_nonzero,
        "n_no_failure_or_zero_bias": n_no_failure,
        "n_timeout": len(timeouts),
        "timeouts": [(r["network"], r["x"], r["y"], r.get("k_x")) for r in timeouts],
        "seconds": round(time.perf_counter() - t_all, 1),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\n[xcomp18] {summary['n_ok']}/{summary['n_skipped_total']} processed; "
          f"{n_finite_nonzero} finite-r/nonzero-bias, {n_no_failure} no-failure-or-zero, "
          f"{len(timeouts)} timeouts; {summary['seconds']}s", flush=True)


if __name__ == "__main__":
    main()
