"""E21: compute the walkthrough toy (and the Figure 2(a) toy) exactly before drawing them.

Every number the walkthrough figure will show is computed here from the submitted code, and
each radius is computed five ways that share as little code as possible; any disagreement is a
consistency failure (exit 2). See PREREG.md for the registered statistic and the predictions.

    python experiments/rebuttal/E21_walkthrough_toy/run.py --smoke --out /tmp/e21_smoke
    python experiments/rebuttal/E21_walkthrough_toy/run.py --out results/rebuttal/E21
"""

from __future__ import annotations

import argparse
import itertools
import json
import multiprocessing as mp
import subprocess
import sys
import threading
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_REPO = HERE.parents[2]

CRITERION = "gac"
UNREACHED = -1  # bkrobust.core.conventions.UNREACHED; re-asserted after import below


# --------------------------------------------------------------------------- the toys

def toys() -> list[dict]:
    """The units, in a fixed order. ``hand`` holds the predictions written in PREREG.md."""
    from bkrobust.demo.example import dag_to_cpdag
    from bkrobust.demo.graph import MPDAG
    from bkrobust.demo.scenario import OUTCOME, TREATMENT, scenarios, true_dag

    walk = MPDAG(["U", "X", "W", "Y"],
                 undirected=[("U", "X"), ("X", "Y"), ("X", "W"), ("W", "Y")])
    fig2 = MPDAG(["A", "W", "X", "Y"],
                 undirected=[("A", "W"), ("W", "X"), ("W", "Y"), ("X", "Y")])
    return [
        {
            "unit": "walkthrough", "smoke": True, "cpdag": walk,
            "knowledge": [("U", "X"), ("X", "Y"), ("X", "W")], "x": "X", "y": "Y",
            "hand": {
                "z": [], "r_val": 2, "r_claim": 2, "n_space": None, "n_upset": 5,
                "upset_valid": {
                    "U->X X->W X->Y": True, "X->W X->Y": True,
                    "X->Y": False, "X->W": False, "": False,
                },
                "witnesses": ["X->W", "X->Y"],
                "n_dags_below_g0": 2, "n_dags_below_g0_valid": 2,
                "non_states": ["U->X X->W"],
            },
        },
        {
            "unit": "fig2a", "smoke": True, "cpdag": fig2,
            "knowledge": [("A", "W"), ("X", "Y")], "x": "X", "y": "Y",
            "hand": {"z": ["W"], "r_val": 1, "printed_in_figure_2b": 3},
        },
        {
            # Added beyond the plan: the drafted reply promises the figure on this example.
            "unit": "running_example_B", "smoke": False, "cpdag": dag_to_cpdag(true_dag()),
            "knowledge": list(scenarios()["B"]["knowledge"]), "x": TREATMENT, "y": OUTCOME,
            "hand": {"z": ["Age", "Smoke"], "r_val": 3, "n_space": 48},
        },
    ]


# --------------------------------------------------------------------------- helpers

def kstr(edges) -> str:
    return " ".join(f"{a}->{b}" for a, b in sorted(edges))


def compute_unit(toy: dict) -> dict:
    """All the exact quantities for one toy, plus the list of consistency failures."""
    from bkrobust.core import conventions
    from bkrobust.core.oracle import extensions
    from bkrobust.demo.example import dag_to_cpdag
    from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag
    from bkrobust.demo.meek import apply_orientations
    from bkrobust.demo.space import bfs_distances, enumerate_space
    from bkrobust.gac.dag_level import is_gac_valid_dag
    from bkrobust.gac.mpdag_level import is_gac_valid_mpdag
    from bkrobust.hybrid import breakdown_radius
    from bkrobust.mpdag_criterion import is_valid_mpdag as backdoor_valid_mpdag
    from bkrobust.robustness.real_survival import commit_z_star
    from bkrobust.search.space_fixed import build_space_correct, knowledge_of

    assert conventions.UNREACHED == UNREACHED
    t0 = time.perf_counter()
    cpdag, k, x, y = toy["cpdag"], [tuple(e) for e in toy["knowledge"]], toy["x"], toy["y"]
    fail: list[str] = []
    row: dict = {"unit": toy["unit"], "criterion": CRITERION, "x": x, "y": y,
                 "cpdag": cpdag.edge_string(), "knowledge": kstr(k), "hand": toy["hand"]}

    # The input must be a genuine CPDAG, or nothing below means what the figure says.
    cp_ext = extensions(cpdag)
    row["n_dags_in_class"] = len(cp_ext)
    if not cp_ext or any(dag_to_cpdag(d) != cpdag for d in cp_ext):
        fail.append("cpdag_is_not_the_cpdag_of_its_extensions")

    g0 = apply_orientations(cpdag, k)
    if g0 is None:
        row["consistency_failures"] = fail + ["knowledge_inconsistent"]
        return row
    kg0 = frozenset(knowledge_of(cpdag, g0))
    row["g0"] = g0.edge_string()
    row["k_g0"] = kstr(kg0)
    row["added_by_meek"] = kstr(kg0 - set(k))

    o_raw = optimal_adjustment_set_mpdag(g0, x, y)
    z, z_status = commit_z_star(g0, x, y)
    row["optimal_set_raw"] = None if o_raw is None else sorted(o_raw)
    row["z_status"] = z_status
    if z is None:
        row["consistency_failures"] = fail + [f"z_not_committed:{z_status}"]
        return row
    row["z"] = sorted(z)

    def gac(g) -> bool:
        return is_gac_valid_mpdag(g, x, y, z)

    def gac_enum(g) -> bool:
        ex = extensions(g)
        return bool(ex) and all(is_gac_valid_dag(d, x, y, z) for d in ex)

    if not gac(g0):
        fail.append("z_invalid_at_g0")

    # ---- the whole space, its order and covering graph (no conjectured shortcut)
    space = build_space_correct(cpdag)
    el = list(space.elements)
    reps, dirs = space.reps, {g: frozenset(g.directed_edges) for g in el}
    ks = {g: frozenset(knowledge_of(cpdag, g)) for g in el}
    row["n_space"] = len(el)
    row["n_space_legacy_enumerate_space"] = len(enumerate_space(cpdag))
    row["n_covers"] = len(space.covers)
    if g0 not in reps:
        fail.append("g0_not_in_space")

    valid = {g: gac(g) for g in el}
    mism = [g.edge_string() for g in el if gac_enum(g) != valid[g]]
    row["oracle_disagreements"] = mism
    if mism:
        fail.append("gac_criterion_vs_enumeration_disagree")
    bd_valid = {g: backdoor_valid_mpdag(g, x, y, z) for g in el}

    lemma_o = sum((reps[g] <= reps[h]) != (dirs[h] <= dirs[g]) for g in el for h in el)
    row["lemma_o_violations"] = lemma_o
    if lemma_o:
        fail.append("lemma_o_violated")

    # Lemma 2: the join of two states is Meek(C, K_G intersect K_H).
    join_bad, n_pairs = [], 0
    for i, g in enumerate(el):
        for h in el[i + 1:]:
            n_pairs += 1
            ub = [u for u in el if reps[g] <= reps[u] and reps[h] <= reps[u]]
            least = [u for u in ub if all(reps[u] <= reps[v] for v in ub)]
            want = apply_orientations(cpdag, sorted(ks[g] & ks[h]))
            if len(least) != 1 or want is None or least[0] != want:
                join_bad.append([g.edge_string(), h.edge_string()])
    row["join_pairs_checked"] = n_pairs
    row["join_violations"] = join_bad[:20]
    if join_bad:
        fail.append("lemma_2_join_violated")

    # ---- radius, way 4: BFS on the covering graph of the full space
    dist = bfs_distances(space.neighbours, g0)
    failing = [g for g in el if not valid[g]]
    reach = [dist[g] for g in failing if g in dist]
    r_bfs = min(reach) if reach else UNREACHED
    wit = sorted(g.edge_string() for g in failing if dist.get(g) == r_bfs)
    row["n_failing_states"] = len(failing)
    row["n_failing_unreachable"] = sum(g not in dist for g in failing)

    # ---- the up-set of G0 and its Hasse diagram (what the figure draws)
    up = [h for h in el if reps[g0] <= reps[h]]
    if {h for h in el if dirs[h] <= dirs[g0]} != set(up):
        fail.append("upset_by_models_vs_by_edges_disagree")
    up_set = set(up)
    up_covers = [(lo, hi) for lo, hi in space.covers if lo in up_set and hi in up_set]
    up_dist = {g0: 0}
    frontier = [g0]
    while frontier:  # BFS upward only: the number of retractions
        nxt = []
        for g in frontier:
            for lo, hi in up_covers:
                if lo == g and hi not in up_dist:
                    up_dist[hi] = up_dist[g] + 1
                    nxt.append(hi)
        frontier = nxt
    up_fail = [up_dist[h] for h in up if not valid[h] and h in up_dist]
    r_up = min(up_fail) if up_fail else UNREACHED
    row["n_upset"] = len(up)
    row["upset"] = [
        {"k": kstr(ks[h]), "graph": h.edge_string(), "valid_gac": valid[h],
         "valid_backdoor": bd_valid[h], "retractions_from_g0": up_dist.get(h),
         "n_dags": len(reps[h])}
        for h in sorted(up, key=lambda h: (-len(ks[h]), kstr(ks[h])))
    ]
    row["upset_covers"] = [
        {"lower": kstr(ks[lo]), "upper": kstr(ks[hi]), "withdrawn": kstr(ks[lo] - ks[hi])}
        for lo, hi in sorted(up_covers, key=lambda p: (kstr(ks[p[0]]), kstr(ks[p[1]])))
    ]
    row["upset_witnesses"] = sorted(kstr(ks[h]) for h in up
                                    if not valid[h] and up_dist.get(h) == r_up)

    # Subsets of K_G0 that are not states: the closure re-adds what was withdrawn.
    if len(kg0) <= 12:
        closed, ghosts = set(), []
        for n in range(len(kg0) + 1):
            for s in itertools.combinations(sorted(kg0), n):
                c = apply_orientations(cpdag, s)
                kc = frozenset(knowledge_of(cpdag, c)) if c is not None else None
                if kc == frozenset(s):
                    closed.add(kc)
                else:
                    ghosts.append({"subset": kstr(s),
                                   "closure_adds": None if kc is None else kstr(kc - set(s))})
        row["non_states"] = ghosts
        if closed != {ks[h] for h in up}:
            fail.append("closed_subsets_of_k_g0_differ_from_upset")

    below = [lo for lo, hi in space.covers if hi == g0]
    row["covered_by_g0"] = [{"graph": d.edge_string(), "is_dag": d.is_dag(),
                             "valid_gac": valid[d]} for d in sorted(below, key=MPDAG_key)]

    # ---- radius, way 5: the closed form over the DAGs of the class (Cor. distance (c))
    bad_d = [d for d in cp_ext if not is_gac_valid_dag(d, x, y, z)]
    costs = [len(kg0 - set(d.directed_edges)) for d in bad_d]
    r_form = min(costs) if costs else UNREACHED
    row["formula_witness_dags"] = sorted(d.edge_string() for d, c in zip(bad_d, costs)
                                         if c == r_form)

    # ---- radius, ways 1-3: the submitted hybrid, three dispatches
    h_crit = breakdown_radius(cpdag, k, x, y, z, g0=g0, criterion=CRITERION)
    h_enum = breakdown_radius(cpdag, k, x, y, z, g0=g0, criterion=CRITERION,
                              use_criterion=False)
    h_sat = breakdown_radius(cpdag, k, x, y, z, g0=g0, criterion=CRITERION, search_budget=0)
    h_bd = breakdown_radius(cpdag, k, x, y, z, g0=g0, criterion="backdoor")
    hy = {}
    for name, h in (("criterion", h_crit), ("enumeration", h_enum), ("sat_forced", h_sat),
                    ("backdoor", h_bd)):
        hy[name] = {"radius": h.radius, "method": h.method, "oracle": h.oracle,
                    "exact": h.exact, "witness": h.witness, "seconds": h.total_seconds}
    row["hybrid"] = hy

    radii = {"hybrid_criterion": h_crit.radius, "hybrid_enumeration": h_enum.radius,
             "hybrid_sat_forced": h_sat.radius, "space_bfs": r_bfs,
             "upset_retractions": r_up, "closed_form": r_form}
    row["radii"] = radii
    row["r_val"] = h_crit.radius
    if len(set(radii.values())) != 1:
        fail.append("radii_disagree")
    if not all(h.exact for h in (h_crit, h_enum, h_sat)):
        fail.append("hybrid_not_exact")
    if h_sat.radius not in (0, UNREACHED) and h_sat.method != "e1_ladder":
        fail.append("sat_leg_not_exercised")
    row["witnesses"] = wit
    row["witness_knowledge"] = sorted(kstr(ks[g]) for g in failing if dist.get(g) == r_bfs)
    for h in (h_crit, h_enum):
        if h.method in ("local_up_fast", "degenerate") and h.witness not in wit:
            fail.append("hybrid_witness_not_a_nearest_failure")

    # ---- claim radius over every subset of the elicited claims
    kk = sorted(set(k))
    r_claim, claim_wit, inconsistent = UNREACHED, [], 0
    for n in range(len(kk) + 1):
        for r in itertools.combinations(kk, n):
            g = apply_orientations(cpdag, [e for e in kk if e not in r])
            if g is None:
                inconsistent += 1
            elif not gac(g):
                claim_wit.append(kstr(r))
        if claim_wit:
            r_claim = n
            break
    singles = [apply_orientations(cpdag, [e for e in kk if e != c]) for c in kk]
    row["r_claim"] = r_claim
    row["claim_witnesses"] = claim_wit
    row["claim_subsets_inconsistent"] = inconsistent
    row["phi_1"] = (sum(g is not None and not gac(g) for g in singles) / len(kk)) if kk else None
    c_costs = [len(set(kk) - set(d.directed_edges)) for d in bad_d]
    row["r_claim_closed_form"] = min(c_costs) if c_costs else UNREACHED
    if row["r_claim_closed_form"] != r_claim:
        fail.append("r_claim_subsets_vs_closed_form_disagree")
    if r_claim != UNREACHED and h_crit.radius != UNREACHED and r_claim > h_crit.radius:
        fail.append("r_claim_exceeds_r_val")

    # ---- the whole space, for drawing (ids index the list)
    idx = {g: i for i, g in enumerate(el)}
    row["space"] = [{"id": idx[g], "k": kstr(ks[g]), "graph": g.edge_string(),
                     "is_dag": g.is_dag(), "valid_gac": valid[g], "valid_backdoor": bd_valid[g],
                     "dist_from_g0": dist.get(g), "in_upset": g in up_set} for g in el]
    row["space_covers"] = sorted([idx[lo], idx[hi]] for lo, hi in space.covers)

    row["consistency_failures"] = fail
    row["seconds"] = time.perf_counter() - t0
    return row


def MPDAG_key(g) -> str:  # noqa: N802 - a sort key, named for what it sorts
    return g.edge_string()


# --------------------------------------------------------------------------- driver

def _worker(args: tuple) -> dict:
    repo, toy_name, smoke, commit = args
    sys.path.insert(0, str(Path(repo) / "src"))
    toy = next(t for t in toys() if t["unit"] == toy_name)
    try:
        row = compute_unit(toy)
        row["status"] = "ok"
    except Exception:  # a crash is recorded as a failed unit, not lost
        row = {"unit": toy_name, "criterion": CRITERION, "status": "error",
               "error": traceback.format_exc(), "consistency_failures": ["crash"]}
    row["commit"] = commit
    row["smoke"] = smoke
    return row


def git(repo: Path, *a: str) -> str:
    try:
        return subprocess.run(["git", "-C", str(repo), *a], capture_output=True, text=True,
                              timeout=30).stdout.strip()
    except Exception:
        return ""


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--smoke", action="store_true", help="the two four-node toys only")
    ap.add_argument("--out", default="results/rebuttal/E21")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--repo", default=str(DEFAULT_REPO), help="repository root (holds src/)")
    ap.add_argument("--units", default="", help="comma-separated unit names to restrict to")
    args = ap.parse_args(argv)
    if not 0 <= args.shard < args.n_shards:
        ap.error("need 0 <= shard < n-shards")

    repo = Path(args.repo).resolve()
    sys.path.insert(0, str(repo / "src"))
    commit = git(repo, "rev-parse", "HEAD") or "unknown"
    src_dirty = bool(git(repo, "status", "--porcelain", "--", "src"))

    all_toys = toys()
    names = [t["unit"] for t in all_toys if t["smoke"] or not args.smoke]
    if args.units:
        keep = set(args.units.split(","))
        names = [n for n in names if n in keep]
    mine = [n for i, n in enumerate(names) if i % args.n_shards == args.shard]

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    tag = "" if args.n_shards == 1 else f".shard{args.shard}of{args.n_shards}"
    rows_path, cfg_path = out / f"rows{tag}.jsonl", out / f"config{tag}.json"
    cfg_path.write_text(json.dumps({
        "experiment": "E21", "criterion": CRITERION, "smoke": args.smoke, "commit": commit,
        "src_dirty": src_dirty, "units_selected": names, "n_units_selected": len(names),
        "units_this_shard": mine, "shard": args.shard, "n_shards": args.n_shards,
        "search_budget_default": 3, "sat_forced_search_budget": 0,
        "methods": ["hybrid_criterion", "hybrid_enumeration", "hybrid_sat_forced",
                    "space_bfs", "upset_retractions", "closed_form"],
    }, indent=1))

    done: dict[str, dict] = {}
    if rows_path.exists():
        for line in rows_path.open():
            if line.strip():
                r = json.loads(line)
                done[r["unit"]] = r
    todo = [n for n in mine if n not in done]
    print(f"E21: {len(mine)} units in shard, {len(done)} already written, {len(todo)} to run; "
          f"commit {commit[:10]} criterion {CRITERION}", flush=True)

    t0 = time.time()
    current = {"unit": None}
    stop = threading.Event()

    def heartbeat() -> None:
        while not stop.wait(30):
            print(f"  alive {time.time() - t0:.0f} s, working on {current['unit']}", flush=True)

    threading.Thread(target=heartbeat, daemon=True).start()
    jobs = [(str(repo), n, args.smoke, commit) for n in todo]
    with rows_path.open("a") as fh:
        if args.workers > 1 and len(jobs) > 1:
            with mp.Pool(min(args.workers, len(jobs))) as pool:
                results = pool.imap_unordered(_worker, jobs)
                for row in results:
                    fh.write(json.dumps(row) + "\n")
                    fh.flush()
                    done[row["unit"]] = row
                    print(f"  {row['unit']}: r_val={row.get('r_val')} "
                          f"failures={row['consistency_failures']}", flush=True)
        else:
            for job in jobs:
                current["unit"] = job[1]
                row = _worker(job)
                fh.write(json.dumps(row) + "\n")
                fh.flush()
                done[row["unit"]] = row
                print(f"  {row['unit']}: r_val={row.get('r_val')} radii={row.get('radii')} "
                      f"failures={row['consistency_failures']} "
                      f"({time.time() - t0:.1f} s)", flush=True)
    stop.set()

    bad = {n: done[n]["consistency_failures"] for n in mine
           if n in done and done[n]["consistency_failures"]}
    if bad:
        print(f"CONSISTENCY FAILURE: {bad}", flush=True)
        return 2
    print(f"E21 shard {args.shard}: all {len(mine)} units consistent", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
