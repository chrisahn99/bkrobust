"""E3S driver: the SAT stage against other exact solvers (PREREGISTRATION.md, E3S).

Every arm solves the same problem from k = 0, with no search stage:
the least number of retracted closure orientations of G0 after which Z fails
(GAC). Arms:
  A1 cpsat_ladder    submitted E1 ladder, CP-SAT feasibility calls k = 0, 1, ...
  A2 cpsat_min       one CP-SAT call minimising the retractions
  A3 cpsat_core      A2 with core-based search only (MaxSAT-style, inside CP-SAT)
  A4 milp_scip       MILP translation, SCIP
  A5 milp_highs      MILP translation, HiGHS
  A6 milp_cbc        MILP translation, CBC
  A7 maxsat_rc2      WCNF translation, RC2 (PySAT)
  A8 sat_cadical     CNF translation, CaDiCaL 1.9.5 ladder via an incremental totalizer
  A9 explicit        retraction search with unbounded depth (radius_local_up_fast)
Each (query, arm) runs in its own forked process, killed at CAP + 30 s; the
build of the E1 model is counted in every arm that uses it. Every returned
witness is decoded and re-checked.
"""

from __future__ import annotations

import argparse
import csv
import json
import multiprocessing as mp
import os
import sys
import time
import traceback

import numpy as np

from common import OUT, UNREACHED, apply_orientations, is_gac_valid_mpdag, k_of

ODIR = OUT / "e3s"
CAP = 300.0
ARMS = ["cpsat_ladder", "cpsat_min", "cpsat_core", "milp_scip", "milp_highs", "milp_cbc",
        "maxsat_rc2", "sat_cadical", "explicit"]


# ------------------------------------------------------------------ populations


def units_e3():
    seen = {}
    for r in csv.DictReader(open(OUT / "e3" / "budget_sat.csv")):
        key = (r["n"], r["p"], r["seed"], r["half"], r["x"], r["y"])
        if key not in seen:
            seen[key] = {"qid": "e3:" + ":".join(key), "pop": "P-E3", "n": int(r["n"]),
                         "p": float(r["p"]), "seed": int(r["seed"]), "half": r["half"] == "True",
                         "x": r["x"], "y": r["y"], "r_ref": int(r["r_default"])}
    return list(seen.values())


def units_satr():
    import e4l_listing as e4l
    want = {}
    for line in open(OUT / "e4l" / "rows.jsonl"):
        r = json.loads(line)
        if r.get("status") == "ok" and r.get("hyb_method") == "e1_ladder":
            want[r["uid"]] = r["hyb_r"]
    out = []
    for gen in (e4l.units_e2, e4l.units_rel, e4l.units_syn):
        for u in gen():
            if u["uid"] in want:
                out.append({"qid": "satr:" + u["uid"], "pop": "P-SATR", "unit": u, "r_ref": want[u["uid"]],
                            "src_pop": u["pop"]})
    return out


def units_scale(per_n=10):
    rows = [r for r in csv.DictReader(open(OUT / "e2" / "timing.csv"))
            if r["r"] not in ("", "inf") and r["method"] not in ("killed", "no_query")]
    out = []
    for n in (24, 32, 48, 64):
        rs = sorted([r for r in rows if int(r["n"]) == n], key=lambda r: (-int(r["r"]), r["family"], r["seed"]))
        for r in rs[:per_n]:
            out.append({"qid": f"scale:{n}:{r['family']}:{r['p']}:{r['seed']}:{r['knowledge']}:{r['x']}:{r['y']}",
                        "pop": "P-SCALE", "n": n, "p": float(r["p"]), "seed": int(r["seed"]),
                        "half": r["knowledge"] == "half", "family": r["family"], "x": r["x"], "y": r["y"],
                        "r_ref": int(r["r"])})
    return out


def build(q):
    """(cpdag, g0, x, y, z) for a query."""
    from oset import optimal_set_closed_form
    if q["pop"] == "P-E3":
        from bkrobust.demo.example import dag_to_cpdag, knowledge_to_recover
        from bkrobust.synth.generators import erdos_renyi_dag
        rng = np.random.default_rng(q["seed"])
        dag = erdos_renyi_dag(q["n"], rng, q["p"])
        cpdag = dag_to_cpdag(dag)
        k = sorted(knowledge_to_recover(dag, cpdag))
        if q["half"]:
            k = k[: len(k) // 2]
        g0 = apply_orientations(cpdag, k)
        x, y = q["x"], q["y"]
    elif q["pop"] == "P-SCALE":
        import e2_timing as e2
        _d, cpdag, g0, _k = e2.instance(q["n"], q["p"], q["seed"], q["half"])
        x, y = q["x"], q["y"]
    else:
        import e4l_listing as e4l
        cpdag, g0, x, y = e4l.build(q["unit"])
    z, st = optimal_set_closed_form(g0, x, y)
    assert st == "ok", st  # Z may be empty (amendment 2, E3S)
    return cpdag, g0, x, y, frozenset(z)


# ------------------------------------------------------------------------ arms


def _check(cpdag, g0, x, y, z, kept, r):
    """Re-check a witness: kept orientations closed, Z invalid there, |K_G0 - kept| = r."""
    k0 = k_of(cpdag, g0)
    kept = frozenset(kept) - set(cpdag.directed_edges)
    g = apply_orientations(cpdag, sorted(kept))
    if g is None or k_of(cpdag, g) != kept or not kept <= k0:
        return False
    return (not is_gac_valid_mpdag(g, x, y, z)) and len(k0 - kept) == r


def run_arm(q, arm, cap):
    from ortools.sat.python import cp_model

    from bkrobust.sat.e1 import _build
    from bkrobust.search.exact_fast import radius_local_up_fast

    res = {"qid": q["qid"], "pop": q["pop"], "arm": arm, "r_ref": q["r_ref"]}
    t_all = time.perf_counter()
    cpdag, g0, x, y, z = build(q)
    res["n_nodes"] = len(cpdag.nodes)
    res["k_g0"] = len(k_of(cpdag, g0))
    t0 = time.perf_counter()
    res["instance_s"] = round(t0 - t_all, 4)

    def left():
        return max(1.0, cap - (time.perf_counter() - t0))

    if arm == "explicit":
        def fails(g):
            return not is_gac_valid_mpdag(g, x, y, z)
        f = radius_local_up_fast(cpdag, g0, fails, max_depth=res["k_g0"])
        r = f.radius
        res.update(r=None if r == UNREACHED else r, exact=bool(f.exact), witness_ok=None,
                   total_s=round(time.perf_counter() - t0, 6), status="ok")
        return res

    model, _ov, retracted, d = _build(cpdag, g0, x, y, z, criterion="gac")
    obj = [v.Index() for v in retracted]
    dvars = {e: v.Index() for e, v in d.items()}
    res["build_s"] = round(time.perf_counter() - t0, 6)

    def decode(val):
        return [e for e, i in dvars.items() if val(i)]

    if arm.startswith("cpsat"):
        proto_vars = len(model.Proto().variables)
        res["n_vars"], res["n_cons"] = proto_vars, len(model.Proto().constraints)
        t1 = time.perf_counter()
        if arm == "cpsat_ladder":
            r, exact, kept = None, False, None
            rungs = []
            res["rungs"] = rungs
            for k in range(0, len(obj) + 1):
                tr0 = time.perf_counter()
                sub = model.Clone()
                if obj:
                    sub.Add(sum(sub.GetBoolVarFromProtoIndex(i) for i in obj) <= k)
                s = cp_model.CpSolver()
                s.parameters.num_workers = 1
                s.parameters.random_seed = 0
                s.parameters.max_time_in_seconds = left()
                st = s.Solve(sub)
                rungs.append((k, round(time.perf_counter() - tr0, 6), int(st)))
                if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
                    r, exact = k, True
                    kept = decode(lambda i: s.Value(sub.GetBoolVarFromProtoIndex(i)))
                    break
                if st != cp_model.INFEASIBLE:
                    break
            else:
                r, exact = UNREACHED, True
        else:
            if obj:
                model.Minimize(sum(model.GetBoolVarFromProtoIndex(i) for i in obj))
            s = cp_model.CpSolver()
            s.parameters.num_workers = 1
            s.parameters.random_seed = 0
            s.parameters.max_time_in_seconds = left()
            if arm == "cpsat_core":
                s.parameters.optimize_with_core = True
            st = s.Solve(model)
            r, exact, kept = None, False, None
            if st == cp_model.OPTIMAL:
                r, exact = int(round(s.ObjectiveValue())), True
                kept = decode(lambda i: s.Value(model.GetBoolVarFromProtoIndex(i)))
            elif st == cp_model.INFEASIBLE:
                r, exact = UNREACHED, True
        res["solve_s"] = round(time.perf_counter() - t1, 6)
    else:
        import e3s_translate as tr
        t1 = time.perf_counter()
        pb = tr.export_proto(model)
        res["n_vars"], res["n_cons"] = len(pb.variables), len(pb.constraints)
        if arm.startswith("milp"):
            backend = {"milp_scip": "SCIP", "milp_highs": "HIGHS", "milp_cbc": "CBC"}[arm]
            mm = tr.to_milp(pb, obj, backend)
            res["trans_s"] = round(time.perf_counter() - t1, 6)
            res["milp_rows"] = mm.n_rows
            s = mm.solver
            s.SetNumThreads(1)
            s.SetTimeLimit(int(left() * 1000))
            t2 = time.perf_counter()
            st = s.Solve()
            res["solve_s"] = round(time.perf_counter() - t2, 6)
            from ortools.linear_solver import pywraplp
            r, exact, kept = None, False, None
            if st == pywraplp.Solver.OPTIMAL:
                r, exact = int(round(s.Objective().Value())), True
                kept = decode(lambda i: mm.x[i].solution_value() > 0.5)
            elif st == pywraplp.Solver.INFEASIBLE:
                r, exact = UNREACHED, True
            res["milp_status"] = int(st)
        else:
            cm = tr.to_cnf(pb, obj)
            res["trans_s"] = round(time.perf_counter() - t1, 6)
            res["clauses"], res["cnf_vars"] = len(cm.hard), cm.top
            t2 = time.perf_counter()
            r, exact, kept = None, False, None
            if arm == "maxsat_rc2":
                from pysat.examples.rc2 import RC2
                from pysat.formula import WCNF
                w = WCNF()
                w.extend(cm.hard)
                for v in cm.soft:
                    w.append([-v], weight=1)
                with RC2(w) as rc2:
                    sol = rc2.compute()
                    if sol is None:
                        r, exact = UNREACHED, True
                    else:
                        r, exact = int(rc2.cost), True
                        val = set(l for l in sol if l > 0)
                        kept = decode(lambda i: cm.var_of[i] in val)
            else:
                from pysat.card import ITotalizer
                from pysat.solvers import Solver
                with Solver(name="cadical195", bootstrap_with=cm.hard) as sv:
                    if cm.soft:
                        tot = ITotalizer(lits=cm.soft, ubound=len(cm.soft), top_id=cm.top)
                        for cl in tot.cnf.clauses:
                            sv.add_clause(cl)
                    rungs = []
                    res["rungs"] = rungs
                    for k in range(0, len(cm.soft) + 1):
                        assum = [-tot.rhs[k]] if (cm.soft and k < len(cm.soft)) else []
                        tr0 = time.perf_counter()
                        ok_k = sv.solve(assumptions=assum)
                        rungs.append((k, round(time.perf_counter() - tr0, 6), int(bool(ok_k))))
                        if ok_k:
                            r, exact = k, True
                            val = set(l for l in sv.get_model() if l > 0)
                            kept = decode(lambda i: cm.var_of[i] in val)
                            break
                    else:
                        r, exact = UNREACHED, True
            res["solve_s"] = round(time.perf_counter() - t2, 6)
    res["r"] = None if r is None or r == UNREACHED else r
    res["unreached"] = r == UNREACHED
    res["exact"] = exact
    res["witness_ok"] = (_check(cpdag, g0, x, y, z, kept, r) if kept is not None else None)
    res["total_s"] = round(time.perf_counter() - t0, 6)
    res["status"] = "ok" if exact else "cap"
    return res


def _child(q, arm, cap, conn):
    try:
        os.environ["OMP_NUM_THREADS"] = "1"
        conn.send(run_arm(q, arm, cap))
    except Exception as e:  # noqa: BLE001
        conn.send({"qid": q["qid"], "pop": q["pop"], "arm": arm, "r_ref": q["r_ref"],
                   "status": f"error:{type(e).__name__}:{e}"[:300],
                   "trace": traceback.format_exc()[-1500:]})
    conn.close()


def schedule(tasks, cap, workers, out_path):
    ctx = mp.get_context("fork")
    running = []
    done = 0
    t_start = time.time()
    with open(out_path, "a") as fout:
        while tasks or running:
            while tasks and len(running) < workers:
                q, arm = tasks.pop(0)
                a, b = ctx.Pipe(duplex=False)
                p = ctx.Process(target=_child, args=(q, arm, cap, b))
                p.start()
                b.close()
                running.append((p, a, q, arm, time.time()))
            time.sleep(0.02)
            still = []
            for p, a, q, arm, t in running:
                res = None
                if a.poll():
                    try:
                        res = a.recv()
                    except EOFError:
                        res = None
                    p.join(5)
                    if res is None:
                        res = {"qid": q["qid"], "pop": q["pop"], "arm": arm, "r_ref": q["r_ref"],
                               "status": "crashed"}
                elif not p.is_alive():
                    res = {"qid": q["qid"], "pop": q["pop"], "arm": arm, "r_ref": q["r_ref"],
                           "status": f"crashed:exit{p.exitcode}"}
                elif time.time() - t > cap + 30:
                    p.kill()
                    p.join(5)
                    res = {"qid": q["qid"], "pop": q["pop"], "arm": arm, "r_ref": q["r_ref"],
                           "status": "killed", "total_s": round(time.time() - t, 3)}
                if res is None:
                    still.append((p, a, q, arm, t))
                    continue
                res["wall_s"] = round(time.time() - t, 3)
                fout.write(json.dumps(res, default=str) + "\n")
                fout.flush()
                done += 1
                if done % 100 == 0:
                    print(f"{done} done, {len(tasks)} queued, {time.time() - t_start:.0f}s", flush=True)
            running = still


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--cap", type=float, default=CAP)
    ap.add_argument("--pops", default="P-E3,P-SATR,P-SCALE")
    ap.add_argument("--arms", default=",".join(ARMS))
    ap.add_argument("--limit", type=int, default=0, help="queries per population (smoke)")
    ap.add_argument("--satr-sample", type=int, default=0)
    ap.add_argument("--out", default="rows.jsonl")
    a = ap.parse_args()
    ODIR.mkdir(parents=True, exist_ok=True)
    pops = a.pops.split(",")
    qs = []
    if "P-E3" in pops:
        qs += units_e3()[: a.limit or None]
    if "P-SCALE" in pops:
        qs += units_scale()[: a.limit or None]
    if "P-SATR" in pops:
        s = units_satr()
        if a.satr_sample:
            rng = np.random.default_rng(0)
            # stratified by source population and radius
            strata = {}
            for q in s:
                strata.setdefault((q["src_pop"], min(q["r_ref"], 7)), []).append(q)
            tot = len(s)
            pick = []
            for key, v in sorted(strata.items()):
                m = max(1, round(a.satr_sample * len(v) / tot))
                idx = rng.permutation(len(v))[:m]
                pick += [v[i] for i in sorted(idx)]
            s = pick
        qs += s[: a.limit or None]
    path = ODIR / a.out
    done = set()
    if path.exists():
        for line in open(path):
            r = json.loads(line)
            done.add((r["qid"], r["arm"]))
    arms = a.arms.split(",")
    tasks = [(q, arm) for q in qs for arm in arms if (q["qid"], arm) not in done]
    print(f"{len(qs)} queries, {len(tasks)} tasks to run, cap {a.cap}s, {a.workers} workers", flush=True)
    schedule(tasks, a.cap, a.workers, path)
    print("done", flush=True)


if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    main()
