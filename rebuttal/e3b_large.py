"""E3b: what the SAT stage costs at large n.

The first design (every b x mode arm on the 12 slowest E2 queries with
n >= 150, 900 s per arm) was stopped: every SAT arm spends its budget building
the encoding (O(n^4) constraints; 3.7 M at n = 150). This version measures, for
the same 12 queries: the depth-3 search time (from E2), the time to build the
E1 encoding (cap 600 s), and, when it builds, one CP-SAT minimise call starting
above depth 3 (cap 600 s). Plus the build time alone on one finite query per n
for n in {24, 48, 100, 150, 200, 300, 400}."""
import csv, multiprocessing as mp, signal, time
import pandas as pd
from ortools.sat.python import cp_model
from common import OUT
import e2_timing as e2
from oset import optimal_set_closed_form
from bkrobust.sat.e1 import _build

class Kill(Exception): pass
def _a(*_): raise Kill

def measure(spec):
    r, label = spec
    signal.signal(signal.SIGALRM, _a)
    dag, cpdag, g0, k = e2.instance(int(r["n"]), float(r["p"]), int(r["seed"]), r["knowledge"] == "half")
    o, _ = optimal_set_closed_form(g0, r["x"], r["y"])
    row = {"label": label, "n": r["n"], "knowledge": r["knowledge"], "k_g0": r["k_g0"], "r_e2": r["r"],
           "search_b3_s_e2": r.get("search_s", ""), "total_e2_s": r["total_s"]}
    signal.alarm(600); t = time.perf_counter()
    try:
        model, _ov, retracted, _d = _build(cpdag, g0, r["x"], r["y"], o, criterion="gac")
        signal.alarm(0)
        row["build_s"] = round(time.perf_counter() - t, 2)
        row["n_constraints"] = len(model.Proto().constraints)
    except Kill:
        row["build_s"] = ">600"
        return row
    if label == "slow":
        model.Add(sum(retracted) >= 4)
        model.Minimize(sum(retracted))
        s = cp_model.CpSolver(); s.parameters.num_workers = 1; s.parameters.max_time_in_seconds = 600
        t = time.perf_counter(); st = s.Solve(model)
        row["minimise_s"] = round(time.perf_counter() - t, 2)
        row["minimise_status"] = s.StatusName(st)
    return row

if __name__ == "__main__":
    d = pd.read_csv(OUT / "e2" / "timing.csv")
    fin = d[(d.query_class == "finite") & d.method.isin(["local_up_fast", "killed", "e1_ladder"])]
    slow = fin[fin.n >= 150].nlargest(12, "total_s")
    scan = pd.concat([fin[fin.n == n].head(1) for n in (24, 48, 100, 150, 200, 300, 400)])
    specs = [(r.to_dict(), "scan") for _, r in scan.iterrows()] + [(r.to_dict(), "slow") for _, r in slow.iterrows()]
    res = []
    with mp.Pool(4) as P:
        for o in P.imap_unordered(measure, specs):
            res.append(o); print(o, flush=True)
    keys = sorted({k for o in res for k in o})
    with open(OUT / "e3" / "large_sat_cost.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(res)
    print("wrote")
