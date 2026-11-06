"""E3b: depth budget and SAT mode on the slowest large-n queries of E2.

The 12 slowest finite-radius queries of E2 with n >= 150 (including the one
killed at 600 s) are rebuilt from their (family, n, p, seed, knowledge) and
rerun under Algorithm 1 with b in {0, 1, 2, 3} and SAT mode ladder / minimise
(same arm() as E3). Per-arm wall limit 900 s."""
import csv, multiprocessing as mp, signal
import pandas as pd
from common import OUT, UNREACHED
import e2_timing as e2
import e3_budget_sat as e3
from oset import optimal_set_closed_form

class Kill(Exception): pass
def _a(*_): raise Kill

def work(r):
    signal.signal(signal.SIGALRM, _a)
    dag, cpdag, g0, k = e2.instance(int(r.n), float(r.p), int(r.seed), r.knowledge == "half")
    o, _ = optimal_set_closed_form(g0, r.x, r.y)
    out = []
    for b in (0, 1, 2, 3):
        for mode in ("ladder", "minimise"):
            signal.alarm(900)
            try:
                rr, ex, s, leg = e3.arm(cpdag, g0, r.x, r.y, o, b, mode)
                signal.alarm(0)
            except Kill:
                rr, ex, s, leg = None, False, 900.0, "killed"
            out.append({"family": r.family, "n": r.n, "knowledge": r.knowledge, "k_g0": r.k_g0,
                        "r_e2": r.r, "b": b, "mode": mode, "r": rr, "exact": ex,
                        "seconds": round(s, 4), "leg": leg})
    return out

if __name__ == "__main__":
    d = pd.read_csv(OUT / "e2" / "timing.csv")
    d = d[(d.query_class == "finite") & (d.n >= 150)].nlargest(12, "total_s")
    res = []
    with mp.Pool(4) as P:
        for o in P.imap_unordered(work, [r for _, r in d.iterrows()]):
            res.extend(o); print(o[0]["n"], o[0]["k_g0"], [(x["b"], x["mode"], x["r"], x["seconds"]) for x in o], flush=True)
    with open(OUT / "e3" / "large_budget.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(res[0])); w.writeheader(); w.writerows(res)
