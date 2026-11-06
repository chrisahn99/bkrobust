"""E1: the release rows whose radius is reported only as a lower bound.

Rows with hop_method == local_up_fast_budget carry r_hop = -1 and a lower
bound. For every such row in an elicited condition (state rebuildable from
the release), run the submitted hybrid with the GAC and with the back-door
predicate, CPDAG check first, SAT ladder beyond depth 3, per-rung limit
--rung-s. Records radius, method, exactness and time per call.
"""

from __future__ import annotations

import argparse
import csv
import multiprocessing as mp
import time

from common import OUT, UNREACHED, hybrid, is_gac_valid_mpdag
from oset import optimal_set_closed_form
from release_rows import knowledge, rows, state


def work(args):
    r, rung_s = args
    dag, cpdag, k, g0, _ = state(r["arm"], r["network"])
    x, y = r["x"], r["y"]
    o, st = optimal_set_closed_form(g0, x, y)
    out = {"arm": r["arm"], "network": r["network"], "x": x, "y": y,
           "lower_bound": r["hop_lower_bound"], "z_size_release": r["z_size"],
           "z_size": len(o) if o is not None else "", "o_status": st,
           "z_valid_at_truth": r["z_valid_at_truth"]}
    if o is None:
        return out
    out["valid_at_chat"] = is_gac_valid_mpdag(cpdag, x, y, o)
    for crit in ("gac", "backdoor"):
        t = time.perf_counter()
        try:
            h = hybrid(cpdag, g0, x, y, o, criterion=crit, time_limit_s=rung_s)
            out[f"{crit}_r"] = "inf" if h.radius == UNREACHED else h.radius
            out[f"{crit}_method"] = h.method
            out[f"{crit}_exact"] = h.exact
            out[f"{crit}_top_s"] = round(h.top_seconds, 4)
            out[f"{crit}_search_s"] = round(h.search_seconds, 4)
            out[f"{crit}_ladder_s"] = round(h.ladder_seconds, 4)
        except Exception as e:  # noqa: BLE001
            out[f"{crit}_r"] = f"error:{type(e).__name__}"
        out[f"{crit}_total_s"] = round(time.perf_counter() - t, 4)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--rung-s", type=float, default=120.0)
    a = ap.parse_args()
    todo = [r for r in rows(set(knowledge())) if r["hop_method"] == "local_up_fast_budget"]
    print(len(todo), "rows", flush=True)
    res = []
    with mp.Pool(a.workers) as p:
        for i, o in enumerate(p.imap_unordered(work, [(r, a.rung_s) for r in todo])):
            res.append(o)
            if i % 25 == 0:
                print(i, o.get("network"), o.get("gac_r"), o.get("gac_method"),
                      o.get("gac_total_s"), flush=True)
    path = OUT / "e1" / "lower_bound_rows.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = sorted({k for o in res for k in o})
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(res)
    print("wrote", path)
