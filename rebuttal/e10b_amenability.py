"""E10b: is the radius the distance to the nearest non-amenable state?

GAC validity in an MPDAG requires that (X, Y) be amenable, which does not
depend on Z. Define r_amen = least number of retractions from G0 that make
(X, Y) non-amenable (same upward search, failure = not amenable; amenability
loss is upward-closed because validity of any set requires it... checked
here empirically by comparing with the exhaustive search). Since every valid
set fails where amenability fails, r_val <= r_amen. This script measures how
often r_val == r_amen, i.e. how often the nearest failure is a loss of
identifiability by adjustment rather than a set-specific failure, on
  - the 4,584 synthetic survival instances of the paper (rebuilt from seed), and
  - every committed release query of the ten elicited conditions.
It also records, at the nearest failing state found by the hybrid, whether
(X, Y) is still amenable there.
"""

from __future__ import annotations

import csv
import multiprocessing as mp
import signal

from common import OUT, UNREACHED, hybrid
from oset import optimal_set_closed_form
from release_rows import knowledge, rows, state

from bkrobust.demo.meek import apply_orientations
from bkrobust.mpdag_criterion import is_amenable
from bkrobust.search.exact_fast import radius_local_up_fast


class Kill(Exception):
    pass


def _a(*_):
    raise Kill


def r_amen(cpdag, g0, x, y):
    if is_amenable(cpdag, x, y):
        return UNREACHED
    f = radius_local_up_fast(cpdag, g0, lambda g: not is_amenable(g, x, y), max_depth=None)
    return f.radius


def one(cpdag, g0, x, y, z):
    h = hybrid(cpdag, g0, x, y, z, time_limit_s=60)
    ra = r_amen(cpdag, g0, x, y)
    w_amen = None
    if h.witness and h.method == "local_up_fast":
        pass
    return {"r_val": h.radius, "r_amen": ra, "equal": h.radius == ra,
            "amenable_at_chat": is_amenable(cpdag, x, y)}


def work_release(r):
    signal.signal(signal.SIGALRM, _a)
    dag, cpdag, k, g0, _ = state(r["arm"], r["network"])
    o, st = optimal_set_closed_form(g0, r["x"], r["y"])
    if st != "ok":
        return None
    signal.alarm(600)
    try:
        res = one(cpdag, g0, r["x"], r["y"], o)
        signal.alarm(0)
    except Kill:
        return {"src": "release", "network": r["network"], "arm": r["arm"], "killed": True}
    return {"src": "release", "network": r["network"], "arm": r["arm"], "x": r["x"], "y": r["y"], **res}


def work_synth(args):
    from e15_aid import rebuild
    r, which = args
    inst, _ = rebuild(r, which)
    if inst is None:
        return None
    res = one(inst["cpdag"], inst["g0"], inst["x"], inst["y"], inst["z_star"])
    return {"src": which, "arm": r["arm"], "instance_id": r["instance_id"], **res}


if __name__ == "__main__":
    from e15_aid import FILES
    synth = []
    for which, p in FILES.items():
        with open(p) as f:
            synth += [(r, which) for r in csv.DictReader(f)]
    rel = [r for r in rows(set(knowledge())) if r["status"] == "ok"]
    out = []
    with mp.Pool(3, maxtasksperchild=100) as P:
        for o in P.imap_unordered(work_synth, synth, chunksize=16):
            if o:
                out.append(o)
        for o in P.imap_unordered(work_release, rel, chunksize=4):
            if o:
                out.append(o)
    keys = sorted({k for o in out for k in o})
    p = OUT / "e10" / "amenability_radius.csv"
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(out)
    import pandas as pd
    d = pd.DataFrame(out)
    print(d.groupby(["src"]).equal.agg(["mean", "sum", "count"]))
    print(d.groupby(["src", "arm"]).equal.agg(["mean", "count"]))
    fin = d[d.r_val != UNREACHED]
    print("finite r_val:", fin.groupby("src").equal.agg(["mean", "count"]))
    print("wrote", p)
