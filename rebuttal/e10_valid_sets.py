"""E10: how much the radius varies across the adjustment sets valid at G0.

Queries: every release row of the ten elicited conditions with status ok
(states rebuilt from the release, K hash checked). Candidate pool C = the
possible ancestors of {X, Y} at G0 (reverse reachability over directed and
undirected edges, a superset) minus forb(X, Y, G0) and {X, Y}. If |C| <= 10,
every subset of C is tested (exhaustive within the ancestor family); otherwise
O, Adjust = C, pa(X) and seeded random subsets of C, up to 20 valid sets.
Each valid set gets its radius (submitted hybrid, GAC) and its validity in the
data-generating DAG. Reported per query: number of valid sets, radius spread,
whether some set has a larger radius than O ("O beaten"), O's rank.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import multiprocessing as mp
import signal

import numpy as np

from common import OUT, UNREACHED, hybrid, is_gac_valid_mpdag
from oset import optimal_set_closed_form
from release_rows import knowledge, rows, state

from bkrobust.gac.mpdag_level import gac_forbidden_set

EXHAUSTIVE = 10
MAX_SAMPLED = 20


class Kill(Exception):
    pass


def _alarm(*_):
    raise Kill


def poss_anc(g, targets):
    seen = set(targets)
    stack = list(targets)
    while stack:
        v = stack.pop()
        for u in g.parents(v) | g.neighbors(v):
            if u not in seen:
                seen.add(u)
                stack.append(u)
    return seen


def rkey(r):
    return float("inf") if r == UNREACHED else r


def work(r):
    signal.signal(signal.SIGALRM, _alarm)
    dag, cpdag, k, g0, _ = state(r["arm"], r["network"])
    x, y = r["x"], r["y"]
    o, st = optimal_set_closed_form(g0, x, y)
    base = {"arm": r["arm"], "network": r["network"], "x": x, "y": y}
    if st != "ok" or not is_gac_valid_mpdag(g0, x, y, o):
        return dict(base, status="no_O")
    forb = set(gac_forbidden_set(g0, x, y)) | {x, y}
    cand = sorted(poss_anc(g0, {x, y}) - forb)
    sets = set()
    mode = "exhaustive" if len(cand) <= EXHAUSTIVE else "sampled"
    if mode == "exhaustive":
        for m in range(len(cand) + 1):
            for c in itertools.combinations(cand, m):
                z = frozenset(c)
                if is_gac_valid_mpdag(g0, x, y, z):
                    sets.add(z)
    else:
        rng = np.random.default_rng(int(hashlib.sha256(f"{r['network']}|{x}|{y}".encode()).hexdigest()[:8], 16))
        pre = [o, frozenset(cand), frozenset(g0.parents(x)) - forb]
        for z in pre:
            if is_gac_valid_mpdag(g0, x, y, z):
                sets.add(z)
        tries = 0
        while len(sets) < MAX_SAMPLED and tries < 400:
            tries += 1
            keep = rng.random(len(cand)) < rng.uniform(0.2, 0.9)
            z = frozenset(c for c, kk in zip(cand, keep) if kk) | o
            if is_gac_valid_mpdag(g0, x, y, z):
                sets.add(z)
    sets.add(o)
    radii = {}
    signal.alarm(900)
    try:
        for z in sorted(sets, key=lambda s: (len(s), sorted(s))):
            h = hybrid(cpdag, g0, x, y, z, time_limit_s=60)
            radii[z] = (h.radius, h.exact, is_gac_valid_mpdag(dag, x, y, z))
        signal.alarm(0)
    except Kill:
        return dict(base, status="killed", n_sets=len(sets), mode=mode)
    ro = rkey(radii[o][0])
    vals = [rkey(v[0]) for v in radii.values()]
    finite = [v for v in vals if v != float("inf")]
    return dict(base, status="ok", mode=mode, n_cand=len(cand), n_sets=len(sets),
                r_O=radii[o][0], r_min=min(vals), r_max=max(vals),
                n_distinct=len(set(vals)), varies=len(set(vals)) > 1,
                O_beaten=any(v > ro for v in vals),
                O_is_max=ro == max(vals), O_rank_desc=1 + sum(v > ro for v in vals),
                n_inf=sum(v == float("inf") for v in vals), n_finite=len(finite),
                all_exact=all(v[1] for v in radii.values()),
                O_valid_truth=radii[o][2],
                n_valid_truth=sum(v[2] for v in radii.values()),
                maxr_set_valid_truth=max(radii.items(), key=lambda kv: rkey(kv[1][0]))[1][2],
                share_valid_truth_by_r=";".join(
                    f"{rk}:{sum(1 for v in radii.values() if rkey(v[0]) == rk and v[2])}/"
                    f"{sum(1 for v in radii.values() if rkey(v[0]) == rk)}"
                    for rk in sorted(set(vals))))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    todo = [r for r in rows(set(knowledge())) if r["status"] == "ok"]
    print(len(todo), "queries", flush=True)
    res = []
    with mp.Pool(a.workers, maxtasksperchild=50) as P:
        for i, o in enumerate(P.imap_unordered(work, todo, chunksize=4)):
            res.append(o)
            if i % 200 == 0:
                print(i, o.get("network"), o.get("status"), o.get("n_sets"), flush=True)
    keys = sorted({k for o in res for k in o})
    path = OUT / "e10" / "valid_sets.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(res)
    print("wrote", path)
