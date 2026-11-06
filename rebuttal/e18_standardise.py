"""E18: Section 7 effect-scale profiles under standardised models.

Instances: the 57 effect-scale queries of the paper (39 from
results/mean_table_relevant_v2/x_component, 18 from .../extended), D_LLM
knowledge, the paper's profile construction (experiments/
mean_table_xcomp_extend_v2.build_record, X-component slice).

Arms, each on coefficient draws 0..D-1 (draw 0 = the paper's seed):
  naive     Sigma of random_sem (|w| ~ U(0.3, 1.5), random sign,
            noise variance ~ U(0.5, 1.5)), as in the paper;
  posthoc   the same model with every variable rescaled to unit variance
            (Sigma' = D^-1/2 Sigma D^-1/2);
  internal  internally standardised SEM (each variable standardised as it is
            generated, before its children: Ormaniec et al. 2025 iSCM).
Recorded per instance/arm/draw: theta_Z, the r_mu signature over the six
tolerances, mu peak, the ceiling B(C_hat)/|theta_Z|, var-sortability.
Checks: naive draw 0 reproduces the stored record (1e-9); posthoc equals naive
on every relative quantity (any difference is a bug). r_val is graphical and
is not recomputed.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import signal
import sys

import numpy as np

from common import OUT, ROOT

sys.path.insert(0, str(ROOT / "experiments"))
import final_table as ft  # noqa: E402
import mean_table_xcomp_extend_v2 as xc  # noqa: E402

from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag, random_sem  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.epsilon.bias import make_context_from_covariance  # noqa: E402

EPS = [0.02, 0.10, 0.25, 0.40, 0.60, 0.90]
SEED = 20260917
HEAVY_S = 30.0


class Kill(Exception):
    pass


def _alarm(*_):
    raise Kill


def topo(dag):
    indeg = {v: len(dag.parents(v)) for v in dag.nodes}
    order, stack = [], sorted(v for v, d in indeg.items() if d == 0)
    while stack:
        v = stack.pop(0)
        order.append(v)
        for c in sorted(dag.children(v)):
            indeg[c] -= 1
            if indeg[c] == 0:
                stack.append(c)
    return order


def cov_from_sem(sem, nodes, internal=False):
    """Covariance of a linear SEM, optionally internally standardised."""
    idx = {v: i for i, v in enumerate(nodes)}
    n = len(nodes)
    S = np.zeros((n, n))
    for v in topo(sem.dag):
        j = idx[v]
        pa = sorted(sem.dag.parents(v))
        w = np.array([sem.weights[(p, v)] for p in pa])
        pi = [idx[p] for p in pa]
        noise = sem.noise_var[v]
        cov_vk = S[pi, :].T @ w if pa else np.zeros(n)   # Cov(v, k) for earlier k
        var_v = (w @ S[np.ix_(pi, pi)] @ w if pa else 0.0) + noise
        if internal:
            s = np.sqrt(var_v)
            cov_vk = cov_vk / s
            var_v = 1.0
        S[j, :] = cov_vk
        S[:, j] = cov_vk
        S[j, j] = var_v
    return S


def var_sortability(S, nodes, dag):
    idx = {v: i for i, v in enumerate(nodes)}
    e = list(dag.directed_edges)
    if not e:
        return None
    up = sum(S[idx[b], idx[b]] > S[idx[a], idx[a]] + 1e-12 for a, b in e)
    eq = sum(abs(S[idx[b], idx[b]] - S[idx[a], idx[a]]) <= 1e-12 for a, b in e)
    return (up + 0.5 * eq) / len(e)


def work(task):
    signal.signal(signal.SIGALRM, _alarm)
    inst, draws, budget = task
    net, x, y = inst["network"], inst["x"], inst["y"]
    parsed = ft.load_networks(names={net}, skip=set(), max_nodes=0)
    knowledge = ft.load_knowledge("D_LLM")
    cpdag, dag = parsed[net]["cpdag"], parsed[net]["dag"]
    g0 = apply_orientations(cpdag, list(knowledge[net]))
    z = frozenset(optimal_adjustment_set_mpdag(g0, x, y))
    nodes = list(cpdag.nodes)
    out = []
    for dr in range(draws):
        seed = ft.derive_seed(SEED if dr == 0 else SEED + 7919 * dr, net, x, y)
        sem = random_sem(dag, np.random.default_rng(seed))
        S_raw = cov_from_sem(sem, nodes)
        sd = np.sqrt(np.diag(S_raw))
        arms = {"naive": S_raw, "posthoc": S_raw / np.outer(sd, sd),
                "internal": cov_from_sem(sem, nodes, internal=True)}
        for arm, S in arms.items():
            ctx = make_context_from_covariance(S, nodes, cpdag, x, y, z)
            scale = abs(ctx.theta_z)
            row = {"network": net, "x": x, "y": y, "r_val": inst["r_val"], "draw": dr,
                   "arm": arm, "theta_z": ctx.theta_z,
                   "var_sortability": var_sortability(S, nodes, dag)}
            if scale < 1e-12:
                out.append(dict(row, status="zero_estimate"))
                continue
            signal.alarm(int(budget))
            try:
                rec = xc.build_record(net, x, y, inst["r_val"], None, cpdag, dag, g0, ctx,
                                      scale, EPS)
                signal.alarm(0)
            except Kill:
                out.append(dict(row, status="timeout"))
                continue
            row.update(status="ok", r_mu={str(e): rec["r_mu_rel"][str(e)] for e in EPS},
                       r_mu_persists={str(e): rec["r_mu_rel_persists"][str(e)] for e in EPS},
                       mu_peak=rec["mu_peak_rel"], b_top=rec["b_top"],
                       max_d=rec["max_distance_rel"],
                       profile=[round(p["mu"], 12) for p in rec["profile_rel"]])
            out.append(row)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--draws", type=int, default=10)
    ap.add_argument("--heavy-draws", type=int, default=3)
    a = ap.parse_args()
    insts = []
    for f in ("x_component", "extended"):
        for line in open(ROOT / f"results/mean_table_relevant_v2/{f}/instances.jsonl"):
            r = json.loads(line)
            if r.get("status") == "ok":
                r["_heavy"] = (r.get("seconds") or 0) > HEAVY_S
                r["_src"] = f
                insts.append(r)
    print(len(insts), "instances", flush=True)
    tasks = [(r, a.heavy_draws if r["_heavy"] else a.draws, 1800 if r["_heavy"] else 300)
             for r in insts]
    tasks.sort(key=lambda t: -int(t[0]["_heavy"]))
    res = []
    with mp.Pool(a.workers) as P:
        for o in P.imap_unordered(work, tasks):
            res.extend(o)
            print(o[0]["network"], o[0]["x"], o[0]["y"], len(o), flush=True)
    path = OUT / "e18" / "standardise.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for r in res:
            f.write(json.dumps(r) + "\n")
    stored = {(r["network"], r["x"], r["y"]): r for r in insts}
    print("wrote", path)
