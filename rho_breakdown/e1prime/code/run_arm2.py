"""
ARM 2 -- the finite-sample analyst simulation.

The pilot's deepest gap (CHRIS-MEETING.md:97): K is TRUE by construction, so
the analyst-facing quantity is never exercised. Here the expert is genuinely
wrong: K_asserted = K_true with r in {0,1,2} statements reversed.

The analyst is handed ONLY (C, K_asserted, Sigma_hat, n). Never tau, never D*.
tau enters exactly once, to adjudicate coverage -- as EVALUATION, not input.

    naive(n)      = est0_hat +- z*SE_n
    robust(rho,n) = [ min_m (est_m - z*SE_m) , max_m (est_m + z*SE_m) ]
                    over consistent+amenable members at distance <= rho, incl. m=0

A wide enough interval always covers, so the deliverable is the coverage-WIDTH
frontier, not coverage (PREREG.md 3d, AUDIT.md A6).

Usage: python run_arm2.py <arm1.json> <out.json> [n_scm] [reps] [n_workers]
"""
import json
import sys
from itertools import combinations
from multiprocessing import Pool

import numpy as np

from adjust import (cov_linear, optimal_adjustment_set, possibly_causal_paths,
                    total_effect_linear)
from graphs import (MeekFail, apply_background_knowledge, dag_to_cpdag,
                    random_dag, undirected_edges)
from scm import make_linear_iscm, sample_linear
from se import Z95, ols_with_se

N_GRID = [200, 1000, 5000, 20000]
R_VALUES = [0, 1, 2]
MAX_RHO = 4
Z = Z95


def regenerate(seed, p, deg):
    """Replay run_arm1.analyse_one's rng stream exactly (E3 gate pattern)."""
    rng = np.random.default_rng(seed)
    D = random_dag(p, deg, rng)
    C = dag_to_cpdag(D)
    U = undirected_edges(C)
    pairs = [(a, b) for a in range(p) for b in range(p)
             if a != b and len(possibly_causal_paths(D, a, b)) > 0]
    x, y = pairs[int(rng.integers(len(pairs)))]
    scm = make_linear_iscm(D, rng)
    Sigma = cov_linear(scm["A"], scm["omega"])
    tau = total_effect_linear(scm["A"], x, y)
    K_all = [(i, j) if D[i, j] == 1 else (j, i) for (i, j) in U]
    order = rng.permutation(len(K_all))
    K_all = [K_all[i] for i in order]
    return D, C, U, x, y, scm, Sigma, tau, K_all


def ball_structure(C, x, y, K, max_rho):
    """Which O* each member of the ball around K yields. Sigma-free, so it is
    computed ONCE and reused across every (rep, n)."""
    k = len(K)
    try:
        G0 = apply_background_knowledge(C, K)
    except MeekFail:
        return dict(caught_free=True)
    O0 = optimal_adjustment_set(G0, x, y)
    if O0 is None:
        return dict(caught_free=False, loud=True)
    out = dict(caught_free=False, loud=False, O0=sorted(O0), members=[])
    for rho in range(1, min(max_rho, k) + 1):
        for flip in combinations(range(k), rho):
            Kp = [(v, u) if i in flip else (u, v) for i, (u, v) in enumerate(K)]
            try:
                G = apply_background_knowledge(C, Kp)
            except MeekFail:
                out["members"].append(dict(rho=rho, ok=False))
                continue
            O = optimal_adjustment_set(G, x, y)
            out["members"].append(dict(rho=rho, ok=O is not None,
                                       O=sorted(O) if O is not None else None))
    return out


def intervals(Sig, x, y, struct, n):
    """naive + robust(rho) from Sigma-hat alone. Returns dict rho -> (lo, hi)."""
    e0, s0, _ = ols_with_se(Sig, x, y, struct["O0"], n)
    if not np.isfinite(e0) or not np.isfinite(s0):
        return None
    los = {0: e0 - Z * s0}
    his = {0: e0 + Z * s0}
    lo, hi = los[0], his[0]
    for rho in range(1, MAX_RHO + 1):
        for m in struct["members"]:
            if m["rho"] != rho or not m["ok"]:
                continue
            e, s, _ = ols_with_se(Sig, x, y, m["O"], n)
            if not np.isfinite(e) or not np.isfinite(s):
                continue
            lo = min(lo, e - Z * s); hi = max(hi, e + Z * s)
        los[rho], his[rho] = lo, hi
    return {r: (los[r], his[r]) for r in los}, e0, s0


def one_scm(args):
    rec, reps, seed_off = args
    seed, p, deg = rec["seed"], rec["p"], rec["deg"]
    D, C, U, x, y, scm, Sigma, tau, K_all = regenerate(seed, p, deg)
    if abs(tau - rec["tau"]) > 1e-9:
        return dict(seed=seed, gate_fail="tau mismatch")
    K_true = K_all[:4]
    if len(K_true) < 4:
        return None

    frng = np.random.default_rng(seed + seed_off)
    structs = {}
    for r in R_VALUES:
        flip = [] if r == 0 else sorted(frng.choice(len(K_true), size=r, replace=False).tolist())
        Ka = [(v, u) if i in flip else (u, v) for i, (u, v) in enumerate(K_true)]
        structs[r] = (flip, ball_structure(C, x, y, Ka, MAX_RHO))

    out = dict(seed=seed, p=p, deg=deg, tau=tau, gate_fail=None,
               caught_free={r: bool(structs[r][1].get("caught_free")) for r in R_VALUES},
               loud={r: bool(structs[r][1].get("loud")) for r in R_VALUES},
               cells=[])
    srng = np.random.default_rng(seed + 7_000_000 + seed_off)
    nmax = max(N_GRID)
    for rep in range(reps):
        Xd = sample_linear(scm, nmax, srng)
        for n in N_GRID:
            Sig = np.cov(Xd[:n], rowvar=False)
            for r in R_VALUES:
                st = structs[r][1]
                if st.get("caught_free") or st.get("loud"):
                    continue
                got = intervals(Sig, x, y, st, n)
                if got is None:
                    continue
                iv, e0, s0 = got
                out["cells"].append(dict(
                    rep=rep, n=n, r=r, est0=float(e0), se0=float(s0),
                    cov={str(k): bool(v[0] <= tau <= v[1]) for k, v in iv.items()},
                    w={str(k): float(v[1] - v[0]) for k, v in iv.items()}))
    return out


def main():
    arm1_path = sys.argv[1]
    out_path = sys.argv[2]
    n_scm = int(sys.argv[3]) if len(sys.argv) > 3 else 400
    reps = int(sys.argv[4]) if len(sys.argv) > 4 else 20
    nw = int(sys.argv[5]) if len(sys.argv) > 5 else 12
    src = json.load(open(arm1_path))[:n_scm]
    jobs = [(r, reps, 991) for r in src]
    with Pool(nw) as pool:
        res = [r for r in pool.imap(one_scm, jobs, chunksize=2) if r]
    fails = [r for r in res if r.get("gate_fail")]
    print(f"[arm2] {len(res)} SCMs, regeneration gate failures: {len(fails)}", flush=True)
    with open(out_path, "w") as f:
        json.dump(dict(n_grid=N_GRID, r_values=R_VALUES, reps=reps, z=Z,
                       gate_regen_failures=len(fails), scms=res), f)
    print(f"[arm2] -> {out_path}", flush=True)


if __name__ == "__main__":
    main()
