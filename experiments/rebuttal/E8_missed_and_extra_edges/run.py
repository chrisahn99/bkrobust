"""E8: missed and extra edges, how local is the damage. See PREREG.md (2026-11-06).

One unit = one SCM draw (ensemble, seed). Each unit writes one JSON line holding every single-edge
deletion and every single-edge addition of that SCM. The SCM stream and the deletion arm follow
the pilot tree, rho_breakdown/x3_skeleton/code/run_x3.py (commit 59066137), with fresh
seeds; the addition arm is new.

    .venv/bin/python experiments/rebuttal/E8_missed_and_extra_edges/run.py --smoke
    .venv/bin/python experiments/rebuttal/E8_missed_and_extra_edges/run.py --workers 8

Exit codes: 0 ok; 2 the two GAC oracles disagreed (run stopped); 3 a theory anomaly was recorded.
"""
import argparse
import json
import os
import subprocess
import sys
import threading
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))
if str(REPO / "src") not in sys.path:
    sys.path.insert(1, str(REPO / "src"))  # fallback when bkrobust is not pip-installed

from vendor_x3 import (MeekFail, PathBlowup, apply_background_knowledge,  # noqa: E402
                       cov_linear, dag_to_cpdag, hop_dist_from, is_amenable,
                       is_undirected, is_valid_adjustment_set, make_linear_iscm,
                       ols_coefficient, optimal_adjustment_set, possibly_causal_paths,
                       random_dag, skeleton, stmt_hop, total_effect_linear,
                       undirected_edges)

from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.gac import is_gac_valid_dag, is_gac_valid_mpdag  # noqa: E402

CRITERION = "gac"
VENDOR_COMMIT = "8486bb3509decf17adf490c078607b225d8df1b0"
MAX_K = 4
SEED_OFFSET = 1_000_000      # X3 used SCM seeds 0..15999; these never overlap
MASTER_SEED = 20261106
# grids unchanged from run_x3.py (itself from run_x1.py), so the population is X3's
ENSEMBLES = {
    "original": dict(ps=range(5, 9), degs=[1.5, 2.0, 2.5], need_u=3, n_draw=12000),
    "licensed": dict(ps=range(5, 11), degs=[1.5, 2.0, 2.5, 4.0, 6.0], need_u=3, n_draw=16000),
    "large": dict(ps=range(15, 26), degs=[1.5, 2.0, 2.5], need_u=3, n_draw=4800),
}
SMOKE_DRAWS = {"original": 24, "licensed": 24, "large": 2}
SMOKE_ADD_CAP = 12


class OracleDisagree(Exception):
    pass


def repo_commit():
    try:
        return subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True,
                              text=True, timeout=10).stdout.strip() or "unknown"
    except Exception:
        return "unknown"


# ------------------------------------------------------------------ two GAC oracles
def to_mpdag(G):
    p = G.shape[0]
    d = [(f"V{i}", f"V{j}") for i in range(p) for j in range(p) if G[i, j] and not G[j, i]]
    u = [(f"V{i}", f"V{j}") for i in range(p) for j in range(i + 1, p) if G[i, j] and G[j, i]]
    return MPDAG([f"V{i}" for i in range(p)], d, u)


def valid_in_dag(D, x, y, Z):
    """GAC validity in a DAG, decided twice; a disagreement invalidates the whole run."""
    a = bool(is_valid_adjustment_set(D, x, y, Z))
    b = bool(is_gac_valid_dag(to_mpdag(D), f"V{x}", f"V{y}", {f"V{z}" for z in Z}))
    if a != b:
        raise OracleDisagree(f"vendored={a} bkrobust={b} Z={sorted(Z)} x={x} y={y}")
    return a


def mpdag_says_valid(G, x, y, Z):
    return bool(is_gac_valid_mpdag(to_mpdag(G), f"V{x}", f"V{y}", {f"V{z}" for z in Z}))


# ------------------------------------------------------------------ one analyst graph
def outcome(G, D, Sigma, x, y, O0, tau):
    """Classify O*(G) against the true DAG D. Returns (record, O or None)."""
    try:
        O = optimal_adjustment_set(G, x, y)
    except PathBlowup:
        return dict(o="blowup"), None
    if O is None:
        return dict(o="not_amenable"), None
    if O == O0:
        return dict(o="unchanged", n_O=len(O)), O
    try:
        valid = valid_in_dag(D, x, y, O)
        est = float(ols_coefficient(Sigma, x, y, O))
    except (PathBlowup, np.linalg.LinAlgError):
        return dict(o="blowup"), None
    return dict(o="changed_valid" if valid else "silent_bias", n_O=len(O),
                rel_bias=float(abs(est - tau) / max(abs(tau), 1e-12))), O


def reachable(D, a):
    seen, stack = {a}, [a]
    while stack:
        v = stack.pop()
        for w in np.flatnonzero(D[v]):
            if int(w) not in seen:
                seen.add(int(w))
                stack.append(int(w))
    return seen


def addition_candidates(D):
    """Every absent pair with every orientation that keeps D acyclic: forced when one node is an
    ancestor of the other, free (both directions tested) otherwise."""
    p = D.shape[0]
    reach = [reachable(D, a) for a in range(p)]
    out = []
    for a in range(p):
        for b in range(a + 1, p):
            if D[a, b] or D[b, a]:
                continue
            if b in reach[a]:
                out.append((a, b, "forced"))
            elif a in reach[b]:
                out.append((b, a, "forced"))
            else:
                out += [(a, b, "free"), (b, a, "free")]
    return out


# ------------------------------------------------------------------ one unit
def analyse_one(job):
    ens, seed, p, deg, add_cap, mpdag_every = job
    cfg = ENSEMBLES[ens]
    base = dict(key=f"{ens}:{seed}", ens=ens, seed=int(seed), p=int(p), deg=float(deg),
                criterion=CRITERION, anomalies=[])
    try:
        return _analyse(base, cfg, seed, p, deg, add_cap, mpdag_every)
    except OracleDisagree as e:
        base.update(status="oracle_disagree", detail=str(e))
        return base
    except PathBlowup:
        # only the unguarded stream steps (query pairs, CPDAG amenability) reach here
        base.update(status="drop_blowup_stream", anomalies=[])
        return base


def _analyse(base, cfg, seed, p, deg, add_cap, mpdag_every):
    # ---- stream as run_x3.py analyse_one, unchanged through K ----------------------------
    rng = np.random.default_rng(seed)
    D = random_dag(p, deg, rng)
    C = dag_to_cpdag(D)
    U = undirected_edges(C)
    if len(U) < cfg["need_u"]:
        return dict(base, status="drop_few_undirected")
    pairs = [(a, b) for a in range(p) for b in range(p)
             if a != b and len(possibly_causal_paths(D, a, b)) > 0]
    if not pairs:
        return dict(base, status="drop_no_causal_pair")
    x, y = pairs[int(rng.integers(len(pairs)))]
    scm = make_linear_iscm(D, rng)
    Sigma = cov_linear(scm["A"], scm["omega"])
    tau = total_effect_linear(scm["A"], x, y)
    if abs(tau) < 1e-6:
        return dict(base, status="drop_tau_zero")
    K_all = [(i, j) if D[i, j] == 1 else (j, i) for (i, j) in U]
    order = rng.permutation(len(K_all))
    K_all = [K_all[i] for i in order]
    if len(K_all) < MAX_K:
        return dict(base, status="drop_few_K")
    K = [(int(a), int(b)) for (a, b) in K_all[:MAX_K]]
    base.update(x=int(x), y=int(y), tau=float(tau), n_edges=int(D.sum()), K=[list(e) for e in K])

    # ---- eligibility (X3 amendment 7) ----------------------------------------------------
    try:
        G0 = apply_background_knowledge(C, K)
    except MeekFail:
        return dict(base, status="drop_G0_meekfail")
    try:
        O0 = optimal_adjustment_set(G0, x, y)
        if O0 is None:
            return dict(base, status="drop_nonamenable")
        if not valid_in_dag(D, x, y, O0):
            return dict(base, status="drop_baseline_invalid")
    except PathBlowup:
        return dict(base, status="drop_blowup")
    O0_S0 = None
    cpdag_amenable = bool(is_amenable(C, x, y))
    if cpdag_amenable:
        try:
            O0_S0 = optimal_adjustment_set(C, x, y)
            if O0_S0 is not None and not valid_in_dag(D, x, y, O0_S0):
                O0_S0 = None
        except PathBlowup:
            O0_S0 = None

    check_mpdag = mpdag_every > 0 and seed % mpdag_every == 0
    n_mpdag = 0
    if check_mpdag:
        n_mpdag += 1
        if not mpdag_says_valid(G0, x, y, O0):
            base["anomalies"].append(dict(kind="mpdag_gac_rejects_O0", O=sorted(O0)))
    distC = hop_dist_from(skeleton(C), [x, y], p)

    def knowledge_arm(Chat, Kr, Dtrue_plus=None):
        """Arm SK on one analyst CPDAG. Dtrue_plus is set in the addition arm: it is in [Ghat]."""
        nonlocal n_mpdag
        try:
            Ghat = apply_background_knowledge(Chat, Kr) if Kr else Chat
        except MeekFail:
            if Dtrue_plus is not None:
                base["anomalies"].append(dict(kind="add_meekfail"))
            return dict(o="meekfail")
        rec, O = outcome(Ghat, D, Sigma, x, y, O0, tau)
        if O is not None and Dtrue_plus is not None:
            try:
                if not valid_in_dag(Dtrue_plus, x, y, O):
                    base["anomalies"].append(dict(kind="add_Ostar_invalid_in_D_plus", O=sorted(O)))
            except PathBlowup:
                pass
        if O is not None and check_mpdag:
            n_mpdag += 1
            if not mpdag_says_valid(Ghat, x, y, O):
                base["anomalies"].append(dict(kind="mpdag_gac_rejects_Ostar", O=sorted(O)))
        return rec

    def no_knowledge_arm(Chat):
        if O0_S0 is None:
            return dict(o="baseline_unusable")
        return outcome(Chat, D, Sigma, x, y, O0_S0, tau)[0]

    # ---- deletion arm: every edge of D ---------------------------------------------------
    dels = []
    for (i, j) in [(int(i), int(j)) for i in range(p) for j in range(p) if D[i, j] == 1]:
        Dp = D.copy()
        Dp[i, j] = 0
        Chat = dag_to_cpdag(Dp)
        Kr = [(a, b) for (a, b) in K if (Chat[a, b] or Chat[b, a]) and is_undirected(Chat, a, b)]
        try:
            sk = knowledge_arm(Chat, Kr)
        except PathBlowup:
            sk = dict(o="blowup")
        dels.append(dict(e=[i, j], hop=int(stmt_hop(distC, i, j)), w=float(abs(scm["A"][i, j])),
                         nK=len(Kr), SK=sk, S0=no_knowledge_arm(Chat)))

    # ---- addition arm: every absent pair, every acyclic orientation ----------------------
    cands = addition_candidates(D)
    n_cand = len(cands)
    if add_cap and n_cand > add_cap:
        pick = np.sort(np.random.default_rng(seed + 7).choice(n_cand, size=add_cap, replace=False))
        cands = [cands[k] for k in pick]
    adds = []
    for (a, b, kind) in cands:
        Dp = D.copy()
        Dp[a, b] = 1
        Chat = dag_to_cpdag(Dp)
        Kr = [(u, v) for (u, v) in K if is_undirected(Chat, u, v)]
        try:
            sk = knowledge_arm(Chat, Kr, Dtrue_plus=Dp)
        except PathBlowup:
            sk = dict(o="blowup")
        try:
            z0_ok = bool(valid_in_dag(Dp, x, y, O0))
        except PathBlowup:
            z0_ok = None
        adds.append(dict(e=[a, b], kind=kind, hop=int(stmt_hop(distC, a, b)), nK=len(Kr),
                         SK=sk, S0=no_knowledge_arm(Chat), Z0_valid_in_D_plus=z0_ok))

    return dict(base, status="ok", n_O0=len(O0), cpdag_amenable=cpdag_amenable,
                n_O0_S0=(len(O0_S0) if O0_S0 is not None else None),
                n_add_candidates=n_cand, n_add_tested=len(adds), mpdag_checked=n_mpdag,
                dels=dels, adds=adds)


# ------------------------------------------------------------------ driver
def build_jobs(smoke, add_cap, mpdag_every, scale):
    rng = np.random.default_rng(MASTER_SEED)
    jobs = []
    for ens, cfg in ENSEMBLES.items():
        ps = list(cfg["ps"])
        n_full = int(round(cfg["n_draw"] * scale))
        # draw (p, deg) for the full stream always, so a smoke unit is the same unit as in the full run
        draws = [(int(rng.choice(ps)), float(rng.choice(cfg["degs"]))) for _ in range(n_full)]
        n = SMOKE_DRAWS[ens] if smoke else n_full
        jobs += [(ens, SEED_OFFSET + s, p, d, add_cap, mpdag_every) for s, (p, d) in enumerate(draws[:n])]
    return jobs


def read_done(path):
    """Keys already written, plus the anomaly and disagreement counts they carry."""
    keys, n_anom, n_dis = set(), 0, 0
    if not path.exists():
        return keys, n_anom, n_dis
    with open(path) as f:
        for line in f:
            try:
                row = json.loads(line)
            except Exception:
                continue  # a line cut by a crash; its unit is simply redone
            keys.add(row["key"])
            n_anom += len(row.get("anomalies", []))
            n_dis += row.get("status") == "oracle_disagree"
    return keys, n_anom, n_dis


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--smoke", action="store_true", help="50 deterministic units, additions capped")
    ap.add_argument("--out", default=str(REPO / "results/rebuttal/E8"))
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--add-cap", type=int, default=None,
                    help="max additions per SCM, 0 = all (default: all; smoke: 12)")
    ap.add_argument("--mpdag-every", type=int, default=None,
                    help="MPDAG-level GAC check on SCMs with seed %% N == 0 (default 25; smoke 1)")
    ap.add_argument("--draw-scale", type=float, default=1.0, help="multiplies the X3 draw counts")
    args = ap.parse_args()
    if args.workers > 8:
        sys.exit("refusing more than 8 workers (15 GB machine)")
    if not (0 <= args.shard < args.n_shards):
        sys.exit("need 0 <= --shard < --n-shards")
    add_cap = args.add_cap if args.add_cap is not None else (SMOKE_ADD_CAP if args.smoke else 0)
    mpdag_every = args.mpdag_every if args.mpdag_every is not None else (1 if args.smoke else 25)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    tag = "smoke" if args.smoke else f"shard{args.shard:03d}of{args.n_shards:03d}"
    path = out / f"E8_{tag}.jsonl"
    commit = repo_commit()
    jobs = build_jobs(args.smoke, add_cap, mpdag_every, args.draw_scale)
    jobs = [j for k, j in enumerate(jobs) if k % args.n_shards == args.shard]
    seen, prior_anom, prior_dis = read_done(path)
    if prior_dis:
        sys.exit(print(f"[E8 {tag}] {path} already holds {prior_dis} oracle disagreement(s); "
                       "fix the harness and start a fresh --out") or 2)
    todo = [j for j in jobs if f"{j[0]}:{j[1]}" not in seen]

    config = dict(experiment="E8", criterion=CRITERION, smoke=bool(args.smoke), commit=commit,
                  vendor_commit=VENDOR_COMMIT, seed_offset=SEED_OFFSET, master_seed=MASTER_SEED,
                  max_k=MAX_K, ensembles=list(ENSEMBLES),
                  n_draw={e: (SMOKE_DRAWS[e] if args.smoke else int(round(c["n_draw"] * args.draw_scale)))
                          for e, c in ENSEMBLES.items()},
                  add_cap=add_cap, add_orientations="all_acyclic", deletions="all_edges",
                  mpdag_every=mpdag_every, shard=args.shard, n_shards=args.n_shards,
                  oracles=["vendor_x3.is_valid_adjustment_set", "bkrobust.gac.is_gac_valid_dag"])
    (out / f"config_{tag}.json").write_text(json.dumps(config, indent=1) + "\n")
    print(f"[E8 {tag}] units={len(jobs)} already_done={len(jobs) - len(todo)} todo={len(todo)} "
          f"workers={args.workers} commit={commit[:10]}", flush=True)

    # a crash can leave a line without its newline; never glue the next record onto it
    if path.exists() and path.stat().st_size and not path.read_bytes().endswith(b"\n"):
        with open(path, "a") as f:
            f.write("\n")

    t0 = time.time()
    state = dict(n=0, ok=0, anomalies=prior_anom)
    stop = threading.Event()

    def heartbeat():
        while not stop.wait(60):
            print(f"[E8 {tag}] {state['n']}/{len(todo)} units ok={state['ok']} "
                  f"anomalies={state['anomalies']} {time.time() - t0:.0f}s", flush=True)

    code = 0
    pool = Pool(args.workers) if args.workers > 1 else None  # fork before the heartbeat thread exists
    threading.Thread(target=heartbeat, daemon=True).start()
    it = pool.imap_unordered(analyse_one, todo, chunksize=1) if pool else map(analyse_one, todo)
    try:
        with open(path, "a") as f:
            for row in it:
                row["commit"] = commit
                f.write(json.dumps(row) + "\n")
                f.flush()
                os.fsync(f.fileno())
                state["n"] += 1
                state["ok"] += row.get("status") == "ok"
                state["anomalies"] += len(row.get("anomalies", []))
                if row.get("status") == "oracle_disagree":
                    print(f"[E8 {tag}] CONSISTENCY FAILURE: two GAC oracles disagree on "
                          f"{row['key']}: {row['detail']}", flush=True)
                    code = 2
                    break
    finally:
        stop.set()
        if pool:
            pool.terminate() if code else pool.close()
            pool.join()
    if code == 0 and state["anomalies"]:
        print(f"[E8 {tag}] THEORY ANOMALIES recorded: {state['anomalies']} (see 'anomalies' in rows)",
              flush=True)
        code = 3
    print(f"[E8 {tag}] finished {state['n']} units in {time.time() - t0:.0f}s, exit {code}", flush=True)
    sys.exit(code)


if __name__ == "__main__":
    main()
