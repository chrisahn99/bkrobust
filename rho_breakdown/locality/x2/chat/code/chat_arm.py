"""The estimated-CPDAG arm: does the firewall survive when C is learned from data?

Every previous experiment used the ORACLE CPDAG, dag_to_cpdag(D). A practitioner runs PC on n
rows and gets Chat != C. This arm re-runs the firewall adjudication on Chat, and adds the
comparison the oracle arm cannot make:

  BASELINE  -- with true BK and NO false statement at all, is the reported O* already invalid
               in D purely because discovery erred?

That baseline is the honest denominator. If discovery error alone invalidates a large share of
reports, a 0% vs 40.4% contrast attributable to background knowledge is swamped and must be
reported as such.

Statement truth is labelled against D, not against Chat: a pair non-adjacent in Chat may be
adjacent in D (a PC false negative), in which case the asserted edge is TRUE and belongs in a
separate cell.
"""
import sys, os, json, time
from collections import defaultdict
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
X2 = os.environ.get("X2CODE", os.path.join(HERE, "..", "..", "code"))
X1 = os.environ.get("X1CODE", os.path.join(HERE, "..", "..", "..", "x1", "code"))
sys.path.insert(0, X2)
import importlib.util
from graphs import (random_dag, dag_to_cpdag, meek_closure, v_structures,
                    has_directed_cycle, undirected_edges, skeleton, is_undirected)
from scm import make_linear_iscm, sample_linear
import adjust
import x2lib as X
from pc import pc, ci_fisher_z, ci_oracle, test_population_recovers_cpdag

_s = importlib.util.spec_from_file_location("x1_ops", os.path.join(X1, "x1_ops.py"))
x1_ops = importlib.util.module_from_spec(_s); _s.loader.exec_module(x1_ops)

POS, ZERO, REFUSE = "POS", "ZERO", "REFUSE"


def report_state(G, x, y):
    O, npaths, amen = X.ostar_and_paths(G, x, y)
    if npaths == 0:
        return ZERO, None
    if not amen:
        return REFUSE, None
    return POS, O


def valid_in_D(D, x, y, Z):
    return adjust.is_valid_adjustment_set(D, x, y, set(Z))


def true_bk_on(Ghat, D, rng):
    """Orient a random subset of Ghat's undirected edges the way D says.
    Edges of Ghat absent from D have no true orientation and are LEFT ALONE (counted)."""
    G = Ghat.copy()
    U = undirected_edges(G)
    n_skipped = 0
    if not U:
        return G, 0
    k = int(rng.integers(0, len(U) + 1))
    for t in rng.permutation(len(U))[:k]:
        u, v = U[t]
        if not is_undirected(G, u, v):
            continue
        if D[u, v] == 1:
            G[v, u] = 0
        elif D[v, u] == 1:
            G[u, v] = 0
        else:
            n_skipped += 1
            continue
        G = meek_closure(G)
    return G, n_skipped


def run(n_scm, ps, degs, ns, alpha, seed0, oracle_mode=False):
    rng = np.random.default_rng(seed0)
    acc = defaultdict(int)
    per_scm = []
    for _ in range(n_scm):
        p = int(rng.choice(ps)); deg = float(rng.choice(degs)); n = int(rng.choice(ns))
        D = random_dag(p, deg, rng)
        scm = make_linear_iscm(D, rng)
        A, om = scm["A"], scm["omega"]
        M = np.linalg.inv(np.eye(p) - A.T)
        Sig = M @ np.diag(om) @ M.T
        if oracle_mode:
            Chat, _ = pc(Sig, p, ci_oracle)
        else:
            Xd = sample_linear(scm, n, rng)
            Shat = np.cov(Xd, rowvar=False)
            Chat, _ = pc(Shat, p, ci_fisher_z, n=n, alpha=alpha)
        C = dag_to_cpdag(D)
        shd = int((skeleton(Chat) != skeleton(C)).sum() // 2)
        loc = defaultdict(int)
        loc["shd_skel"] = shd
        loc["n_scm"] = 1
        G0, nskip = true_bk_on(Chat, D, rng)
        loc["bk_edges_absent_in_D"] = nskip
        ref = v_structures(Chat)
        NA = [(a, b) for a in range(p) for b in range(a + 1, p) if skeleton(Chat)[a, b] == 0]
        for x in range(p):
            for y in range(p):
                if x == y:
                    continue
                s0, O0 = report_state(G0, x, y)
                if s0 != POS:
                    continue
                # ---- BASELINE: no false statement at all, only discovery error + true BK
                loc["base", "n"] += 1
                base_ok = valid_in_D(D, x, y, O0)
                if not base_ok:
                    loc["base", "invalid"] += 1
                bs = "BV" if base_ok else "BI"
                for (a0, b0) in NA:
                    for (a, b) in ((a0, b0), (b0, a0)):
                        true_stmt = bool(D[a, b] == 1)
                        H, info = x1_ops.bk_assert(G0, [(a, b)])
                        ok = ((not info["conflict"]) and (not has_directed_cycle(H))
                              and x1_ops.pdag_extendable(H))
                        tag = "S1" if ok else "S1rej"
                        tt = "T" if true_stmt else "F"
                        s1, O1 = report_state(H, x, y)
                        if s1 != POS:
                            loc[tag, tt, bs, "exit"] += 1
                            continue
                        mv = (O1 != O0)
                        sd = valid_in_D(D, x, y, O1)
                        loc[tag, tt, bs, "mv" if mv else "sm", "V" if sd else "I"] += 1
        per_scm.append(dict(loc))
        for k2, v in loc.items():
            acc[k2] += v
    return acc, per_scm


def boot(per, tag, tt, bs="BV", B=2000, seed=0):
    rng = np.random.default_rng(seed)
    def rate(sample):
        mv = sum(d.get((tag, tt, bs, "mv", "V"), 0) + d.get((tag, tt, bs, "mv", "I"), 0) for d in sample)
        iv = sum(d.get((tag, tt, bs, "mv", "I"), 0) for d in sample)
        return iv / mv if mv else float("nan")
    pt = rate(per)
    draws = []
    for _ in range(B):
        idx = rng.integers(0, len(per), len(per))
        r = rate([per[i] for i in idx])
        if r == r:
            draws.append(r)
    if not draws:
        return pt, float("nan"), float("nan")
    return pt, float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))


if __name__ == "__main__":
    bad, tot = test_population_recovers_cpdag(60)
    print(f"[gate] population-PC == dag_to_cpdag: {tot-bad}/{tot} -> {'PASS' if bad==0 else 'FAIL'}")
    if bad:
        sys.exit(1)
    n_scm = int(sys.argv[1]); ps = [int(z) for z in sys.argv[2].split(",")]
    ns = [int(z) for z in sys.argv[3].split(",")]; seed = int(sys.argv[4])
    alpha = float(sys.argv[5]) if len(sys.argv) > 5 else 0.01
    degs = [float(z) for z in sys.argv[6].split(",")] if len(sys.argv) > 6 else [1.5, 2.0, 2.5]
    oracle = len(sys.argv) > 7 and sys.argv[7] == "oracle"
    t = time.time()
    acc, per = run(n_scm, ps, degs, ns, alpha, seed, oracle_mode=oracle)
    el = round(time.time() - t, 1)
    print(f"elapsed {el}s  n_scm={len(per)}  mode={'ORACLE-CI' if oracle else f'PC alpha={alpha}'}  n={ns}")
    nq = sum(d.get(("base", "n"), 0) for d in per)
    bi = sum(d.get(("base", "invalid"), 0) for d in per)
    msh = np.mean([d.get("shd_skel", 0) for d in per])
    print(f"skeleton SHD(Chat,C) mean = {msh:.3f}")
    print(f"BASELINE (true BK, no false statement): {bi}/{nq} = {bi/max(nq,1):.4f} invalid")
    print("\n--- CONDITIONED ON THE BASELINE REPORT BEING VALID (the firewall test) ---")
    for tag in ("S1", "S1rej"):
        for tt in ("F", "T"):
            mv = sum(d.get((tag, tt, "BV", "mv", "V"), 0) + d.get((tag, tt, "BV", "mv", "I"), 0) for d in per)
            iv = sum(d.get((tag, tt, "BV", "mv", "I"), 0) for d in per)
            ex = sum(d.get((tag, tt, "BV", "exit"), 0) for d in per)
            if mv == 0 and ex == 0:
                continue
            pt, lo, hi = boot(per, tag, tt, "BV")
            lbl = "FALSE stmt" if tt == "F" else "TRUE stmt"
            print(f"{tag:6s} {lbl}: moved={mv:7d} newly_invalid={iv:6d} rate={pt:.4f} "
                  f"CI95=[{lo:.4f},{hi:.4f}]  visible_abort={ex}")
    print("--- baseline ALREADY invalid (discovery error; BK cannot be blamed) ---")
    for tag in ("S1", "S1rej"):
        mv = sum(d.get((tag, "F", "BI", "mv", "V"), 0) + d.get((tag, "F", "BI", "mv", "I"), 0) for d in per)
        iv = sum(d.get((tag, "F", "BI", "mv", "I"), 0) for d in per)
        if mv:
            print(f"{tag:6s} FALSE stmt: moved={mv:7d} still_invalid={iv:6d} rate={iv/mv:.4f}")
    out = {"|".join(map(str, k)) if isinstance(k, tuple) else k: v for k, v in acc.items()}
    os.makedirs(os.path.join(HERE, "..", "results"), exist_ok=True)
    fn = os.path.join(HERE, "..", "results",
                      f"chat_{'oracle' if oracle else 'pc'}_p{'_'.join(map(str,ps))}_n{'_'.join(map(str,ns))}_s{seed}.json")
    json.dump({"acc": out, "n_scm": len(per), "alpha": alpha, "ns": ns, "ps": ps,
               "oracle_mode": oracle, "elapsed_s": el}, open(fn, "w"), indent=1)
    print("->", fn)
