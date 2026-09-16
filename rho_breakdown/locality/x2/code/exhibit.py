"""
X2 COUNTEREXAMPLE DOSSIER -- PREREG 4.5.

A counterexample is not a row in a table.  This script emits, for every
rho >= 2 ALL-FLIPS-FAR O*-change found in the E1' archive:

  1. (seed, p, deg, ensemble) + the one-line regeneration command
  2. C and D as edge lists
  3. (X, Y), delta(.) for every vertex, the flipped statements and their dmin/dmax
  4. O*(G0), O*(G'), and beta for each
  5. the MEEK RULE CHAIN: the ordered (rule, a, b) orientations fired by the
     closure -- traced by an instrumented copy of graphs.meek_closure that is
     CHECKED against the original on every call
  6. validity of O*(G') in D and |beta(O') - tau| / |tau|
  7. minimality: the smallest p at which the phenomenon occurs in the searched
     ensembles, plus an exhaustive rho=2 far-pair search at p <= 5.

Usage: python exhibit.py <arch_dir> <out.json> [n_exhibits]
"""
import json
import sys
from collections import defaultdict
from itertools import combinations, product

import numpy as np

import x2lib as X
from adjust import cov_linear, is_valid_adjustment_set, total_effect_linear
from graphs import (MeekFail, adjacent, apply_background_knowledge,
                    dag_to_cpdag, directed_edges, has_directed_cycle,
                    is_directed, is_undirected, meek_closure, random_dag,
                    skeleton, undirected_edges, v_structures)
from scm import make_linear_iscm
from se import ols_with_se

ARCH = sys.argv[1] if len(sys.argv) > 1 else "/home/costaj/latent-causal/e1prime-se/results"
OUT = sys.argv[2] if len(sys.argv) > 2 else "../results/exhibit.json"
NEX = int(sys.argv[3]) if len(sys.argv) > 3 else 3
TAU_FLOOR = 1e-3
CELLS = [("original", "K4", 4), ("licensed", "K4", 4), ("large", "K4", 4),
         ("k8", "K4", 3), ("k8", "K6", 3), ("k8", "K8", 3)]


# ------------------------------------------------------------- traced closure
def meek_closure_traced(G):
    """Byte-for-byte the logic of graphs.meek_closure, plus a trace of which
    rule oriented which edge.  Verified equal to meek_closure on every call."""
    G = G.copy()
    p = G.shape[0]
    trace = []
    changed = True
    while changed:
        changed = False
        for a, b in [(a, b) for a in range(p) for b in range(p) if is_undirected(G, a, b)]:
            if not is_undirected(G, a, b):
                continue
            rule = None
            for c in range(p):
                if c in (a, b):
                    continue
                if is_directed(G, c, a) and not adjacent(G, c, b):
                    rule = ("R1", int(c))
                    break
            if rule is None:
                for c in range(p):
                    if c in (a, b):
                        continue
                    if is_directed(G, a, c) and is_directed(G, c, b):
                        rule = ("R2", int(c))
                        break
            if rule is None:
                cand = [c for c in range(p) if c not in (a, b)
                        and is_undirected(G, a, c) and is_directed(G, c, b)]
                for c, d in combinations(cand, 2):
                    if not adjacent(G, c, d):
                        rule = ("R3", int(c), int(d))
                        break
            if rule is None:
                for d in range(p):
                    if d in (a, b) or not is_undirected(G, a, d) or adjacent(G, b, d):
                        continue
                    for c in range(p):
                        if c in (a, b, d):
                            continue
                        if is_directed(G, d, c) and is_directed(G, c, b) and adjacent(G, a, c):
                            rule = ("R4", int(c), int(d))
                            break
                    if rule is not None:
                        break
            if rule is not None:
                G[b, a] = 0
                changed = True
                trace.append(dict(rule=rule[0], oriented=[int(a), int(b)],
                                  witness=[int(z) for z in rule[1:]]))
    return G, trace


def apply_bk_traced(C, K):
    """apply_background_knowledge with a per-statement Meek trace."""
    G = C.copy()
    chain = []
    for (i, j) in K:
        if is_directed(G, j, i):
            raise MeekFail("conflict")
        if G[i, j] == 0 and G[j, i] == 0:
            raise MeekFail("not an edge")
        G[j, i] = 0
        G, tr = meek_closure_traced(G)
        chain.append(dict(asserted=[int(i), int(j)], propagated=tr))
    if has_directed_cycle(G) or v_structures(G) != v_structures(C):
        raise MeekFail("guard")
    return G, chain


def edge_list(G):
    p = G.shape[0]
    out = []
    for i in range(p):
        for j in range(i + 1, p):
            if is_directed(G, i, j):
                out.append(f"{i}->{j}")
            elif is_directed(G, j, i):
                out.append(f"{j}->{i}")
            elif is_undirected(G, i, j):
                out.append(f"{i}--{j}")
    return out


def flips_of(k, max_rho):
    return [fl for rho in range(1, min(max_rho, k) + 1)
            for fl in combinations(range(k), rho)]


# ------------------------------------------------------------- find them all
def find_all():
    found = []
    for ens, key, max_rho in CELLS:
        recs = json.load(open(f"{ARCH}/arm1_{ens}.json"))
        for rec in recs:
            if abs(rec["tau"]) < TAU_FLOOR:
                continue
            arm = rec["arms"].get(key)
            if arm is None or not arm.get("mpdag_amenable"):
                continue
            hop = [float(h) if h < X.SENTINEL else np.inf for h in arm["stmt_hopdist"]]
            for flip, m in zip(flips_of(arm["n_K"], max_rho), arm["members"]):
                if len(flip) < 2 or not m.get("consistent") or not m.get("amenable"):
                    continue
                if not m.get("ostar_changed"):
                    continue
                hs = [hop[i] for i in flip]
                if not all(np.isfinite(h) and h >= 1 for h in hs):
                    continue
                found.append(dict(ens=ens, arm=key, max_rho=max_rho,
                                  seed=rec["seed"], p=rec["p"], deg=rec["deg"],
                                  x=rec["x"], y=rec["y"], tau=rec["tau"],
                                  flip=list(flip), rho=len(flip),
                                  dmins=[float(h) for h in hs],
                                  ostar_valid=bool(m.get("ostar_valid")),
                                  est=m.get("est"), est0=arm["est0"],
                                  K=[list(e) for e in arm["K"]]))
    return found


def dossier(f):
    """Regenerate the SCM exactly and produce the full PREREG 4.5 record."""
    rng = np.random.default_rng(f["seed"])
    D = random_dag(f["p"], f["deg"], rng)
    C = dag_to_cpdag(D)
    U = undirected_edges(C)
    from adjust import possibly_causal_paths
    pairs = [(a, b) for a in range(f["p"]) for b in range(f["p"])
             if a != b and len(possibly_causal_paths(D, a, b)) > 0]
    x, y = pairs[int(rng.integers(len(pairs)))]
    scm = make_linear_iscm(D, rng)
    Sigma = cov_linear(scm["A"], scm["omega"])
    tau = total_effect_linear(scm["A"], x, y)
    K_all = [(i, j) if D[i, j] == 1 else (j, i) for (i, j) in U]
    order = rng.permutation(len(K_all))
    K_all = [K_all[i] for i in order]
    K = K_all[:len(f["K"])]
    ok = dict(xy=(x, y) == (f["x"], f["y"]),
              tau=abs(tau - f["tau"]) < 1e-9,
              K=[list(e) for e in K] == f["K"])
    dist = X.hop_dist_inf(C, [x, y])
    G0, chain0 = apply_bk_traced(C, K)
    assert np.array_equal(G0, apply_background_knowledge(C, K)), "traced != library"
    O0, _, a0 = X.ostar_and_paths(G0, x, y)
    Kp = [(v, u) if i in f["flip"] else (u, v) for i, (u, v) in enumerate(K)]
    G1, chain1 = apply_bk_traced(C, Kp)
    assert np.array_equal(G1, apply_background_knowledge(C, Kp)), "traced != library"
    O1, _, a1 = X.ostar_and_paths(G1, x, y)
    b0 = ols_with_se(Sigma, x, y, O0, np.inf)[0]
    b1 = ols_with_se(Sigma, x, y, O1, np.inf)[0]
    return dict(
        regen_ok=ok,
        cmd=(f"python -c \"import numpy as np;from graphs import random_dag,dag_to_cpdag;"
             f"rng=np.random.default_rng({f['seed']});D=random_dag({f['p']},{f['deg']},rng)\""),
        ens=f["ens"], arm=f["arm"], seed=f["seed"], p=f["p"], deg=f["deg"],
        D_edges=edge_list(D), C_edges=edge_list(C),
        x=int(x), y=int(y), tau=float(tau),
        delta={int(v): (float(dist[v]) if np.isfinite(dist[v]) else None)
               for v in range(f["p"])},
        K=[list(e) for e in K],
        K_perturbed=[list(e) for e in Kp],
        flipped=[list(K[i]) for i in f["flip"]],
        flipped_dmin=[float(X.stmt_dist(dist, K[i])[0]) for i in f["flip"]],
        flipped_dmax=[float(X.stmt_dist(dist, K[i])[1]) for i in f["flip"]],
        G0_edges=edge_list(G0), G1_edges=edge_list(G1),
        meek_chain_base=chain0, meek_chain_perturbed=chain1,
        O0=sorted(int(z) for z in O0), O1=sorted(int(z) for z in O1),
        beta0=float(b0), beta1=float(b1),
        O0_valid_in_D=bool(is_valid_adjustment_set(D, x, y, O0)),
        O1_valid_in_D=bool(is_valid_adjustment_set(D, x, y, O1)),
        rel_bias=float(abs(b1 - tau) / abs(tau)),
        rel_bias_base=float(abs(b0 - tau) / abs(tau)))


# ------------------------------------------------------------- minimality
def exhaustive_rho2_far(p):
    """Exhaustive: does ANY pair of far statements move O* at p nodes?
    Reuses run_lemma's reachable-MPDAG enumeration."""
    from run_lemma import all_dags, reachable_mpdags
    cp = {}
    for D in all_dags(p):
        C = dag_to_cpdag(D)
        cp.setdefault(C.tobytes(), C)
    n_trials = n_hits = n_am = 0
    hits = []
    for C in cp.values():
        ref = v_structures(C)
        for G0 in reachable_mpdags(C):
            U0 = undirected_edges(G0)
            if len(U0) < 2:
                continue
            for x in range(p):
                for y in range(p):
                    if x == y:
                        continue
                    O0, _, a0 = X.ostar_and_paths(G0, x, y)
                    if not a0:
                        continue
                    dist = X.hop_dist_inf(G0, [x, y])
                    for e1, e2 in combinations(U0, 2):
                        d1 = X.stmt_dist(dist, e1)[0]
                        d2 = X.stmt_dist(dist, e2)[0]
                        if not (np.isfinite(d1) and d1 >= 1 and np.isfinite(d2) and d2 >= 1):
                            continue
                        for o1 in (e1, (e1[1], e1[0])):
                            for o2 in (e2, (e2[1], e2[0])):
                                H = G0.copy()
                                H[o1[1], o1[0]] = 0
                                H = meek_closure(H)
                                if H[o2[0], o2[1]] == 0:
                                    continue
                                H[o2[1], o2[0]] = 0
                                H = meek_closure(H)
                                n_trials += 1
                                if has_directed_cycle(H) or v_structures(H) != ref:
                                    continue
                                O1, _, a1 = X.ostar_and_paths(H, x, y)
                                if not a1:
                                    continue
                                n_am += 1
                                if O1 != O0:
                                    n_hits += 1
                                    if len(hits) < 5:
                                        hits.append(dict(p=p, x=x, y=y,
                                                         s1=list(map(int, o1)),
                                                         s2=list(map(int, o2)),
                                                         G0=edge_list(G0),
                                                         O0=sorted(int(z) for z in O0),
                                                         O1=sorted(int(z) for z in O1)))
    return dict(p=p, n_cpdags=len(cp), n_trials=n_trials, n_amenable=n_am,
                n_hits=n_hits, hits=hits)


def main():
    found = find_all()
    byp = defaultdict(int)
    byens = defaultdict(int)
    for f in found:
        byp[f["p"]] += 1
        byens[f["ens"] + "/" + f["arm"]] += 1
    print(f"[exhibit] {len(found)} rho>=2 ALL-FAR O*-changes in the archive")
    print(f"[exhibit] by p: {dict(sorted(byp.items()))}")
    print(f"[exhibit] by cell: {dict(byens)}")
    pmin = min(byp) if byp else None
    print(f"[exhibit] minimum p in the searched ensembles: {pmin}")

    found.sort(key=lambda f: (f["p"], f["rho"], f["seed"]))
    dos = []
    for f in found[:NEX]:
        try:
            dos.append(dossier(f))
        except Exception as e:
            dos.append(dict(error=str(e), f=f))
    mini = {}
    for p in (4, 5):
        if pmin is not None and p >= pmin:
            break
        mini[str(p)] = exhaustive_rho2_far(p)
        print(f"[exhibit] exhaustive rho=2 far-pair search p={p}: "
              f"{mini[str(p)]['n_hits']} hits of {mini[str(p)]['n_amenable']} amenable trials")
    out = dict(n_found=len(found), by_p=dict(sorted(byp.items())),
               by_cell=dict(byens), p_min=pmin,
               n_bias=sum(1 for f in found if not f["ostar_valid"]),
               dossiers=dos, minimality=mini,
               all_found=[{k: v for k, v in f.items() if k != "K"} for f in found])
    with open(OUT, "w") as f2:
        json.dump(out, f2, indent=1, default=float)
    print(f"[exhibit] -> {OUT}")
    for d in dos:
        if "error" in d:
            print("  ERROR", d["error"])
            continue
        print("\n" + "-" * 70)
        print(f"  {d['ens']}/{d['arm']} seed={d['seed']} p={d['p']} deg={d['deg']} "
              f"regen_ok={d['regen_ok']}")
        print(f"  D: {' '.join(d['D_edges'])}")
        print(f"  C: {' '.join(d['C_edges'])}")
        print(f"  X={d['x']} Y={d['y']} tau={d['tau']:.6f}  delta={d['delta']}")
        print(f"  K={d['K']}  flipped={d['flipped']} dmin={d['flipped_dmin']} dmax={d['flipped_dmax']}")
        print(f"  O*(G0)={d['O0']} beta={d['beta0']:.6f} valid={d['O0_valid_in_D']}")
        print(f"  O*(G')={d['O1']} beta={d['beta1']:.6f} valid={d['O1_valid_in_D']}  "
              f"|bias|/|tau|={d['rel_bias']:.4f}")
        for st in d["meek_chain_perturbed"]:
            pr = "; ".join(f"{t['rule']}({t['witness']}) => {t['oriented'][0]}->{t['oriented'][1]}"
                           for t in st["propagated"])
            print(f"    assert {st['asserted'][0]}->{st['asserted'][1]}  propagated: {pr or '(none)'}")


if __name__ == "__main__":
    main()
