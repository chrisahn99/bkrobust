"""X2 REFUTER: search for a dmin>=2 counterexample in families X2 did not search.

Design rationale (from the X2 archive re-analysis):
  * X2's dmin>=2 denominator (6,726) is 83% `large` (p 15-25, deg 1.5-2.5), whose
    OWN dmin=1 rate is 0.0002.  Its within-SCM matched near rate is 1/5,253.
    So the "0 of 6,726" has ~6 expected events pooled, ~0.5 matched.
  * The events live in DENSE, SMALL graphs (mean skeleton density 0.74; 41/70
    have eccentricity 1, i.e. NO far edge can exist there at all).
  * => search graphs that are dense LOCALLY (so undirected edges survive Meek and
    O* is movable) but have real DISTANCE structure (so dmin>=2 edges exist).
Families: clique-chain, clique-chain+ER noise, BA scale-free, ER at p 7-14 with
argmax-delta query, and "two blobs + bridge".
"""
import sys, json, time, itertools, collections
import numpy as np
sys.path.insert(0, "/home/costaj/latent-causal/x2-locality/code")
from graphs import (random_dag, dag_to_cpdag, undirected_edges, directed_edges,
                    apply_background_knowledge, MeekFail, skeleton, is_directed)
from adjust import (optimal_adjustment_set, is_valid_adjustment_set,
                    possibly_causal_paths, causal_nodes)
from multiprocessing import Pool

INF = 1 << 20

def hop(C, srcs, p):
    S = skeleton(C); d = np.full(p, INF, dtype=int)
    fr = list(srcs)
    for v in fr: d[v] = 0
    while fr:
        nxt = []
        for v in fr:
            for w in np.flatnonzero(S[v]):
                if d[w] > d[v] + 1:
                    d[w] = d[v] + 1; nxt.append(int(w))
        fr = nxt
    return d

# ---------------------------------------------------------------- skeleton families
def skel_clique_chain(rng, n_cl, s, p_extra=0.0):
    p = n_cl * s
    S = np.zeros((p, p), dtype=np.int8)
    for c in range(n_cl):
        idx = list(range(c*s, (c+1)*s))
        for a, b in itertools.combinations(idx, 2): S[a, b] = S[b, a] = 1
        if c:                                   # single bridge to previous clique
            u = rng.integers((c-1)*s, c*s); v = rng.integers(c*s, (c+1)*s)
            S[u, v] = S[v, u] = 1
    if p_extra > 0:
        for a, b in itertools.combinations(range(p), 2):
            if S[a, b] == 0 and rng.random() < p_extra: S[a, b] = S[b, a] = 1
    return S

def skel_ba(rng, p, m):
    S = np.zeros((p, p), dtype=np.int8)
    for a, b in itertools.combinations(range(m+1), 2): S[a, b] = S[b, a] = 1
    deg = S.sum(0).astype(float)
    for v in range(m+1, p):
        w = deg[:v] + 0.5; w = w / w.sum()
        tg = rng.choice(v, size=min(m, v), replace=False, p=w)
        for t in tg: S[v, t] = S[t, v] = 1
        deg = S.sum(0).astype(float)
    return S

def skel_two_blob(rng, s, bridge_len):
    p = 2*s + bridge_len
    S = np.zeros((p, p), dtype=np.int8)
    A = list(range(s)); B = list(range(s+bridge_len, p)); br = list(range(s, s+bridge_len))
    for grp in (A, B):
        for a, b in itertools.combinations(grp, 2): S[a, b] = S[b, a] = 1
    path = [int(rng.integers(0, s))] + br + [int(rng.choice(B))]
    for u, v in zip(path, path[1:]): S[u, v] = S[v, u] = 1
    return S

def dag_from_skel(S, rng):
    p = S.shape[0]; order = rng.permutation(p); pos = np.empty(p, int)
    for i, v in enumerate(order): pos[v] = i
    D = np.zeros((p, p), dtype=np.int8)
    for a in range(p):
        for b in range(a+1, p):
            if S[a, b]:
                if pos[a] < pos[b]: D[a, b] = 1
                else: D[b, a] = 1
    return D

def draw_family(fam, rng):
    if fam == "clique_chain":
        return skel_clique_chain(rng, int(rng.integers(3, 5)), int(rng.integers(3, 5)))
    if fam == "clique_chain_noise":
        return skel_clique_chain(rng, int(rng.integers(3, 5)), int(rng.integers(3, 5)), 0.05)
    if fam == "ba":
        return skel_ba(rng, int(rng.integers(7, 15)), int(rng.integers(2, 4)))
    if fam == "two_blob":
        return skel_two_blob(rng, int(rng.integers(3, 6)), int(rng.integers(1, 4)))
    if fam == "er":
        p = int(rng.integers(7, 15)); d = float(rng.choice([3.0, 4.0, 5.0, 6.0]))
        return skeleton(random_dag(p, d, rng))
    raise ValueError(fam)

# ---------------------------------------------------------------- one graph
def one(args):
    seed, fam, qrule = args
    rng = np.random.default_rng(seed)
    S = draw_family(fam, rng)
    p = S.shape[0]
    D = dag_from_skel(S, rng)
    C = dag_to_cpdag(D)
    U = undirected_edges(C)
    if len(U) < 2: return None
    pairs = [(a, b) for a in range(p) for b in range(p)
             if a != b and len(possibly_causal_paths(D, a, b)) > 0]
    if not pairs: return None
    if qrule == "argmax_delta":
        best, bd = [], -1
        for (a, b) in pairs:
            d = hop(C, [a], p)[b]
            if d < INF and d > bd: best, bd = [(a, b)], d
            elif d == bd: best.append((a, b))
        x, y = best[int(rng.integers(len(best)))]
    else:
        x, y = pairs[int(rng.integers(len(pairs)))]
    dist = hop(C, [x, y], p)
    # base knowledge: a random subset of TRUE statements on undirected edges of C
    K_all = [(i, j) if D[i, j] == 1 else (j, i) for (i, j) in U]
    idx = rng.permutation(len(K_all))
    kb = int(rng.integers(0, len(K_all)+1))
    Kb = [K_all[i] for i in idx[:kb]]
    try:
        G0 = apply_background_knowledge(C, Kb)
    except MeekFail:
        return None
    U0 = undirected_edges(G0)          # still-open edges in the base MPDAG
    if not U0: return None
    O0 = optimal_adjustment_set(G0, x, y)
    out = collections.Counter(); hits = []
    # ---- rho=1 sweep: reversal (both orientations of one open edge)
    for (u, v) in U0:
        dm = int(min(dist[u], dist[v])); dM = int(max(dist[u], dist[v]))
        if dm >= INF: bucket = "disc"
        else: bucket = str(min(dm, 4))
        Os = {}
        for (a, b) in ((u, v), (v, u)):
            try: G = apply_background_knowledge(G0, [(a, b)])
            except MeekFail: Os[(a, b)] = "incons"; continue
            Os[(a, b)] = optimal_adjustment_set(G, x, y)
        oA, oB = Os[(u, v)], Os[(v, u)]
        out[f"rev/{bucket}/D_all"] += 1
        if oA == "incons" or oB == "incons": continue
        out[f"rev/{bucket}/D_con"] += 1
        if oA is None or oB is None: continue
        out[f"rev/{bucket}/D_am"] += 1
        if oA != oB:
            out[f"rev/{bucket}/N_set"] += 1
            if dm >= 2 and dm < INF:
                # which orientation is the misstatement (disagrees with D)?
                trueor = (u, v) if D[u, v] == 1 else (v, u)
                wrong = (v, u) if trueor == (u, v) else (u, v)
                Ow = Os[wrong]
                hits.append(dict(kind="rev", seed=int(seed), fam=fam, p=int(p),
                                 x=int(x), y=int(y), edge=[int(u), int(v)],
                                 dmin=dm, dmax=dM, kb=kb,
                                 O_true=sorted(Os[trueor]), O_wrong=sorted(Ow),
                                 valid_wrong=bool(is_valid_adjustment_set(D, x, y, Ow)),
                                 valid_true=bool(is_valid_adjustment_set(D, x, y, Os[trueor]))))
    # ---- rho=1 "add from the amenable base" (X2 claims a structural zero)
    if O0 is not None:
        for (u, v) in U0:
            dm = int(min(dist[u], dist[v]))
            bucket = "disc" if dm >= INF else str(min(dm, 4))
            for (a, b) in ((u, v), (v, u)):
                try: G = apply_background_knowledge(G0, [(a, b)])
                except MeekFail: continue
                O = optimal_adjustment_set(G, x, y)
                out[f"add/{bucket}/D_con"] += 1
                if O is None: continue
                out[f"add/{bucket}/D_am"] += 1
                if O != O0:
                    out[f"add/{bucket}/N_set"] += 1
                    if 2 <= dm < INF:
                        hits.append(dict(kind="add", seed=int(seed), fam=fam, p=int(p),
                                         x=int(x), y=int(y), stmt=[int(a), int(b)],
                                         dmin=dm, kb=kb, O0=sorted(O0), O1=sorted(O),
                                         true_in_D=bool(D[a, b] == 1),
                                         valid_O1=bool(is_valid_adjustment_set(D, x, y, O)),
                                         valid_O0=bool(is_valid_adjustment_set(D, x, y, O0))))
    return dict(counts=dict(out), hits=hits, p=int(p), ecc=int(max(d for d in dist if d < INF)))

def main():
    n = int(sys.argv[1]); fam = sys.argv[2]; qrule = sys.argv[3]; nw = int(sys.argv[4])
    s0 = int(sys.argv[5]); tag = sys.argv[6]
    jobs = [(s0 + i, fam, qrule) for i in range(n)]
    tot = collections.Counter(); hits = []; nkeep = 0; eccs = collections.Counter()
    t0 = time.time()
    with Pool(nw) as pool:
        for r in pool.imap_unordered(one, jobs, chunksize=8):
            if not r: continue
            nkeep += 1; tot.update(r["counts"]); hits.extend(r["hits"]); eccs[r["ecc"]] += 1
    res = dict(fam=fam, qrule=qrule, n_draw=n, n_keep=nkeep, secs=time.time()-t0,
               ecc=dict(eccs), counts=dict(tot), n_hits=len(hits), hits=hits[:200])
    json.dump(res, open(f"/home/costaj/latent-causal/x2-locality-refute/results/{tag}.json", "w"))
    def g(k): return tot.get(k, 0)
    print(f"[{tag}] fam={fam} q={qrule} kept {nkeep}/{n} in {res['secs']:.0f}s  ecc={dict(eccs)}")
    for op in ("rev", "add"):
        for b in ("0", "1", "2", "3", "4", "disc"):
            da = g(f"{op}/{b}/D_am"); ns = g(f"{op}/{b}/N_set")
            if da: print(f"   {op} dmin={b:4s} N_set={ns:7d} / D_am={da:8d}  = {ns/da:.5f}")
    print(f"   HITS at dmin>=2: {len(hits)}")

if __name__ == "__main__":
    main()
