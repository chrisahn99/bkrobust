"""Arm 2: dense / scale-free / long-chordless, FULL class enumeration when
feasible, plus instrumentation of the two halves of the adjustment criterion.

Instruments (necessary conditions for a violation):
  cn_shrink   : cn_H  NOT >= cn_G0      (kills the synthesiser's proof handle)
  forb_shrink : forb_H NOT >= forb_G0   (necessary for the FORBIDDEN half)
  O_lost      : O*_G0 \\ O*_H != {}      (necessary-ish for the D-SEP half)
"""
import sys, time, json, random
sys.path.insert(0, "/private/tmp/claude-501/-Users-josecosta-mugango/bceac0a4-93a7-4757-b3b3-a3d4e6c42d55/scratchpad/cfw_ref_A")
from cfwlib import *
import cfwlib as _CF
from rand_arm import extend_biased, flip_pref


def scale_free_dag(p, m, rng):
    D = np.zeros((p, p), dtype=np.int8)
    order = rng.permutation(p)
    deg = np.ones(p)
    for i in range(1, p):
        v = order[i]
        prev = order[:i]
        w = deg[prev] / deg[prev].sum()
        k = min(m, i)
        ch = rng.choice(i, size=k, replace=False, p=w)
        for t in ch:
            u = prev[t]
            D[u, v] = 1
            deg[u] += 1; deg[v] += 1
    return D


def chain_dag(p, rng, extra=0):
    """Long chordless path skeleton (+ a few random extras)."""
    D = np.zeros((p, p), dtype=np.int8)
    order = rng.permutation(p)
    for i in range(p - 1):
        a, b = order[i], order[i + 1]
        if rng.random() < 0.5: D[a, b] = 1
        else: D[b, a] = 1
    # fix acyclicity by re-deriving from a topological relabel
    A = D.astype(bool); R = A.copy()
    for _ in range(p): R = R | (R @ A)
    if np.any(np.diag(R)):
        D = np.zeros((p, p), dtype=np.int8)
        for i in range(p - 1):
            D[order[i], order[i + 1]] = 1
    for _ in range(extra):
        i, j = rng.integers(0, p, 2)
        if i == j: continue
        a, b = order[min(i, j)], order[max(i, j)]
        D[a, b] = 1
    A = D.astype(bool); R = A.copy()
    for _ in range(p): R = R | (R @ A)
    if np.any(np.diag(R)): return None
    return D


def gen(rng, mode):
    if mode == "er":
        p = int(rng.choice([6, 7, 8])); deg = float(rng.choice([2.0, 2.5, 3.0, 3.5, 4.0, 5.0]))
        return random_dag(p, deg, rng)
    if mode == "sf":
        p = int(rng.choice([7, 8, 9])); m = int(rng.choice([1, 2, 3]))
        return scale_free_dag(p, m, rng)
    if mode == "chain":
        p = int(rng.choice([7, 8, 9, 10]))
        return chain_dag(p, rng, extra=int(rng.integers(0, 3)))
    if mode == "dense":
        p = int(rng.choice([6, 7]))
        return random_dag(p, float(rng.choice([4.0, 5.0, 6.0])), rng)


def run(nscm, seed0, modes, tlimit, ufull=14, ndraw=25):
    rng = np.random.default_rng(seed0); pyr = random.Random(seed0 + 7)
    st = new_stats(); hits = []; t0 = time.time()
    it = 0
    while it < nscm and time.time() - t0 < tlimit:
        it += 1
        mode = modes[it % len(modes)]
        D0 = gen(rng, mode)
        if D0 is None: continue
        p = D0.shape[0]
        C = dag_to_cpdag(D0); ref = v_structures(C)
        U = undirected_edges(C)
        G0 = C.copy()
        r = rng.random()
        if len(U) and r < 0.55:                      # partial knowledge
            k = int(rng.integers(0, len(U)))
            for t in rng.permutation(len(U))[:k]:
                u, v = U[t]
                if not is_undirected(G0, u, v): continue
                a, b = (u, v) if rng.random() < 0.5 else (v, u)
                Hc = G0.copy(); Hc[b, a] = 0; Hc = meek_closure(Hc)
                if has_directed_cycle(Hc) or v_structures(Hc) != ref: continue
                G0 = Hc
        U0 = undirected_edges(G0)
        S = skeleton(G0)
        NA = [(a, b) for a in range(p) for b in range(a + 1, p) if S[a, b] == 0]
        if not NA: continue
        st["n_scm"] += 1
        st["n_U0"] += len(U0)
        full = consistent_dag_extensions(G0, ref_vstructs=ref) if len(U0) <= ufull else None
        if full is not None: st["n_full_enum"] += 1
        for x in range(p):
            for y in range(p):
                if x == y: continue
                O0, np0, amen0 = X.ostar_and_paths(G0, x, y)
                if not amen0: continue
                cn0 = set()
                for pth in X.pcp_capped(G0, x, y): cn0.update(pth[1:])
                fb0 = adjust.poss_de(G0, cn0) | {x}
                st["n_query_pos"] += 1
                for (a0, b0) in NA:
                    for (a, b) in ((a0, b0), (b0, a0)):
                        H, info = x1_ops.bk_assert(G0, [(a, b)])
                        st["n_trial"] += 1
                        if info["conflict"] or has_directed_cycle(H) or not x1_ops.pdag_extendable(H):
                            st["n_rejected"] += 1; continue
                        st["n_stamped"] += 1
                        O1, np1, amen1 = X.ostar_and_paths(H, x, y)
                        if not amen1:
                            st["n_abort"] += 1
                            st["n_abort_zero" if np1 == 0 else "n_abort_unamen"] += 1
                            continue
                        st["n_pos"] += 1
                        cn1 = set()
                        for pth in X.pcp_capped(H, x, y): cn1.update(pth[1:])
                        fb1 = adjust.poss_de(H, cn1) | {x}
                        if not (cn0 <= cn1): st["cn_shrink"] += 1
                        if not (fb0 <= fb1):
                            st["forb_shrink"] += 1
                            if len(hits) < 6:
                                hits.append(dict(kind="forb_shrink", p=p, G0=G0.tolist(),
                                                 H=H.tolist(), x=int(x), y=int(y),
                                                 edge=[int(a), int(b)],
                                                 fb0=sorted(map(int, fb0)), fb1=sorted(map(int, fb1))))
                        if O0 - O1: st["O_lost"] += 1
                        moved = (O1 != O0)
                        if moved: st["n_moved"] += 1
                        if full is not None:
                            cands = full
                        else:
                            pref = flip_pref(G0, H); seen = set(); cands = []
                            for _ in range(ndraw):
                                Dc = extend_biased(G0, ref, pref, pyr)
                                if Dc is not None and Dc.tobytes() not in seen:
                                    seen.add(Dc.tobytes()); cands.append(Dc)
                            for _ in range(ndraw):
                                Dc = extend_biased(G0, ref, {}, pyr)
                                if Dc is not None and Dc.tobytes() not in seen:
                                    seen.add(Dc.tobytes()); cands.append(Dc)
                        st["n_Dcheck"] += len(cands)
                        bad = [Dc for Dc in cands if not valid_in(Dc, x, y, O1)]
                        if bad:
                            st["n_violation"] += 1
                            if moved: st["n_violation_moved"] += 1
                            hits.insert(0, dict(kind="VIOLATION", p=p, G0=G0.tolist(),
                                                H=H.tolist(), x=int(x), y=int(y),
                                                edge=[int(a), int(b)],
                                                O0=sorted(map(int, O0)), O1=sorted(map(int, O1)),
                                                Dbad=bad[0].tolist(), moved=bool(moved)))
        if len(_CF._vcache) > 800000: _CF._vcache.clear()
        if it % 20 == 0:
            print(it, dict(st), round(time.time() - t0, 1), flush=True)
    return st, hits


if __name__ == "__main__":
    n = int(sys.argv[1]); s = int(sys.argv[2]); tl = float(sys.argv[3])
    modes = sys.argv[4].split(",") if len(sys.argv) > 4 else ["er", "sf", "chain", "dense"]
    st, hits = run(n, s, modes, tl)
    print("FINAL", json.dumps(dict(st), indent=1))
    print("n_hits", len(hits), "n_violation", st.get("n_violation", 0))
    json.dump({"stats": dict(st), "hits": hits[:6]}, open(f"arm2_{s}.json", "w"))
