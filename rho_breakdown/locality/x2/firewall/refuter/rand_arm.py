"""Random + ADVERSARIALLY-BIASED strong-form search at larger p.

For each stamped H that survives the three checks and stays amenable, we look
for a bad D in [G0] two ways:
  (1) full enumeration of [G0]  when |U(G0)| is small;
  (2) biased completions: orient every undirected edge of G0 AGAINST the
      orientation H forced on it (the mechanism the theorem must survive),
      then randomly for the rest.
Also instruments the synthesiser's proof handle: cn_H >= cn_G0 ?  forb_H >= forb_G0 ?
"""
import sys, time, json, random
sys.path.insert(0, "/private/tmp/claude-501/-Users-josecosta-mugango/bceac0a4-93a7-4757-b3b3-a3d4e6c42d55/scratchpad/cfw_ref_A")
from cfwlib import *
import cfwlib as _CF
from graphs import adjacent

def extend_biased(G0, ref, pref, rng, tries=6):
    """One consistent DAG extension of G0, preferring pref[(u,v)] orientation."""
    for _ in range(tries):
        G = G0.copy()
        U = undirected_edges(G)
        rng.shuffle(U)
        ok = True
        for (u, v) in U:
            if not is_undirected(G, u, v):
                continue
            want = pref.get((u, v))
            if want is None:
                cand = [(u, v), (v, u)]
                rng.shuffle(cand)
            else:
                cand = [want, (want[1], want[0])]
            placed = False
            for (a, b) in cand:
                Hc = G.copy(); Hc[b, a] = 0; Hc = meek_closure(Hc)
                if has_directed_cycle(Hc) or v_structures(Hc) != ref:
                    continue
                if not x1_ops.pdag_extendable(Hc):
                    continue
                G = Hc; placed = True; break
            if not placed:
                ok = False; break
        if ok and not undirected_edges(G):
            return G
    return None


def flip_pref(G0, H):
    """For every edge undirected in G0 and directed in H, prefer the OPPOSITE."""
    pref = {}
    for (u, v) in undirected_edges(G0):
        if is_directed(H, u, v):
            pref[(u, v)] = (v, u)
        elif is_directed(H, v, u):
            pref[(u, v)] = (u, v)
    return pref


def run(nscm, ps, degs, seed0, ndraw=12, ufull=13, tlimit=1e9, verbose=True):
    rng = np.random.default_rng(seed0)
    pyr = random.Random(seed0 + 1)
    st = new_stats(); hits = []
    t0 = time.time()
    for it in range(nscm):
        if time.time() - t0 > tlimit: break
        p = int(rng.choice(ps)); deg = float(rng.choice(degs))
        D0 = random_dag(p, deg, rng)
        C = dag_to_cpdag(D0); ref = v_structures(C)
        # G0: C itself, or a random partial orientation
        U = undirected_edges(C)
        G0 = C.copy()
        if len(U) and rng.random() < 0.6:
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
        full = consistent_dag_extensions(G0, ref_vstructs=ref) if len(U0) <= ufull else None
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
                        if not (cn0 <= cn1): st["cn_shrank"] += 1
                        if not (fb0 <= fb1): st["forb_shrank"] += 1
                        moved = (O1 != O0)
                        if moved: st["n_moved"] += 1
                        # --- hunt for a bad D
                        cands = []
                        if full is not None:
                            cands = full
                        else:
                            pref = flip_pref(G0, H)
                            seen = set()
                            for _ in range(ndraw):
                                Dc = extend_biased(G0, ref, pref, pyr)
                                if Dc is not None and Dc.tobytes() not in seen:
                                    seen.add(Dc.tobytes()); cands.append(Dc)
                            for _ in range(ndraw // 2):
                                Dc = extend_biased(G0, ref, {}, pyr)
                                if Dc is not None and Dc.tobytes() not in seen:
                                    seen.add(Dc.tobytes()); cands.append(Dc)
                        st["n_Dcheck"] += len(cands)
                        bad = [Dc for Dc in cands if not valid_in(Dc, x, y, O1)]
                        if bad:
                            st["n_violation"] += 1
                            if moved: st["n_violation_moved"] += 1
                            if len(hits) < 10:
                                hits.append(dict(p=p, G0=G0.tolist(), H=H.tolist(),
                                                 x=int(x), y=int(y), edge=[int(a), int(b)],
                                                 O0=sorted(int(z) for z in O0),
                                                 O1=sorted(int(z) for z in O1),
                                                 Dbad=bad[0].tolist(), moved=bool(moved)))
        if len(_CF._vcache) > 600000: _CF._vcache.clear()
        if verbose and it % 25 == 0:
            print(it, dict(st), round(time.time() - t0, 1), flush=True)
    return st, hits


if __name__ == "__main__":
    n = int(sys.argv[1]); ps = [int(z) for z in sys.argv[2].split(",")]
    s = int(sys.argv[3]); degs = [float(z) for z in sys.argv[4].split(",")]
    tl = float(sys.argv[5]) if len(sys.argv) > 5 else 1e9
    st, hits = run(n, ps, degs, s, tlimit=tl)
    print("FINAL", json.dumps(dict(st), indent=1))
    print("n_hits", len(hits))
    json.dump({"stats": dict(st), "hits": hits}, open(f"rand_{'_'.join(map(str,ps))}_{s}.json", "w"))
