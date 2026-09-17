"""Arm 3: hunt the EXACT preconditions of each half of the adjustment criterion.

forb_D subset forb_G0 for every D in [G0] (checked, not assumed).  So:
  FORBIDDEN half can only break if   O*_H  intersects  forb_G0     -> `pre_forb`
  DSEP half is checked directly against every D in [G0].
Whenever `pre_forb` fires we enumerate the WHOLE class no matter its size.
"""
import sys, time, json, random
sys.path.insert(0, "/private/tmp/claude-501/-Users-josecosta-mugango/bceac0a4-93a7-4757-b3b3-a3d4e6c42d55/scratchpad/cfw_ref_A")
from cfwlib import *
import cfwlib as _CF
from rand_arm import extend_biased, flip_pref
from arm2 import gen

def run(nscm, seed0, modes, tlimit, ufull=15, ndraw=30, check_forbD=False):
    rng = np.random.default_rng(seed0); pyr = random.Random(seed0 + 3)
    st = new_stats(); hits = []; t0 = time.time(); it = 0
    while it < nscm and time.time() - t0 < tlimit:
        it += 1
        D0 = gen(rng, modes[it % len(modes)])
        if D0 is None: continue
        p = D0.shape[0]
        C = dag_to_cpdag(D0); ref = v_structures(C)
        U = undirected_edges(C); G0 = C.copy()
        if len(U) and rng.random() < 0.5:
            k = int(rng.integers(0, len(U)))
            for t in rng.permutation(len(U))[:k]:
                u, v = U[t]
                if not is_undirected(G0, u, v): continue
                a, b = (u, v) if rng.random() < 0.5 else (v, u)
                Hc = G0.copy(); Hc[b, a] = 0; Hc = meek_closure(Hc)
                if has_directed_cycle(Hc) or v_structures(Hc) != ref: continue
                G0 = Hc
        U0 = undirected_edges(G0); S = skeleton(G0)
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
                if check_forbD and full is not None:
                    for Dc in full:
                        if not (adjust.forb(Dc, x, y) <= fb0):
                            st["forbD_not_subset"] += 1
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
                        inter = set(O1) & fb0
                        if inter:
                            st["pre_forb"] += 1
                            if len(hits) < 8:
                                hits.append(dict(kind="pre_forb", p=p, G0=G0.tolist(),
                                                 H=H.tolist(), x=int(x), y=int(y),
                                                 edge=[int(a), int(b)],
                                                 O1=sorted(map(int, O1)),
                                                 fb0=sorted(map(int, fb0)),
                                                 inter=sorted(map(int, inter))))
                        moved = (O1 != O0)
                        if moved: st["n_moved"] += 1
                        if full is not None:
                            cands = full
                        elif inter:
                            cands = consistent_dag_extensions(G0, ref_vstructs=ref)
                        else:
                            pref = flip_pref(G0, H); seen = set(); cands = []
                            for _ in range(ndraw):
                                for pf in (pref, {}):
                                    Dc = extend_biased(G0, ref, pf, pyr)
                                    if Dc is not None and Dc.tobytes() not in seen:
                                        seen.add(Dc.tobytes()); cands.append(Dc)
                        st["n_Dcheck"] += len(cands)
                        bad = [Dc for Dc in cands if not valid_in(Dc, x, y, O1)]
                        if bad:
                            st["n_violation"] += 1
                            hits.insert(0, dict(kind="VIOLATION", p=p, G0=G0.tolist(),
                                                H=H.tolist(), x=int(x), y=int(y),
                                                edge=[int(a), int(b)],
                                                O0=sorted(map(int, O0)), O1=sorted(map(int, O1)),
                                                Dbad=bad[0].tolist(), moved=bool(moved)))
        if len(_CF._vcache) > 800000: _CF._vcache.clear()
        if it % 25 == 0:
            print(it, dict(st), round(time.time() - t0, 1), flush=True)
    return st, hits


if __name__ == "__main__":
    n = int(sys.argv[1]); s = int(sys.argv[2]); tl = float(sys.argv[3])
    modes = sys.argv[4].split(",")
    cf = len(sys.argv) > 5 and sys.argv[5] == "checkforbD"
    st, hits = run(n, s, modes, tl, check_forbD=cf)
    print("FINAL", json.dumps(dict(st), indent=1))
    json.dump({"stats": dict(st), "hits": hits[:8]}, open(f"arm3_{s}.json", "w"))
