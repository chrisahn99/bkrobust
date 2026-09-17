"""Arm 5: MOVED-only, full-class strong-form check, parallel, dense + varied.

Only trials where O*_H != O*_G0 can break the theorem (if O*_H == O*_G0 the set
is HPM-valid in every D in [G0] already), so all D-enumeration budget goes there.
Also carries the two preconditions:
  pre_forb : O*_H  meets forb_G0   (necessary for the FORBIDDEN half)
  O_lost   : O*_G0 \\ O*_H nonempty (route to the DSEP half)
"""
import sys, time, json, random
sys.path.insert(0, "/private/tmp/claude-501/-Users-josecosta-mugango/bceac0a4-93a7-4757-b3b3-a3d4e6c42d55/scratchpad/cfw_ref_A")
from cfwlib import *
import cfwlib as _CF
from arm2 import gen, scale_free_dag, chain_dag
from multiprocessing import Pool


def one(args):
    seed, mode = args
    rng = np.random.default_rng(seed)
    st = new_stats(); hits = []
    D0 = gen(rng, mode)
    if D0 is None: return dict(st), hits
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
    if not NA or len(U0) > 17 or len(U0) < 9: return dict(st), hits
    st["n_scm"] += 1; st["n_U0"] += len(U0)
    full = consistent_dag_extensions(G0, ref_vstructs=ref)
    st["class_size"] += len(full)
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
                    if O1 == O0: continue
                    st["n_moved"] += 1
                    if set(O1) & fb0: st["pre_forb"] += 1
                    if O0 - O1: st["O_lost"] += 1
                    if O1 - O0: st["O_gained"] += 1
                    bad = [Dc for Dc in full if not valid_in(Dc, x, y, O1)]
                    st["n_Dcheck"] += len(full)
                    if bad:
                        st["n_violation"] += 1
                        hits.append(dict(kind="VIOLATION", p=p, seed=int(seed), mode=mode,
                                         G0=G0.tolist(), H=H.tolist(), x=int(x), y=int(y),
                                         edge=[int(a), int(b)],
                                         O0=sorted(map(int, O0)), O1=sorted(map(int, O1)),
                                         Dbad=bad[0].tolist()))
    _CF._vcache.clear()
    return dict(st), hits


if __name__ == "__main__":
    n = int(sys.argv[1]); s0 = int(sys.argv[2]); nw = int(sys.argv[3])
    modes = sys.argv[4].split(",")
    t0 = time.time(); tot = new_stats(); HITS = []
    args = [(s0 + i, modes[i % len(modes)]) for i in range(n)]
    with Pool(nw) as pool:
        for k, (st, hits) in enumerate(pool.imap_unordered(one, args, chunksize=4)):
            for kk, v in st.items(): tot[kk] += v
            HITS.extend(hits)
            if k % 2000 == 0:
                print(k, dict(tot), round(time.time()-t0,1), flush=True)
                if HITS: print("!!!! VIOLATION", flush=True)
    print("FINAL", json.dumps(dict(tot), indent=1))
    print("n_violation_hits", len(HITS))
    json.dump({"stats": dict(tot), "hits": HITS[:8]}, open(f"arm8_{s0}.json", "w"))
