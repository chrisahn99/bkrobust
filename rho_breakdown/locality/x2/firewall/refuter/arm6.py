"""Arm 6: detector-power control + hypothesis ablations + the k=2 boundary.

variants
  main    : the theorem exactly (1 non-adjacent stamp, all three checks)
  control : same but the Dor-Tarsi check is IGNORED (stamp anyway)  -> must fire
  nocyc   : ignore only the directed-cycle check
  k2      : TWO non-adjacent stamps, all three checks on the final H
  k3      : THREE non-adjacent stamps, all three checks
  nomeek  : G0 is a PDAG that is NOT Meek-closed (drop the closure after the
            partial orientation) -- tests whether Meek-closedness of G0 is
            load-bearing
"""
import sys, time, json, random
from itertools import combinations
sys.path.insert(0, "/private/tmp/claude-501/-Users-josecosta-mugango/bceac0a4-93a7-4757-b3b3-a3d4e6c42d55/scratchpad/cfw_ref_A")
from cfwlib import *
import cfwlib as _CF
from arm2 import gen
from multiprocessing import Pool


def one(args):
    seed, mode, variant = args
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
            if variant == "nomeek":
                Hc = G0.copy(); Hc[b, a] = 0
                if has_directed_cycle(Hc): continue
                G0 = Hc
            else:
                Hc = G0.copy(); Hc[b, a] = 0; Hc = meek_closure(Hc)
                if has_directed_cycle(Hc) or v_structures(Hc) != ref: continue
                G0 = Hc
    U0 = undirected_edges(G0); S = skeleton(G0)
    NA = [(a, b) for a in range(p) for b in range(a + 1, p) if S[a, b] == 0]
    if not NA or len(U0) > 15: return dict(st), hits
    full = consistent_dag_extensions(G0, ref_vstructs=ref)
    if not full: return dict(st), hits
    st["n_scm"] += 1; st["class_size"] += len(full)
    kk = {"k2": 2, "k3": 3}.get(variant, 1)
    if kk > 1 and len(NA) < kk: return dict(st), hits
    # candidate statement sets
    if kk == 1:
        Ks = [[(a, b)] for (a0, b0) in NA for (a, b) in ((a0, b0), (b0, a0))]
    else:
        Ks = []
        for _ in range(60):
            idx = rng.choice(len(NA), size=kk, replace=False)
            K = []
            for t in idx:
                a, b = NA[int(t)]
                K.append((int(a), int(b)) if rng.random() < 0.5 else (int(b), int(a)))
            Ks.append(K)
    for x in range(p):
        for y in range(p):
            if x == y: continue
            O0, np0, amen0 = X.ostar_and_paths(G0, x, y)
            if not amen0: continue
            st["n_query_pos"] += 1
            for K in Ks:
                H, info = x1_ops.bk_assert(G0, K)
                st["n_trial"] += 1
                cyc = has_directed_cycle(H)
                ext = x1_ops.pdag_extendable(H)
                if variant == "control":
                    ok = (not info["conflict"]) and (not ext)      # DT FAILS (cycles included)
                elif variant == "control_acyc":
                    ok = (not info["conflict"]) and (not cyc) and (not ext)
                elif variant == "nocyc":
                    ok = (not info["conflict"]) and cyc and ext
                else:
                    ok = (not info["conflict"]) and (not cyc) and ext
                if not ok:
                    st["n_rejected"] += 1; continue
                st["n_stamped"] += 1
                O1, np1, amen1 = X.ostar_and_paths(H, x, y)
                if not amen1:
                    st["n_abort"] += 1; continue
                st["n_pos"] += 1
                if O1 == O0:
                    st["n_same"] += 1; continue
                st["n_moved"] += 1
                bad = [Dc for Dc in full if not valid_in(Dc, x, y, O1)]
                st["n_Dcheck"] += len(full)
                if bad:
                    st["n_violation"] += 1
                    if len(hits) < 3:
                        hits.append(dict(variant=variant, p=p, seed=int(seed),
                                         G0=G0.tolist(), H=H.tolist(), K=[list(map(int,e)) for e in K],
                                         x=int(x), y=int(y),
                                         O0=sorted(map(int, O0)), O1=sorted(map(int, O1)),
                                         Dbad=bad[0].tolist()))
    _CF._vcache.clear()
    return dict(st), hits


if __name__ == "__main__":
    n = int(sys.argv[1]); s0 = int(sys.argv[2]); nw = int(sys.argv[3])
    variant = sys.argv[4]; modes = sys.argv[5].split(",")
    t0 = time.time(); tot = new_stats(); HITS = []
    args = [(s0 + i, modes[i % len(modes)], variant) for i in range(n)]
    with Pool(nw) as pool:
        for k, (st, h) in enumerate(pool.imap_unordered(one, args, chunksize=4)):
            for kk2, v in st.items(): tot[kk2] += v
            HITS.extend(h)
            if k % 3000 == 0: print(k, variant, dict(tot), round(time.time()-t0,1), flush=True)
    print("FINAL", variant, json.dumps(dict(tot), indent=1))
    print("n_violation", tot.get("n_violation", 0), "of moved", tot.get("n_moved", 0))
    json.dump({"variant": variant, "stats": dict(tot), "hits": HITS[:4]},
              open(f"arm6_{variant}_{s0}.json", "w"))
