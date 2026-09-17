"""Arm 9: is the AMENABILITY ABORT (disjunct (i)) load-bearing?

Take exactly the trials the theorem lets escape through (i) -- H passes the
three coherence checks but is NOT amenable rel (X,Y) -- and report
pa_H(cn_H) \\ forb_H anyway.  Count how often that set is invalid in some
D in [G0].  A high rate means clause (i) is doing real work, not decoration."""
import sys, time, json
sys.path.insert(0, "/private/tmp/claude-501/-Users-josecosta-mugango/bceac0a4-93a7-4757-b3b3-a3d4e6c42d55/scratchpad/cfw_ref_A")
from cfwlib import *
import cfwlib as _CF
from arm2 import gen
from multiprocessing import Pool

def ostar_ignore_amen(G, x, y):
    paths = X.pcp_capped(G, x, y)
    if not paths: return None
    cn = set()
    for pth in paths: cn.update(pth[1:])
    fb = adjust.poss_de(G, cn) | {x}
    return frozenset(adjust.parents_of_set(G, cn) - fb)

def one(args):
    seed, mode = args
    rng = np.random.default_rng(seed); st = new_stats()
    D0 = gen(rng, mode)
    if D0 is None: return dict(st)
    p = D0.shape[0]
    C = dag_to_cpdag(D0); ref = v_structures(C); G0 = C.copy()
    U = undirected_edges(C)
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
    if not NA or len(U0) > 15: return dict(st)
    full = consistent_dag_extensions(G0, ref_vstructs=ref)
    if not full: return dict(st)
    st["n_scm"] += 1
    for x in range(p):
        for y in range(p):
            if x == y: continue
            O0, np0, amen0 = X.ostar_and_paths(G0, x, y)
            if not amen0: continue
            for (a0, b0) in NA:
                for (a, b) in ((a0, b0), (b0, a0)):
                    H, info = x1_ops.bk_assert(G0, [(a, b)])
                    if info["conflict"] or has_directed_cycle(H) or not x1_ops.pdag_extendable(H):
                        continue
                    O1, np1, amen1 = X.ostar_and_paths(H, x, y)
                    if amen1: continue
                    st["n_abort"] += 1
                    if np1 == 0:
                        st["n_abort_zero"] += 1
                        # ZERO report: claims no effect. Sound iff cn_D empty.
                        for Dc in full:
                            st["n_zero_check"] += 1
                            if adjust.causal_nodes(Dc, x, y):
                                st["n_zero_unsound"] += 1; break
                        continue
                    Oz = ostar_ignore_amen(H, x, y)
                    if Oz is None: continue
                    st["n_abort_unamen"] += 1
                    if Oz == O0: st["n_same"] += 1; continue
                    st["n_moved"] += 1
                    st["n_Dcheck"] += len(full)
                    if any(not valid_in(Dc, x, y, Oz) for Dc in full):
                        st["n_violation"] += 1
    _CF._vcache.clear()
    return dict(st)

if __name__ == "__main__":
    n = int(sys.argv[1]); s0 = int(sys.argv[2]); nw = int(sys.argv[3]); modes = sys.argv[4].split(",")
    t0 = time.time(); tot = new_stats()
    args = [(s0 + i, modes[i % len(modes)]) for i in range(n)]
    with Pool(nw) as pool:
        for k, st in enumerate(pool.imap_unordered(one, args, chunksize=4)):
            for kk, v in st.items(): tot[kk] += v
            if k % 2000 == 0: print(k, dict(tot), round(time.time()-t0,1), flush=True)
    print("FINAL", json.dumps(dict(tot), indent=1))
