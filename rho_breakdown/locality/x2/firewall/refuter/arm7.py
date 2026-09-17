"""Arm 7: LARGE p (10-16). [G0] is sampled, not enumerated -- every draw is
biased to orient the undirected edges AGAINST whatever H forced on them, which
is the exact mechanism a counterexample would need."""
import sys, time, json, random
sys.path.insert(0, "/private/tmp/claude-501/-Users-josecosta-mugango/bceac0a4-93a7-4757-b3b3-a3d4e6c42d55/scratchpad/cfw_ref_A")
from cfwlib import *
import cfwlib as _CF
from rand_arm import extend_biased, flip_pref
from arm2 import scale_free_dag, chain_dag
from multiprocessing import Pool


def one(args):
    seed, mode = args
    rng = np.random.default_rng(seed); pyr = random.Random(seed + 5)
    st = new_stats(); hits = []
    p = int(rng.choice([10, 12, 14, 16]))
    if mode == "er":
        D0 = random_dag(p, float(rng.choice([1.5, 2.0, 2.5, 3.0])), rng)
    elif mode == "sf":
        D0 = scale_free_dag(p, int(rng.choice([1, 2])), rng)
    else:
        D0 = chain_dag(p, rng, extra=int(rng.integers(0, 4)))
    if D0 is None: return dict(st), hits
    C = dag_to_cpdag(D0); ref = v_structures(C)
    U = undirected_edges(C); G0 = C.copy()
    if len(U) and rng.random() < 0.6:
        k = int(rng.integers(0, len(U)))
        for t in rng.permutation(len(U))[:k]:
            u, v = U[t]
            if not is_undirected(G0, u, v): continue
            a, b = (u, v) if rng.random() < 0.5 else (v, u)
            Hc = G0.copy(); Hc[b, a] = 0; Hc = meek_closure(Hc)
            if has_directed_cycle(Hc) or v_structures(Hc) != ref: continue
            G0 = Hc
    S = skeleton(G0)
    NA = [(a, b) for a in range(p) for b in range(a + 1, p) if S[a, b] == 0]
    if not NA: return dict(st), hits
    st["n_scm"] += 1; st["p_sum"] += p
    # sample queries and stamps to keep the per-graph cost bounded
    QS = [(int(x), int(y)) for x in range(p) for y in range(p) if x != y]
    pyr.shuffle(QS); QS = QS[:40]
    pyr.shuffle(NA); NAs = NA[:25]
    for (x, y) in QS:
        try:
            O0, np0, amen0 = X.ostar_and_paths(G0, x, y, cap=20000)
        except Exception:
            st["n_pathcap"] += 1; continue
        if not amen0: continue
        st["n_query_pos"] += 1
        for (a0, b0) in NAs:
            for (a, b) in ((a0, b0), (b0, a0)):
                H, info = x1_ops.bk_assert(G0, [(a, b)])
                st["n_trial"] += 1
                if info["conflict"] or has_directed_cycle(H) or not x1_ops.pdag_extendable(H):
                    st["n_rejected"] += 1; continue
                st["n_stamped"] += 1
                try:
                    O1, np1, amen1 = X.ostar_and_paths(H, x, y, cap=20000)
                except Exception:
                    st["n_pathcap"] += 1; continue
                if not amen1:
                    st["n_abort"] += 1; continue
                st["n_pos"] += 1
                if O1 == O0:
                    st["n_same"] += 1; continue
                st["n_moved"] += 1
                pref = flip_pref(G0, H); seen = set(); cands = []
                for _ in range(6):
                    for pf in (pref, pref, {}):
                        Dc = extend_biased(G0, ref, pf, pyr, tries=3)
                        if Dc is not None and Dc.tobytes() not in seen:
                            seen.add(Dc.tobytes()); cands.append(Dc)
                st["n_Dcheck"] += len(cands)
                bad = [Dc for Dc in cands if not valid_in(Dc, x, y, O1)]
                if bad:
                    st["n_violation"] += 1
                    hits.append(dict(p=p, seed=int(seed), G0=G0.tolist(), H=H.tolist(),
                                     x=int(x), y=int(y), edge=[int(a), int(b)],
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
        for k, (st, h) in enumerate(pool.imap_unordered(one, args, chunksize=2)):
            for kk, v in st.items(): tot[kk] += v
            HITS.extend(h)
            if k % 200 == 0:
                print(k, dict(tot), round(time.time()-t0,1), flush=True)
                if HITS: print("!!!! VIOLATION", flush=True)
    print("FINAL", json.dumps(dict(tot), indent=1))
    print("n_violation_hits", len(HITS))
    json.dump({"stats": dict(tot), "hits": HITS[:6]}, open(f"arm7_{s0}.json","w"))
