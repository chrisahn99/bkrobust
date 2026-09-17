import sys, time, json
sys.path.insert(0, "/private/tmp/claude-501/-Users-josecosta-mugango/bceac0a4-93a7-4757-b3b3-a3d4e6c42d55/scratchpad/cfw_ref_A")
from cfwlib import *
from multiprocessing import Pool

p = int(sys.argv[1]); nw = int(sys.argv[2]) if len(sys.argv) > 2 else 8

def work(Cb):
    C = np.frombuffer(Cb, dtype=np.int8).reshape(p, p).copy()
    st = new_stats(); hits = []
    ref = v_structures(C)
    for G0 in reachable_mpdags(C):
        scan_G0(G0, st, hits, ref=ref)
    import cfwlib as _c
    if len(_c._vcache) > 400000: _c._vcache.clear()
    return dict(st), hits

if __name__ == "__main__":
    t0 = time.time()
    cp = {}
    for D in enum_dags(p):
        C = dag_to_cpdag(D); cp[C.tobytes()] = C
    keys = list(cp.keys())
    print("p", p, "n_cpdag", len(keys), "t", round(time.time()-t0,1), flush=True)
    tot = new_stats(); HITS = []
    with Pool(nw) as pool:
        for n, (st, hits) in enumerate(pool.imap_unordered(work, keys, chunksize=8)):
            for k, v in st.items(): tot[k] += v
            HITS.extend(hits)
            if n % 500 == 0:
                print(n, dict(tot), round(time.time()-t0,1), flush=True)
                if HITS: print("!!! HIT", flush=True)
    print("FINAL p=%d" % p, json.dumps(dict(tot), indent=1))
    print("elapsed", round(time.time()-t0,1), "n_hits", len(HITS))
    json.dump({"stats": dict(tot), "hits": HITS[:12]}, open(f"exh_p{p}.json","w"))
