import sys, time, json
sys.path.insert(0, "/private/tmp/claude-501/-Users-josecosta-mugango/bceac0a4-93a7-4757-b3b3-a3d4e6c42d55/scratchpad/cfw_ref_A")
from cfwlib import *

p = int(sys.argv[1])
t0 = time.time()
cp = {}
for D in enum_dags(p):
    C = dag_to_cpdag(D); cp[C.tobytes()] = C
print("p", p, "n_dag_enum_done n_cpdag", len(cp), "t", round(time.time()-t0,1), flush=True)
st = new_stats(); hits = []
for i, C in enumerate(cp.values()):
    for G0 in reachable_mpdags(C):
        scan_G0(G0, st, hits, ref=v_structures(C))
    if i % 25 == 0:
        print(i, dict(st), round(time.time()-t0,1), flush=True)
print("FINAL p=%d" % p, json.dumps(dict(st), indent=1))
print("elapsed", round(time.time()-t0,1), "n_hits", len(hits))
json.dump({"stats": dict(st), "hits": hits}, open(f"exh_p{p}.json", "w"))
