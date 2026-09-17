import sys
from collections import defaultdict
import numpy as np
sys.path.insert(0,".")
from lib import *
def ostar_raw(G,x,y):
    cn=adjust.causal_nodes(G,x,y)
    if not cn: return None
    fb=adjust.poss_de(G,cn)|{x}
    return frozenset(adjust.parents_of_set(G,cn)-fb)
def job(args):
    kb,p=args
    C=np.frombuffer(kb,dtype=np.int8).reshape(p,p).copy(); ref=v_structures(C); acc=defaultdict(int)
    for G0 in reachable_mpdags(C):
        NA=nonadjacent_pairs(G0)
        if not NA: continue
        ext=consistent_dag_extensions(G0,ref_vstructs=ref)
        if not ext: continue
        Q=[(x,y) for x in range(p) for y in range(p) if x!=y and X.ostar_and_paths(G0,x,y)[2]]
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info=bk_assert(G0,[(a,b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                for (x,y) in Q:
                    _,npx,am1=X.ostar_and_paths(H,x,y)
                    if am1: continue
                    acc["abort"]+=1
                    acc["abort_nopath"]+=int(npx==0)
                    Oraw=ostar_raw(H,x,y)
                    if Oraw is None: acc["raw_none"]+=1; continue
                    for D in ext:
                        acc["ungated_chk"]+=1
                        if not adjust.is_valid_adjustment_set(D,x,y,set(Oraw)): acc["ungated_INVALID"]+=1
    return dict(acc)
if __name__=="__main__":
    p=int(sys.argv[1]); cp={}
    for D in all_dags(p):
        C=dag_to_cpdag(D); cp.setdefault(C.tobytes(),C)
    from multiprocessing import Pool
    tot=defaultdict(int)
    with Pool(10) as pool:
        for a in pool.imap_unordered(job,[(k,p) for k in cp],chunksize=4):
            for k,v in a.items(): tot[k]+=v
    print("p=",p,dict(sorted(tot.items())),flush=True)
