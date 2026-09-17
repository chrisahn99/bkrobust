"""HPM instantiated on H itself: O*(x,y,H) valid in EVERY DAG of [H]. Census p=3,4."""
import sys
from collections import defaultdict
import numpy as np
sys.path.insert(0,".")
from lib import *
def job(args):
    kb,p=args
    C=np.frombuffer(kb,dtype=np.int8).reshape(p,p).copy(); ref=v_structures(C)
    acc=defaultdict(int)
    for G0 in reachable_mpdags(C):
        NA=nonadjacent_pairs(G0)
        if not NA: continue
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info=bk_assert(G0,[(a,b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                extH=consistent_dag_extensions(H,ref_vstructs=v_structures(H))
                acc["extH_nonempty"]+=int(len(extH)>0)
                for x in range(p):
                    for y in range(p):
                        if x==y: continue
                        O1,_,am1=X.ostar_and_paths(H,x,y)
                        if not am1: continue
                        for Dh in extH:
                            acc["chk"]+=1
                            if not adjust.is_valid_adjustment_set(Dh,x,y,set(O1)): acc["HPM_FAIL"]+=1
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
