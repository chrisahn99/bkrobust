"""Sharpness of the G0-amenability hypothesis, split by npaths(G0)==0."""
import sys
from collections import defaultdict
import numpy as np
BASE="/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious"
sys.path.insert(0,BASE+"/x2/code"); sys.path.insert(0,BASE+"/x2/cfirewall/code")
from graphs import consistent_dag_extensions, dag_to_cpdag, has_directed_cycle, v_structures
import adjust, x2lib as X
from run_lemma import all_dags, reachable_mpdags
import cfirewall as CF
x1=CF.x1_ops
def scan(C):
    ref=v_structures(C); acc=defaultdict(int)
    for G0 in reachable_mpdags(C):
        p=G0.shape[0]; NA=x1.nonadjacent_pairs(G0)
        if not NA: continue
        ext=consistent_dag_extensions(G0,ref_vstructs=ref)
        if not ext: continue
        QN=[]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,n0,am0=X.ostar_and_paths(G0,x,y)
                if not am0: QN.append((x,y,n0))
        if not QN: continue
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info=x1.bk_assert(G0,[(a,b)])
                if info["conflict"] or has_directed_cycle(H) or not x1.pdag_extendable(H): continue
                for (x,y,n0) in QN:
                    Z,n1,am1=X.ostar_and_paths(H,x,y)
                    if not am1: continue
                    lbl = "nopath" if n0==0 else "unamen"
                    nbad=sum(1 for D in ext if not adjust.is_valid_adjustment_set(D,x,y,set(Z)))
                    acc[(lbl,"rep")]+=1
                    acc[(lbl,"allvalid" if nbad==0 else ("allbad" if nbad==len(ext) else "somebad"))]+=1
    return {str(k):v for k,v in acc.items()}
def _job(arg):
    kb,p=arg
    return scan(np.frombuffer(kb,dtype=np.int8).reshape(p,p).copy())
if __name__=="__main__":
    p=int(sys.argv[1]); cp={}
    for D in all_dags(p): C=dag_to_cpdag(D); cp.setdefault(C.tobytes(),C)
    from multiprocessing import Pool
    tot=defaultdict(int)
    with Pool(4) as pool:
        for acc in pool.imap_unordered(_job,[(k,p) for k in cp],chunksize=4):
            for k,v in acc.items(): tot[k]+=v
    for k in sorted(tot,key=str): print(f"{k:30s} {tot[k]}")
