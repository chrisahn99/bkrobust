"""Sharpness of the 'visible abort' clause: on the aborts (H coherent but NOT
amenable rel (X,Y)), what if the tool reported pa_H(cn_H) \ forb_H anyway?"""
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
        Q=[(x,y) for x in range(p) for y in range(p)
           if x!=y and X.ostar_and_paths(G0,x,y)[2]]
        if not Q: continue
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info=x1.bk_assert(G0,[(a,b)])
                if info["conflict"] or has_directed_cycle(H) or not x1.pdag_extendable(H): continue
                for (x,y) in Q:
                    Z,n1,am1=X.ostar_and_paths(H,x,y)
                    if am1: continue
                    acc["abort"]+=1
                    if n1==0: acc["abort_nopath"]+=1; continue
                    cn=adjust.causal_nodes(H,x,y)
                    S=adjust.parents_of_set(H,cn)-adjust.forb(H,x,y)
                    nbad=sum(1 for D in ext if not adjust.is_valid_adjustment_set(D,x,y,S))
                    acc["would_report"]+=1
                    acc["would_be_wrong"]+= (nbad>0)
                    acc["would_be_wrong_allD"]+= (nbad==len(ext))
    return {k:v for k,v in acc.items()}
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
    for k in sorted(tot): print(f"{k:22s} {tot[k]}")
