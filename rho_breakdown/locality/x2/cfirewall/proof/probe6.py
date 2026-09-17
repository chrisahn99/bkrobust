"""On the NON-LIFTABLE (trial,D) pairs: what certificate works?"""
import sys; sys.path.insert(0,'/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof')
from common import *
from collections import defaultdict
import numpy as np, json
from multiprocessing import Pool

def scan_cpdag(kb_p):
    kb,p = kb_p
    C = np.frombuffer(kb,dtype=np.int8).reshape(p,p).copy()
    ref = v_structures(C); acc=defaultdict(int); wit=[]
    for G0 in reachable_mpdags(C):
        NA = nonadjacent_pairs(G0)
        if not NA: continue
        Q=[]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,np0,am0 = ostar(G0,x,y)
                if am0: Q.append((x,y,O0))
        if not Q: continue
        ext = consistent_dag_extensions(G0, ref_vstructs=ref)
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info = bk_assert(G0,[(a,b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                vH=v_structures(H)
                lift={}
                for k,D in enumerate(ext):
                    Dp=D.copy(); Dp[a,b]=1
                    lift[k] = (not has_directed_cycle(Dp)) and all(is_directed(Dp,u,v) for (u,v) in directed_edges(H)) and v_structures(Dp)==vH
                for (x,y,O0) in Q:
                    OH,nph,amH = ostar(H,x,y)
                    if not amH: continue
                    for k,D in enumerate(ext):
                        if lift[k]: acc["lift_pair"]+=1; continue
                        acc["nolift_pair"]+=1
                        cnD=cn(D,x,y); fbD=forb(D,x,y); OD=pa_set(D,cnD)-fbD
                        if OH==O0: acc["nl_OH_eq_O0"]+=1
                        if OD<=OH: acc["nl_OD_sub_OH"]+=1
                        if OH<=O0: acc["nl_OH_sub_O0"]+=1
                        if not cnD: acc["nl_cnD_empty"]+=1
                        if OH & fbD: acc["nl_BAD_forb"]+=1
                        if not valid(D,x,y,OH): acc["nl_INVALID"]+=1
                        # does some OTHER member of [G0] lift and give the same info?  no. record shape
                        if OH!=O0 and not (OD<=OH):
                            acc["nl_moved_and_not_super"]+=1
                            if len(wit)<4: wit.append((C.tolist(),G0.tolist(),D.tolist(),x,y,a,b,sorted(O0),sorted(OH),sorted(OD)))
    return dict(acc), wit

def main():
    p=int(sys.argv[1]); nw=int(sys.argv[2]) if len(sys.argv)>2 else 10
    cps={}
    for D in all_dags(p):
        C=dag_to_cpdag(D); cps.setdefault(C.tobytes(),C)
    jobs=[(k,p) for k in cps]
    tot=defaultdict(int); W=[]
    with Pool(nw) as pool:
        for acc,w in pool.imap_unordered(scan_cpdag,jobs,chunksize=4):
            for k,v in acc.items(): tot[k]+=v
            W.extend(w[:max(0,4-len(W))])
    print("p",p,json.dumps(dict(sorted(tot.items())),indent=1))
    for it in W[:4]: print("WIT",json.dumps(it))
if __name__=="__main__": main()
