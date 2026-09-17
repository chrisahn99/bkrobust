"""Why is D not liftable?  cycle / contradiction / new-v-structure."""
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
        anyQ=False
        Q=[]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,np0,am0 = ostar(G0,x,y)
                if am0: Q.append((x,y))
        if not Q: continue
        ext = consistent_dag_extensions(G0, ref_vstructs=ref)
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info = bk_assert(G0,[(a,b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                vH = v_structures(H)
                for D in ext:
                    acc["D_total"]+=1
                    Dp = D.copy(); Dp[a,b]=1
                    cyc = has_directed_cycle(Dp)
                    contr = not all(is_directed(Dp,u,v) for (u,v) in directed_edges(H))
                    vs = (v_structures(Dp)!=vH)
                    lift = (not cyc) and (not contr) and (not vs)
                    if lift: acc["lift"]+=1; continue
                    acc["nolift"]+=1
                    key = ("cyc" if cyc else "")+("|contr" if contr else "")+("|vs" if vs else "")
                    acc["nolift_"+key]+=1
                    if key=="|vs" and len(wit)<4:
                        wit.append((C.tolist(),G0.tolist(),D.tolist(),a,b,sorted(map(list,vH)),sorted(map(list,v_structures(Dp)))))
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
