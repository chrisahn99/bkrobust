"""Is H a genuine MPDAG?  Is O*_H valid in every D+ in [H]?  How do [H] and [G0] relate?"""
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
        if not ext: continue
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info = bk_assert(G0,[(a,b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                acc["stmt"]+=1
                extH = consistent_dag_extensions(H, ref_vstructs=v_structures(H))
                acc["extH_empty"] += (len(extH)==0)
                if extH:
                    cog = common_orientation_graph(extH, skeleton(H))
                    if not np.array_equal(cog,H): acc["H_not_MPDAG"]+=1
                # does H contradict D on an S-edge?
                for D in ext:
                    contr = any(is_directed(H,u,v) and is_directed(D,v,u) for (u,v) in directed_edges(H))
                    acc["D_contradicted_by_H"] += contr
                    acc["D_total"]+=1
                    # is there D+ in [H] whose S-restriction equals D?
                    ok=False
                    for Dp in extH:
                        Dm=Dp.copy(); Dm[a,b]=0; Dm[b,a]=0
                        if np.array_equal(Dm,D): ok=True; break
                    acc["D_liftable_to_H"] += ok
                for (x,y,O0) in Q:
                    OH,nph,amH = ostar(H,x,y)
                    if not amH: continue
                    acc["trial"]+=1
                    for Dp in extH:
                        acc["hcheck"]+=1
                        if not valid(Dp,x,y,OH): acc["OH_invalid_in_extH"]+=1
    return dict(acc), wit

def main():
    p=int(sys.argv[1]); nw=int(sys.argv[2]) if len(sys.argv)>2 else 10
    cps={}
    for D in all_dags(p):
        C=dag_to_cpdag(D); cps.setdefault(C.tobytes(),C)
    jobs=[(k,p) for k in cps]
    tot=defaultdict(int)
    with Pool(nw) as pool:
        for acc,w in pool.imap_unordered(scan_cpdag,jobs,chunksize=4):
            for k,v in acc.items(): tot[k]+=v
    print("p",p,json.dumps(dict(sorted(tot.items())),indent=1))
if __name__=="__main__": main()
