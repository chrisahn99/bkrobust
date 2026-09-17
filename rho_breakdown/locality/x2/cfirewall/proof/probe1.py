"""Census probe: test R1 (F1), R2 (O*_D subseteq O*_H), R3 (cn_D nonempty)."""
import sys; sys.path.insert(0,'/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof')
from common import *
from collections import defaultdict
import numpy as np, json
from multiprocessing import Pool

def scan_cpdag(kb_p):
    kb,p = kb_p
    C = np.frombuffer(kb,dtype=np.int8).reshape(p,p).copy()
    ref = v_structures(C); acc=defaultdict(int); wit=defaultdict(list)
    for G0 in reachable_mpdags(C):
        NA = nonadjacent_pairs(G0)
        if not NA: continue
        Q=[]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,np0,am0 = ostar(G0,x,y)
                if am0: Q.append((x,y,O0,cn(G0,x,y),forb(G0,x,y)))
        if not Q: continue
        ext = consistent_dag_extensions(G0, ref_vstructs=ref)
        if not ext: continue
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info = bk_assert(G0,[(a,b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                for (x,y,O0,cn0,fb0) in Q:
                    OH,nph,amH = ostar(H,x,y)
                    if not amH:
                        acc["abort"]+=1; continue
                    acc["trial"]+=1
                    if OH!=O0: acc["moved"]+=1
                    # F1
                    if OH & fb0:
                        acc["F1_viol"]+=1
                        if len(wit["F1"])<3: wit["F1"].append((C.tolist(),G0.tolist(),x,y,a,b,sorted(OH),sorted(fb0)))
                    for D in ext:
                        cnD = cn(D,x,y); fbD = forb(D,x,y)
                        acc["dcheck"]+=1
                        if not cnD:
                            acc["cnD_empty"]+=1
                            if len(wit["cnD"])<3: wit["cnD"].append((C.tolist(),G0.tolist(),D.tolist(),x,y,a,b))
                            continue
                        OD = pa_set(D,cnD)-fbD
                        if not (OD <= OH):
                            acc["R2_viol"]+=1
                            if len(wit["R2"])<3: wit["R2"].append((C.tolist(),G0.tolist(),D.tolist(),x,y,a,b,sorted(OD),sorted(OH)))
                        if OH & fbD:
                            acc["R1_viol"]+=1
                        if not (fbD <= fb0):
                            acc["fbD_not_sub_fb0"]+=1
    return dict(acc), {k:v for k,v in wit.items()}

def main():
    p=int(sys.argv[1]); nw=int(sys.argv[2]) if len(sys.argv)>2 else 10
    cps={}
    for D in all_dags(p):
        C=dag_to_cpdag(D); cps.setdefault(C.tobytes(),C)
    jobs=[(k,p) for k in cps]
    tot=defaultdict(int); W=defaultdict(list)
    with Pool(nw) as pool:
        for acc,w in pool.imap_unordered(scan_cpdag,jobs,chunksize=4):
            for k,v in acc.items(): tot[k]+=v
            for k,v in w.items():
                if len(W[k])<3: W[k].extend(v[:3-len(W[k])])
    print("p",p,json.dumps(dict(sorted(tot.items())),indent=1))
    for k,v in W.items():
        print("### WITNESS",k); 
        for it in v[:2]: print(json.dumps(it))
if __name__=="__main__": main()
