"""Finer containments between D-world and H-world."""
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
                    if not amH: continue
                    cnH=cn(H,x,y); fbH=forb(H,x,y); paH=pa_set(H,cnH)
                    for D in ext:
                        cnD = cn(D,x,y); fbD = forb(D,x,y)
                        if not cnD: continue
                        OD = pa_set(D,cnD)-fbD
                        acc["n"]+=1
                        if not (cnD<=cnH): acc["cnD_notsub_cnH"]+=1
                        if not (cnD<=cnH|fbH): acc["cnD_notsub_cnH_u_fbH"]+=1
                        if OD & fbH: acc["OD_meets_fbH"]+=1
                        if not (OD<=paH): acc["OD_notsub_paH"]+=1
                        # the *needed* statement N: entry parents
                        # v in pa_D(cnD)\fbD  that has a child in cnD, and is not in cnD
                        for v in OD:
                            if v in OH: continue
                            acc["OD_minus_OH"]+=1
                            if v in fbH: acc["OD_minus_OH_inFbH"]+=1
                            if v not in paH: acc["OD_minus_OH_notPaH"]+=1
                            # is v connected to X at all in D^pbd?
    return dict(acc), {}

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
