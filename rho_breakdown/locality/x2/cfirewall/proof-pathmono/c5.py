import sys, json
from collections import defaultdict
import numpy as np
sys.path.insert(0,"/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof-pathmono")
from lib import *

def scan_cpdag(C):
    ref=v_structures(C); acc=defaultdict(int); wit=[]
    for G0 in reachable_mpdags(C):
        NA=nonadjacent_pairs(G0)
        if not NA: continue
        p=G0.shape[0]
        Q=[]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,_,am0=X.ostar_and_paths(G0,x,y)
                if am0: Q.append((x,y,O0,adjust.forb(G0,x,y)))
        if not Q: continue
        ext=consistent_dag_extensions(G0,ref_vstructs=ref)
        if not ext: continue
        ODc={}
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info=bk_assert(G0,[(a,b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                for (x,y,O0,fb0) in Q:
                    O1,_,am1=X.ostar_and_paths(H,x,y)
                    if not am1: continue
                    moved=(O1!=O0)
                    for iD,D in enumerate(ext):
                        acc["trial"]+=1
                        if (iD,x,y) not in ODc:
                            ODc[(iD,x,y)]=(adjust.optimal_adjustment_set(D,x,y),adjust.forb(D,x,y))
                        OD,fbD=ODc[(iD,x,y)]
                        Dp=D.copy(); Dp[a,b]=1; Dp[b,a]=0
                        A=(not has_directed_cycle(Dp)) and in_class(Dp,H)
                        S=(OD is not None) and set(OD)<=set(O1)
                        if set(O1)&fb0: acc["OH_meets_fb0"]+=1
                        if set(O1)&fbD: acc["OH_meets_fbD"]+=1
                        acc["A"]+=int(A); acc["S"]+=int(S); acc["AorS"]+=int(A or S)
                        if not (A or S):
                            acc["resid"]+=1
                            acc["resid_moved"]+=int(moved)
                            acc["resid_a_is_x"]+=int(a==x); acc["resid_b_is_x"]+=int(b==x)
                            acc["resid_a_is_y"]+=int(a==y); acc["resid_b_is_y"]+=int(b==y)
                            acc["resid_ab_touches_xy"]+=int(a in (x,y) or b in (x,y))
                            acc["resid_valid"]+=int(adjust.is_valid_adjustment_set(D,x,y,set(O1)))
                            acc["resid_OH_sub_OD"]+=int(set(O1)<=set(OD))
                            acc["resid_Dp_cyclic"]+=int(bool(has_directed_cycle(Dp)))
                            if len(wit)<200:
                                wit.append(dict(G0=G0.tolist(),H=H.tolist(),D=D.tolist(),x=int(x),y=int(y),
                                  a=int(a),b=int(b),O0=sorted(map(int,O0)),O1=sorted(map(int,O1)),
                                  OD=sorted(map(int,OD)),moved=bool(moved)))
    return dict(acc),wit

def _job(args):
    kb,p=args
    return scan_cpdag(np.frombuffer(kb,dtype=np.int8).reshape(p,p).copy())

def main():
    p=int(sys.argv[1]); nw=int(sys.argv[2]) if len(sys.argv)>2 else 10
    cp={}
    for D in all_dags(p):
        C=dag_to_cpdag(D); cp.setdefault(C.tobytes(),C)
    print(len(cp),"cpdags",flush=True)
    from multiprocessing import Pool
    tot=defaultdict(int); W=[]; n=0
    with Pool(nw) as pool:
        for acc,w in pool.imap_unordered(_job,[(k,p) for k in cp],chunksize=4):
            for k,v in acc.items(): tot[k]+=v
            if len(W)<200: W.extend(w[:200-len(W)])
            n+=1
            if n%1000==0: print(n,dict(tot).get("trial"),dict(tot).get("resid"),flush=True)
    print("p=",p,json.dumps(dict(sorted(tot.items())),indent=1))
    json.dump(dict(counts={k:int(v) for k,v in tot.items()},wit=W),open("c5_%d.json"%p,"w"),indent=1)
if __name__=="__main__": main()
