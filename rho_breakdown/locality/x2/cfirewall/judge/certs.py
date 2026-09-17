"""Coverage of certA (lift) / certB (superset) / certC (back-door parents) over the census."""
import sys, json, time
from collections import defaultdict
from multiprocessing import Pool
import numpy as np
from jlib import *

def scan(C):
    p=C.shape[0]; ref=v_structures(C); acc=defaultdict(int); wit=[]
    for G0 in reachable_mpdags(C):
        NA=nonadjacent_pairs(G0)
        if not NA: continue
        Q=[]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,n0,am0,cn0,fb0=ostar(G0,x,y)
                if am0: Q.append((x,y,O0))
        if not Q: continue
        ext=consistent_dag_extensions(G0,ref_vstructs=ref)
        if not ext: continue
        # per-D precompute
        PD=[]
        for D in ext:
            paX={x:{u for u in range(p) if is_directed(D,u,x)} for x in range(p)}
            deX={x:poss_de(D,{x}) for x in range(p)}
            PD.append((paX,deX))
        DQ={}
        for di,D in enumerate(ext):
            for (x,y,_) in Q:
                if (di,x,y) in DQ: continue
                cnD=causal_nodes(D,x,y); fbD=(poss_de(D,cnD)|{x}) if cnD else {x}
                DQ[(di,x,y)]=(cnD,fbD,(parents_of_set(D,cnD)-fbD) if cnD else set())
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info=bk_assert(G0,[(a,b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                liftc={}
                for (x,y,O0) in Q:
                    O1,n1,am1,cn1,fb1=ostar(H,x,y)
                    if not am1: continue
                    acc["report"]+=1
                    S1=set(O1)
                    for di,D in enumerate(ext):
                        if di not in liftc: liftc[di]=lifts(D,a,b,H)[0]
                        A=liftc[di]
                        cnD,fbD,paD=DQ[(di,x,y)]
                        B=(not (S1&fbD)) and bool(cnD) and (paD<=S1)
                        paX,deX=PD[di]
                        Cc=(paX[x]<=S1) and not (S1 & deX[x])
                        acc["n"]+=1; acc["A"]+=A; acc["B"]+=B; acc["C"]+=Cc
                        acc["AB"]+= (A or B); acc["AC"]+=(A or Cc); acc["BC"]+=(B or Cc)
                        acc["ABC"]+=(A or B or Cc)
                        if not (A or B or Cc):
                            acc["UNCOV"]+=1
                            v=is_valid_adjustment_set(D,x,y,S1)
                            acc["UNCOV_invalid"]+= (not v)
                            if len(wit)<30: wit.append(dict(G0=G0.tolist(),H=H.tolist(),D=D.tolist(),x=x,y=y,a=a,b=b,O1=sorted(O1),O0=sorted(O0),valid=bool(v)))
                        # certC soundness: whenever C holds, validity must hold
                        if Cc:
                            v=is_valid_adjustment_set(D,x,y,S1)
                            if not v: acc["certC_UNSOUND"]+=1
    return dict(acc),wit

def main():
    p=int(sys.argv[1]); nw=int(sys.argv[2])
    Cs=cpdags(p); acc=defaultdict(int); W=[]; t0=time.time()
    with Pool(nw) as pool:
        for a,w in pool.imap_unordered(scan,Cs,chunksize=2):
            for k,v in a.items(): acc[k]+=v
            W.extend(w[:3])
    acc["elapsed"]=round(time.time()-t0,1)
    print(json.dumps(dict(sorted(acc.items())),indent=1))
    if W: json.dump(W[:30],open("uncov_p%d.json"%p,"w"),indent=1); print("uncovered witnesses saved",len(W))
if __name__=="__main__": main()
