"""LEMMA S (proved below) instantiated with F = de_D(cn_H).
   S:  D DAG; F with Y in F, de_D(F)=F, X notin F, ch_D(X) n F subset cn_D;
       Z n F = {} and pa_D(F)\F subset Z  ==>  Z blocks every X..Y path in D^pbd.
   Test the five conditions on the census, and certify soundness against the oracle."""
import sys, json
from collections import defaultdict
from multiprocessing import Pool
import numpy as np
from jlib import *

def scan(C):
    p=C.shape[0]; ref=v_structures(C); acc=defaultdict(int)
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
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info=bk_assert(G0,[(a,b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                for (x,y,O0) in Q:
                    O1,n1,am1,cn1,fb1=ostar(H,x,y)
                    if not am1: continue
                    Z=set(O1)
                    for D in ext:
                        F=poss_de(D,cn1)               # DAG: de_D(cn_H), closed by construction
                        cnD=causal_nodes(D,x,y)
                        fbD=(poss_de(D,cnD)|{x}) if cnD else {x}
                        chX={w for w in range(p) if is_directed(D,x,w)}
                        c2 = x not in F
                        c3 = (chX & F) <= cnD
                        c4 = not (Z & F)
                        c5 = (parents_of_set(D,F)-F) <= Z
                        c0 = not (Z & fbD)
                        acc["n"]+=1
                        acc["c2"]+=c2; acc["c3"]+=c3; acc["c4"]+=c4; acc["c5"]+=c5
                        allc = c0 and c2 and c3 and c4 and c5
                        acc["S_fires"]+=allc
                        if allc and not is_valid_adjustment_set(D,x,y,Z): acc["S_UNSOUND"]+=1
                        # cross with the earlier certificates
                        ok=lifts(D,a,b,H)[0]
                        R=reach_outside_forb(D,x,y,cnD,fbD) if cnD else set()
                        need=((parents_of_set(D,cnD)-fbD)&R) if cnD else None
                        L=(need is not None) and c0 and (need<=Z)
                        paX={u for u in range(p) if is_directed(D,u,x)}
                        Cc=(paX<=Z) and not (Z & poss_de(D,{x}))
                        acc["ALC"]+= (ok or L or Cc)
                        acc["ALCS"]+= (ok or L or Cc or allc)
                        if not (ok or L or Cc or allc): acc["OPEN"]+=1
    return dict(acc)

if __name__=="__main__":
    p=int(sys.argv[1]); acc=defaultdict(int)
    with Pool(10) as pool:
        for a in pool.imap_unordered(scan,cpdags(p),chunksize=2):
            for k,v in a.items(): acc[k]+=v
    print("p=%d"%p, json.dumps(dict(sorted(acc.items())),indent=1))
