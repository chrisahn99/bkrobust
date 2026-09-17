"""Enumerate the ENTIRE p=5 open set (certA v certL v certC all fail) and profile it."""
import sys, json
from collections import defaultdict, Counter
from multiprocessing import Pool
import numpy as np
from jlib import *

def scan(C):
    p=C.shape[0]; ref=v_structures(C); out=[]; acc=defaultdict(int)
    for G0 in reachable_mpdags(C):
        NA=nonadjacent_pairs(G0)
        if not NA: continue
        Q=[]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,n0,am0,cn0,fb0=ostar(G0,x,y)
                if am0: Q.append((x,y,O0,cn0))
        if not Q: continue
        ext=consistent_dag_extensions(G0,ref_vstructs=ref)
        if not ext: continue
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info=bk_assert(G0,[(a,b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                for (x,y,O0,cn0) in Q:
                    O1,n1,am1,cn1,fb1=ostar(H,x,y)
                    if not am1: continue
                    S1=set(O1)
                    for D in ext:
                        ok,Dp,why=lifts(D,a,b,H)
                        if ok: continue
                        cnD=causal_nodes(D,x,y); fbD=(poss_de(D,cnD)|{x}) if cnD else {x}
                        R=reach_outside_forb(D,x,y,cnD,fbD) if cnD else set()
                        need=((parents_of_set(D,cnD)-fbD)&R) if cnD else None
                        if (need is not None) and (not(S1&fbD)) and (need<=S1): continue
                        paX={u for u in range(p) if is_directed(D,u,x)}
                        if (paX<=S1) and not (S1&poss_de(D,{x})): continue
                        acc["open"]+=1
                        acc["why_"+why]+=1
                        acc["cnD_sub_cnH"]+= (cnD<=cn1)
                        acc["cnH_sub_cnD"]+= (cn1<=cnD)
                        acc["fbD_sub_fbH"]+= (fbD<=fb1)
                        acc["O1_sub_O0"]+= (S1<set(O0))
                        acc["misorients"]+= any(not is_directed(D,i,j) for (i,j) in directed_edges(H) if (i,j)!=(a,b))
                        acc["dropped_in_cnH"]+= ((set(O0)-S1)<=cn1)
                        acc["valid"]+= bool(is_valid_adjustment_set(D,x,y,S1))
                        if len(out)<12:
                            out.append(dict(G0=G0.tolist(),H=H.tolist(),D=D.tolist(),x=x,y=y,a=a,b=b,O0=sorted(O0),O1=sorted(O1),why=why))
    return dict(acc),out

if __name__=="__main__":
    Cs=cpdags(5); acc=defaultdict(int); OUT=[]
    with Pool(10) as pool:
        for a,o in pool.imap_unordered(scan,Cs,chunksize=2):
            for k,v in a.items(): acc[k]+=v
            OUT.extend(o[:2])
    print(json.dumps(dict(sorted(acc.items())),indent=1))
    json.dump(OUT[:200],open("OPEN_SET_p5.json","w"))
