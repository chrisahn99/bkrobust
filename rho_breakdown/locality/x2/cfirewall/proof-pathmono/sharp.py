"""Sharpness witnesses at p=4 (smallest found)."""
import sys, json
from collections import defaultdict
import numpy as np
sys.path.insert(0,".")
from lib import *

def es(G):
    G=np.array(G,dtype=np.int8); p=G.shape[0]; o=[]
    for i in range(p):
        for j in range(p):
            if is_directed(G,i,j): o.append("%d->%d"%(i,j))
    for i in range(p):
        for j in range(i+1,p):
            if is_undirected(G,i,j): o.append("%d--%d"%(i,j))
    return " ".join(o)

def ostar_raw(G,x,y):
    """pa(cn)\\forb with NO amenability gate."""
    cn=adjust.causal_nodes(G,x,y)
    if not cn: return None
    fb=adjust.poss_de(G,cn)|{x}
    return frozenset(adjust.parents_of_set(G,cn)-fb)

W={"S1_amenability":[], "S2_coherence":[]}
for p in (3,4):
    cp={}
    for D in all_dags(p):
        C=dag_to_cpdag(D); cp.setdefault(C.tobytes(),C)
    for kb,C in cp.items():
        ref=v_structures(C)
        for G0 in reachable_mpdags(C):
            NA=nonadjacent_pairs(G0)
            if not NA: continue
            ext=consistent_dag_extensions(G0,ref_vstructs=ref)
            if not ext: continue
            Q=[]
            for x in range(p):
                for y in range(p):
                    if x!=y:
                        O0,_,am0=X.ostar_and_paths(G0,x,y)
                        if am0: Q.append((x,y,O0))
            for (a0,b0) in NA:
                for (a,b) in ((a0,b0),(b0,a0)):
                    H,info=bk_assert(G0,[(a,b)])
                    coherent = (not info["conflict"]) and (not has_directed_cycle(H)) and pdag_extendable(H)
                    for (x,y,O0) in Q:
                        O1,_,am1=X.ostar_and_paths(H,x,y)
                        Oraw=ostar_raw(H,x,y)
                        if coherent and not am1 and Oraw is not None:
                            bad=[D for D in ext if not adjust.is_valid_adjustment_set(D,x,y,set(Oraw))]
                            if bad and len(W["S1_amenability"])<3:
                                W["S1_amenability"].append(dict(p=p,G0=es(G0),H=es(H),stmt="%d->%d"%(a,b),
                                    x=x,y=y,Ostar_G0=sorted(map(int,O0)),Ostar_H_nogate=sorted(map(int,Oraw)),
                                    D_where_invalid=es(bad[0])))
                        if (not coherent) and am1:
                            bad=[D for D in ext if not adjust.is_valid_adjustment_set(D,x,y,set(O1))]
                            if bad and len(W["S2_coherence"])<3:
                                why=[]
                                if info["conflict"]: why.append("conflict")
                                if has_directed_cycle(H): why.append("cycle")
                                if not pdag_extendable(H): why.append("not-extendable")
                                W["S2_coherence"].append(dict(p=p,G0=es(G0),H=es(H),stmt="%d->%d"%(a,b),
                                    x=x,y=y,fails=why,Ostar_G0=sorted(map(int,O0)),
                                    Ostar_H=sorted(map(int,O1)),D_where_invalid=es(bad[0])))
    if all(len(v)>=1 for v in W.values()): break
print(json.dumps(W,indent=1))
