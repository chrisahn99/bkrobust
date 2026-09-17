"""Which component of the S1 guard is doing the work?
   Decompose rejected class-S statements into: conflict / directed cycle / Dor-Tarsi-only.
   For each, does an ALREADY-IDENTIFIED sound report get corrupted?"""
import sys
from collections import defaultdict
sys.path.insert(0,"/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/falsegain/code")
from fglib import *
import fglib, numpy as np, importlib.util
from graphs import random_dag, is_undirected, dag_agrees_with
spec=importlib.util.spec_from_file_location("x1_ops","/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x1/code/x1_ops.py")
x1_ops=importlib.util.module_from_spec(spec); spec.loader.exec_module(x1_ops)
def true_bk(C,D,rng):
    G=C.copy(); U=undirected_edges(G)
    if not U: return G
    for t in rng.permutation(len(U))[:int(rng.integers(0,len(U)+1))]:
        u,v=U[t]
        if not is_undirected(G,u,v): continue
        if D[u,v]==1: G[v,u]=0
        else: G[u,v]=0
        G=meek_closure(G)
    return G
rng=np.random.default_rng(int(sys.argv[2])); acc=defaultdict(int); wit=[]
for it in range(int(sys.argv[1])):
    p=int(rng.choice([5,6,7,8])); deg=float(rng.choice([1.5,2.0,2.5,3.0]))
    D=random_dag(p,deg,rng); C=dag_to_cpdag(D); G0=true_bk(C,D,rng)
    for (a0,b0) in x1_ops.nonadjacent_pairs(C):
        for (a,b) in ((a0,b0),(b0,a0)):
            H,info=x1_ops.bk_assert(G0,[(a,b)])
            cyc=bool(has_directed_cycle(H)); con=bool(info["conflict"]); ext=x1_ops.pdag_extendable(H)
            if con: g="conflict"
            elif cyc: g="cycle"
            elif not ext: g="dortarsi_only"
            else: g="ACCEPT"
            for x in range(p):
                for y in range(p):
                    if x==y: continue
                    s0,O0=report_state(G0,x,y)
                    if s0==REFUSE: continue
                    if not sound(D,x,y,s0,O0): continue
                    s1,O1=report_state(H,x,y); n1=sound(D,x,y,s1,O1)
                    acc[(g,s0,s1,str(n1))]+=1
                    if g=="dortarsi_only" and n1 is False and len(wit)<3:
                        wit.append((p,D.tolist(),G0.tolist(),H.tolist(),x,y,(a,b),s0,s1))
    if len(fglib._valid_cache)>150000: fglib._valid_cache.clear(); fglib._cn_cache.clear()
for k,v in sorted(acc.items()): print("|".join(k),v)
print("dortarsi_only corruption witnesses:", len(wit))
import json; json.dump(wit,open("guard_wit.json","w"),indent=1,default=str)
