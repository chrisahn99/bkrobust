import sys
from collections import defaultdict
sys.path.insert(0,"/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/falsegain/code")
from fglib import *
import fglib, numpy as np, importlib.util
from graphs import random_dag, is_undirected
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
    NA=x1_ops.nonadjacent_pairs(C)
    for rho in (2,3,4):
        if len(NA)<rho: continue
        for rep in range(6):
            idx=rng.choice(len(NA),size=rho,replace=False); K=[]
            for t in idx:
                a,b=NA[int(t)]; K.append((int(a),int(b)) if rng.random()<0.5 else (int(b),int(a)))
            H,info=x1_ops.bk_assert(G0,K)
            cyc=bool(has_directed_cycle(H)); con=bool(info["conflict"]); ext=x1_ops.pdag_extendable(H)
            g = "conflict" if con else ("cycle" if cyc else ("dortarsi_only" if not ext else "ACCEPT"))
            for x in range(p):
                for y in range(p):
                    if x==y: continue
                    s0,O0=report_state(G0,x,y)
                    if s0==REFUSE or not sound(D,x,y,s0,O0): continue
                    s1,O1=report_state(H,x,y); n1=sound(D,x,y,s1,O1)
                    acc[(g,rho,s0,s1,str(n1))]+=1
                    if g=="dortarsi_only" and n1 is False and len(wit)<3:
                        wit.append((p,D.tolist(),G0.tolist(),H.tolist(),x,y,K,s0,s1))
    if len(fglib._valid_cache)>150000: fglib._valid_cache.clear(); fglib._cn_cache.clear()
for k,v in sorted(acc.items(),key=lambda z:str(z[0])): print("|".join(map(str,k)),v)
print("dortarsi_only corruptions:", len(wit))
