"""Is a false statement ever CAUGHT by Meek's consistency check while the edge is
still UNDIRECTED in the analyst's current MPDAG?"""
import sys
sys.path.insert(0,"/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/falsegain/code")
from fglib import *
import numpy as np
from collections import defaultdict
from graphs import random_dag, is_undirected, apply_background_knowledge, MeekFail
rng=np.random.default_rng(4242); acc=defaultdict(int)
for it in range(600):
    p=int(rng.choice([5,6,7,8])); deg=float(rng.choice([1.5,2.0,2.5,3.0]))
    D=random_dag(p,deg,rng); C=dag_to_cpdag(D); U=undirected_edges(C)
    if not U: continue
    Ktrue=[(u,v) if D[u,v]==1 else (v,u) for (u,v) in U]
    for rho in (1,2,3):
        for rep in range(4):
            nt=int(rng.integers(0,max(len(U)-rho,0)+1)); perm=rng.permutation(len(U))
            bi=list(perm[:nt]); fi=list(perm[nt:nt+rho])
            if len(fi)<rho: continue
            Kb=[Ktrue[i] for i in bi]; Kf=[(Ktrue[i][1],Ktrue[i][0]) for i in fi]
            try: G0=apply_background_knowledge(C,Kb)
            except MeekFail: continue
            fresh=all(is_undirected(G0,a,b) for (a,b) in Kf)
            try:
                apply_background_knowledge(C,Kb+Kf); ok=True
            except MeekFail: ok=False
            acc[(rho,"fresh" if fresh else "stale","accepted" if ok else "CAUGHT")]+=1
for k,v in sorted(acc.items()): print("|".join(map(str,k)),v)
