"""Lemma 2 (Lift Lemma) as a standalone claim, sampled at p=5,6:
   D in [G0], (alpha) dir(H)|_S subseteq dir(D), (beta) no b~>a in D   ==>  D+(a->b) in [H]."""
import sys; sys.path.insert(0,'.')
from common import *
import numpy as np
p=int(sys.argv[1]); n=int(sys.argv[2]); rng=np.random.default_rng(int(sys.argv[3]))
tested=0; bad=0
for _ in range(n):
    deg=float(rng.choice([1.5,2.0,2.5,3.0,3.5,4.0]))
    D0=random_dag(p,deg,rng); C=dag_to_cpdag(D0); ref=v_structures(C)
    G0=C.copy(); t=int(rng.integers(0,len(undirected_edges(C))+1))
    for _ in range(t):
        U=undirected_edges(G0)
        if not U: break
        uu,vv=U[int(rng.integers(len(U)))]
        aa,bb=(uu,vv) if rng.random()<0.5 else (vv,uu)
        Hh=G0.copy(); Hh[bb,aa]=0; Hh=meek_closure(Hh)
        if has_directed_cycle(Hh) or v_structures(Hh)!=ref: continue
        G0=Hh
    if len(undirected_edges(G0))>12: continue
    ext=consistent_dag_extensions(G0,ref_vstructs=ref)
    if not ext: continue
    for (u,v) in nonadjacent_pairs(G0):
        for (a,b) in ((u,v),(v,u)):
            H,info=bk_assert(G0,[(a,b)])
            if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
            dirHS=[(s,t) for (s,t) in directed_edges(H) if not (s==a and t==b)]
            for D in ext:
                if not all(is_directed(D,s,t) for (s,t) in dirHS): continue      # alpha
                Dp=D.copy(); Dp[a,b]=1
                if has_directed_cycle(Dp): continue                              # beta
                tested+=1
                if v_structures(Dp)!=v_structures(H) or not all(is_directed(Dp,s,t) for (s,t) in directed_edges(H)):
                    bad+=1
                    if bad==1: print("WITNESS",G0.tolist(),D.tolist(),a,b)
print("LIFT LEMMA p=%d: tested=%d violations=%d"%(p,tested,bad))
