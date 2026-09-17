import sys, time
sys.path.insert(0,"/Users/josecosta/bkrobust/src"); sys.path.insert(0,"/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag, is_valid_adjustment_set_dag,
                                    is_valid_adjustment_set_mpdag)
from harness import random_dag
def build(rng):
    dag=random_dag(rng,6,0.35)
    if not dag.directed_edges: return None
    cp=dag_to_cpdag(dag); und=sorted(tuple(sorted(e)) for e in cp.undirected_edges)
    if not (2<=len(und)<=6): return None
    nodes=sorted(dag.nodes); cand=[(x,y) for x in nodes for y in nodes if x!=y]; rng.shuffle(cand)
    for x,y in cand:
        idx=rng.permutation(len(und))[:min(3,len(und))]
        K=[]
        for t in idx:
            a,b=und[t]; K.append((a,b) if (a,b) in dag.directed_edges else (b,a))
        g0=apply_orientations(cp,K)
        if g0 is None: continue
        z=optimal_adjustment_set_mpdag(g0,x,y)
        if z is None or not is_valid_adjustment_set_mpdag(g0,x,y,z): continue
        exts=enumerate_dag_extensions(cp)
        br={frozenset(i for i,(a,b) in enumerate(K) if (b,a) in d.directed_edges)
            for d in exts if not is_valid_adjustment_set_dag(d,x,y,z)}
        if not br: return None
        return dict(K=K,br=br)
    return None
rng=np.random.default_rng(20260907); P=[]; t0=time.time()
while len(P)<400 and time.time()-t0<300:
    p=build(rng)
    if p: P.append(p)
rK=np.array([min(len(C) for C in p['br']) for p in P])
s1=[p for p,r in zip(P,rK) if r==1]
n1=np.array([sum(1 for C in p['br'] if len(C)==1) for p in s1])
print(f"clean problems={len(P)}; r_K=1 stratum = {100*np.mean(rK==1):.1f}%")
print(f"within it: exactly ONE single-claim breaking pattern in {100*np.mean(n1==1):.1f}% of problems")
print(f"=> for those, r_w = -log(1-p) of a single elicited claim: one human number, zero graph content.")
print(f"mean |K| = {np.mean([len(p['K']) for p in P]):.2f}")
