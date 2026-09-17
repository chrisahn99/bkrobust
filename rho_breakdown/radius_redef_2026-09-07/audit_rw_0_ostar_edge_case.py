import sys
sys.path.insert(0,"/Users/josecosta/bkrobust/src")
sys.path.insert(0,"/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions, is_consistent_extension
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag, optimal_adjustment_set_dag,
                                    is_valid_adjustment_set_dag, is_valid_adjustment_set_mpdag)
from harness import random_dag

rng=np.random.default_rng(20260907)
found=0
for trial in range(3000):
    dag=random_dag(rng,6,0.35)
    if not dag.directed_edges: continue
    cp=dag_to_cpdag(dag)
    und=sorted(tuple(sorted(e)) for e in cp.undirected_edges)
    if not (2<=len(und)<=6): continue
    nodes=sorted(dag.nodes); cand=[(x,y) for x in nodes for y in nodes if x!=y]; rng.shuffle(cand)
    for x,y in cand:
        km=min(3,len(und)); idx=rng.permutation(len(und))[:km]
        K=[]
        for t in idx:
            a,b=und[t]; K.append((a,b) if (a,b) in dag.directed_edges else (b,a))
        g0=apply_orientations(cp,K)
        if g0 is None: continue
        z=optimal_adjustment_set_mpdag(g0,x,y)
        if z is None: continue
        exts=enumerate_dag_extensions(cp)
        zero=[d for d in exts if all((b,a) not in d.directed_edges for a,b in K)
              and not is_valid_adjustment_set_dag(d,x,y,z)]
        if zero:
            d=zero[0]
            print("X,Y =",x,y,"  Z =",sorted(z))
            print("K =",K)
            print("g0 directed:",sorted(g0.directed_edges),"  g0 undirected:",sorted(g0.undirected_edges))
            print("witness d directed:",sorted(d.directed_edges))
            print("is_consistent_extension(d,g0) =",is_consistent_extension(d,g0))
            print("d in enumerate_dag_extensions(g0)?",any(dd.directed_edges==d.directed_edges for dd in enumerate_dag_extensions(g0)))
            print("O*(d) =",sorted(optimal_adjustment_set_dag(d,x,y)) if optimal_adjustment_set_dag(d,x,y) is not None else None)
            print("is_valid_adjustment_set_mpdag(g0,x,y,z) =",is_valid_adjustment_set_mpdag(g0,x,y,z))
            print("true dag directed:",sorted(dag.directed_edges))
            print("valid in true dag?",is_valid_adjustment_set_dag(dag,x,y,z))
            found+=1
        break
    if found>=3: break
