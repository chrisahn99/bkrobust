"""Larger elicitations: does higher-order (size>=2) unsafe structure grow with |K*|?"""
from __future__ import annotations
import os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import random_dag, true_orientations, canonicalise, Problem, minimal_unsafe, min_hitting_sets, powerset
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag, is_valid_adjustment_set_mpdag

kmin,kmax=int(sys.argv[1]),int(sys.argv[2]); target=int(sys.argv[3]); seed=int(sys.argv[4])
rng=np.random.default_rng(seed); rows=[]; t0=time.time(); tries=0
while len(rows)<target and time.time()-t0<600:
    tries+=1
    n=int(rng.integers(9,13)); p=float(rng.uniform(0.14,0.26))
    dag=random_dag(rng,n,p)
    try: cp=dag_to_cpdag(dag)
    except ValueError: continue
    A=true_orientations(dag,cp)
    if len(A)<kmin: continue
    size=int(rng.integers(kmin,len(A)+1))
    K=[A[i] for i in sorted(rng.choice(len(A),size=size,replace=False))]
    G0=apply_orientations(cp,K)
    if G0 is None: continue
    nodes=list(dag.nodes)
    pairs=[(a,b) for a in nodes for b in nodes if a!=b and b in dag.descendants(a)]
    if not pairs: continue
    x,y=pairs[int(rng.integers(len(pairs)))]
    Z=optimal_adjustment_set_mpdag(G0,x,y)
    if Z is None or not is_valid_adjustment_set_mpdag(G0,x,y,Z): continue
    Kc=canonicalise(cp,K,G0)
    if not (kmin<=len(Kc)<=kmax): continue
    m=len(Kc); pr=Problem(cp,Kc,x,y,Z)
    mins=minimal_unsafe(pr,m); c,hs=min_hitting_sets(mins,m)
    phi1=sum(1 for i in range(m) if not pr.safe((i,)))
    nW=sum(1 for J in powerset(m) if pr.world(J) is not None)
    core=set(hs[0])
    nM=sum(1 for J in powerset(m) if pr.world(J) is not None and ((not J) or set(J)&core))
    rows.append(dict(m=m,c=c,phi1=phi1,mm=max((len(u) for u in mins),default=0),
                     checks=pr.n_checks,nhs=len(hs),nW=nW,nM=nM))
g=lambda k: np.array([r[k] for r in rows])
from collections import Counter
print(f"|K*| in [{kmin},{kmax}]  n={len(rows)} tries={tries} secs={time.time()-t0:.1f}")
if rows:
    c=g('c'); m=g('m')
    print(f"  |K*| mean={m.mean():.2f}; c distinct={len(set(c.tolist()))} dist={sorted(Counter(c.tolist()).items())} mean={c.mean():.3f}")
    print(f"  P(c=0)={np.mean(c==0):.4f} P(c=|K*|)={np.mean(c==m):.4f} modal mass={max(Counter(c.tolist()).values())/len(rows):.4f}")
    print(f"  c==phi1 on {np.mean(c==g('phi1')):.4f}; P(minimal unsafe set of size>=2)={np.mean(g('mm')>=2):.4f}")
    print(f"  core unique on {np.mean(g('nhs')==1):.4f}")
    print(f"  validity checks/problem mean={g('checks').mean():.2f} max={g('checks').max()} vs 2^m={np.mean(2.0**m):.1f} 3^m={np.mean(3.0**m):.1f}")
    print(f"  hedge support |W|={g('nW').mean():.2f} |M|={g('nM').mean():.2f} ratio={np.mean(g('nW')/g('nM')):.3f}")
