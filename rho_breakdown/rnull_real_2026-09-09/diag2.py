import sys; sys.path.insert(0,"/Users/josecosta/bkrobust/src")
import numpy as np
from rnull_real import *
dag,W,var = load("ecoli70"); sem=sem_of(dag,W,var); cpdag=dag_to_cpdag(dag)
sig=sem.covariance(); idx={v:i for i,v in enumerate(dag.nodes)}
K=sorted(select_knowledge(dag,cpdag,1.0)); g0=apply_orientations(cpdag,K)
rng=np.random.default_rng(20260909)
prs=pairs_for(dag,cpdag,60,rng=rng)
X,ix=sample(dag,W,var,20000,rng); sh=np.cov(X,rowvar=False)
print(f"{'x':>7} {'y':>7} {'casco j=1 (oraculo)':>28} {'0?':>3} {'casco j=1 (n=20000)':>28} {'0?':>3}")
diff=0; tot=0
for x,y in prs:
    Z=optimal_adjustment_set_mpdag(g0,x,y)
    if Z is None: continue
    comp=component_of(cpdag,x); Kl=[e for e in K if component_of(cpdag,e[0])==comp]
    if not Kl: continue
    a=radii(cpdag,Kl,x,y,Z,sig,idx,1); b=radii(cpdag,Kl,x,y,Z,sh,ix,1)
    if 1 not in a or 1 not in b: continue
    tot+=1
    za,zb=a[1]["has_zero"],b[1]["has_zero"]
    if za!=zb:
        diff+=1
        print(f"{x:>7} {y:>7} [{a[1]['lo']:+.4f},{a[1]['hi']:+.4f}] {str(za):>5} "
              f"  [{b[1]['lo']:+.4f},{b[1]['hi']:+.4f}] {str(zb):>5}   <-- DIVERGE")
print(f"\n{diff} divergencias em {tot} instancias  ({100*diff/max(1,tot):.1f}%)")
