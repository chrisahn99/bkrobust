import sys, time
sys.path.insert(0,"/Users/josecosta/bkrobust/src")
import numpy as np
from rnull_real import *
t=time.time(); dag,W,var = load("ecoli70"); print("load", round(time.time()-t,2), flush=True)
t=time.time(); cpdag = dag_to_cpdag(dag); print("cpdag", round(time.time()-t,2), flush=True)
K = sorted(select_knowledge(dag,cpdag,1.0)); print("K", len(K), flush=True)
print("comps", sorted(len(c) for c in undirected_components(cpdag)), flush=True)
t=time.time(); prs = pairs_for(dag,cpdag,50,rng=np.random.default_rng(1))
print("pairs", round(time.time()-t,2), len(prs), prs[:4], flush=True)
sem=sem_of(dag,W,var); sig=sem.covariance(); idx={v:i for i,v in enumerate(dag.nodes)}
x,y = prs[0]; comp=component_of(cpdag,x)
Kloc=[e for e in K if component_of(cpdag,e[0])==comp]
print("x,y",x,y,"|C|",len(comp),"|Kloc|",len(Kloc), flush=True)
g0=apply_orientations(cpdag,K); Z=optimal_adjustment_set_mpdag(g0,x,y)
print("Z", sorted(Z) if Z else None, flush=True)
t=time.time(); c={}; th=theta(cpdag,x,y,sig,idx,c); print("theta(cpdag)", round(time.time()-t,2), "vals",len(th) if th else 0, flush=True)
t=time.time(); v=is_valid_adjustment_set_mpdag(cpdag,x,y,Z); print("valid?", round(time.time()-t,2), v, flush=True)
