"""Two checks: (1) is exp(-r_w) the probability the plain reading says it is?
   (2) does the log-odds cost ever reorder patterns against plain counting?"""
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
    if not (2<=len(und)<=5): return None
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
        pats=[(frozenset(i for i,(a,b) in enumerate(K) if (b,a) in d.directed_edges),
               not is_valid_adjustment_set_dag(d,x,y,z)) for d in exts]
        if not any(b for _,b in pats): continue
        return dict(K=K,pats=pats)
    return None

rng=np.random.default_rng(777); P=[]; t0=time.time()
while len(P)<400 and time.time()-t0<300:
    p=build(rng)
    if p: P.append(p)
print(f"problems={len(P)}")

rr=np.random.default_rng(5); ratios=[]; claimed=[]; actual=[]; reorder=0; tot=0
for p in P:
    m=len(p['K']); conf=rr.uniform(0.6,0.95,m); w=-np.log(1-conf)
    breaking={C for C,b in p['pats'] if b}
    rw=min(sum(w[i] for i in C) for C in breaking)
    # the plain reading's number
    claimed_prob=np.exp(-rw)
    # the actual prior probability that SOME feasible error pattern breaks Z,
    # under the same independence prior, conditioned on the feasible set
    def patprob(C): return np.prod([1-conf[i] for i in C])*np.prod([conf[i] for i in range(m) if i not in C])
    seen={}
    for C,b in p['pats']: seen[C]=b
    Zn=sum(patprob(C) for C in seen)
    actual_prob=sum(patprob(C) for C,b in seen.items() if b)/Zn
    claimed.append(claimed_prob); actual.append(actual_prob); ratios.append(actual_prob/claimed_prob)
    # does log-odds cost ever order two patterns against their sizes?
    Cs=list(seen)
    for i in range(len(Cs)):
        for j in range(len(Cs)):
            if len(Cs[i])<len(Cs[j]):
                tot+=1
                if sum(w[k] for k in Cs[i])>sum(w[k] for k in Cs[j]): reorder+=1
claimed=np.array(claimed); actual=np.array(actual); ratios=np.array(ratios)
print(f"\nexp(-r_w), the number the plain reading calls 'the prior probability of breakage':")
print(f"   mean {claimed.mean():.4f}   median {np.median(claimed):.4f}")
print(f"actual P(some feasible pattern breaks Z) under the SAME independence prior:")
print(f"   mean {actual.mean():.4f}   median {np.median(actual):.4f}")
print(f"ratio actual/claimed: median {np.median(ratios):.2f}x   90th pct {np.percentile(ratios,90):.2f}x   max {ratios.max():.2f}x")
print(f"   the reported number UNDERSTATES the true breakage probability in {100*np.mean(ratios>1):.1f}% of problems")
print(f"\nlog-odds cost reorders a smaller pattern above a larger one in {100*reorder/max(tot,1):.2f}% of comparable pairs "
      f"(n={tot}) -- p in (0.6,0.95) gives w in ({-np.log(1-0.6):.2f},{-np.log(1-0.95):.2f}), ratio {np.log(1-0.6)/np.log(1-0.95):.2f}x")
