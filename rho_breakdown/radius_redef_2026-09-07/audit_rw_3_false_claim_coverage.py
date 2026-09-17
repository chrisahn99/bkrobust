"""One elicited claim is FALSE. Does the exp(-cost) prior sub-level set still cover the truth,
and at what width, against a uniform-count budget of the same size?"""
import sys, time
from collections import Counter
sys.path.insert(0,"/Users/josecosta/bkrobust/src")
sys.path.insert(0,"/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag, optimal_adjustment_set_dag,
                                    is_valid_adjustment_set_dag, random_sem, adjusted_estimand)
from harness import random_dag

def build(rng,n=6,p=0.35,kmax=3):
    dag=random_dag(rng,n,p)
    if not dag.directed_edges: return None
    cp=dag_to_cpdag(dag)
    und=sorted(tuple(sorted(e)) for e in cp.undirected_edges)
    if not (2<=len(und)<=5): return None
    nodes=sorted(dag.nodes); cand=[(x,y) for x in nodes for y in nodes if x!=y]; rng.shuffle(cand)
    for x,y in cand:
        km=min(kmax,len(und)); idx=rng.permutation(len(und))[:km]
        Ktrue=[]
        for t in idx:
            a,b=und[t]; Ktrue.append((a,b) if (a,b) in dag.directed_edges else (b,a))
        jbad=int(rng.integers(len(Ktrue)))              # exactly one claim is REVERSED
        K=[(b,a) if i==jbad else (a,b) for i,(a,b) in enumerate(Ktrue)]
        g0=apply_orientations(cp,K)
        if g0 is None: continue
        z=optimal_adjustment_set_mpdag(g0,x,y)
        if z is None: continue
        exts=enumerate_dag_extensions(cp)
        if len(exts)<3: continue
        dtrue=next((d for d in exts if d.directed_edges==dag.directed_edges),None)
        if dtrue is None: continue
        sem=random_sem(dag,rng); tau=sem.true_total_effect(x,y)
        rows=[]
        for d in exts:
            C=frozenset(i for i,(a,b) in enumerate(K) if (b,a) in d.directed_edges)
            zd=optimal_adjustment_set_dag(d,x,y)
            rows.append((C,adjusted_estimand(sem,x,y,zd),d.directed_edges==dag.directed_edges))
        # analyst's own point estimate under the false knowledge
        est0=adjusted_estimand(sem,x,y,z)
        return dict(K=K,jbad=jbad,rows=rows,tau=tau,est0=est0,
                    broken=not is_valid_adjustment_set_dag(dag,x,y,z))
    return None

def sublevel(costs, prior_mass=0.95):
    """smallest cost sub-level set carrying >=95% of the exp(-cost) prior"""
    c=np.array(costs); wgt=np.exp(-c); wgt=wgt/wgt.sum()
    o=np.argsort(c); cum=0.0
    for i in o:
        cum+=wgt[i]
        if cum>=prior_mass: return c[i]
    return c.max()

rng=np.random.default_rng(4242); P=[]; t0=time.time()
while len(P)<400 and time.time()-t0<400:
    pr=build(rng)
    if pr: P.append(pr)
print(f"problems={len(P)}  ({time.time()-t0:.0f}s)   Z actually invalid at the truth in "
      f"{100*np.mean([p['broken'] for p in P]):.1f}% of them")

regimes={"calibrated (false claim least confident)":"cal",
         "random p ~ U(0.6,0.95)":"rand",
         "anti-calibrated (false claim MOST confident)":"anti"}
rr=np.random.default_rng(11)
print(f"\n{'regime':46s} {'cover':>7} {'width':>8} {'|S|/|C|':>8}   || uniform budget, matched size: cover / width")
for label,mode in regimes.items():
    cov=[];wid=[];frac=[];covU=[];widU=[]
    for p in P:
        m=len(p['K']); pe=np.sort(rr.uniform(0.6,0.95,m))[::-1]  # descending
        conf=np.empty(m); order=list(rr.permutation(m))
        if mode=="cal":   order=[i for i in order if i!=p['jbad']]+[p['jbad']]   # false claim last = lowest p
        elif mode=="anti":order=[p['jbad']]+[i for i in order if i!=p['jbad']]   # false claim first = highest p
        for rank,i in enumerate(order): conf[i]=pe[rank]
        w=-np.log(1-conf)
        costs=[sum(w[i] for i in C) for C,_,_ in p['rows']]
        t=sublevel(costs)
        S=[k for k,c in enumerate(costs) if c<=t+1e-12]
        vals=[p['rows'][k][1] for k in S]
        cov.append(min(vals)-1e-9<=p['tau']<=max(vals)+1e-9); wid.append(max(vals)-min(vals))
        frac.append(len(S)/len(costs))
        # uniform-count budget truncated to the SAME number of DAGs (ties broken by count then index)
        cnt=[len(C) for C,_,_ in p['rows']]
        ordU=sorted(range(len(cnt)),key=lambda k:(cnt[k],k))[:len(S)]
        vU=[p['rows'][k][1] for k in ordU]
        covU.append(min(vU)-1e-9<=p['tau']<=max(vU)+1e-9); widU.append(max(vU)-min(vU))
    print(f"{label:46s} {np.mean(cov):7.3f} {np.mean(wid):8.4f} {np.mean(frac):8.3f}   ||  {np.mean(covU):.3f} / {np.mean(widU):.4f}")

allv=[[r[1] for r in p['rows']] for p in P]
print(f"\nblanket hedge over all of [Chat]: cover {np.mean([min(v)<=p['tau']<=max(v) for v,p in zip(allv,P)]):.3f}  width {np.mean([max(v)-min(v) for v in allv]):.4f}")
print(f"analyst point estimate (no hedge): |bias| mean {np.mean([abs(p['est0']-p['tau']) for p in P]):.4f}, "
      f"median {np.median([abs(p['est0']-p['tau']) for p in P]):.4f}")
