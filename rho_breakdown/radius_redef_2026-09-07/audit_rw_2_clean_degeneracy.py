"""Clean audit of r_w. Filters out the degenerate-query stratum (Z invalid at cost 0)."""
import sys, itertools, time
from collections import Counter
sys.path.insert(0,"/Users/josecosta/bkrobust/src")
sys.path.insert(0,"/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from scipy.stats import spearmanr
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag, optimal_adjustment_set_dag,
                                    is_valid_adjustment_set_dag, is_valid_adjustment_set_mpdag)
from harness import random_dag

def build(rng,n=6,p=0.35,kmax=3):
    dag=random_dag(rng,n,p)
    if not dag.directed_edges: return None
    cp=dag_to_cpdag(dag)
    und=sorted(tuple(sorted(e)) for e in cp.undirected_edges)
    if not (2<=len(und)<=6): return None
    nodes=sorted(dag.nodes); cand=[(x,y) for x in nodes for y in nodes if x!=y]; rng.shuffle(cand)
    for x,y in cand:
        km=min(kmax,len(und)); idx=rng.permutation(len(und))[:km]
        K=[]
        for t in idx:
            a,b=und[t]; K.append((a,b) if (a,b) in dag.directed_edges else (b,a))
        g0=apply_orientations(cp,K)
        if g0 is None: continue
        z=optimal_adjustment_set_mpdag(g0,x,y)
        if z is None: continue
        # THE FILTER the existing script lacks: Z must be valid throughout [G0],
        # i.e. no zero-cost error pattern may break it.
        if not is_valid_adjustment_set_mpdag(g0,x,y,z): continue
        exts=enumerate_dag_extensions(cp)
        if not exts: continue
        pats=[]
        for d in exts:
            C=frozenset(i for i,(a,b) in enumerate(K) if (b,a) in d.directed_edges)
            pats.append((C, not is_valid_adjustment_set_dag(d,x,y,z), d))
        breaking=sorted({C for C,bad,_ in pats if bad}, key=len)
        if not breaking: return dict(status="UNREACHED",K=K)
        return dict(status="ok",dag=dag,cp=cp,x=x,y=y,z=z,K=K,exts=exts,pats=pats,
                    breaking=breaking,nund=len(und))
    return None

def r_w(breaking,w): return min(sum(w[i] for i in C) for C in breaking)

rng=np.random.default_rng(20260907); probs=[]; unreached=0; tried=0; t0=time.time()
while len(probs)<400 and tried<60000 and time.time()-t0<420:
    tried+=1; pr=build(rng)
    if pr is None: continue
    if pr["status"]=="UNREACHED": unreached+=1; continue
    probs.append(pr)
print(f"problems={len(probs)}  UNREACHED(discarded)={unreached}  tried={tried}  {time.time()-t0:.1f}s")
print(f"mean |K|={np.mean([len(p['K']) for p in probs]):.2f}  mean |[Chat]|={np.mean([len(p['exts']) for p in probs]):.1f}")

rK=np.array([min(len(C) for C in p['breaking']) for p in probs])
c=Counter(rK.tolist())
print("\n=== r_K (uniform weights) on the CLEAN sample ===")
for k in sorted(c): print(f"  r_K={k}: {c[k]:4d}  {100*c[k]/len(probs):5.1f}%")
print(f"  ties on the modal value: {100*max(c.values())/len(probs):.1f}%")

print("\n=== r_w distinctness vs elicitation regime, AND the cosmetic control ===")
def draw(name,rng,m):
    return {"U(0.6,0.95)":lambda: rng.uniform(0.6,0.95,m),
            "U(0.85,0.95)":lambda: rng.uniform(0.85,0.95,m),
            "5-pt scale":lambda: rng.choice([0.6,0.7,0.8,0.9,0.95],m),
            "3-pt scale":lambda: rng.choice([0.7,0.85,0.95],m),
            "all at cap 0.99":lambda: np.full(m,0.99),
            "all 0.90 flat":lambda: np.full(m,0.90)}[name]()
rng2=np.random.default_rng(1); store={}
for name in ["U(0.6,0.95)","U(0.85,0.95)","5-pt scale","3-pt scale","all 0.90 flat","all at cap 0.99"]:
    vals=np.array([round(r_w(p['breaking'],-np.log(1-draw(name,rng2,len(p['K'])))),9) for p in probs])
    store[name]=vals; cc=Counter(vals.tolist())
    print(f"  r_w {name:16s} distinct={len(cc):4d}/{len(vals)}  largest tie={100*max(cc.values())/len(vals):5.1f}%  rho(r_w,r_K)={spearmanr(vals,rK).statistic:.3f}")
rng9=np.random.default_rng(99)
ctrl=rK+rng9.uniform(0,1,len(rK))
cc=Counter(np.round(ctrl,9).tolist())
print(f"  CONTROL r_K+U(0,1)      distinct={len(cc):4d}/{len(ctrl)}  largest tie={100*max(cc.values())/len(ctrl):5.1f}%  rho(ctrl,r_K)={spearmanr(ctrl,rK).statistic:.3f}")

print("\n=== variance decomposition of r_w (300 redraws of p on a FIXED graph) ===")
rng3=np.random.default_rng(2); R=np.zeros((len(probs),300))
for i,p in enumerate(probs):
    for j in range(300):
        R[i,j]=r_w(p['breaking'],-np.log(1-rng3.uniform(0.6,0.95,len(p['K']))))
tot=R.var(); btw=R.mean(axis=1).var(); wth=R.var(axis=1).mean()
print(f"  total={tot:.4f}  between-problem(graph)={btw:.4f} ({100*btw/tot:.1f}%)  within-problem(dice)={wth:.4f} ({100*wth/tot:.1f}%)")

print("\n=== permutation attack: same confidence multiset, reassigned across claims ===")
rng5=np.random.default_rng(4); sp=[]; ra=[]
for p in probs:
    m=len(p['K'])
    if m<2: continue
    w0=-np.log(1-np.sort(rng5.uniform(0.6,0.95,m)))
    v=[r_w(p['breaking'],np.array(pm)) for pm in itertools.permutations(w0)]
    sp.append(max(v)-min(v)); ra.append(max(v)/min(v) if min(v)>0 else np.nan)
sp=np.array(sp); ra=np.array(ra)
print(f"  n={len(sp)}  changes r_w at all in {100*np.mean(sp>1e-9):.1f}% of problems; mean max/min ratio {np.nanmean(ra):.2f}x  max {np.nanmax(ra):.2f}x")
print(f"  mean spread {sp.mean():.3f} vs mean r_w {store['U(0.6,0.95)'].mean():.3f}")

print("\n=== monotone inflation: r_w is nondecreasing in every self-reported p ===")
for lo,hi in [(0.6,0.95),(0.8,0.98),(0.9,0.99)]:
    rr=np.random.default_rng(7)
    v=np.array([r_w(p['breaking'],-np.log(1-rr.uniform(lo,hi,len(p['K'])))) for p in probs])
    print(f"  p ~ U({lo},{hi}): mean r_w = {v.mean():.3f}   (same graphs, same damage, {v.mean()/store['U(0.6,0.95)'].mean():.2f}x the certificate)")
