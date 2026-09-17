import sys, itertools, time
from collections import Counter
sys.path.insert(0,"/Users/josecosta/bkrobust/src")
sys.path.insert(0,"/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag, is_valid_adjustment_set_dag)
from harness import random_dag

def build(rng, n=6, p=0.35, kmax=3):
    dag = random_dag(rng,n,p)
    if not dag.directed_edges: return None
    cp = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cp.undirected_edges)
    if not (2 <= len(und) <= 6): return None
    nodes=sorted(dag.nodes); cand=[(x,y) for x in nodes for y in nodes if x!=y]; rng.shuffle(cand)
    for x,y in cand:
        km = min(kmax, len(und))
        idx = rng.permutation(len(und))[:km]
        K=[]
        for t in idx:
            a,b=und[t]; K.append((a,b) if (a,b) in dag.directed_edges else (b,a))
        g0 = apply_orientations(cp,K)
        if g0 is None: continue
        z = optimal_adjustment_set_mpdag(g0,x,y)
        if z is None: continue
        exts = enumerate_dag_extensions(cp)
        if not exts: continue
        # C(d) = indices of claims REVERSED in d ; breaks(d)
        pats=[]
        for d in exts:
            C = frozenset(i for i,(a,b) in enumerate(K) if (b,a) in d.directed_edges)
            bad = not is_valid_adjustment_set_dag(d,x,y,z)
            pats.append((C,bad))
        breaking = sorted({C for C,bad in pats if bad}, key=len)
        if not breaking: continue          # UNREACHED
        # closure restatement of K (granularity attack): every und edge oriented in g0
        undset=set(und)
        Kclos=[e for e in g0.directed_edges if tuple(sorted(e)) in undset]
        patsC=[]
        for d in exts:
            C = frozenset(i for i,(a,b) in enumerate(Kclos) if (b,a) in d.directed_edges)
            bad = not is_valid_adjustment_set_dag(d,x,y,z)
            patsC.append((C,bad))
        breakC = sorted({C for C,bad in patsC if bad}, key=len)
        return dict(K=K,Kclos=Kclos,breaking=breaking,breakC=breakC,
                    nund=len(und), nexts=len(exts), x=x,y=y,z=z)
    return None

def r_w(breaking, w):
    return min(sum(w[i] for i in C) for C in breaking)

rng=np.random.default_rng(20260907)
probs=[]; tried=0; t0=time.time()
while len(probs)<400 and tried<40000 and time.time()-t0<300:
    tried+=1
    pr=build(rng)
    if pr: probs.append(pr)
print(f"problems={len(probs)}  tried={tried}  {time.time()-t0:.1f}s")
print(f"mean |K|={np.mean([len(p['K']) for p in probs]):.2f}  mean |Kclos|={np.mean([len(p['Kclos']) for p in probs]):.2f}  mean |[Chat]|={np.mean([p['nexts'] for p in probs]):.1f}")

# ---- r_K (uniform weights) ----
rK=np.array([min(len(C) for C in p['breaking']) for p in probs])
cK=Counter(rK.tolist())
print("\n=== r_K (all w=1) ===")
for k in sorted(cK): print(f"  r_K={k}: {cK[k]:4d}  {100*cK[k]/len(probs):5.1f}%")
# how often is the cheapest damaging pattern a SINGLE claim, and is it UNIQUE?
sing=[p for p,r in zip(probs,rK) if r==1]
nmin=[sum(1 for C in p['breaking'] if len(C)==1) for p in sing]
print(f"  among r_K=1 problems: mean #single-claim breaks = {np.mean(nmin):.2f}, exactly one = {100*np.mean(np.array(nmin)==1):.1f}%")

# ---- r_w under different elicitation regimes ----
def regime_draw(name, rng, m):
    if name=="U(0.6,0.95)": return rng.uniform(0.6,0.95,m)
    if name=="5-point scale": return rng.choice([0.6,0.7,0.8,0.9,0.95],m)
    if name=="3-point scale": return rng.choice([0.7,0.85,0.95],m)
    if name=="all 0.90 (flat)": return np.full(m,0.90)
    if name=="U(0.85,0.95) tight": return rng.uniform(0.85,0.95,m)
    raise KeyError(name)

print("\n=== r_w degeneracy vs elicitation regime (the claim: 101 distinct values, largest tie 30.6%) ===")
rng2=np.random.default_rng(1)
store={}
for name in ["U(0.6,0.95)","U(0.85,0.95) tight","5-point scale","3-point scale","all 0.90 (flat)"]:
    vals=[]
    for p in probs:
        w=-np.log(1-regime_draw(name,rng2,len(p['K'])))
        vals.append(round(r_w(p['breaking'],w),9))
    store[name]=np.array(vals)
    c=Counter(vals)
    top=max(c.values())
    print(f"  {name:22s} distinct={len(c):4d}/{len(vals)}  largest tie={100*top/len(vals):5.1f}%  mean={np.mean(vals):.3f}")

# ---- variance decomposition: is r_w measuring the graph or the dice? ----
print("\n=== where does r_w's spread come from? (200 independent redraws of p per problem) ===")
rng3=np.random.default_rng(2)
R=np.zeros((len(probs),200))
for i,p in enumerate(probs):
    for j in range(200):
        w=-np.log(1-rng3.uniform(0.6,0.95,len(p['K'])))
        R[i,j]=r_w(p['breaking'],w)
within=R.var(axis=1).mean(); between=R.mean(axis=1).var(); total=R.var()
print(f"  Var(r_w) total={total:.4f}  between-problem(graph)={between:.4f} ({100*between/total:.1f}%)  within-problem(weight draw)={within:.4f} ({100*within/total:.1f}%)")
print(f"  => {100*within/total:.1f}% of r_w's variation is the analyst's dice on a FIXED graph.")
# how much does r_w add over r_K?
from scipy.stats import spearmanr
v=store["U(0.6,0.95)"]
print(f"  Spearman(r_w, r_K) = {spearmanr(v,rK).statistic:.3f}")
# within the r_K==1 stratum, does r_w order anything about the graph?
m1=rK==1
print(f"  within r_K=1 stratum (n={m1.sum()}): corr(r_w, #single-claim-breaks) = {spearmanr(v[m1],np.array([sum(1 for C in p['breaking'] if len(C)==1) for p in probs])[m1]).statistic:.3f}")

# ---- granularity / gaming: restate K as its Meek closure, SAME G0, SAME Z ----
print("\n=== GAMING 1: restate K as its own Meek closure (same G0, same Z, same beliefs) ===")
rng4=np.random.default_rng(3); infl=[]; same=[]
for p in probs:
    if not p['breakC']: continue
    pe=rng4.uniform(0.6,0.95,len(p['Kclos']))
    # elicited claims keep their own p; cascaded ones get p drawn the same way
    a=r_w(p['breaking'],-np.log(1-pe[:len(p['K'])]))
    b=r_w(p['breakC'],-np.log(1-pe))
    infl.append(b-a); same.append((len(p['Kclos'])-len(p['K'])))
infl=np.array(infl); same=np.array(same)
print(f"  mean cascade size (|Kclos|-|K|) = {same.mean():.2f}, max {same.max()}")
print(f"  r_w after restatement - before: mean {infl.mean():+.3f}  inflated {100*np.mean(infl>1e-9):.1f}%  deflated {100*np.mean(infl<-1e-9):.1f}%")
sub=same>0
if sub.sum(): print(f"  restricted to a real cascade (n={sub.sum()}): mean {infl[sub].mean():+.3f}  inflated {100*np.mean(infl[sub]>1e-9):.1f}%")

# ---- adversarial permutation of the SAME confidence multiset ----
print("\n=== GAMING 2: analyst keeps the same confidence multiset, permutes which claim gets which ===")
rng5=np.random.default_rng(4); ratios=[]; spreads=[]
for p in probs:
    m=len(p['K'])
    if m<2: continue
    pe=np.sort(rng5.uniform(0.6,0.95,m))
    w0=-np.log(1-pe)
    vals=[r_w(p['breaking'],np.array(perm)) for perm in itertools.permutations(w0)]
    lo,hi=min(vals),max(vals)
    spreads.append(hi-lo); ratios.append(hi/lo if lo>0 else np.nan)
spreads=np.array(spreads); ratios=np.array(ratios)
print(f"  n={len(spreads)}  max/min r_w over permutations: mean ratio {np.nanmean(ratios):.2f}x  median {np.nanmedian(ratios):.2f}x  max {np.nanmax(ratios):.2f}x")
print(f"  absolute spread: mean {spreads.mean():.3f}  (compare to mean r_w {store['U(0.6,0.95)'].mean():.3f})")
print(f"  fraction where permuting alone changes r_w at all: {100*np.mean(spreads>1e-9):.1f}%")

# ---- the cap ----
print("\n=== the p ceiling ===")
for cap in [0.95,0.99,0.999]:
    w=-np.log(1-cap); print(f"  p cap {cap} -> w={w:.2f}; an analyst declaring the cap on every claim gets r_K*{w:.2f}")
