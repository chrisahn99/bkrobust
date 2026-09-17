#!/usr/bin/env python3
"""(a) frontier shape on LONG ladders (m*>=3); (b) the query-local row CRH^loc_j;
(c) corrected RHIG_1/blanket ratio with the (x,y)-keyed tau cache."""
import sys, itertools, time
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from bkrobust.demo.graph import canon
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import optimal_adjustment_set_dag, random_sem, adjusted_estimand
from harness import random_dag

def tau_range(g, sem, x, y, cache):
    key = (g.edge_string(), x, y)
    if key in cache: return cache[key]
    vals = []
    for D in enumerate_dag_extensions(g):
        try: vals.append(adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(D, x, y)))
        except Exception: pass
    out = (min(vals), max(vals)) if vals else None
    cache[key] = out; return out

def oriented_beyond(cpdag, g0):
    und = {canon(a,b) for a,b in cpdag.undirected_edges}
    return sorted((a,b) for (a,b) in g0.directed_edges if canon(a,b) in und)

def minimal_generating_sets(cpdag, g0, Or):
    gens=[]
    for r in range(0,len(Or)+1):
        for S in itertools.combinations(Or,r):
            if any(set(G)<=set(S) for G in gens): continue
            h=apply_orientations(cpdag,S)
            if h is not None and h==g0: gens.append(S)
    return gens

def hull(cpdag, subsets, sem, x, y, cache):
    lo,hi,ok=np.inf,-np.inf,False
    for T in subsets:
        g=apply_orientations(cpdag,T)
        if g is None: continue
        r=tau_range(g,sem,x,y,cache)
        if r is None: continue
        ok=True; lo=min(lo,r[0]); hi=max(hi,r[1])
    return (lo,hi) if ok else None

def crh_sets(gens,j):
    out=[]
    for S in gens: out += list(itertools.combinations(S, max(len(S)-j,0)))
    return out

def local_nodes(cpdag,x,y):
    L={x,y}
    for a,b in list(cpdag.directed_edges)+[tuple(e) for e in cpdag.undirected_edges]:
        if a in (x,y): L.add(b)
        if b in (x,y): L.add(a)
    return L

def crhloc_sets(gens,j,L):
    out=[]
    for S in gens:
        loc=[e for e in S if e[0] in L or e[1] in L]
        far=[e for e in S if not (e[0] in L or e[1] in L)]
        for d in range(0,min(j,len(loc))+1):
            for drop in itertools.combinations(loc,d):
                out.append(tuple(far)+tuple(e for e in loc if e not in drop))
    return out

def sweep(seed, want, mstar_min, n_nodes=7, p=0.32, kmax=7):
    rng=np.random.default_rng(seed); rows=[]; t0=time.time(); tried=0
    while len(rows)<want and time.time()-t0<420 and tried<400000:
        tried+=1
        dag=random_dag(rng,n_nodes,p)
        if not dag.directed_edges: continue
        cpdag=dag_to_cpdag(dag)
        und=sorted(canon(*e) for e in cpdag.undirected_edges)
        if not (2<=len(und)<=kmax): continue
        nodes=sorted(dag.nodes); cand=[(a,b) for a in nodes for b in nodes if a!=b]
        rng.shuffle(cand); sem=random_sem(dag,rng); cache={}
        for x,y in cand[:6]:
            blank=tau_range(cpdag,sem,x,y,cache)
            if blank is None or blank[1]-blank[0]<1e-8: continue
            km=int(rng.integers(1,min(5,len(und))+1))
            idx=rng.permutation(len(und))[:km]
            K=[]
            for t in idx:
                a,b=und[t]; K.append((a,b) if (a,b) in dag.directed_edges else (b,a))
            g0=apply_orientations(cpdag,K)
            if g0 is None: continue
            Or=oriented_beyond(cpdag,g0); gens=minimal_generating_sets(cpdag,g0,Or)
            if not gens: continue
            mstar=min(len(S) for S in gens)
            if mstar<mstar_min: continue
            try: tau_true=adjusted_estimand(sem,x,y,optimal_adjustment_set_dag(dag,x,y))
            except Exception: continue
            bw=blank[1]-blank[0]
            W=[]; 
            for j in range(0,mstar+1):
                r=hull(cpdag,crh_sets(gens,j),sem,x,y,cache)
                W.append(None if r is None else r[1]-r[0])
            L=local_nodes(cpdag,x,y)
            rl=hull(cpdag,crhloc_sets(gens,1,L),sem,x,y,cache)
            rl2=hull(cpdag,crhloc_sets(gens,99,L),sem,x,y,cache)
            r1=hull(cpdag,crh_sets(gens,1),sem,x,y,cache)
            rows.append(dict(mstar=mstar,k=len(und),blanket=bw,W=W,
                w_crh1=None if r1 is None else r1[1]-r1[0],
                cov_crh1=None if r1 is None else r1[0]-1e-9<=tau_true<=r1[1]+1e-9,
                w_loc1=None if rl is None else rl[1]-rl[0],
                cov_loc1=None if rl is None else rl[0]-1e-9<=tau_true<=rl[1]+1e-9,
                w_locall=None if rl2 is None else rl2[1]-rl2[0],
                cov_locall=None if rl2 is None else rl2[0]-1e-9<=tau_true<=rl2[1]+1e-9,
                nloc=len([e for e in gens[0] if e[0] in L or e[1] in L])))
            break
    return rows,tried,time.time()-t0

for mm in (3,1):
    rows,tried,el=sweep(23+mm, 80 if mm==3 else 150, mm)
    n=len(rows); B=np.array([r['blanket'] for r in rows])
    print(f"\n===== stratum m* >= {mm} : n={n} (tried {tried}, {el:.0f}s) =====")
    print("  m*:", dict(sorted(Counter(r['mstar'] for r in rows).items())))
    for j in range(0,5):
        v=[(r['W'][j]/r['blanket']) for r in rows if j<len(r['W'])]
        if v: print(f"   j={j}: mean width/blanket = {np.mean(v):.3f}  (n={len(v)})  P(=blanket)={np.mean([abs(r['W'][j]-r['blanket'])<1e-9 for r in rows if j<len(r['W'])]):.3f}")
    sat=Counter()
    for r in rows:
        for j,w in enumerate(r['W']):
            if abs(w-r['blanket'])<1e-9: sat[j]+=1; break
        else: sat['never']+=1
    print("   j_sat:",dict(sorted(sat.items(),key=lambda t:str(t[0]))),
          " P(>=2 strict rises)=%.3f"%np.mean([sum(1 for a,b in zip(r['W'],r['W'][1:]) if b-a>1e-9)>=2 for r in rows]))
    W1=np.array([r['w_crh1'] for r in rows]); WL=np.array([r['w_loc1'] for r in rows]); WLA=np.array([r['w_locall'] for r in rows])
    print(f"   CRH_1      : sum-ratio {W1.sum()/B.sum():.3f}  coverage {np.mean([r['cov_crh1'] for r in rows]):.3f}")
    print(f"   CRH^loc_1  : sum-ratio {WL.sum()/B.sum():.3f}  coverage {np.mean([r['cov_loc1'] for r in rows]):.3f}  (local facts/problem {np.mean([r['nloc'] for r in rows]):.2f})")
    print(f"   CRH^loc_all: sum-ratio {WLA.sum()/B.sum():.3f}  coverage {np.mean([r['cov_locall'] for r in rows]):.3f}")

print("\n===== atom masses of the reported width (non-trivial stratum) =====")
for mm in (1,3):
    rows,_,_ = sweep(23+mm, 80 if mm==3 else 150, mm)
    for j in (0,1):
        v=np.array([r['W'][j]/r['blanket'] for r in rows if j<len(r['W'])])
        print(f"  m*>={mm} j={j}: P(=0)={np.mean(v<1e-9):.3f} P(=blanket)={np.mean(np.abs(v-1)<1e-9):.3f} "
              f"P(interior)={np.mean((v>1e-9)&(v<1-1e-9)):.3f} distinct values={len(np.unique(np.round(v,6)))} n={len(v)}")
