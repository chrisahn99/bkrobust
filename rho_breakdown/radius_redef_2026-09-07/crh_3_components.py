#!/usr/bin/env python3
"""Is the saturation of CRH_1 governed by the CHORDAL-COMPONENT structure of Chat?
Component = connected component of the undirected subgraph of the CPDAG."""
import sys, itertools, time
from collections import Counter, defaultdict
sys.path.insert(0,"/Users/josecosta/bkrobust/src"); sys.path.insert(0,"/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np, networkx as nx
from bkrobust.demo.graph import canon
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import optimal_adjustment_set_dag, random_sem, adjusted_estimand
from harness import random_dag

def tau_range(g,sem,x,y,c):
    k=(g.edge_string(),x,y)
    if k in c: return c[k]
    v=[]
    for D in enumerate_dag_extensions(g):
        try: v.append(adjusted_estimand(sem,x,y,optimal_adjustment_set_dag(D,x,y)))
        except Exception: pass
    c[k]=(min(v),max(v)) if v else None; return c[k]
def orb(cp,g0):
    u={canon(a,b) for a,b in cp.undirected_edges}
    return sorted((a,b) for a,b in g0.directed_edges if canon(a,b) in u)
def mingen(cp,g0,Or):
    gs=[]
    for r in range(len(Or)+1):
        for S in itertools.combinations(Or,r):
            if any(set(G)<=set(S) for G in gs): continue
            h=apply_orientations(cp,S)
            if h is not None and h==g0: gs.append(S)
    return gs
def hull(cp,subs,sem,x,y,c):
    lo,hi,ok=np.inf,-np.inf,False
    for T in subs:
        g=apply_orientations(cp,T)
        if g is None: continue
        r=tau_range(g,sem,x,y,c)
        if r is None: continue
        ok=True; lo=min(lo,r[0]); hi=max(hi,r[1])
    return (lo,hi) if ok else None
def comps(cp):
    G=nx.Graph(); G.add_nodes_from(cp.nodes); G.add_edges_from([tuple(e) for e in cp.undirected_edges])
    return [set(c) for c in nx.connected_components(G) if len(c)>1]

for (NN,PP,KMAX) in ((7,0.32,7),(12,0.16,9)):
    rng=np.random.default_rng(5); rows=[]; t0=time.time(); tried=0
    while len(rows)<120 and time.time()-t0<300 and tried<400000:
        tried+=1
        dag=random_dag(rng,NN,PP)
        if not dag.directed_edges: continue
        cp=dag_to_cpdag(dag); und=sorted(canon(*e) for e in cp.undirected_edges)
        if not (2<=len(und)<=KMAX): continue
        CS=comps(cp)
        nodes=sorted(dag.nodes); cand=[(a,b) for a in nodes for b in nodes if a!=b]
        rng.shuffle(cand); sem=random_sem(dag,rng); c={}
        for x,y in cand[:6]:
            bl=tau_range(cp,sem,x,y,c)
            if bl is None or bl[1]-bl[0]<1e-8: continue
            km=int(rng.integers(1,min(5,len(und))+1)); idx=rng.permutation(len(und))[:km]
            K=[]
            for t in idx:
                a,b=und[t]; K.append((a,b) if (a,b) in dag.directed_edges else (b,a))
            g0=apply_orientations(cp,K)
            if g0 is None: continue
            Or=orb(cp,g0); gs=mingen(cp,g0,Or)
            if not gs: continue
            S=min(gs,key=len); mstar=len(S)
            touched={i for i,cc in enumerate(CS) for (a,b) in S if a in cc or b in cc}
            r1=hull(cp,[T for G in gs for T in itertools.combinations(G,max(len(G)-1,0))],sem,x,y,c)
            if r1 is None: continue
            rows.append(dict(ncomp=len(CS),ntouch=len(touched),mstar=mstar,k=len(und),
                             blanket=bl[1]-bl[0],w1=r1[1]-r1[0]))
            break
    n=len(rows); B=np.array([r['blanket'] for r in rows]); W=np.array([r['w1'] for r in rows])
    print(f"\n=== n_nodes={NN} p={PP}: {n} non-trivial problems ({tried} draws, {time.time()-t0:.0f}s) ===")
    print("  #chordal components of Chat:", dict(sorted(Counter(r['ncomp'] for r in rows).items())),
          " #components carrying a generating fact:", dict(sorted(Counter(r['ntouch'] for r in rows).items())))
    print(f"  overall: CRH_1 sum-ratio {W.sum()/B.sum():.3f}  P(CRH_1=blanket) {np.mean(np.abs(W-B)<1e-9):.3f}")
    for t in sorted(set(r['ntouch'] for r in rows)):
        m=np.array([r['ntouch']==t for r in rows])
        if m.sum()>=5:
            print(f"   facts spread over {t} component(s): n={m.sum():3d}  sum-ratio {W[m].sum()/B[m].sum():.3f}  P(=blanket) {np.mean(np.abs(W[m]-B[m])<1e-9):.3f}")
