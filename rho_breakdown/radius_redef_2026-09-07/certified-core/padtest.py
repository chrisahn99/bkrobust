"""Padding immunity: does canonicalisation make c invariant to Meek-implied claims?"""
from __future__ import annotations
import os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import make_problem, canonicalise, Problem, minimal_unsafe, min_hitting_sets
from bkrobust.demo.graph import canon
from bkrobust.demo.meek import apply_orientations

def core_of(cpdag,K,x,y,Z):
    m=len(K); p=Problem(cpdag,K,x,y,Z)
    c,_=min_hitting_sets(minimal_unsafe(p,m),m)
    return c,m,p.n_checks

rng=np.random.default_rng(2718281)
rows=[]; t0=time.time(); tries=0
while len(rows)<400 and tries<200000 and time.time()-t0<600:
    tries+=1
    P=make_problem(rng,1,7)
    if P is None: continue
    cpdag,Kc,G0,x,y,Z=P["cpdag"],P["Kc"],P["G0"],P["x"],P["y"],P["Z"]
    pads=[(t,h) for (t,h) in sorted(G0.directed_edges)
          if canon(t,h) in cpdag.undirected_edges and (t,h) not in Kc]
    if not pads: continue
    Kp=list(Kc)+pads
    if len(Kp)>10: continue
    gp=apply_orientations(cpdag,Kp)
    if gp is None or gp.directed_edges!=G0.directed_edges: continue
    c0,m0,_=core_of(cpdag,Kc,x,y,Z)
    cR,mR,_=core_of(cpdag,Kp,x,y,Z)                 # raw padded, no canonicalisation
    Kpc=canonicalise(cpdag,Kp,G0)
    cC,mC,_=core_of(cpdag,Kpc,x,y,Z)                # padded then canonicalised
    rows.append((m0,c0,mR,cR,mC,cC))
a=np.array(rows)
m0,c0,mR,cR,mC,cC=[a[:,i] for i in range(6)]
k0=1-c0/m0; kR=1-cR/mR
print(f"n={len(rows)} tries={tries} secs={time.time()-t0:.1f}")
print(f"  mean |K*|={m0.mean():.2f} -> padded raw |K|={mR.mean():.2f} -> padded+canonicalised |K*|={mC.mean():.2f}")
print(f"  RAW padding: kappa rose on {np.mean(kR>k0):.4f}, fell on {np.mean(kR<k0):.4f}, mean delta={np.mean(kR-k0):+.4f}, max={np.max(kR-k0):+.4f}")
print(f"  RAW padding: c changed on {np.mean(cR!=c0):.4f} (up {np.mean(cR>c0):.4f}, down {np.mean(cR<c0):.4f})")
print(f"  CANONICALISED padding: |K*| changed on {np.mean(mC!=m0):.4f}; c changed on {np.mean(cC!=c0):.4f}; kappa changed on {np.mean(np.abs((1-cC/mC)-k0)>1e-12):.4f}")
