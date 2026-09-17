"""Does the tau(D)=adjust-on-O*(D) convention change SEL_A's own numbers?"""
import sys, itertools, time
sys.path.insert(0,"/Users/josecosta/bkrobust/src")
sys.path.insert(0,"/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_dag, is_valid_adjustment_set_dag,
                                    adjusted_estimand)
from audit_sel_2_coverage import build, Z196

def hull(cp, claims, sem, x, y, filt):
    g = apply_orientations(cp, list(claims))
    if g is None: return None
    vals=[]
    for D in enumerate_dag_extensions(g):
        try:
            z = optimal_adjustment_set_dag(D,x,y)
            if filt and not is_valid_adjustment_set_dag(D,x,y,z):  continue
            vals.append(adjusted_estimand(sem,x,y,z))
        except Exception: pass
    return (min(vals),max(vals)) if vals else None

rng=np.random.default_rng(2024); rows=[]; t0=time.time()
while len(rows)<300 and time.time()-t0<180:
    pr=build(rng,0)
    if pr is None: continue
    out={}
    for filt in (False,True):
        m=len(pr['K']); th0=pr['theta0']; d=Z196*np.sqrt(pr['avar']/10000)
        D=[]
        for i in range(m):
            h=hull(pr['cp'],[pr['K'][j] for j in range(m) if j!=i],pr['sem'],pr['x'],pr['y'],filt)
            D.append(0.0 if h is None else max(abs(h[0]-th0),abs(h[1]-th0)))
        b=hull(pr['cp'],pr['K'],pr['sem'],pr['x'],pr['y'],filt)  # sanity: should be a point
        bl=hull(pr['cp'],[],pr['sem'],pr['x'],pr['y'],filt)
        out[filt]=dict(D=np.array(D), fire=any(x>d for x in D),
                       base_w=0 if b is None else b[1]-b[0],
                       blanket=0 if bl is None else bl[1]-bl[0])
    rows.append(out)
print(f"problems={len(rows)}")
for filt in (False,True):
    f=np.mean([r[filt]['fire'] for r in rows]); bw=np.mean([r[filt]['blanket'] for r in rows])
    base=np.mean([r[filt]['base_w'] for r in rows])
    mx=np.mean([r[filt]['D'].max() for r in rows])
    print(f"  tau over {'AMENABLE extensions only' if filt else 'ALL extensions (as coded)'}: "
          f"screen fires {100*f:.1f}%  mean blanket width {bw:.4f}  mean max Delta {mx:.4f}  "
          f"mean width of hull(K) [must be 0] {base:.2e}")
chg=np.mean([r[False]['fire']!=r[True]['fire'] for r in rows])
print(f"  the two conventions disagree on whether the screen fires in {100*chg:.1f}% of problems")
