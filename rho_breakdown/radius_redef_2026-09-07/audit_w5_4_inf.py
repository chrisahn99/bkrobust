"""Where do the INF certificates go wrong? Is it just non-amenability?"""
import sys, math
sys.path.insert(0, "/Users/josecosta/bkrobust/src"); sys.path.insert(0, ".")
import numpy as np
from bkrobust.demo.evaluate import (optimal_adjustment_set_dag, is_valid_adjustment_set_dag,
                                    adjusted_estimand, random_sem)
from audit_w0_rk import m_of, gen
INF = math.inf
def parents(d,x): return {a for (a,b) in d.directed_edges if b==x}
for nf in (1,2):
    rng=np.random.default_rng(4242+nf); probs,_=gen(rng,500,n=7,p=0.40,k_claims=3,n_false=nf)
    rng2=np.random.default_rng(31337+nf); L=[]
    for pr in probs:
        dag,x,y,Z,K=pr['dag'],pr['x'],pr['y'],pr['z'],pr['K']
        sem=random_sem(dag,rng2); bZ=adjusted_estimand(sem,x,y,Z); tau=sem.true_total_effect(x,y)
        rec=[(m_of(d,K), adjusted_estimand(sem,x,y,optimal_adjustment_set_dag(d,x,y)),
              adjusted_estimand(sem,x,y,parents(d,x))) for d in pr['fibre']]
        sp=max(abs(b-bZ) for _,b,_ in rec)
        if sp<=1e-9: continue
        amen=is_valid_adjustment_set_dag(dag,x,y,optimal_adjustment_set_dag(dag,x,y))
        L.append(dict(rec=rec,bZ=bZ,tau=tau,amen=amen))
    n=len(L); print(f"\nf={nf}: {n} live problems, non-amenable {100*np.mean([not p['amen'] for p in L]):.1f}%")
    for dl in (.25,.5,1.0):
        for j,lab in ((1,'O*(d)  '),(2,'pa_d(X)')):
            inf=[p for p in L if not any(abs(r[j]-p['bZ'])>dl for r in p['rec'])]
            if not inf: continue
            bad=[p for p in inf if abs(p['bZ']-p['tau'])>dl]
            na=[p for p in inf if not p['amen']]
            badna=[p for p in bad if not p['amen']]
            print(f"  delta={dl:<4.2g} beta={lab} INF on {100*len(inf)/n:5.1f}% | of those "
                  f"|bZ-tau|>delta on {100*len(bad)/len(inf):5.1f}% | non-amenable share of INF "
                  f"{100*len(na)/len(inf):5.1f}% | of the BAD ones, non-amenable {100*len(badna)/max(len(bad),1):5.1f}%")
