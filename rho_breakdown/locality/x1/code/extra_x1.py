import json, sys, numpy as np
sys.path.insert(0,"/home/costaj/latent-causal/x1-spurious/code")
from stats_x1 import wilson, cluster_bootstrap_rr
from x1_ops import UNREACH
def silent(m): return bool(m.get("consistent") and m.get("amenable") and m.get("ostar_valid") is False)
def boot(num,den,B=10000,seed=20260819):
    num=np.asarray(num,float); den=np.asarray(den,float); m=len(num)
    pt=num.sum()/den.sum(); rng=np.random.default_rng(seed); o=[];d=0
    while d<B:
        b=min(500,B-d); idx=rng.integers(0,m,size=(b,m))
        with np.errstate(divide="ignore",invalid="ignore"): o.append(num[idx].sum(1)/den[idx].sum(1))
        d+=b
    o=np.concatenate(o); o=o[np.isfinite(o)]
    return pt,int(num.sum()),int(den.sum()),np.quantile(o,0.025),np.quantile(o,0.975)
out={}
for e in ("original","licensed","large"):
    raw=json.load(open(f"/home/costaj/latent-causal/x1-spurious/results/x1_{e}.json"))["scms"]
    n=len(raw)
    # P3 primary: S-loc pure_spurious non-query, member, rho<=4
    a=np.zeros(n); b=np.zeros(n)
    for i,s in enumerate(raw):
        for m in s["Sloc"]:
            if not all(s["lab_loc"][t][0]=="pure_spurious" and not s["lab_loc"][t][1] for t in m["flip"]): continue
            b[i]+=1
            if silent(m): a[i]+=1
    pt,k,nn,lo,hi = boot(a,b)
    # S1 (careful tool) vs R conditional RR, rho=1 reachable
    ns=np.zeros(n); ds=np.zeros(n); nr=np.zeros(n); dr=np.zeros(n)
    for i,s in enumerate(raw):
        for m in s["Suni"]:
            if m["rho"]!=1 or m["hopC"]>=UNREACH or not m["consistent_S1"]: continue
            ds[i]+=1;  ns[i]+= silent(m)
        for m in s["R"]:
            if m["rho"]!=1 or m["hopC"]>=UNREACH or not m["consistent"]: continue
            dr[i]+=1;  nr[i]+= silent(m)
    rr,rlo,rhi = cluster_bootstrap_rr(ns,ds,nr,dr,B=10000)
    # S-loc pooled reversal-stratum share of the pooled damage
    print(f"[{e}] P3 primary (S-loc pure_spurious NON-query, member, rho<=4): "
          f"{pt:.4f} k={k} n={nn} wilson{[round(v,4) for v in wilson(k,nn)]} CLUSTER[{lo:.4f},{hi:.4f}]")
    print(f"      A26 window [0.087,0.100] overlap with CLUSTER CI: {not (hi<0.087 or lo>0.100)}")
    print(f"[{e}] S1(careful tool) vs R conditional RR rho=1 reachable: RR={rr:.3f} [{rlo:.3f},{rhi:.3f}]  "
          f"S1 rate={ns.sum()/ds.sum():.4f} (k={int(ns.sum())} n={int(ds.sum())})  R rate={nr.sum()/dr.sum():.4f} (k={int(nr.sum())} n={int(dr.sum())})")
