import json, sys, math
import numpy as np
UNREACH = 1 << 20
def wil(k,n):
    if n==0: return (float('nan'),)*3
    z=1.959963984540054; p=k/n; d=1+z*z/n
    c=(p+z*z/(2*n))/d; h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return p,max(0,c-h),min(1,c+h)
def f(k,n):
    p,lo,hi=wil(k,n); return f"{p:.4f} [{lo:.4f},{hi:.4f}] k={k} n={n}"
def silent(m): return bool(m.get("consistent") and m.get("amenable") and m.get("ostar_valid") is False)
def harm(m):   return bool(silent(m) or (m.get("consistent") and not m.get("amenable")))

tot_hop2_coh_k=tot_hop2_coh_n=0
tot_hop2_all_k=tot_hop2_all_n=0
for ens in ("original","licensed","large"):
    raw=json.load(open(f"/home/costaj/latent-causal/x1-spurious/results/x1_{ens}.json"))
    scms=raw["scms"]
    print(f"== {ens}")
    # (a) unconditional total harm, statement unit, reachable, rho=1, paired cluster bootstrap
    a=[];b=[];c=[];d=[]
    for s in scms:
        na=nb=nc=nd=0
        for m in s["Suni"]:
            if m["rho"]!=1 or m["hopC"]>=UNREACH: continue
            nb+=1; na+=harm(m)
        for m in s["R"]:
            if m["rho"]!=1 or m["hopC"]>=UNREACH: continue
            nd+=1; nc+=harm(m)
        a.append(na);b.append(nb);c.append(nc);d.append(nd)
    a=np.array(a,float);b=np.array(b,float);c=np.array(c,float);d=np.array(d,float)
    ps=a.sum()/b.sum(); pr=c.sum()/d.sum()
    rng=np.random.default_rng(5); n=len(a); dif=[]
    for _ in range(10000):
        i=rng.integers(0,n,n); dif.append(a[i].sum()/b[i].sum()-c[i].sum()/d[i].sum())
    dif=np.array(dif)
    print(f"  UNCOND TOTAL HARM (silent OR loud; MeekFail = caught, NOT harm), rho=1, reachable")
    print(f"    S-uni {f(int(a.sum()),int(b.sum()))}")
    print(f"    R     {f(int(c.sum()),int(d.sum()))}")
    print(f"    diff S-R = {ps-pr:+.4f}  cluster CI [{np.percentile(dif,2.5):+.4f},{np.percentile(dif,97.5):+.4f}]")
    # (b) hop>=2 zero, coherent and not
    for lab,filt in (("all (b-LOAD, no check)", lambda m: True),
                     ("consistent_S1", lambda m: m.get("consistent_S1"))):
        k=n2=0
        for s in scms:
            for m in s["Suni"]:
                if m["rho"]!=1 or m["hopC"]>=UNREACH or m["hopC"]<2: continue
                if not filt(m): continue
                n2+=1; k+=silent(m)
        print(f"    hop>=2 silent, {lab}: {f(k,n2)}")
        if lab=="consistent_S1": tot_hop2_coh_k+=k; tot_hop2_coh_n+=n2
        else: tot_hop2_all_k+=k; tot_hop2_all_n+=n2
    # (c) hop share of the uniform spurious draw (pruning cost)
    from collections import Counter
    cn=Counter()
    for s in scms:
        for m in s["Suni"]:
            if m["rho"]!=1: continue
            h=m["hopC"]; cn["unreach" if h>=UNREACH else (str(h) if h<2 else ">=2")]+=1
    T=sum(cn.values())
    print("    hop share of a uniform N(C) draw:", {k:f"{v/T:.3f}" for k,v in sorted(cn.items())})
print(f"\n POOLED across 3 ensembles, hop>=2, rho=1:")
print(f"   b-LOAD semantics : {f(tot_hop2_all_k,tot_hop2_all_n)}")
print(f"   consistent_S1    : {f(tot_hop2_coh_k,tot_hop2_coh_n)}")
