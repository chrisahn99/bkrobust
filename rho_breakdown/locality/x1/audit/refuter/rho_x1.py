import json, numpy as np
UNREACH=1<<20
def silent(m): return bool(m.get("consistent") and m.get("amenable") and m.get("ostar_valid") is False)
def boot(scms,rho,sfilt,B=4000,seed=3):
    a=[];b=[];c=[];d=[]
    for s in scms:
        na=nb=nc=nd=0
        for m in s["Suni"]:
            if m["rho"]!=rho or m["hopC"]>=UNREACH or not sfilt(m): continue
            nb+=1; na+=silent(m)
        for m in s["R"]:
            if m["rho"]!=rho or m["hopC"]>=UNREACH or not m.get("consistent"): continue
            nd+=1; nc+=silent(m)
        a.append(na);b.append(nb);c.append(nc);d.append(nd)
    a=np.array(a,float);b=np.array(b,float);c=np.array(c,float);d=np.array(d,float)
    if b.sum()==0 or c.sum()==0: return None
    pt=(a.sum()/b.sum())/(c.sum()/d.sum())
    rng=np.random.default_rng(seed);n=len(a);o=[]
    for _ in range(B):
        i=rng.integers(0,n,n)
        if b[i].sum()==0 or c[i].sum()==0: continue
        o.append((a[i].sum()/b[i].sum())/(c[i].sum()/d[i].sum()))
    o=np.array(o)
    return pt,float(np.percentile(o,2.5)),float(np.percentile(o,97.5)),int(a.sum()),int(b.sum()),int(c.sum()),int(d.sum())
for ens in ("original","licensed","large"):
    scms=json.load(open(f"/home/costaj/latent-causal/x1-spurious/results/x1_{ens}.json"))["scms"]
    print(f"== {ens}")
    for lab,fl in (("S2 (b-LOAD, registered)",lambda m:True),
                   ("S1 (coherence-checked)",lambda m:m.get("consistent_S1")),
                   ("mpdag_valid only",lambda m:m.get("ext",True) and not m.get("cycle",False))):
        row=[]
        for rho in (1,2,3,4):
            r=boot(scms,rho,fl)
            row.append(f"rho{rho} RR={r[0]:.3f}[{r[1]:.3f},{r[2]:.3f}] ({r[3]}/{r[4]} vs {r[5]}/{r[6]})" if r else f"rho{rho} n/a")
        print(f"  {lab}:")
        for x in row: print("     ",x)
