"""Geometry of the counterexamples + spot checks on X1-P1 and the eligibility filter."""
import json, sys, math
import numpy as np
from collections import Counter
sys.path.insert(0,"/home/costaj/latent-causal/x1-spurious/code")
from graphs import random_dag, dag_to_cpdag, skeleton, undirected_edges
from adjust import possibly_causal_paths, cov_linear, total_effect_linear
from scm import make_linear_iscm
from x1_ops import draw_suni, hop_dist_from
UNREACH=1<<20; MAX_K=4
def silent(m): return bool(m.get("consistent") and m.get("amenable") and m.get("ostar_valid") is False)
def regen_dist(seed,p,deg):
    rng=np.random.default_rng(seed); D=random_dag(p,deg,rng); C=dag_to_cpdag(D)
    U=undirected_edges(C)
    pairs=[(a,b) for a in range(p) for b in range(p) if a!=b and len(possibly_causal_paths(D,a,b))>0]
    x,y=pairs[int(rng.integers(len(pairs)))]
    make_linear_iscm(D,rng)
    K_all=[(i,j) if D[i,j]==1 else (j,i) for (i,j) in U]
    order=rng.permutation(len(K_all)); K_all=[K_all[i] for i in order]
    srng=np.random.default_rng([int(seed),0xA55E27]); reps=draw_suni(C,MAX_K,srng)
    dist=hop_dist_from(skeleton(C),[x,y],p)
    return dist,reps,x,y
for ens in ("original","licensed","large"):
    raw=json.load(open(f"/home/costaj/latent-causal/x1-spurious/results/x1_{ens}.json"))
    scms=raw["scms"]
    print(f"== {ens}")
    # X1-P1 spot check: non-extendable but acyclic at rho=1
    na=nc=nx=ncoh=n=0
    for s in scms:
        for m in s["Suni"]:
            if m["rho"]!=1: continue
            n+=1
            if not m["ext"] and not m["cycle"]: na+=1
            if m["cycle"]: nc+=1
            if m["conflict"]: nx+=1
            if m["consistent_S1"]: ncoh+=1
    print(f"  rho=1 n={n}: conflict {nx/n:.4f}  cycle {nc/n:.4f}  "
          f"non-ext&acyclic {na} (report claims 0)  C_S1(=1-coh) {1-ncoh/n:.4f}")
    # geometry of the coherent counterexamples
    mx=Counter(); pairsd=Counter(); nev=0
    for s in scms:
        hits=[m for m in s["Suni"] if m["rho"]==1 and 1<=m["hopC"]<UNREACH and silent(m) and m.get("consistent_S1")]
        if not hits: continue
        dist,reps,x,y=regen_dist(s["seed"],s["p"],s["deg"])
        for m in hits:
            a,b=reps[m["flip"][0]]
            da,db=int(dist[a]),int(dist[b])
            lo,hi=min(da,db),max(da,db)
            nev+=1; mx[hi if hi<UNREACH else "unreach"]+=1; pairsd[(lo,hi if hi<UNREACH else "U")]+=1
    print(f"  coherent hop>=1 counterexamples n={nev}: max-endpoint-dist {dict(sorted(mx.items(),key=str))}")
    print(f"     (min,max) endpoint distances: {dict(sorted(pairsd.items(),key=str))}")
    # eligibility: is arm-R damage correlated with |N(C)|?
    bins={}
    for s in scms:
        bkt = "nN<=4" if s["n_N"]<=4 else ("nN 5-8" if s["n_N"]<=8 else "nN>8")
        for m in s["R"]:
            if m["rho"]!=1 or m["hopC"]>=UNREACH or not m.get("consistent"): continue
            c=bins.setdefault(bkt,[0,0]); c[1]+=1; c[0]+=silent(m)
    print("  arm R|cons silent by |N(C)| bucket:", {k:f"{v[0]}/{v[1]}={v[0]/v[1]:.4f}" for k,v in sorted(bins.items())})
