"""REFUTER: (i) cluster-bootstrap the rho-trend of the headline RR;
(ii) discriminate the two hop>=1 damage mechanisms the implementer declared NOT RUN."""
import json, sys
import numpy as np
from itertools import combinations
sys.path.insert(0, "/home/costaj/latent-causal/x1-spurious/code")
from graphs import (apply_background_knowledge, dag_to_cpdag, random_dag,
                    skeleton, undirected_edges, meek_closure, is_directed)
from adjust import (cov_linear, is_valid_adjustment_set, optimal_adjustment_set,
                    possibly_causal_paths, total_effect_linear, causal_nodes, forb)
from scm import make_linear_iscm
from x1_ops import bk_assert, draw_suni, hop_dist_from, stmt_hop, pdag_extendable
MAX_K = 4; UNREACH = 1 << 20

def silent(m): return bool(m.get("consistent") and m.get("amenable") and m.get("ostar_valid") is False)

def boot_rr(scms, rho, B=4000, seed=11):
    """paired SCM cluster bootstrap of RR = P(silent|Suni)/P(silent|cons,R), reachable."""
    a=[];b=[];c=[];d=[]
    for s in scms:
        na=nb=nc=nd=0
        for m in s["Suni"]:
            if m["rho"]!=rho or m["hopC"]>=UNREACH: continue
            nb+=1; na+=silent(m)
        for m in s["R"]:
            if m["rho"]!=rho or m["hopC"]>=UNREACH or not m.get("consistent"): continue
            nd+=1; nc+=silent(m)
        a.append(na);b.append(nb);c.append(nc);d.append(nd)
    a=np.array(a,float);b=np.array(b,float);c=np.array(c,float);d=np.array(d,float)
    pt=(a.sum()/b.sum())/(c.sum()/d.sum())
    rng=np.random.default_rng(seed); n=len(a); out=[]
    for _ in range(B):
        i=rng.integers(0,n,n)
        if b[i].sum()==0 or d[i].sum()==0 or c[i].sum()==0: continue
        out.append((a[i].sum()/b[i].sum())/(c[i].sum()/d[i].sum()))
    out=np.array(out)
    return pt, float(np.percentile(out,2.5)), float(np.percentile(out,97.5)), int(a.sum()), int(b.sum()), int(c.sum()), int(d.sum())

def boot_harm(scms, rho=1, B=4000, seed=13):
    """UNCONDITIONAL total harm (silent OR non-amenable) RR, reachable."""
    a=[];b=[];c=[];d=[]
    for s in scms:
        na=nb=nc=nd=0
        for m in s["Suni"]:
            if m["rho"]!=rho or m["hopC"]>=UNREACH: continue
            nb+=1; na += (silent(m) or not m.get("amenable"))
        for m in s["R"]:
            if m["rho"]!=rho or m["hopC"]>=UNREACH: continue
            nd+=1; nc += (silent(m) or not m.get("amenable"))
        a.append(na);b.append(nb);c.append(nc);d.append(nd)
    a=np.array(a,float);b=np.array(b,float);c=np.array(c,float);d=np.array(d,float)
    pt=(a.sum()/b.sum())/(c.sum()/d.sum())
    rng=np.random.default_rng(seed); n=len(a); out=[]
    for _ in range(B):
        i=rng.integers(0,n,n)
        out.append((a[i].sum()/b[i].sum())/(c[i].sum()/d[i].sum()))
    out=np.array(out)
    return pt, float(np.percentile(out,2.5)), float(np.percentile(out,97.5)), a.sum()/b.sum(), c.sum()/d.sum()

def regen(seed, p, deg):
    rng = np.random.default_rng(seed)
    D = random_dag(p, deg, rng); C = dag_to_cpdag(D); U = undirected_edges(C)
    pairs = [(a,b) for a in range(p) for b in range(p)
             if a != b and len(possibly_causal_paths(D,a,b)) > 0]
    x, y = pairs[int(rng.integers(len(pairs)))]
    scm = make_linear_iscm(D, rng); Sigma = cov_linear(scm["A"], scm["omega"])
    tau = total_effect_linear(scm["A"], x, y)
    K_all = [(i,j) if D[i,j]==1 else (j,i) for (i,j) in U]
    order = rng.permutation(len(K_all)); K_all=[K_all[i] for i in order]
    K=[(int(a),int(b)) for (a,b) in K_all[:MAX_K]]
    srng = np.random.default_rng([int(seed),0xA55E27])
    reps_uni = draw_suni(C, MAX_K, srng)
    return D,C,Sigma,x,y,tau,K,reps_uni

def mechanism(ens, only_coherent):
    raw=json.load(open(f"/home/costaj/latent-causal/x1-spurious/results/x1_{ens}.json"))
    scms=raw["scms"]
    n_on_path=n_off_path=0; n_cases=0
    examples=[]
    for s in scms:
        hits=[m for m in s["Suni"] if m["rho"]==1 and 1<=m["hopC"]<UNREACH and silent(m)
              and ((not only_coherent) or m.get("consistent_S1"))]
        if not hits: continue
        D,C,Sigma,x,y,tau,K,reps=regen(s["seed"],s["p"],s["deg"])
        for m in hits:
            t=m["flip"][0]; a,b=reps[t]
            Kp=[reps[i] if i==t else K[i] for i in range(4)]
            G,info=bk_assert(C,Kp)
            paths=possibly_causal_paths(G,x,y)
            on=False
            for pth in paths:
                for u,v in zip(pth,pth[1:]):
                    if (u,v)==(a,b) or (u,v)==(b,a): on=True;break
                if on: break
            n_cases+=1
            if on: n_on_path+=1
            else:
                n_off_path+=1
                if len(examples)<6: examples.append((s["seed"],s["p"],(x,y),(a,b),m["hopC"]))
    tag = "coherent only (consistent_S1)" if only_coherent else "all"
    print(f"  [{ens}] hop>=1 silent events, {tag}: n={n_cases}  "
          f"added edge ON a proper possibly-causal path X..Y in G: {n_on_path} "
          f"({n_on_path/max(n_cases,1):.3f});  OFF-path (=> Meek propagation): {n_off_path}")
    if examples: print("    off-path examples (seed,p,(x,y),(a,b),hop):", examples)

if __name__=="__main__":
    for ens in ("original","licensed","large"):
        raw=json.load(open(f"/home/costaj/latent-causal/x1-spurious/results/x1_{ens}.json"))
        scms=raw["scms"]
        print(f"== {ens}: headline RR by rho, paired SCM cluster bootstrap B=4000")
        for rho in (1,2,3,4):
            pt,lo,hi,ka,na,kc,nc=boot_rr(scms,rho)
            print(f"   rho={rho}  RR={pt:.3f} [{lo:.3f},{hi:.3f}]   S {ka}/{na}   R|cons {kc}/{nc}")
        pt,lo,hi,ps,pr=boot_harm(scms)
        print(f"   UNCONDITIONAL TOTAL HARM (silent or non-amenable), rho=1: "
              f"S {ps:.4f}  R {pr:.4f}  RR={pt:.3f} [{lo:.3f},{hi:.3f}]")
    print("\n== mechanism discrimination for hop>=1 silent events (rho=1)")
    for ens in ("original","licensed","large"):
        mechanism(ens, False)
        mechanism(ens, True)
