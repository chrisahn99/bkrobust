"""Bigger + hop-stratified + SCM-level clustering for a bootstrap."""
import sys, json, time
from collections import defaultdict
sys.path.insert(0,"/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/falsegain/code")
from fglib import *
import fglib
import x2lib as X
import numpy as np, importlib.util
from graphs import is_undirected, random_dag
spec = importlib.util.spec_from_file_location("x1_ops","/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x1/code/x1_ops.py")
x1_ops = importlib.util.module_from_spec(spec); spec.loader.exec_module(x1_ops)

def true_bk_mpdag(C,D,rng):
    G=C.copy(); U=undirected_edges(G)
    if not U: return G
    k=int(rng.integers(0,len(U)+1))
    for t in rng.permutation(len(U))[:k]:
        u,v=U[t]
        if not is_undirected(G,u,v): continue
        if D[u,v]==1: G[v,u]=0
        else: G[u,v]=0
        G=meek_closure(G)
    return G

def run(nseed,ps,degs,seed0):
    rng=np.random.default_rng(seed0); acc=defaultdict(int); per_scm=[]
    for it in range(nseed):
        p=int(rng.choice(ps)); deg=float(rng.choice(degs))
        D=random_dag(p,deg,rng); C=dag_to_cpdag(D); ref=v_structures(C)
        G0=true_bk_mpdag(C,D,rng)
        U0=undirected_edges(G0); NA=x1_ops.nonadjacent_pairs(C)
        loc=defaultdict(int)
        for x in range(p):
            for y in range(p):
                if x==y: continue
                s0,O0=report_state(G0,x,y)
                if s0!=POS: continue
                dist=X.hop_dist_inf(G0,[x,y])
                for (u,v) in U0:
                    dmin,_=X.stmt_dist(dist,(u,v))
                    for (a,b) in ((u,v),(v,u)):
                        if D[a,b]==1: continue
                        H=G0.copy(); H[b,a]=0; H=meek_closure(H)
                        if has_directed_cycle(H) or v_structures(H)!=ref: continue
                        s1,O1=report_state(H,x,y)
                        if s1!=POS: loc[("R","exit")]+=1; continue
                        mv=(O1!=O0); sd=sound(D,x,y,s1,O1)
                        loc[("R","mv" if mv else "sm","V" if sd else "I")]+=1
                        if mv: loc[("R","mv",X.dbucket(dmin),"V" if sd else "I")]+=1
                for (a0,b0) in NA:
                    dmin,_=X.stmt_dist(dist,(a0,b0))
                    for (a,b) in ((a0,b0),(b0,a0)):
                        H,info=x1_ops.bk_assert(G0,[(a,b)])
                        ok=(not info["conflict"]) and (not has_directed_cycle(H)) and x1_ops.pdag_extendable(H)
                        tag="S1" if ok else "S1rej"
                        s1,O1=report_state(H,x,y)
                        if s1!=POS: loc[(tag,"exit")]+=1; continue
                        mv=(O1!=O0); sd=sound(D,x,y,s1,O1)
                        loc[(tag,"mv" if mv else "sm","V" if sd else "I")]+=1
                        if mv: loc[(tag,"mv",X.dbucket(dmin),"V" if sd else "I")]+=1
        per_scm.append(dict(loc))
        for k,v in loc.items(): acc[k]+=v
        if len(fglib._valid_cache)>500000: fglib._valid_cache.clear(); fglib._cn_cache.clear()
    return acc, per_scm

if __name__=="__main__":
    n=int(sys.argv[1]); ps=[int(z) for z in sys.argv[2].split(",")]; s=int(sys.argv[3])
    degs=[float(z) for z in sys.argv[4].split(",")] if len(sys.argv)>4 else [1.5,2.0,2.5,3.0]
    t=time.time(); acc,per=run(n,ps,degs,s); print("elapsed",round(time.time()-t,1),"n_scm",len(per))
    for k,v in sorted(acc.items(),key=lambda z:str(z[0])): print("|".join(map(str,k)),v)
    # SCM-level cluster bootstrap on the S1rej invalid-move rate and the S1 invalid-move rate
    rng=np.random.default_rng(0)
    def rate(sample,tag):
        mv=sum(d.get((tag,"mv","V"),0)+d.get((tag,"mv","I"),0) for d in sample)
        iv=sum(d.get((tag,"mv","I"),0) for d in sample)
        return iv/mv if mv else float('nan')
    for tag in ["S1","S1rej"]:
        bs=[]
        for _ in range(2000):
            idx=rng.integers(0,len(per),len(per)); bs.append(rate([per[i] for i in idx],tag))
        bs=[b for b in bs if b==b]
        print(f"BOOT {tag} invalid|moved: point={rate(per,tag):.4f} CI95=[{np.percentile(bs,2.5):.4f},{np.percentile(bs,97.5):.4f}]")
    json.dump({"acc":{"|".join(map(str,k)):v for k,v in acc.items()},"n_scm":len(per)},
              open(f"move2_{'_'.join(map(str,ps))}_{s}.json","w"),indent=1)
