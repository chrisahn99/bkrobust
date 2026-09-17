"""Targeted random search, p=5..9:
   (A) can class R (`add`, Meek-consistent) make a POS report UNSOUND?      [expect: no]
   (B) can class S (spurious, S1-accepted) make a POS report UNSOUND?       [expect: ?]
   (C) how far (dmin) can a FALSE GAIN reach?
"""
import sys, json, time
from collections import defaultdict
sys.path.insert(0,"/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/falsegain/code")
from fglib import *
import fglib, numpy as np, importlib.util
from graphs import is_undirected
from graphs import random_dag
spec = importlib.util.spec_from_file_location("x1_ops","/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x1/code/x1_ops.py")
x1_ops = importlib.util.module_from_spec(spec); spec.loader.exec_module(x1_ops)

def true_bk_mpdag(C, D, rng):
    """M(C,K) for a random TRUE K: orient a random subset of undirected edges as in D."""
    G = C.copy()
    U = undirected_edges(G)
    if not U: return G
    k = int(rng.integers(0, len(U)+1))
    order = rng.permutation(len(U))
    for t in order[:k]:
        u,v = U[t]
        if not is_undirected(G,u,v): continue
        if D[u,v]==1: G[v,u]=0
        else: G[u,v]=0
        G = meek_closure(G)
    return G

def run(nseed, ps, degs, seed0=0):
    rng = np.random.default_rng(seed0)
    acc = defaultdict(int); wits=defaultdict(list); maxd=defaultdict(float)
    for it in range(nseed):
        p = int(rng.choice(ps)); deg = float(rng.choice(degs))
        D = random_dag(p, deg, rng)
        C = dag_to_cpdag(D)
        ref = v_structures(C)
        G0 = true_bk_mpdag(C, D, rng)
        assert dag_agrees_with(D, G0)
        U0 = undirected_edges(G0); NA = x1_ops.nonadjacent_pairs(C)
        for x in range(p):
            for y in range(p):
                if x==y: continue
                s0,O0 = report_state(G0,x,y)
                snd0 = sound(D,x,y,s0,O0)
                acc[("base",s0,snd0)]+=1
                if snd0 is False: acc[("BASE_UNSOUND",s0)]+=1
                dist = X.hop_dist_inf(G0,[x,y])
                # ---- class R
                for (u,v) in U0:
                    dmin,_ = X.stmt_dist(dist,(u,v))
                    for (a,b) in ((u,v),(v,u)):
                        if D[a,b]==1: continue        # only FALSE statements
                        H = G0.copy(); H[b,a]=0; H=meek_closure(H)
                        if has_directed_cycle(H) or v_structures(H)!=ref: continue
                        s1,O1 = report_state(H,x,y); snd1=sound(D,x,y,s1,O1)
                        acc[("R",s0,s1,snd1)]+=1
                        if s0==POS and snd1 is False and len(wits["R_POS_unsound"])<3:
                            wits["R_POS_unsound"].append(dict(p=p,D=D.tolist(),G0=G0.tolist(),H=H.tolist(),x=x,y=y,s=[int(a),int(b)]))
                        if s0==REFUSE and s1==POS and snd1 is False:
                            maxd["R_falsegain"]=max(maxd["R_falsegain"],dmin)
                            acc[("R_fg_dmin",X.dbucket(dmin))]+=1
                # ---- class S
                for (a0,b0) in NA:
                    dmin,_ = X.stmt_dist(dist,(a0,b0))
                    for (a,b) in ((a0,b0),(b0,a0)):
                        H,info = x1_ops.bk_assert(G0,[(a,b)])
                        ok1 = (not info["conflict"]) and (not has_directed_cycle(H)) and x1_ops.pdag_extendable(H)
                        s1,O1 = report_state(H,x,y); snd1=sound(D,x,y,s1,O1)
                        tag = "S1" if ok1 else "S1rej"
                        acc[(tag,s0,s1,snd1)]+=1
                        if ok1 and s0==POS and snd1 is False and len(wits["S1_POS_unsound"])<3:
                            wits["S1_POS_unsound"].append(dict(p=p,D=D.tolist(),G0=G0.tolist(),H=H.tolist(),x=x,y=y,s=[int(a),int(b)]))
                        if ok1 and s0==REFUSE and s1==POS and snd1 is False:
                            maxd["S1_falsegain"]=max(maxd["S1_falsegain"],dmin)
                            acc[("S1_fg_dmin",X.dbucket(dmin))]+=1
                        if ok1 and s0==ZERO and s1==POS and snd1 is False:
                            acc[("S1_fabricate_dmin",X.dbucket(dmin))]+=1
                            if len(wits["S1_ZERO_to_POS_unsound"])<3:
                                wits["S1_ZERO_to_POS_unsound"].append(dict(p=p,D=D.tolist(),G0=G0.tolist(),H=H.tolist(),x=x,y=y,s=[int(a),int(b)]))
        if len(fglib._valid_cache)>400000: fglib._valid_cache.clear(); fglib._cn_cache.clear()
        if it%200==0: print(f"  it={it}",flush=True)
    return acc, wits, maxd

if __name__=="__main__":
    n=int(sys.argv[1]); ps=[int(z) for z in sys.argv[2].split(",")]; s0=int(sys.argv[3])
    t=time.time()
    acc,wits,maxd = run(n, ps, [1.5,2.0,2.5,3.0], s0)
    print("elapsed", time.time()-t)
    for k,v in sorted(acc.items(), key=lambda z: str(z[0])): print("|".join(map(str,k)), v)
    print("MAXD", dict(maxd))
    json.dump({"counts":{"|".join(map(str,k)):v for k,v in acc.items()},
               "wits":{k:v for k,v in wits.items()},"maxd":dict(maxd)},
              open(f"/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/falsegain/code/search_{'_'.join(map(str,ps))}_{s0}.json","w"),indent=1)
