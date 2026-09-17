"""EXPERIMENT FG-S: false identifiability under class-S (spurious required edge).

Base G0 = M(C,K), K TRUE (D in [G0]).  Statement s = a->b on a pair NON-ADJACENT
in skeleton(C) = skeleton(D): guaranteed false, and inexpressible in Meek's
consistency check (apply_background_knowledge raises).  Operator = x1_ops.bk_assert.
Two acceptance views: S1 = Dor-Tarsi extendable & acyclic & no conflict; S2 = b-LOAD (accept all).
"""
import sys, json
from collections import defaultdict
from multiprocessing import Pool
sys.path.insert(0,"/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/falsegain/code")
sys.path.insert(0,"/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x1/code")
from fglib import *
import fglib
import numpy as np
import importlib.util
spec = importlib.util.spec_from_file_location("x1_ops","/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x1/code/x1_ops.py")
x1_ops = importlib.util.module_from_spec(spec); spec.loader.exec_module(x1_ops)

def job(arg):
    Cb, p = arg
    C = np.frombuffer(Cb,dtype=np.int8).reshape(p,p).copy()
    mp = reachable_mpdags(C)
    mec = [G for G in mp if len(undirected_edges(G))==0]
    NA = x1_ops.nonadjacent_pairs(C)
    acc = defaultdict(int); wit=[]
    for D in mec:
        for G0 in mp:
            if not dag_agrees_with(D,G0): continue
            for x in range(p):
                for y in range(p):
                    if x==y: continue
                    s0,O0 = report_state(G0,x,y)
                    dist = X.hop_dist_inf(G0,[x,y])
                    for (a0,b0) in NA:
                        dmin,_ = X.stmt_dist(dist,(a0,b0))
                        bkt = X.dbucket(dmin)
                        for (a,b) in ((a0,b0),(b0,a0)):
                            H, info = x1_ops.bk_assert(G0, [(a,b)])
                            ok1 = (not info["conflict"]) and (not has_directed_cycle(H)) \
                                  and x1_ops.pdag_extendable(H)
                            s1,O1 = report_state(H,x,y)
                            snd1 = sound(D,x,y,s1,O1)
                            acc[("cell",bkt,"S1" if ok1 else "S1rej",s0,s1,snd1)]+=1
                            if ok1 and s0==REFUSE and s1 in (POS,ZERO) and snd1 is False and len(wit)<25:
                                wit.append(dict(p=int(p),D=D.tolist(),C=C.tolist(),G0=G0.tolist(),
                                  H=H.tolist(),x=int(x),y=int(y),s=[int(a),int(b)],dmin=float(dmin),
                                  state1=s1,O1=sorted(int(z) for z in O1) if O1 is not None else None))
    fglib._valid_cache.clear(); fglib._cn_cache.clear()
    return dict(acc), wit

def main():
    p=int(sys.argv[1]); nw=int(sys.argv[2]) if len(sys.argv)>2 else 10
    ncap=int(sys.argv[3]) if len(sys.argv)>3 else 0
    from run_lemma import all_dags
    cp={}
    for D in all_dags(p):
        C=dag_to_cpdag(D); cp.setdefault(C.tobytes(),C)
    keys=list(cp.keys())
    if ncap and ncap<len(keys):
        rng=np.random.default_rng(7); keys=[keys[i] for i in rng.choice(len(keys),ncap,replace=False)]
    print(f"p={p}: {len(keys)} CPDAGs",flush=True)
    tot=defaultdict(int); wits=[]; n=0
    with Pool(nw) as pool:
        for acc,w in pool.imap_unordered(job,[(k,p) for k in keys],chunksize=2):
            n+=1
            for k,v in acc.items(): tot[k]+=v
            for z in w:
                if len(wits)<25: wits.append(z)
            if n%500==0: print(f"  {n}/{len(keys)}",flush=True)
    out={"|".join(map(str,k)):v for k,v in tot.items()}
    json.dump(dict(p=p,n_cpdags=len(keys),counts=out,witnesses=wits),
        open(f"/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/falsegain/code/fgS_p{p}.json","w"),indent=1)
    for k,v in sorted(out.items()): print(k,v)

if __name__=="__main__": main()
