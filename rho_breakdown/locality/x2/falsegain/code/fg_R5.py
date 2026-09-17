"""FG-R, parallel, grouped by CPDAG (its MEC = the DAGs mapping to it)."""
import sys, json
from collections import defaultdict
from multiprocessing import Pool
sys.path.insert(0,"/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/falsegain/code")
from fglib import *
import fglib
import numpy as np

def job(arg):
    Cb, p = arg
    C = np.frombuffer(Cb, dtype=np.int8).reshape(p,p).copy()
    ref = v_structures(C)
    mp = reachable_mpdags(C)
    mec = [G for G in mp if len(undirected_edges(G))==0]
    acc = defaultdict(int); wit=[]
    for D in mec:
        for G0 in mp:
            if not dag_agrees_with(D, G0): continue
            U0 = undirected_edges(G0)
            if not U0: continue
            for x in range(p):
                for y in range(p):
                    if x==y: continue
                    s0,O0 = report_state(G0,x,y)
                    acc[("base",s0,sound(D,x,y,s0,O0))]+=1
                    dist = X.hop_dist_inf(G0,[x,y])
                    for (u,v) in U0:
                        dmin,_ = X.stmt_dist(dist,(u,v))
                        b = X.dbucket(dmin)
                        for (a,bb) in ((u,v),(v,u)):
                            H = G0.copy(); H[bb,a]=0; H=meek_closure(H)
                            if has_directed_cycle(H) or v_structures(H)!=ref:
                                acc[("meekfail",b)]+=1; continue
                            truth = "T" if D[a,bb]==1 else "F"
                            s1,O1 = report_state(H,x,y)
                            snd1 = sound(D,x,y,s1,O1)
                            acc[("cell",b,truth,s0,s1,snd1)]+=1
                            if s0==REFUSE and s1 in (POS,ZERO) and snd1 is False and dmin>=1 and len(wit)<25:
                                wit.append(dict(p=int(p),D=D.tolist(),C=C.tolist(),G0=G0.tolist(),
                                    H=H.tolist(),x=int(x),y=int(y),s=[int(a),int(bb)],
                                    dmin=float(dmin),truth=truth,state1=s1,
                                    O1=sorted(int(z) for z in O1) if O1 is not None else None))
    fglib._valid_cache.clear(); fglib._cn_cache.clear()
    return dict(acc), wit

def main():
    p = int(sys.argv[1]); nw = int(sys.argv[2]) if len(sys.argv)>2 else 10
    cp = {}
    from run_lemma import all_dags
    for D in all_dags(p):
        C = dag_to_cpdag(D); cp.setdefault(C.tobytes(), C)
    keys = list(cp.keys())
    print(f"p={p}: {len(keys)} CPDAGs", flush=True)
    tot=defaultdict(int); wits=[]; n=0
    with Pool(nw) as pool:
        for acc,w in pool.imap_unordered(job, [(k,p) for k in keys], chunksize=4):
            n+=1
            for k,v in acc.items(): tot[k]+=v
            for z in w:
                if len(wits)<25: wits.append(z)
            if n%500==0: print(f"  {n}/{len(keys)}", flush=True)
    out = {"|".join(map(str,k)):v for k,v in tot.items()}
    json.dump(dict(p=p,n_cpdags=len(keys),counts=out,witnesses=wits),
              open(f"/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/falsegain/code/fgR_p{p}.json","w"),indent=1)
    for k,v in sorted(out.items()): print(k,v)

if __name__=="__main__": main()
