"""EXPERIMENT FG-R: false identifiability under class-R (within-skeleton) BK.

Base G0 = M(C,K) for a TRUE K (so D in [G0]) -- the analyst's honest state.
One further statement s: orient an undirected edge of G0.  Exactly one of the
two orientations agrees with D (TRUE), the other does not (FALSE-but-consistent).
Classify (state0,state1,soundness1) x dmin x truth.
"""
import sys, json
from collections import defaultdict
sys.path.insert(0,"/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/falsegain/code")
from fglib import *

def scan_dag(D, C, mp):
    p = D.shape[0]
    ref = v_structures(C)
    acc = defaultdict(int)
    wit = []
    for G0 in mp:
        if not dag_agrees_with(D, G0):
            continue
        U0 = undirected_edges(G0)
        if not U0: continue
        dist_all = {}
        for x in range(p):
            for y in range(p):
                if x == y: continue
                s0, O0 = report_state(G0, x, y)
                snd0 = sound(D, x, y, s0, O0)
                acc[("base", s0, snd0)] += 1
                dist = X.hop_dist_inf(G0, [x, y])
                for (u,v) in U0:
                    dmin,_ = X.stmt_dist(dist, (u,v))
                    b = X.dbucket(dmin)
                    for (a,bb) in ((u,v),(v,u)):
                        H = G0.copy(); H[bb,a]=0; H = meek_closure(H)
                        if has_directed_cycle(H) or v_structures(H)!=ref:
                            acc[("meekfail",b)] += 1; continue
                        truth = "T" if D[a,bb]==1 else "F"
                        s1, O1 = report_state(H, x, y)
                        snd1 = sound(D, x, y, s1, O1)
                        acc[("cell", b, truth, s0, s1, snd1)] += 1
                        if s0==REFUSE and s1 in (POS,ZERO) and snd1 is False:
                            if len(wit) < 40:
                                wit.append(dict(p=int(p), D=D.tolist(), C=C.tolist(),
                                    G0=G0.tolist(), H=H.tolist(), x=int(x), y=int(y),
                                    s=[int(a),int(bb)], dmin=float(dmin), truth=truth,
                                    state1=s1, O1=sorted(int(z) for z in O1) if O1 is not None else None))
    return acc, wit

def main(p, nd=None, seed=0):
    tot = defaultdict(int); wits=[]
    cpc = {}
    import numpy as np
    dags = list(all_dags(p))
    if nd is not None and nd < len(dags):
        rng = np.random.default_rng(seed)
        idx = rng.choice(len(dags), size=nd, replace=False)
        dags = [dags[i] for i in idx]
    for n,D in enumerate(dags):
        C = dag_to_cpdag(D)
        kb = C.tobytes()
        if kb not in cpc:
            cpc[kb] = reachable_mpdags(C)
        acc, w = scan_dag(D, C, cpc[kb])
        for k,v in acc.items(): tot[k]+=v
        for x in w:
            if len(wits)<40: wits.append(x)
        if n % 500 == 0: print(f"  {n}/{len(dags)}", flush=True)
    return tot, wits, len(dags)

if __name__=="__main__":
    p = int(sys.argv[1]); nd = int(sys.argv[2]) if len(sys.argv)>2 and sys.argv[2]!="all" else None
    tot,wits,ndags = main(p, nd)
    out = {"|".join(map(str,k)): v for k,v in sorted(tot.items(), key=lambda t: str(t[0]))}
    json.dump(dict(p=p,ndags=ndags,counts=out,witnesses=wits),
              open(f"/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/falsegain/code/fgR_p{p}.json","w"),indent=1)
    for k,v in sorted(out.items()):
        print(k, v)
