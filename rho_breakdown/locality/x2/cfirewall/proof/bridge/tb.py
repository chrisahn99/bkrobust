"""
Diagnostics on the bridge.
For each (C, G0, stmt (a,b), H coherent, query (x,y) with G0 & H amenable, D in [G0]):
  - D+ = D + a->b : acyclic? preserves H? uc(D+) subset uc(P) where P=G0+a->b?
  - if D+ in ext(H): is O*(X,Y,H) == O(X,Y,D+) (DAG optimal set)?
  - validity of O*_H in D+ ; validity in D
"""
from common import *
from collections import defaultdict
import sys, json
from itertools import combinations

def ucolliders(W):
    """unshielded colliders u->v<-w with u,w non-adjacent, using DIRECTED edges only"""
    p = W.shape[0]; out=set()
    for v in range(p):
        pars=[u for u in range(p) if is_directed(W,u,v)]
        for u,w in combinations(sorted(pars),2):
            if not adjacent(W,u,w): out.add((u,v,w))
    return out

def dag_optimal(W,x,y):
    return adjust.optimal_adjustment_set(W,x,y)

def scan_cpdag(C):
    p = C.shape[0]; ref=v_structures(C)
    acc=defaultdict(int); wit=[]
    for G0 in reachable_mpdags(C):
        NA=nonadjacent_pairs(G0)
        if not NA: continue
        Q=[]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,n0,am0 = X.ostar_and_paths(G0,x,y)
                if am0: Q.append((x,y,O0))
        if not Q: continue
        ext=consistent_dag_extensions(G0,ref_vstructs=ref)
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info=bk_assert(G0,[(a,b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                P=G0.copy(); P[a,b]=1; P[b,a]=0
                ucP=ucolliders(P); ucH=ucolliders(H)
                for (x,y,O0) in Q:
                    O1,n1,am1 = X.ostar_and_paths(H,x,y)
                    if not am1: continue
                    for D in ext:
                        acc["trial_D"]+=1
                        Dp=D.copy(); Dp[a,b]=1
                        cyc=has_directed_cycle(Dp)
                        newuc = (not cyc) and (not (ucolliders(Dp) <= ucP))
                        pres = (not cyc) and preserves(Dp,H)
                        inext = (not cyc) and pres and (ucolliders(Dp)<=ucH)
                        if cyc: acc["cyc"]+=1
                        elif not pres: acc["nopres"]+=1
                        else: acc["bridge_ok"]+=1
                        if (not cyc) and newuc: acc["newuc"]+=1
                        if (not cyc) and (not newuc) and (not pres): acc["nopres_but_nonewuc"]+=1
                        if cyc and (not newuc): pass
                        if inext:
                            acc["inextH"]+=1
                            Od = dag_optimal(Dp,x,y)
                            if Od is None: acc["Dp_nonamen"]+=1
                            elif Od==O1: acc["E7_eq"]+=1
                            else:
                                acc["E7_ne"]+=1
                                if len(wit)<8: wit.append(dict(tag="E7",C=C.tolist(),G0=G0.tolist(),H=H.tolist(),D=D.tolist(),a=int(a),b=int(b),x=int(x),y=int(y),O1=sorted(map(int,O1)),Od=sorted(map(int,Od))))
                        # validity
                        vD = adjust.is_valid_adjustment_set(D,x,y,set(O1))
                        if not vD:
                            acc["VIOLATION"]+=1
                            if len(wit)<8: wit.append(dict(tag="VIOL",C=C.tolist(),G0=G0.tolist(),H=H.tolist(),D=D.tolist(),a=int(a),b=int(b),x=int(x),y=int(y),O1=sorted(map(int,O1))))
                        if not cyc:
                            vDp = adjust.is_valid_adjustment_set(Dp,x,y,set(O1))
                            acc["validDp" if vDp else "invalidDp"]+=1
                            if (not vDp) and inext:
                                acc["invalidDp_inext"]+=1
    return dict(acc),wit

def main():
    p=int(sys.argv[1])
    cp={}
    for D in all_dags(p):
        C=dag_to_cpdag(D); cp.setdefault(C.tobytes(),C)
    print(f"p={p}: {len(cp)} CPDAGs",flush=True)
    tot=defaultdict(int); W=[]
    for k in cp:
        acc,w=scan_cpdag(cp[k])
        for kk,vv in acc.items(): tot[kk]+=vv
        if len(W)<8: W.extend(w[:8-len(W)])
    print(json.dumps(dict(sorted(tot.items())),indent=1))
    json.dump(W,open(f"tb_wit_p{p}.json","w"),indent=1)
main()
