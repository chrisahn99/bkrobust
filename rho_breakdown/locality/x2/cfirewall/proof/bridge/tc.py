"""Characterise the residual (bridge-failing) cases."""
from common import *
from collections import defaultdict
import sys, json
from itertools import combinations

def ucolliders(W):
    p=W.shape[0]; out=set()
    for v in range(p):
        pars=[u for u in range(p) if is_directed(W,u,v)]
        for u,w in combinations(sorted(pars),2):
            if not adjacent(W,u,w): out.add((u,v,w))
    return out

def scan_cpdag(C):
    p=C.shape[0]; ref=v_structures(C)
    acc=defaultdict(int); wit=[]
    for G0 in reachable_mpdags(C):
        NA=nonadjacent_pairs(G0)
        if not NA: continue
        Q=[]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,n0,am0=X.ostar_and_paths(G0,x,y)
                if am0: Q.append((x,y,O0))
        if not Q: continue
        ext=consistent_dag_extensions(G0,ref_vstructs=ref)
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info=bk_assert(G0,[(a,b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                P=G0.copy(); P[a,b]=1; P[b,a]=0
                ucP=ucolliders(P)
                for (x,y,O0) in Q:
                    O1,n1,am1=X.ostar_and_paths(H,x,y)
                    if not am1: continue
                    for D in ext:
                        Dp=D.copy(); Dp[a,b]=1
                        cyc=has_directed_cycle(Dp)
                        ok = (not cyc) and (ucolliders(Dp)<=ucP)
                        moved = (O1!=O0)
                        tag = "bridge" if ok else ("cyc" if cyc else "newuc")
                        acc[tag]+=1
                        acc[tag+"_moved" if moved else tag+"_same"]+=1
                        if not ok:
                            # does H disagree with D on any edge?
                            dis = not preserves(D,H) # H's directed edges vs D (ignoring a->b which is absent in D)
                            dis2 = any(not is_directed(D,i,j) for (i,j) in directed_edges(H) if (i,j)!=(a,b))
                            acc[tag+"_disagree"] += int(dis2)
                            if moved and dis2:
                                acc[tag+"_moved_disagree"]+=1
                                if len(wit)<6:
                                    wit.append(dict(tag=tag,C=C.tolist(),G0=G0.tolist(),H=H.tolist(),
                                                    D=D.tolist(),a=int(a),b=int(b),x=int(x),y=int(y),
                                                    O0=sorted(map(int,O0)),O1=sorted(map(int,O1))))
    return dict(acc),wit

def main():
    p=int(sys.argv[1])
    cp={}
    for D in all_dags(p):
        C=dag_to_cpdag(D); cp.setdefault(C.tobytes(),C)
    tot=defaultdict(int); W=[]
    for k in cp:
        acc,w=scan_cpdag(cp[k])
        for kk,vv in acc.items(): tot[kk]+=vv
        if len(W)<6: W.extend(w[:6-len(W)])
    print(json.dumps(dict(sorted(tot.items())),indent=1))
    json.dump(W,open(f"tc_wit_p{p}.json","w"),indent=1)
main()
