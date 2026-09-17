"""Test candidate invariants:  cn_D <= cn_H ;  forb_D <= forb_H ;  O_H cap forb_D = 0 ;
   O_H cap forb_G0 = 0 (F1) ; blocking clause separately."""
from common import *
from collections import defaultdict
import sys, json, networkx as nx
from itertools import combinations

def cnset(G,x,y):
    return adjust.causal_nodes(G,x,y)
def forbset(G,x,y):
    return adjust.forb(G,x,y)

def blocks_noncausal(D,x,y,Z):
    """clause (2) only: Z d-separates x,y in proper back-door graph of D"""
    cn=adjust.causal_nodes(D,x,y)
    Dp=D.copy()
    for w in cn:
        if is_directed(Dp,x,w): Dp[x,w]=0
    g=nx.DiGraph(); g.add_nodes_from(range(D.shape[0]))
    for i in range(D.shape[0]):
        for j in range(D.shape[0]):
            if Dp[i,j]==1: g.add_edge(i,j)
    return nx.is_d_separator(g,{x},{y},set(Z))

def ucolliders(W):
    p=W.shape[0]; out=set()
    for v in range(p):
        pars=[u for u in range(p) if is_directed(W,u,v)]
        for u,w in combinations(sorted(pars),2):
            if not adjacent(W,u,w): out.add((u,v,w))
    return out

def scan_cpdag(C):
    p=C.shape[0]; ref=v_structures(C)
    acc=defaultdict(int)
    for G0 in reachable_mpdags(C):
        NA=nonadjacent_pairs(G0)
        if not NA: continue
        Q=[]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,n0,am0=X.ostar_and_paths(G0,x,y)
                if am0: Q.append((x,y,O0,cnset(G0,x,y),forbset(G0,x,y)))
        if not Q: continue
        ext=consistent_dag_extensions(G0,ref_vstructs=ref)
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info=bk_assert(G0,[(a,b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                P=G0.copy(); P[a,b]=1; P[b,a]=0
                ucP=ucolliders(P)
                for (x,y,O0,cn0,fb0) in Q:
                    O1,n1,am1=X.ostar_and_paths(H,x,y)
                    if not am1: continue
                    cnH=cnset(H,x,y); fbH=forbset(H,x,y)
                    if O1 & fb0: acc["F1_FAIL"]+=1
                    for D in ext:
                        Dp=D.copy(); Dp[a,b]=1
                        cyc=has_directed_cycle(Dp)
                        ok=(not cyc) and (ucolliders(Dp)<=ucP)
                        tag="bridge" if ok else ("cyc" if cyc else "newuc")
                        cnD=cnset(D,x,y); fbD=forbset(D,x,y)
                        acc[tag]+=1
                        if not (cnD<=cnH): acc[tag+"_cnD_notsub"]+=1
                        if not (fbD<=fbH): acc[tag+"_fbD_notsub"]+=1
                        if not (fbD<=fb0):  acc[tag+"_fbD_notsub_G0"]+=1
                        if O1 & fbD: acc[tag+"_clause1_FAIL"]+=1
                        if not blocks_noncausal(D,x,y,O1): acc[tag+"_clause2_FAIL"]+=1
                        if not cnD: acc[tag+"_cnD_EMPTY"]+=1
    return dict(acc)

def main():
    p=int(sys.argv[1]); cp={}
    for D in all_dags(p):
        C=dag_to_cpdag(D); cp.setdefault(C.tobytes(),C)
    tot=defaultdict(int)
    for k in cp:
        for kk,vv in scan_cpdag(cp[k]).items(): tot[kk]+=vv
    print(json.dumps(dict(sorted(tot.items())),indent=1))
main()
