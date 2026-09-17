"""Exhaustive scan of Conjecture A over EVERY reachable MPDAG of a few
maximal-chain-component CPDAGs: the complete graph (fully undirected CPDAG)
and the chordless-cycle CPDAG, at p = 6,7,8."""
import sys, json
sys.path.insert(0,"/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/code")
import numpy as np, x2lib as X
from graphs import *
from adjust import causal_nodes, poss_de, parents_of_set

def scan(C, name, cap=400000):
    p=C.shape[0]; ref=v_structures(C)
    seen={C.tobytes():C}; frontier=[C]
    while frontier:
        nxt=[]
        for G in frontier:
            for (u,v) in undirected_edges(G):
                for (a,b) in ((u,v),(v,u)):
                    H=G.copy(); H[b,a]=0; H=meek_closure(H)
                    if has_directed_cycle(H) or v_structures(H)!=ref: continue
                    kb=H.tobytes()
                    if kb not in seen:
                        seen[kb]=H; nxt.append(H)
        frontier=nxt
        if len(seen)>cap: return dict(name=name,p=p,ABORTED=len(seen))
    s=dict(name=name,p=p,n_mpdag=len(seen),D_all=0,N_cycle=0,N_vstruct=0,D_con=0,
           D_am=0,N_ident=0,forb_shrink=0,cn_shrink=0,N_set=0)
    viol=[]
    for G in seen.values():
        U=undirected_edges(G)
        if not U: continue
        succ=[]
        for (u,v) in U:
            for (a,b) in ((u,v),(v,u)):
                H=G.copy(); H[b,a]=0; H=meek_closure(H)
                bad=1 if has_directed_cycle(H) else (2 if v_structures(H)!=ref else 0)
                succ.append(((a,b),H,bad))
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,_,a0=X.ostar_and_paths(G,x,y)
                if not a0: continue
                cn0=causal_nodes(G,x,y); fb0=poss_de(G,cn0)|{x}
                for (ab,H,bad) in succ:
                    s["D_all"]+=1
                    if bad==1: s["N_cycle"]+=1; continue
                    if bad==2: s["N_vstruct"]+=1; continue
                    s["D_con"]+=1
                    O1,_,a1=X.ostar_and_paths(H,x,y)
                    if not a1: s["N_ident"]+=1; continue
                    s["D_am"]+=1
                    cn1=causal_nodes(H,x,y); fb1=poss_de(H,cn1)|{x}
                    if cn1!=cn0: s["cn_shrink"]+=1
                    if fb1!=fb0: s["forb_shrink"]+=1
                    if O1!=O0:
                        s["N_set"]+=1
                        if len(viol)<2: viol.append(dict(x=x,y=y,e=list(ab),G0=G.tolist(),H=H.tolist()))
    s["viol"]=viol
    return s

out=[]
for p in (5,6,7):
    K=np.ones((p,p),dtype=np.int8); np.fill_diagonal(K,0)
    out.append(scan(K,f"complete-K{p}"))
    print(json.dumps(out[-1],default=float)[:600],flush=True)
for p in (6,7,8):
    Cy=np.zeros((p,p),dtype=np.int8)
    for i in range(p): Cy[i,(i+1)%p]=Cy[(i+1)%p,i]=1
    # chordless cycle skeleton with NO v-structures -> fully undirected? not a CPDAG in general;
    # take the CPDAG of a DAG on that skeleton with a single sink
    D=np.zeros((p,p),dtype=np.int8)
    for i in range(p-1): D[i,i+1]=1
    D[0,p-1]=1
    out.append(scan(dag_to_cpdag(D),f"cycle-{p}"))
    print(json.dumps(out[-1],default=float)[:600],flush=True)
