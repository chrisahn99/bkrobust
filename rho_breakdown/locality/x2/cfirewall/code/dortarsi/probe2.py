"""Contingency table: what separates valid from invalid, treatment vs control."""
import sys, json
from collections import defaultdict
import numpy as np
BASE = "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious"
sys.path.insert(0, BASE + "/x2/code"); sys.path.insert(0, BASE + "/x2/cfirewall/code")
from graphs import (consistent_dag_extensions, dag_to_cpdag, dag_agrees_with,
                    has_directed_cycle, meek_closure, random_dag, skeleton,
                    undirected_edges, v_structures, directed_edges, is_directed)
import adjust, x2lib as X
from run_lemma import all_dags, reachable_mpdags
import cfirewall as CF
x1 = CF.x1_ops

def ancestors(D, v):
    p=D.shape[0]; seen={v}; st=[v]
    while st:
        u=st.pop()
        for w in range(p):
            if D[w,u]==1 and D[u,w]==0 and w not in seen: seen.add(w); st.append(w)
    return seen

def poss_desc(G,v): return adjust.poss_de(G,{v})

def partially_directed_cycle(G):
    """Is there a cycle v1 ~> v2 ~> ... ~> v1 with each step directed-or-undirected
    and at least one directed?  = a directed edge u->v with a possibly-directed
    path v ~> u."""
    p=G.shape[0]
    for u in range(p):
        for v in range(p):
            if is_directed(G,u,v) and u in adjust.poss_de(G,{v}): return True
    return False

def scan_cpdag(C):
    ref=v_structures(C); acc=defaultdict(int); wit=defaultdict(list)
    for G0 in reachable_mpdags(C):
        p=G0.shape[0]
        NA=x1.nonadjacent_pairs(G0)
        if not NA: continue
        Q=[]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,n0,am0=X.ostar_and_paths(G0,x,y)
                if am0: Q.append((x,y,O0,adjust.forb(G0,x,y),adjust.causal_nodes(G0,x,y)))
        if not Q: continue
        ext=consistent_dag_extensions(G0,ref_vstructs=ref)
        if not ext: continue
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info=x1.bk_assert(G0,[(a,b)])
                cyc=has_directed_cycle(H)
                coh=(not info["conflict"]) and (not cyc) and x1.pdag_extendable(H)
                arm="T" if coh else "K"
                acc[arm+"_pdc"] += partially_directed_cycle(H)
                # does G0 have a possibly-directed path b ~> a ?
                pdba = a in adjust.poss_de(G0,{b})
                vH=v_structures(H)
                for (x,y,O0,F0,CN0) in Q:
                    Z,n1,am1=X.ostar_and_paths(H,x,y)
                    if not am1: acc[arm+"_abort"]+=1; continue
                    moved=(Z!=O0)
                    for D in ext:
                        ok=adjust.is_valid_adjustment_set(D,x,y,set(Z))
                        acyc = b not in ancestors(D,a)
                        key=(arm,int(moved),int(acyc),int(ok))
                        acc[key]+=1
                        if not ok:
                            # diagnostics on the failure
                            f1 = bool(set(Z) & adjust.forb(D,x,y))
                            acc[(arm,"bad","clause1" if f1 else "clause2")]+=1
                            acc[(arm,"bad","hitsforbG0")] += bool(set(Z)&F0)
                            acc[(arm,"bad","pdba")] += pdba
                            if len(wit[arm])<6:
                                wit[arm].append(dict(G0=G0.tolist(),H=H.tolist(),a=int(a),b=int(b),
                                    x=int(x),y=int(y),Z=sorted(map(int,Z)),O0=sorted(map(int,O0)),
                                    D=D.tolist(),clause1=f1,cyc=bool(cyc)))
                        if acyc:
                            Dp=D.copy(); Dp[a,b]=1; Dp[b,a]=0
                            inH=dag_agrees_with(Dp,H) and v_structures(Dp)==vH
                            okp=adjust.is_valid_adjustment_set(Dp,x,y,set(Z))
                            acc[(arm,"Dp",int(inH),int(okp),int(ok))]+=1
                        else:
                            Dm=D.copy(); Dm[b,a]=1; Dm[a,b]=0
                            okm=adjust.is_valid_adjustment_set(Dm,x,y,set(Z))
                            acc[(arm,"Dm",int(okm),int(ok))]+=1
                            acc[(arm,"cycbranch_pdba")]+=pdba
    return {str(k):v for k,v in acc.items()}, {k:v for k,v in wit.items()}

def _job(arg):
    kb,p=arg
    C=np.frombuffer(kb,dtype=np.int8).reshape(p,p).copy()
    return scan_cpdag(C)

if __name__=="__main__":
    p=int(sys.argv[1]); nw=int(sys.argv[2]) if len(sys.argv)>2 else 8
    cp={}
    for D in all_dags(p):
        C=dag_to_cpdag(D); cp.setdefault(C.tobytes(),C)
    from multiprocessing import Pool
    tot=defaultdict(int); W=defaultdict(list)
    with Pool(nw) as pool:
        for acc,w in pool.imap_unordered(_job,[(k,p) for k in cp],chunksize=4):
            for k,v in acc.items(): tot[k]+=v
            for k,v in w.items():
                if len(W[k])<6: W[k].extend(v[:6-len(W[k])])
    for k in sorted(tot): print(f"{k:60s} {tot[k]}")
    json.dump(dict(counts={k:int(v) for k,v in tot.items()},wit={k:v for k,v in W.items()}),
              open(f"/private/tmp/claude-501/-Users-josecosta-mugango/bceac0a4-93a7-4757-b3b3-a3d4e6c42d55/scratchpad/dt/probe2_p{p}.json","w"))
