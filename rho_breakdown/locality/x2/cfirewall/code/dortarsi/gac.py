"""Generalized adjustment criterion ON THE MPDAG: is it exactly equivalent to
'valid in every D in [G0]'?  Tested over ALL subsets Z, full census."""
import sys, itertools
from collections import defaultdict
import numpy as np
BASE = "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious"
sys.path.insert(0, BASE + "/x2/code"); sys.path.insert(0, BASE + "/x2/cfirewall/code")
from graphs import (consistent_dag_extensions, dag_to_cpdag, is_directed,
                    is_undirected, adjacent, v_structures)
import adjust
from run_lemma import all_dags, reachable_mpdags

def all_paths(G, x, y):
    """all simple paths x..y in skeleton(G)"""
    p=G.shape[0]; out=[]
    def dfs(path, vis):
        v=path[-1]
        if v==y: out.append(list(path)); return
        for w in range(p):
            if w in vis or not adjacent(G,v,w): continue
            vis.add(w); path.append(w); dfs(path,vis); path.pop(); vis.remove(w)
    dfs([x],{x}); return out

def node_status(G,u,w,v):
    if is_directed(G,u,w) and is_directed(G,v,w): return "col"
    if is_directed(G,w,u) or is_directed(G,w,v): return "dnc"
    if is_undirected(G,u,w) and is_undirected(G,w,v) and not adjacent(G,u,v): return "dnc"
    return "undef"

def is_causal(G,path):
    return all(is_directed(G,path[i],path[i+1]) for i in range(len(path)-1))

def gac(G,x,y,Z,variant):
    Z=set(Z)
    if Z & adjust.forb(G,x,y): return False
    for path in all_paths(G,x,y):
        if len(path)<2: continue
        st=[node_status(G,path[i-1],path[i],path[i+1]) for i in range(1,len(path)-1)]
        if any(s=="undef" for s in st): continue          # not definite status
        if is_causal(G,path): continue                    # causal
        blocked=False
        for i,s in enumerate(st,start=1):
            w=path[i]
            if s=="dnc" and w in Z: blocked=True; break
            if s=="col":
                dd = adjust.poss_de(G,{w}) if variant=="poss" else _de(G,w)
                if not (dd & Z): blocked=True; break
        if not blocked: return False
    return True

def _de(G,w):
    p=G.shape[0]; seen={w}; st=[w]
    while st:
        u=st.pop()
        for v in range(p):
            if is_directed(G,u,v) and v not in seen: seen.add(v); st.append(v)
    return seen

def truth(G0,ext,x,y,Z):
    return all(adjust.is_valid_adjustment_set(D,x,y,set(Z)) for D in ext)

def run(p):
    cp={}
    for D in all_dags(p): C=dag_to_cpdag(D); cp.setdefault(C.tobytes(),C)
    acc=defaultdict(int); wit=defaultdict(list)
    for C in cp.values():
        ref=v_structures(C)
        for G0 in reachable_mpdags(C):
            ext=consistent_dag_extensions(G0,ref_vstructs=ref)
            if not ext: continue
            for x in range(p):
                for y in range(p):
                    if x==y: continue
                    O0,n0,am0 = __import__("x2lib").ostar_and_paths(G0,x,y)
                    if not am0: continue
                    others=[v for v in range(p) if v not in (x,y)]
                    for r in range(len(others)+1):
                        for Z in itertools.combinations(others,r):
                            t=truth(G0,ext,x,y,Z)
                            for var in ("poss","dir"):
                                g=gac(G0,x,y,Z,var)
                                acc[(var,int(t),int(g))]+=1
                                if t!=g and len(wit[var])<4:
                                    wit[var].append((G0.tolist(),x,y,list(Z),t,g))
    return acc,wit

if __name__=="__main__":
    for p in (3,4):
        acc,wit=run(p)
        print(f"=== p={p}")
        for k in sorted(acc, key=str): print("   ",k,acc[k])
        for var in wit:
            for w in wit[var][:2]: print("   MISMATCH",var,w)
