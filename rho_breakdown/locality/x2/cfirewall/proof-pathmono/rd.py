import sys
import numpy as np
sys.path.insert(0,"/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof-pathmono")
from lib import *

def all_simple_paths(G,x,y):
    p=G.shape[0]; S=skeleton(G); out=[]
    def dfs(path,vis):
        v=path[-1]
        if v==y: out.append(list(path)); return
        for w in np.flatnonzero(S[v]):
            w=int(w)
            if w in vis: continue
            vis.add(w); path.append(w); dfs(path,vis); path.pop(); vis.remove(w)
    dfs([x],{x}); return out

def required_parents(D,x,y):
    """R_D: nodes v_j whose membership in Z the blocking argument needs."""
    fb=adjust.forb(D,x,y); cn=adjust.causal_nodes(D,x,y)
    R=set()
    for pth in all_simple_paths(D,x,y):
        k=len(pth)-1
        causal=all(is_directed(D,pth[i],pth[i+1]) for i in range(k))
        if causal: continue
        idx=[i for i in range(k+1) if pth[i] not in fb]
        if not idx: continue
        j=idx[-1]
        # by the lemma the edge must be v_j -> v_{j+1}
        assert is_directed(D,pth[j],pth[j+1]), "orientation lemma broken"
        if all(is_directed(D,pth[i],pth[i+1]) for i in range(j+1,k)):
            R.add(pth[j])
    return R,fb,cn
