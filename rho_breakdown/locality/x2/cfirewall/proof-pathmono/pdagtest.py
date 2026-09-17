"""Is every acyclic Meek-closed PDAG extendable?  Exhaustive p=3,4; sampled p=5,6."""
import sys, itertools
import numpy as np
sys.path.insert(0,".")
from lib import *
def check(G):
    if not np.array_equal(meek_closure(G),G): return None
    if has_directed_cycle(G): return None
    return pdag_extendable(G)
bad=[]; n=0
for p in (3,4):
    idx=[(i,j) for i in range(p) for j in range(i+1,p)]
    for st in itertools.product([0,1,2,3],repeat=len(idx)):
        G=np.zeros((p,p),dtype=np.int8)
        for (i,j),s in zip(idx,st):
            if s==1: G[i,j]=1
            elif s==2: G[j,i]=1
            elif s==3: G[i,j]=G[j,i]=1
        r=check(G)
        if r is None: continue
        n+=1
        if not r: bad.append((p,directed_edges(G),undirected_edges(G)))
print("closed+acyclic:",n," non-extendable:",len(bad)); print(bad[:3])
rng=np.random.default_rng(3)
n2=0;b2=0;w=[]
for _ in range(300000):
    p=int(rng.integers(5,7))
    idx=[(i,j) for i in range(p) for j in range(i+1,p)]
    G=np.zeros((p,p),dtype=np.int8)
    for (i,j) in idx:
        s=int(rng.integers(0,4))
        if s==1: G[i,j]=1
        elif s==2: G[j,i]=1
        elif s==3: G[i,j]=G[j,i]=1
    r=check(G)
    if r is None: continue
    n2+=1
    if not r:
        b2+=1
        if len(w)<2: w.append((p,directed_edges(G),undirected_edges(G)))
print("p5-6 sample closed+acyclic:",n2,"non-extendable:",b2); print(w)
