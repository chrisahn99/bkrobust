"""Verify the three proof-carrying lemmas exhaustively at small p."""
import sys, itertools
from collections import defaultdict
import numpy as np
sys.path.insert(0,".")
from lib import *
from rd import all_simple_paths

def blocks(D,pth,Z):
    """standard d-sep blocking of a fixed path in a DAG"""
    k=len(pth)-1
    for i in range(1,k):
        u,v,w=pth[i-1],pth[i],pth[i+1]
        coll = is_directed(D,u,v) and is_directed(D,w,v)
        if coll:
            de=adjust.poss_de(D,{v})   # in a DAG poss_de = de
            if not (de & set(Z)): return True
        else:
            if v in Z: return True
    return False

def crit(D,x,y,Z):
    """Z n forb = {} and Z blocks every proper non-causal path x..y"""
    if set(Z) & adjust.forb(D,x,y): return False
    for pth in all_simple_paths(D,x,y):
        k=len(pth)-1
        if all(is_directed(D,pth[i],pth[i+1]) for i in range(k)): continue  # causal
        if not blocks(D,pth,Z): return False
    return True

acc=defaultdict(int); w1=[];w2=[]
for p in (3,4):
    for D in all_dags(p):
        for x in range(p):
            for y in range(p):
                if x==y: continue
                V=[v for v in range(p) if v not in (x,y)]
                for r in range(len(V)+1):
                    for Z in itertools.combinations(V,r):
                        Z=set(Z)
                        a=adjust.is_valid_adjustment_set(D,x,y,Z); b=crit(D,x,y,Z)
                        acc["cmp"]+=1
                        if a!=b:
                            acc["crit_mismatch"]+=1
                            if len(w1)<3: w1.append((p,D.tolist(),x,y,sorted(Z),a,b))
                        if a:
                            # L3: edge deletion monotonicity
                            for (i,j) in directed_edges(D):
                                Dm=D.copy(); Dm[i,j]=0
                                acc["del"]+=1
                                if not adjust.is_valid_adjustment_set(Dm,x,y,Z):
                                    acc["L3_fail"]+=1
                                    if len(w2)<3: w2.append((p,D.tolist(),x,y,sorted(Z),(i,j)))
print(dict(acc)); print("crit mismatches:",w1); print("L3 fails:",w2)
