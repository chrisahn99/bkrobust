# Test: for a DAG D and query (x,y), is every Z with O*_D <= Z <= V\forb_D valid?
import sys, itertools
from collections import defaultdict
import numpy as np
sys.path.insert(0,".")
from lib import *
acc=defaultdict(int); wit=[]
for p in (4,5):
    for D in all_dags(p):
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O=adjust.optimal_adjustment_set(D,x,y)
                if O is None: continue
                fb=adjust.forb(D,x,y)
                rest=[v for v in range(p) if v not in fb and v not in O]
                for r in range(len(rest)+1):
                    for extra in itertools.combinations(rest,r):
                        Z=set(O)|set(extra)
                        v=adjust.is_valid_adjustment_set(D,x,y,Z)
                        acc["n"]+=1; acc["bad"]+=int(not v)
                        if not v and len(wit)<5:
                            wit.append((p,D.tolist(),x,y,sorted(O),sorted(Z)))
print(dict(acc)); print(wit[:3])
