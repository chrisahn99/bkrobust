"""LEM-SUP: in a DAG D, does  Z cap forb_D = {} and Z supseteq O*_D  imply Z valid?"""
import sys; sys.path.insert(0,'.')
from common import *
import itertools, numpy as np
p=int(sys.argv[1]) if len(sys.argv)>1 else 4
n=0; bad=0; wit=None
for D in all_dags(p):
    for x in range(p):
        for y in range(p):
            if x==y: continue
            F = forb(D,x,y); C_=cn(D,x,y)
            O = pa_set(D,C_)-F
            if not valid(D,x,y,O):      # no valid set exists / O* itself invalid
                continue
            rest=[v for v in range(p) if v not in F and v!=y]
            for r in range(len(rest)+1):
                for extra in itertools.combinations(rest,r):
                    Z=set(O)|set(extra)
                    if Z & F: continue
                    if y in Z: continue
                    n+=1
                    if not valid(D,x,y,Z):
                        bad+=1
                        if wit is None: wit=(D.tolist(),x,y,sorted(O),sorted(Z))
print("p",p,"tested",n,"violations",bad)
if wit: print("WITNESS D,x,y,O*,Z =",wit)
