"""DELETION LEMMA, exhaustive: Z valid for (X,Y) in DAG D  =>  Z valid in D minus any edge."""
import sys; sys.path.insert(0,'/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof')
from common import *
import itertools
p=int(sys.argv[1]); n=0; bad=0; wit=None
for D in all_dags(p):
    E=[(i,j) for i in range(p) for j in range(p) if D[i,j]==1]
    for x in range(p):
        for y in range(p):
            if x==y: continue
            others=[v for v in range(p) if v not in (x,y)]
            for r in range(len(others)+1):
                for Z in itertools.combinations(others,r):
                    if not adjust.is_valid_adjustment_set(D,x,y,set(Z)): continue
                    for (i,j) in E:
                        Dm=D.copy(); Dm[i,j]=0; n+=1
                        if not adjust.is_valid_adjustment_set(Dm,x,y,set(Z)):
                            bad+=1
                            if wit is None: wit=(D.tolist(),x,y,list(Z),(i,j))
print("DELETION p=%d tested=%d violations=%d"%(p,n,bad));  print(" WIT",wit) if wit else None
