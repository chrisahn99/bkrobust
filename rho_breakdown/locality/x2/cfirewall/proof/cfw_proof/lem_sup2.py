"""LEM-SUP restricted: cn_D nonempty."""
import sys; sys.path.insert(0,'.')
from common import *
import itertools
p=int(sys.argv[1])
n=0; bad=0; wit=[]
for D in all_dags(p):
    for x in range(p):
        for y in range(p):
            if x==y: continue
            C_=cn(D,x,y)
            if not C_: continue
            F = forb(D,x,y)
            O = pa_set(D,C_)-F
            if not valid(D,x,y,O): print("O* INVALID!",D.tolist(),x,y); continue
            rest=[v for v in range(p) if v not in F and v!=y]
            for r in range(len(rest)+1):
                for extra in itertools.combinations(rest,r):
                    Z=set(O)|set(extra)
                    if (Z & F) or y in Z: continue
                    n+=1
                    if not valid(D,x,y,Z):
                        bad+=1
                        if len(wit)<3: wit.append((D.tolist(),x,y,sorted(O),sorted(Z)))
print("p",p,"tested",n,"violations",bad)
for w in wit: print("W",w)
