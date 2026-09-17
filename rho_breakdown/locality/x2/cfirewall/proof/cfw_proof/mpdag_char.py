"""Over ALL PDAGs on p nodes:  Meek-closed + acyclic + extendable  ==>  maximally oriented?
   and  (+ amenable rel (x,y))  ==>  O*(x,y,H) valid in every DAG in [H]?   (= HPM)"""
import sys; sys.path.insert(0,'.')
from common import *
from itertools import product, combinations
import numpy as np
p=int(sys.argv[1])
idx=[(i,j) for i in range(p) for j in range(i+1,p)]
n=0; notmax=0; hpmbad=0; w1=None; w2=None
for state in product([0,1,2,3],repeat=len(idx)):   # 0 none 1 i->j 2 j->i 3 undirected
    H=np.zeros((p,p),dtype=np.int8)
    for (i,j),s in zip(idx,state):
        if s==1: H[i,j]=1
        elif s==2: H[j,i]=1
        elif s==3: H[i,j]=H[j,i]=1
    if has_directed_cycle(H): continue
    if not np.array_equal(meek_closure(H),H): continue
    if not pdag_extendable(H): continue
    ext=consistent_dag_extensions(H,ref_vstructs=v_structures(H))
    if not ext: continue
    n+=1
    cog=common_orientation_graph(ext,skeleton(H))
    if not np.array_equal(cog,H):
        notmax+=1
        if w1 is None: w1=(H.tolist(),cog.tolist())
    for x in range(p):
        for y in range(p):
            if x==y: continue
            O=adjust.optimal_adjustment_set(H,x,y)
            if O is None: continue
            for D in ext:
                if not adjust.is_valid_adjustment_set(D,x,y,set(O)):
                    hpmbad+=1
                    if w2 is None: w2=(H.tolist(),D.tolist(),x,y,sorted(O))
print("p=%d  Meek-closed+acyclic+extendable PDAGs=%d   NOT maximally oriented=%d   HPM failures=%d"%(p,n,notmax,hpmbad))
if w1: print(" W1",w1)
if w2: print(" W2",w2)
