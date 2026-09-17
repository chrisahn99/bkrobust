# Collect EVERY residual (trial x D) instance at p=4 and characterise it.
import json
from collections import Counter, defaultdict
import numpy as np
from jlib import *

R=[]
for C in cpdags(4):
    ref=v_structures(C)
    for G0 in reachable_mpdags(C):
        NA=nonadjacent_pairs(G0)
        if not NA: continue
        Q=[]
        for x in range(4):
            for y in range(4):
                if x==y: continue
                O0,n0,am0,cn0,fb0=ostar(G0,x,y)
                if am0: Q.append((x,y,O0,cn0,fb0))
        if not Q: continue
        ext=consistent_dag_extensions(G0,ref_vstructs=ref)
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info=bk_assert(G0,[(a,b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                for (x,y,O0,cn0,fb0) in Q:
                    O1,n1,am1,cn1,fb1=ostar(H,x,y)
                    if not am1: continue
                    for D in ext:
                        ok,Dp,why=lifts(D,a,b,H)
                        cnD=causal_nodes(D,x,y); fbD=(poss_de(D,cnD)|{x}) if cnD else {x}
                        paD=parents_of_set(D,cnD)-fbD if cnD else set()
                        certB=(not (O1&fbD)) and bool(cnD) and (paD<=set(O1))
                        if ok or certB: continue
                        R.append(dict(G0=G0.tolist(),H=H.tolist(),D=D.tolist(),x=x,y=y,a=a,b=b,
                                      O1=sorted(O1),O0=sorted(O0),why=why,
                                      valid=bool(is_valid_adjustment_set(D,x,y,set(O1))),
                                      cnD=sorted(cnD),fbD=sorted(fbD),paD=sorted(paD),
                                      cnH=sorted(cn1),fbH=sorted(fb1),cn0=sorted(cn0),fb0=sorted(fb0),
                                      moved=bool(O1!=O0)))
print("residual instances:",len(R))
print("invalid among them:",sum(1 for r in R if not r["valid"]))
print("why:",Counter(r["why"] for r in R))
print("moved:",Counter(r["moved"] for r in R))
print("O1 vs O0:",Counter(("O1==O0" if r["O1"]==r["O0"] else "diff") for r in R))
# is O1 a superset/subset of O0?
c=Counter()
for r in R:
    s1,s0=set(r["O1"]),set(r["O0"])
    c["sub" if s1<s0 else "sup" if s1>s0 else "eq" if s1==s0 else "inc"]+=1
print("O1 rel O0:",c)
# how many distinct (G0,x,y,a,b) trials
print("distinct trials:",len({(str(r['G0']),r['x'],r['y'],r['a'],r['b']) for r in R}))
json.dump(R,open("residual_p4.json","w"))
for r in R[:6]:
    G0=np.array(r["G0"],dtype=np.int8); H=np.array(r["H"],dtype=np.int8); D=np.array(r["D"],dtype=np.int8)
    def sh(G):
        p=G.shape[0]; out=[]
        for i in range(p):
            for j in range(p):
                if i<j:
                    if is_undirected(G,i,j): out.append("%d--%d"%(i,j))
                    elif is_directed(G,i,j): out.append("%d->%d"%(i,j))
                    elif is_directed(G,j,i): out.append("%d<-%d"%(i,j))
        return " ".join(out)
    print("---- x=%d y=%d stmt %d->%d why=%s valid=%s"%(r["x"],r["y"],r["a"],r["b"],r["why"],r["valid"]))
    print("   G0:",sh(G0)); print("   H :",sh(H)); print("   D :",sh(D))
    print("   O0=%s O1=%s cnD=%s fbD=%s paD=%s cnH=%s fbH=%s"%(r["O0"],r["O1"],r["cnD"],r["fbD"],r["paD"],r["cnH"],r["fbH"]))
