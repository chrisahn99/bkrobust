"""Minimal explicit witnesses for the sharpness section."""
import sys
sys.path.insert(0,"/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/falsegain/code")
from fglib import *
import numpy as np, importlib.util
from graphs import directed_edges, undirected_edges as UE, dag_agrees_with
from run_lemma import all_dags
spec=importlib.util.spec_from_file_location("x1_ops","/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x1/code/x1_ops.py")
x1_ops=importlib.util.module_from_spec(spec); spec.loader.exec_module(x1_ops)
from scm import make_linear_iscm
from adjust import cov_linear, total_effect_linear, ols_coefficient

def show(tag,D,G0,H,x,y,st,s0,O0,s1,O1):
    rng=np.random.default_rng(3); bias=[]
    if s1==POS:
        for _ in range(3):
            sc=make_linear_iscm(D,rng); Sig=cov_linear(sc["A"],sc["omega"])
            bias.append((total_effect_linear(sc["A"],x,y), ols_coefficient(Sig,x,y,sorted(O1))))
    print(f"--- {tag}: p={D.shape[0]} x={x} y={y} statement={st[0]}->{st[1]}")
    print(f"    D  edges {directed_edges(D)}")
    print(f"    G0 dir {directed_edges(G0)}  und {UE(G0)}   D in [G0]? {dag_agrees_with(D,G0)}")
    print(f"    H  dir {directed_edges(H)}  und {UE(H)}   extendable={x1_ops.pdag_extendable(H)}")
    print(f"    R(G0)={s0} O*={sorted(O0) if O0 else O0}  ->  R(H)={s1} O*={sorted(O1) if O1 else O1}")
    print(f"    (tau,beta) over 3 SCMs: {[(round(a,4),round(b,4)) for a,b in bias]}")

found={}
for p in (4,5):
    for D in all_dags(p):
        if len(found)>=3: break
        C=dag_to_cpdag(D)
        for G0 in reachable_mpdags(C):
            if not dag_agrees_with(D,G0): continue
            NA=x1_ops.nonadjacent_pairs(C)
            for (a0,b0) in NA:
                for (a,b) in ((a0,b0),(b0,a0)):
                    H,info=x1_ops.bk_assert(G0,[(a,b)])
                    ok=(not info["conflict"]) and (not has_directed_cycle(H)) and x1_ops.pdag_extendable(H)
                    if ok: continue
                    for x in range(p):
                        for y in range(p):
                            if x==y: continue
                            s0,O0=report_state(G0,x,y)
                            if s0==REFUSE: continue
                            if not sound(D,x,y,s0,O0): continue
                            s1,O1=report_state(H,x,y)
                            if sound(D,x,y,s1,O1) is not False: continue
                            key=("POS_corrupt" if s0==POS else "ZERO_to_"+s1)
                            if key not in found:
                                found[key]=1; show(key,D,G0,H,x,y,(a,b),s0,O0,s1,O1)
    if len(found)>=3: break
