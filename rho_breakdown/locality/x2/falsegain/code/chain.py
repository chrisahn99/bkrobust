"""WITNESS FAMILY W_k: false gain at hop distance k, for every k>=0, on p=k+3 nodes.
Skeleton = induced path  c_0 - c_1 - ... - c_k - x - y.
True DAG D : source at y  ->  y->x->c_k->...->c_1->c_0   (so tau_D(x,y)=0).
Statement  : c_0 -> c_1   (FALSE in D; Meek-consistent with C=cpdag(D)).
"""
import sys
sys.path.insert(0,"/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/falsegain/code")
from fglib import *
import numpy as np
from graphs import apply_background_knowledge, MeekFail, directed_edges
from scm import make_linear_iscm
from adjust import cov_linear, total_effect_linear, ols_coefficient

def build(k):
    p = k+3
    # nodes: c_0..c_k = 0..k ; x = k+1 ; y = k+2
    D = np.zeros((p,p),dtype=np.int8)
    x, y = k+1, k+2
    D[y,x]=1
    D[x,k]=1
    for i in range(k,0,-1): D[i,i-1]=1
    return D, x, y

for k in range(0,9):
    D,x,y = build(k); p=k+3
    C = dag_to_cpdag(D)
    nund = len(undirected_edges(C))
    s = (0,1)                       # c_0 -> c_1 ; TRUE D has 1->0
    false_stmt = (D[0,1]!=1)
    try:
        H = apply_background_knowledge(C, [s]); consistent=True
    except MeekFail:
        H=None; consistent=False
    dist = X.hop_dist_inf(C,[x,y]); dmin,_ = X.stmt_dist(dist,s)
    s0,O0 = report_state(C,x,y)
    s1,O1 = (report_state(H,x,y) if H is not None else (None,None))
    snd1 = sound(D,x,y,s1,O1) if H is not None else None
    # numeric
    rng=np.random.default_rng(100+k); bias=[]
    if s1==POS:
        for _ in range(5):
            sc=make_linear_iscm(D,rng); Sig=cov_linear(sc["A"],sc["omega"])
            tau=total_effect_linear(sc["A"],x,y)
            b=ols_coefficient(Sig,x,y,sorted(O1))
            bias.append(abs(b-tau))
    print(f"k={k} p={p} und(C)={nund} consistent={consistent} dmin={dmin} false={false_stmt} "
          f"R(C)={s0} -> R(H)={s1} O*={sorted(O1) if O1 else O1} sound={snd1} "
          f"H_dir={directed_edges(H) if H is not None else None} maxbias={max(bias) if bias else None}")
