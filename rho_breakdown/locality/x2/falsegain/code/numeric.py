"""NUMERIC ARM: does the graphical definition of (un)soundness bind numerically?
For random (D,G0,s) with a class-R or class-S false statement, draw M linear
iSCMs on D and compare beta(y~x+O*(H)) against tau_D(x,y).
SOUND  -> expect |beta-tau| = 0 to machine precision, ALWAYS.
UNSOUND-> expect |beta-tau| > 0 for a.e. draw (genericity)."""
import sys, json
from collections import defaultdict
sys.path.insert(0,"/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/falsegain/code")
from fglib import *
import fglib, numpy as np, importlib.util
from graphs import random_dag, is_undirected
from scm import make_linear_iscm
from adjust import cov_linear, total_effect_linear, ols_coefficient
spec = importlib.util.spec_from_file_location("x1_ops","/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x1/code/x1_ops.py")
x1_ops = importlib.util.module_from_spec(spec); spec.loader.exec_module(x1_ops)

def true_bk_mpdag(C,D,rng):
    G=C.copy(); U=undirected_edges(G)
    if not U: return G
    for t in rng.permutation(len(U))[:int(rng.integers(0,len(U)+1))]:
        u,v=U[t]
        if not is_undirected(G,u,v): continue
        if D[u,v]==1: G[v,u]=0
        else: G[u,v]=0
        G=meek_closure(G)
    return G

def main(nseed, M, seed0):
    rng=np.random.default_rng(seed0)
    stat=defaultdict(lambda: [0,0,0.0])   # [n, n_biased, max|d|]
    for it in range(nseed):
        p=int(rng.choice([5,6,7])); deg=float(rng.choice([1.5,2.0,2.5]))
        D=random_dag(p,deg,rng); C=dag_to_cpdag(D); ref=v_structures(C)
        G0=true_bk_mpdag(C,D,rng)
        scms=[]
        for _ in range(M):
            s=make_linear_iscm(D,rng); scms.append((cov_linear(s["A"],s["omega"]), s["A"]))
        U0=undirected_edges(G0); NA=x1_ops.nonadjacent_pairs(C)
        cand=[]
        for (u,v) in U0:
            for (a,b) in ((u,v),(v,u)):
                if D[a,b]==1: continue
                H=G0.copy(); H[b,a]=0; H=meek_closure(H)
                if has_directed_cycle(H) or v_structures(H)!=ref: continue
                cand.append(("R",H,(a,b)))
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info=x1_ops.bk_assert(G0,[(a,b)])
                ok=(not info["conflict"]) and (not has_directed_cycle(H)) and x1_ops.pdag_extendable(H)
                cand.append(("S1" if ok else "S1rej",H,(a,b)))
        for x in range(p):
            for y in range(p):
                if x==y: continue
                s0,O0=report_state(G0,x,y)
                for (tag,H,st) in cand:
                    s1,O1=report_state(H,x,y)
                    if s1!=POS: continue
                    snd=sound(D,x,y,s1,O1)
                    cell=(tag, "gain" if s0==REFUSE else ("fab" if s0==ZERO else "keep"),
                          "sound" if snd else "unsound")
                    Z=sorted(O1)
                    for (Sig,A) in scms:
                        tau=total_effect_linear(A,x,y)
                        try: beta=ols_coefficient(Sig,x,y,Z)
                        except np.linalg.LinAlgError: continue
                        d=abs(beta-tau)
                        e=stat[cell]; e[0]+=1; e[1]+= int(d>1e-9); e[2]=max(e[2],d)
        if len(fglib._valid_cache)>200000: fglib._valid_cache.clear(); fglib._cn_cache.clear()
        if it%50==0: print("  it",it,flush=True)
    return stat

if __name__=="__main__":
    st=main(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))
    for k,v in sorted(st.items()):
        print("|".join(k), f"n={v[0]} n_biased={v[1]} rate={v[1]/max(v[0],1):.6f} max|b-t|={v[2]:.4g}")
