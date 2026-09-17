import sys
sys.argv=[sys.argv[0]]
exec(open("numeric.py").read().split('if __name__')[0])
import numpy as np
from collections import defaultdict
def run(nseed,M,seed0):
    rng=np.random.default_rng(seed0); vals=defaultdict(list)
    for it in range(nseed):
        p=int(rng.choice([5,6,7,8])); deg=float(rng.choice([1.5,2.0,2.5]))
        D=random_dag(p,deg,rng); C=dag_to_cpdag(D); ref=v_structures(C); G0=true_bk_mpdag(C,D,rng)
        scms=[]
        for _ in range(M):
            s=make_linear_iscm(D,rng); scms.append((cov_linear(s["A"],s["omega"]),s["A"]))
        cand=[]
        for (u,v) in undirected_edges(G0):
            for (a,b) in ((u,v),(v,u)):
                if D[a,b]==1: continue
                H=G0.copy(); H[b,a]=0; H=meek_closure(H)
                if has_directed_cycle(H) or v_structures(H)!=ref: continue
                cand.append(("R",H))
        for (a0,b0) in x1_ops.nonadjacent_pairs(C):
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info=x1_ops.bk_assert(G0,[(a,b)])
                if (not info["conflict"]) and (not has_directed_cycle(H)) and x1_ops.pdag_extendable(H):
                    cand.append(("S1",H))
        for x in range(p):
            for y in range(p):
                if x==y: continue
                s0,O0=report_state(G0,x,y)
                if s0!=REFUSE: continue
                for tag,H in cand:
                    s1,O1=report_state(H,x,y)
                    if s1!=POS or sound(D,x,y,s1,O1): continue
                    Z=sorted(O1)
                    for Sig,A in scms:
                        tau=total_effect_linear(A,x,y)
                        try: b=ols_coefficient(Sig,x,y,Z)
                        except Exception: continue
                        vals[tag].append((abs(b-tau), abs(tau), abs(b)))
        if len(fglib._valid_cache)>150000: fglib._valid_cache.clear(); fglib._cn_cache.clear()
    return vals
v=run(150,4,505)
for k,L in v.items():
    a=np.array(L)
    print(k,"n=",len(a))
    for nm,col in [("|beta-tau|",0),("|tau|",1),("|beta|",2)]:
        q=np.quantile(a[:,col],[0.05,0.25,0.5,0.75,0.95])
        print("   ",nm,"q05/25/50/75/95 =", " ".join(f"{z:.3f}" for z in q))
    print("    P(tau==0 but beta!=0) =", float(np.mean((a[:,1]<1e-12)&(a[:,2]>1e-9))))
    print("    P(sign flip vs tau)   =", "n/a (abs only)")
    print("    P(|beta-tau| > |tau|) =", float(np.mean(a[:,0]>a[:,1])))
