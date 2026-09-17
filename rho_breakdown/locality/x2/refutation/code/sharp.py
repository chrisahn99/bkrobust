import json, glob, itertools, collections, math
BASE="/home/costaj/latent-causal/e1prime-se/results/"
CELLS=[("original","K4"),("licensed","K4"),("large","K4"),("k8","K4"),("k8","K6"),("k8","K8")]
INF=1<<20
def cp_upper(n):  # one-sided 95% upper for 0/n
    return 1-0.05**(1.0/n) if n>0 else float('nan')
out={}
glob_tot=collections.defaultdict(lambda: collections.Counter())
for ens,arm in CELLS:
    recs=json.load(open(BASE+f"arm1_{ens}.json"))
    st=collections.defaultdict(lambda: collections.Counter())
    per_scm=[]  # (has_near_far_trial, near N/D, far N/D)
    for r in recs:
        a=r["arms"].get(arm)
        if not a or not a.get("mpdag_amenable"): continue
        hd=a["stmt_hopdist"]; k=a["n_K"]; top=min(r["max_rho"],k)
        idx=0
        loc=collections.defaultdict(lambda: collections.Counter())
        for rho in range(1,top+1):
            for flip in itertools.combinations(range(k),rho):
                m=a["members"][idx]; idx+=1
                ds=[hd[i] for i in flip]
                if any(d==0 for d in ds): continue          # not "all far"
                if any(d>=INF for d in ds): continue        # drop disconnected entirely
                key=min(ds)
                kk=1 if key==1 else 2  # 1 vs >=2
                loc[kk]["D_all"]+=1; st[key]["D_all"]+=1
                if not m.get("consistent"): continue
                loc[kk]["D_con"]+=1; st[key]["D_con"]+=1
                if not m.get("amenable"): continue
                loc[kk]["D_am"]+=1; st[key]["D_am"]+=1
                if m.get("ostar_changed"):
                    loc[kk]["N_set"]+=1; st[key]["N_set"]+=1
        if loc[1]["D_am"] and loc[2]["D_am"]:
            per_scm.append((loc[1]["N_set"],loc[1]["D_am"],loc[2]["N_set"],loc[2]["D_am"]))
    near=st[1]; far=sum((st[k] for k in st if k>=2), collections.Counter())
    rn=near["N_set"]/near["D_am"] if near["D_am"] else 0
    exp=rn*far["D_am"]
    out[f"{ens}/{arm}"]=dict(
        near=dict(near), far=dict(far), near_rate=rn,
        far_D_am=far["D_am"], expected_if_same_rate=exp,
        p_zero_given_same_rate=math.exp(-exp) if exp<700 else 0.0,
        by_dmin={str(k):dict(v) for k,v in sorted(st.items())},
        matched_n_scm=len(per_scm),
        matched_near=(sum(a for a,_,_,_ in per_scm), sum(b for _,b,_,_ in per_scm)),
        matched_far=(sum(c for _,_,c,_ in per_scm), sum(d for _,_,_,d in per_scm)),
    )
    for k,v in st.items(): glob_tot[1 if k==1 else 2].update(v)
pool_near=glob_tot[1]; pool_far=glob_tot[2]
rn=pool_near["N_set"]/pool_near["D_am"]
out["POOLED"]=dict(near=dict(pool_near),far=dict(pool_far),near_rate=rn,
                   expected_if_same_rate=rn*pool_far["D_am"],
                   p_zero=math.exp(-rn*pool_far["D_am"]),
                   rule_of_three_upper=3/pool_far["D_am"],
                   cp_upper=cp_upper(pool_far["D_am"]))
print(json.dumps(out,indent=1))
