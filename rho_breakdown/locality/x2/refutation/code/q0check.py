"""Is S8's Q_0 refutation the SAME events as S5's rho>=2 refutation?
Derivation: est0 is in both balls; E_est(inf) == E_set (X2 sec 8, 0 discrepancies);
so the pruned (near-only) interval differs from the full interval ONLY IF some
member with >=1 flip at dmin>=1 changed O*.  Check numerically."""
import json, itertools, collections
BASE="/home/costaj/latent-causal/e1prime-se/results/"
CELLS=[("original","K4"),("licensed","K4"),("large","K4"),("k8","K4"),("k8","K6"),("k8","K8")]
INF=1<<20
for ens,arm in CELLS:
    recs=json.load(open(BASE+f"arm1_{ens}.json"))
    n=0; q0=0; anyfar=0; both=0
    for r in recs:
        a=r["arms"].get(arm)
        if not a or not a.get("mpdag_amenable"): continue
        hd=a["stmt_hopdist"]; k=a["n_K"]; top=min(r["max_rho"],k)
        ests_full=[a["est0"]]; ests_near=[a["est0"]]; far_changed=False
        idx=0
        for rho in range(1,top+1):
            for flip in itertools.combinations(range(k),rho):
                m=a["members"][idx]; idx+=1
                if not (m.get("consistent") and m.get("amenable")): continue
                e=m["est"]; ests_full.append(e)
                allnear = all(hd[i]==0 for i in flip)
                if allnear: ests_near.append(e)
                else:
                    if m.get("ostar_changed"): far_changed=True
        n+=1
        d = (abs(min(ests_full)-min(ests_near))>1e-12) or (abs(max(ests_full)-max(ests_near))>1e-12)
        q0+=d; anyfar+=far_changed; both+= (d and far_changed)
    print(f"{ens}/{arm:3s} n={n:5d}  Q0(interval differs)={q0/n:.4f}  "
          f"P(some far-flip member changed O*)={anyfar/n:.4f}  "
          f"P(both)={both/n:.4f}  identical={q0==anyfar==both}")
