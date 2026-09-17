"""The decomposition the report never runs: RR(silent) = RR(O* moved) x RR(invalid | moved)."""
import json, math, numpy as np
UNREACH=1<<20
def wil(k,n):
    if n==0: return (float('nan'),)*3
    z=1.959963984540054;p=k/n;d=1+z*z/n
    c=(p+z*z/(2*n))/d;h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return p,max(0,c-h),min(1,c+h)
def f(k,n):
    p,lo,hi=wil(k,n);return f"{p:.4f} [{lo:.4f},{hi:.4f}] k={k} n={n}"
def silent(m): return bool(m.get("consistent") and m.get("amenable") and m.get("ostar_valid") is False)
for ens in ("original","licensed","large"):
    raw=json.load(open(f"/home/costaj/latent-causal/x1-spurious/results/x1_{ens}.json"))
    scms=raw["scms"]; print(f"== {ens}  (rho=1, reachable)")
    for arm,cond in (("R", lambda m:m.get("consistent")), ("Suni", lambda m:True)):
        n=chg=sil=am=0
        for s in scms:
            for m in s[arm]:
                if m["rho"]!=1 or m["hopC"]>=UNREACH or not cond(m): continue
                n+=1
                if m.get("amenable"):
                    am+=1
                    if m.get("ostar_changed"): chg+=1
                if silent(m): sil+=1
        print(f"  {arm:5s} n={n:5d}  P(O* moved to a DIFFERENT SET) {f(chg,n)}")
        print(f"        P(silent) {f(sil,n)}    P(invalid | moved) {f(sil,chg)}")
    # hop>=1 ostar_changed contrast, the clean-channel locality result
    print("  --- hop>=1 reachable, O* moved to a different set:")
    for arm,cond in (("R", lambda m:m.get("consistent")), ("Suni", lambda m:True),
                     ("Suni|S1", lambda m:m.get("consistent_S1"))):
        a = "Suni" if arm.startswith("Suni") else "R"
        n=chg=0
        for s in scms:
            for m in s[a]:
                if m["rho"]!=1 or not (1<=m["hopC"]<UNREACH) or not cond(m): continue
                n+=1
                if m.get("amenable") and m.get("ostar_changed"): chg+=1
        print(f"      {arm:9s} {f(chg,n)}")
    # arm D: does withdrawal EVER move O* to a different set?
    n=chg=nonam=0
    for s in scms:
        for m in s["Drop"]:
            if m["rho"]!=1 or m["hopC"]>=UNREACH: continue
            n+=1
            if m.get("amenable"):
                if m.get("ostar_changed"): chg+=1
            else: nonam+=1
    print(f"  ARM D withdrawal: O* moved to a DIFFERENT SET {f(chg,n)} ; O* became UNDEFINED {f(nonam,n)}")
