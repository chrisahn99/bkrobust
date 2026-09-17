"""REFUTER pass on X1 raw per-SCM records. Independent of analyse_x1.py."""
import json, sys, math
import numpy as np
UNREACH = 1 << 20

def wil(k, n):
    if n == 0: return (float('nan'), float('nan'), float('nan'))
    z = 1.959963984540054
    p = k / n
    d = 1 + z*z/n
    c = (p + z*z/(2*n)) / d
    h = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / d
    return (p, max(0.0, c-h), min(1.0, c+h))

def fmt(k, n):
    p, lo, hi = wil(k, n)
    return f"{p:.4f} [{lo:.4f},{hi:.4f}] k={k} n={n}"

def silent(m):
    return bool(m.get("consistent") and m.get("amenable") and m.get("ostar_valid") is False)
def loud(m):
    return bool(m.get("consistent") and not m.get("amenable"))
def mpdag_valid(m):
    return bool(m.get("ext", True) and not m.get("cycle", False))

def hopband(h):
    if h >= UNREACH: return "unreach"
    h = int(h)
    if h == 0: return "0"
    if h == 1: return "1"
    if h == 2: return "2"
    return ">=3"

def main(ens):
    raw = json.load(open(f"/home/costaj/latent-causal/x1-spurious/results/x1_{ens}.json"))
    scms = raw["scms"]
    print(f"===================== {ens}  n_scm={len(scms)}")

    # ---------- 1. full outcome accounting per arm, rho=1, reachable
    print("\n-- [1] FULL OUTCOME ACCOUNTING, rho=1, reachable statements only")
    for arm in ("R", "Suni", "Sloc", "Drop"):
        n = ns = nl = nch = nam = ncons = 0
        ncons_silent = ncons_loud = ncons_n = 0
        for s in scms:
            for m in s[arm]:
                if m["rho"] != 1: continue
                if m["hopC"] >= UNREACH: continue
                n += 1
                if m.get("consistent"): ncons += 1
                if silent(m): ns += 1
                if loud(m): nl += 1
                if m.get("amenable"):
                    nam += 1
                    if m.get("ostar_changed"): nch += 1
        print(f"  {arm:5s} n={n:6d} consistent={ncons:6d} "
              f"amenable={nam:6d} ostar_changed={nch:6d} "
              f"silent={ns:5d} loud(non-amenable&cons)={nl:5d} "
              f"harm(silent+loud)={fmt(ns+nl, n)}")
    # arm R conditional on consistent
    kR = nR = lR = 0
    for s in scms:
        for m in s["R"]:
            if m["rho"] != 1 or m["hopC"] >= UNREACH: continue
            if not m.get("consistent"): continue
            nR += 1
            if silent(m): kR += 1
            if loud(m): lR += 1
    print(f"  R|consistent  silent={fmt(kR,nR)}  loud={fmt(lR,nR)}  harm={fmt(kR+lR,nR)}")
    kS = nS = lS = 0
    for s in scms:
        for m in s["Suni"]:
            if m["rho"] != 1 or m["hopC"] >= UNREACH: continue
            nS += 1
            if silent(m): kS += 1
            if loud(m): lS += 1
    print(f"  Suni (all)    silent={fmt(kS,nS)}  loud={fmt(lS,nS)}  harm={fmt(kS+lS,nS)}")
    print(f"  ==> HARM RR (S/R|cons) = {((kS+lS)/nS)/((kR+lR)/nR):.3f}")
    print(f"  ==> LOUD RR (S/R|cons) = {(lS/nS)/((lR/nR) if lR else float('nan')):.3f}")

    # ---------- 2. where do Suni's silent events live?
    print("\n-- [2] Suni rho=1 silent events by coherence status (reachable)")
    cells = {}
    for s in scms:
        for m in s["Suni"]:
            if m["rho"] != 1 or m["hopC"] >= UNREACH: continue
            key = (bool(m.get("consistent_S1")), mpdag_valid(m))
            c = cells.setdefault(key, [0, 0])
            c[1] += 1
            if silent(m): c[0] += 1
    for key in sorted(cells):
        k, n = cells[key]
        print(f"  consistent_S1={key[0]!s:5s} mpdag_valid={key[1]!s:5s}: silent={fmt(k,n)}")

    # ---------- 3. hop-stratified silent, several conditionings
    print("\n-- [3] hop-stratified silent rate, rho=1")
    def strat(arm, filt, label):
        d = {}
        for s in scms:
            for m in s[arm]:
                if m["rho"] != 1: continue
                if not filt(m): continue
                b = hopband(m["hopC"])
                c = d.setdefault(b, [0, 0]); c[1] += 1
                if silent(m): c[0] += 1
        order = ["0", "1", "2", ">=3", "unreach"]
        print(f"  {label}")
        tot_k = tot_n = 0
        for b in order:
            if b not in d: continue
            k, n = d[b]
            print(f"      hop {b:7s} {fmt(k,n)}")
            if b not in ("0", "unreach"): tot_k += k; tot_n += n
        h0 = d.get("0", [0, 0])
        r0 = h0[0]/h0[1] if h0[1] else float('nan')
        r1 = tot_k/tot_n if tot_n else float('nan')
        print(f"      hop>=1 reachable {fmt(tot_k, tot_n)}   materiality(hop>=1 / hop0) = {r1/r0 if r0 else float('nan'):.3f}")
        return d
    strat("R", lambda m: m.get("consistent"), "ARM R | consistent")
    strat("Suni", lambda m: True, "ARM S-uni (b-LOAD semantics, no check)")
    strat("Suni", lambda m: m.get("consistent_S1"), "ARM S-uni | consistent_S1 (careful tool)")
    strat("Suni", lambda m: mpdag_valid(m), "ARM S-uni | mpdag_valid")
    strat("Suni", lambda m: not m.get("conflict"), "ARM S-uni | no conflict")

    # ---------- 4. ostar_changed by hop (the deeper locality property)
    print("\n-- [4] P(ostar_changed | amenable) by hop, rho=1")
    for arm, filt in (("R", lambda m: m.get("consistent") and m.get("amenable")),
                      ("Suni", lambda m: m.get("amenable")),
                      ("Suni_S1", None)):
        d = {}
        a = "Suni" if arm == "Suni_S1" else arm
        f = filt if filt else (lambda m: m.get("amenable") and m.get("consistent_S1"))
        for s in scms:
            for m in s[a]:
                if m["rho"] != 1 or not f(m): continue
                b = hopband(m["hopC"])
                c = d.setdefault(b, [0, 0]); c[1] += 1
                if m.get("ostar_changed"): c[0] += 1
        print(f"  {arm}: " + "  ".join(f"hop{b}={d[b][0]}/{d[b][1]}" for b in ["0","1","2",">=3","unreach"] if b in d))

    # ---------- 5. |O0| stratification
    print("\n-- [5] silent rate by |O0| stratum, rho=1, reachable")
    for lab, pred in (("|O0|=0", lambda s: len(s["O0"]) == 0), ("|O0|>0", lambda s: len(s["O0"]) > 0)):
        kR = nR = kS = nS = 0
        for s in scms:
            if not pred(s): continue
            for m in s["R"]:
                if m["rho"] != 1 or m["hopC"] >= UNREACH or not m.get("consistent"): continue
                nR += 1; kR += silent(m)
            for m in s["Suni"]:
                if m["rho"] != 1 or m["hopC"] >= UNREACH: continue
                nS += 1; kS += silent(m)
        rr = (kS/nS)/(kR/nR) if nR and kR else float('nan')
        print(f"  {lab}: R|cons {fmt(kR,nR)}   S-uni {fmt(kS,nS)}   RR={rr:.3f}")

    # ---------- 6. hop>=1 silent events: how many at hop exactly 1
    print("\n-- [6] Suni hop>=1 silent event census, rho=1")
    from collections import Counter
    c = Counter()
    for s in scms:
        for m in s["Suni"]:
            if m["rho"] != 1: continue
            if m["hopC"] >= UNREACH or m["hopC"] < 1: continue
            if silent(m): c[int(m["hopC"])] += 1
    print("  ", dict(sorted(c.items())))

    # ---------- 7. rho sweep of the headline RR
    print("\n-- [7] silent rate by rho, reachable")
    for rho in (1, 2, 3, 4):
        kR = nR = kS = nS = 0
        for s in scms:
            for m in s["R"]:
                if m["rho"] != rho or m["hopC"] >= UNREACH or not m.get("consistent"): continue
                nR += 1; kR += silent(m)
            for m in s["Suni"]:
                if m["rho"] != rho or m["hopC"] >= UNREACH: continue
                nS += 1; kS += silent(m)
        rr = (kS/nS)/(kR/nR) if nR and kR else float('nan')
        print(f"  rho={rho}: R|cons {fmt(kR,nR)}   S-uni {fmt(kS,nS)}   RR={rr:.3f}")

if __name__ == "__main__":
    for e in sys.argv[1:]:
        main(e)
