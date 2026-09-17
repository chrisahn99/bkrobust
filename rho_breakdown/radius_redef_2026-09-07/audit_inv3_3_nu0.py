#!/usr/bin/env python3
"""inverse-3 audit 3: the definition asserts 'a Bernstein polynomial with nu(0) = 0'.
nu(0) = nu_0 = 1[ Z = O*(G0) is NOT certified valid in G0 ].  Is it really 0?
"""
import sys, itertools
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import make_problem, ball, r_val
from bkrobust.demo.evaluate import is_valid_adjustment_set_mpdag
from audit_inv3_1_core import corruption_table, nu_profile, nu

if __name__ == "__main__":
    rng = np.random.default_rng(20260907)
    got, tried = 0, 0
    rec = []
    while got < 500 and tried < 60000:
        tried += 1
        pr = make_problem(rng)
        if pr is None:
            continue
        rb = ball(pr)
        if not rb:
            continue
        m, rows = corruption_table(pr)
        if m == 0:
            continue
        got += 1
        nu_j, a_j = nu_profile(m, rows)
        v0 = bool(is_valid_adjustment_set_mpdag(pr["g0"], pr["x"], pr["y"], pr["z"]))
        rec.append(dict(m=m, nu_j=nu_j, nu0=nu_j[0], zvalid=v0, rval=r_val(rb),
                        z=pr["z"], nu01=nu(0.01, m, nu_j), nu10=nu(0.10, m, nu_j)))

    n = len(rec)
    bad = [r for r in rec if r["nu0"] != 0.0]
    print(f"problemas: {n}")
    print(f"\nnu(0) = nu_0 = 1  (Z = O*(G0) NAO certificado valido no PROPRIO G0):")
    print(f"   {len(bad)}/{n} = {100*len(bad)/n:.1f}%   <-- a definicao afirma 'nu(0) = 0 by construction'")
    print(f"   is_valid_adjustment_set_mpdag(G0,x,y,Z) falso em {100*np.mean([not r['zvalid'] for r in rec]):.1f}%")
    print(f"   Z vazio nesses casos: {100*np.mean([len(r['z'])==0 for r in bad]) if bad else 0:.1f}%")
    print(f"   Z vazio no geral:     {100*np.mean([len(r['z'])==0 for r in rec]):.1f}%")

    S1 = [r for r in rec if r["rval"] == 1]
    b1 = [r for r in S1 if r["nu0"] != 0]
    print(f"\ndentro do estrato r_val = 1 (n={len(S1)}): nu_0 = 1 em {len(b1)} = {100*len(b1)/len(S1):.1f}%")
    v = np.array([r["nu10"] for r in S1])
    vb = np.array([r["nu10"] for r in b1]); vg = np.array([r["nu10"] for r in S1 if r["nu0"]==0])
    print(f"   nu(0.10) mediana no estrato: {np.median(v):.4f}")
    print(f"      onde nu_0 = 1 (artefato):   mediana {np.median(vb):.4f}  n={len(vb)}")
    print(f"      onde nu_0 = 0 (legitimo):   mediana {np.median(vg):.4f}  n={len(vg)}")
    print(f"   distintos em q=0.10, so nos legitimos: {len(set(np.round(vg,12)))}")
    print(f"   nu(0.01) mediana no estrato: {np.median([r['nu01'] for r in S1]):.4f} "
          f"| so legitimos: {np.median([r['nu01'] for r in S1 if r['nu0']==0]):.4f}")

    # sanity ceiling: with nu_0 = 0 the polynomial cannot exceed P(at least one wrong)
    print(f"\nteto sao: nu(q) <= 1-(1-q)^m quando nu_0 = 0")
    for q in (0.01, 0.10):
        viol = [r for r in rec if r["nu_j"][0]==0 and nu(q, r["m"], r["nu_j"]) > 1-(1-q)**r["m"] + 1e-12]
        print(f"   q={q}: violacoes {len(viol)} (esperado 0);  "
              f"max nu(q) entre legitimos = {max((nu(q,r['m'],r['nu_j']) for r in rec if r['nu_j'][0]==0), default=0):.4f}"
              f" vs teto 1-(1-q)^3 = {1-(1-q)**3:.4f}")
