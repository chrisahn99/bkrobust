#!/usr/bin/env python3
"""inverse-3 audit 7: gaming, on the CLEAN generator (nu(0)=0 guaranteed).
Pad K with claims Meek-IMPLIED by K.  G0 is bit-identical, Z is identical, the
analyst's epistemic state is identical.  Does nu(q) move?
"""
import sys
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from bkrobust.demo.meek import apply_orientations
from audit_inv3_4_clean import clean_problem
from audit_inv3_1_core import corruption_table, nu_profile, nu, nu_cond

QS = [0.01, 0.05, 0.10, 0.20]
if __name__ == "__main__":
    rng = np.random.default_rng(20260907)
    got = tried = 0
    A, B, AC, BC, npad = {q: [] for q in QS}, {q: [] for q in QS}, {q: [] for q in QS}, {q: [] for q in QS}, []
    while got < 250 and tried < 200000:
        tried += 1
        pr = clean_problem(rng)
        if pr is None:
            continue
        g0, Kset = pr["g0"], {tuple(e) for e in pr["K"]}
        implied = [tuple(e) for e in g0.directed_edges
                   if tuple(sorted(e)) in set(pr["undirected"]) and tuple(e) not in Kset]
        if not implied:
            continue
        Kpad = list(pr["K"]) + implied
        gpad = apply_orientations(pr["cpdag"], Kpad)
        if gpad is None or set(gpad.directed_edges) != set(g0.directed_edges) \
           or set(gpad.undirected_edges) != set(g0.undirected_edges):
            continue
        got += 1
        npad.append(len(implied))
        m0, r0 = corruption_table(pr)
        n0, _ = nu_profile(m0, r0)
        p2 = dict(pr); p2["K"] = Kpad
        m1, r1 = corruption_table(p2)
        n1, _ = nu_profile(m1, r1)
        for q in QS:
            A[q].append(nu(q, m0, n0)); B[q].append(nu(q, m1, n1))
            AC[q].append(nu_cond(q, m0, r0)); BC[q].append(nu_cond(q, m1, r1))

    print(f"pares (K, K+implicados) LIMPOS: {got};  claims implicados adicionados: "
          f"mediana {np.median(npad):.0f} max {max(npad)}")
    print("\nG0 identico, Z identico, conhecimento identico -- e mesmo assim:")
    for q in QS:
        a, b = np.array(A[q]), np.array(B[q])
        ac, bc = np.array(AC[q]), np.array(BC[q])
        nz = a > 1e-15
        r = b[nz] / a[nz]
        print(f"  q={q:<5} nu muda em {100*np.mean(np.abs(a-b)>1e-12):5.1f}%  | "
              f"razao nu_pad/nu (onde nu>0): mediana {np.median(r):.4f} min {r.min():.4f} max {r.max():.4f} "
              f"| queda >=2x em {100*np.mean(r<=0.5):.1f}%   || nu_declarada muda em "
              f"{100*np.mean(np.abs(ac-bc)>1e-12):.1f}%")
