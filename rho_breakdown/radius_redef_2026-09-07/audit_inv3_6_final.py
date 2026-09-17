#!/usr/bin/env python3
"""inverse-3 audit 6: the three remaining questions.
(a) Is supp(nu) exactly {r_val finite}?  (exact, clean generator)
(b) Do the corruption-reachable DAGs already span the WHOLE equivalence class?
    If so nu(q) can only be narrower than the blanket hedge by TRIMMING, and any
    trimming a defensible q produces is the J<=1 truncation.
(c) Soundness: does nu(q) <= alpha bound miscoverage?  And can nu be small while
    the damage is catastrophic (the heavy tail the constraints forbid bounding)?
"""
import sys, itertools
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_dag, random_sem,
                                    adjusted_estimand)
from audit_inv3_4_clean import clean_problem
from audit_inv3_1_core import corruption_table, nu_profile, nu
from harness import ball, r_val

if __name__ == "__main__":
    rng = np.random.default_rng(4242)
    got, tried = 0, 0
    agree = 0
    span_eq = span_tot = 0
    small_nu_big_dmg = []
    rec = []
    while got < 300 and tried < 200000:
        tried += 1
        pr = clean_problem(rng)
        if pr is None:
            continue
        rb = ball(pr)
        if not rb:
            continue
        m, rows = corruption_table(pr)
        if m == 0:
            continue
        got += 1
        nu_j, _ = nu_profile(m, rows)
        rv = r_val(rb)
        v10 = nu(0.10, m, nu_j)
        agree += ((v10 == 0.0) == (rv is None))
        # (b) span
        full = {frozenset(d.directed_edges) for d in enumerate_dag_extensions(pr["cpdag"], limit=500)}
        reach = set()
        for j in range(m + 1):
            for J in itertools.combinations(range(m), j):
                KJ = [(b, a) if i in J else (a, b) for i, (a, b) in enumerate(pr["K"])]
                g = apply_orientations(pr["cpdag"], KJ)
                if g is None:
                    continue
                for d in enumerate_dag_extensions(g, limit=200):
                    reach.add(frozenset(d.directed_edges))
        span_tot += 1
        span_eq += (reach == full)
        # (c) worst damage vs nu
        sem_worst = 0.0
        for j, adm, brk, J in rows:
            if not (adm and brk):
                continue
            KJ = [(b, a) if i in J else (a, b) for i, (a, b) in enumerate(pr["K"])]
            g = apply_orientations(pr["cpdag"], KJ)
            for d in enumerate_dag_extensions(g, limit=30):
                s2 = random_sem(d, rng)
                t = adjusted_estimand(s2, pr["x"], pr["y"], optimal_adjustment_set_dag(d, pr["x"], pr["y"]))
                e = adjusted_estimand(s2, pr["x"], pr["y"], pr["z"])
                if abs(t) > 1e-9:
                    sem_worst = max(sem_worst, abs(e - t) / abs(t))
        rec.append((v10, sem_worst, rv, m, len(reach), len(full)))

    print(f"problemas limpos: {got}")
    print(f"\n[a] supp(nu) == {{r_val finito}}: concordancia {100*agree/got:.2f}%  "
          f"({got-agree} discordancias em {got})")
    print(f"[b] DAGs alcancaveis por corrupcao == classe de equivalencia inteira: "
          f"{100*span_eq/span_tot:.1f}%  ({span_eq}/{span_tot})")
    frac = np.mean([r[4]/r[5] for r in rec])
    print(f"    fracao media da classe alcancada: {frac:.4f}")

    v = np.array([r[0] for r in rec]); w = np.array([r[1] for r in rec])
    print(f"\n[c] SOUNDNESS / cauda pesada")
    for thr in (0.05, 0.10, 0.20):
        sel = v <= thr
        if sel.sum():
            print(f"    nu(0.10) <= {thr}:  n={sel.sum():3d}  pior vies relativo "
                  f"max {w[sel].max():7.2f}  p90 {np.percentile(w[sel],90):6.2f}  "
                  f"|  fracao com vies > 1 (100%): {100*np.mean(w[sel]>1):.1f}%")
    sel = v == 0.0
    print(f"    nu == 0 exatamente (a definicao diz 'risco nulo'): n={sel.sum()}  "
          f"pior vies max {w[sel].max():.4f}")
    hi = v >= 0.30
    if hi.sum():
        print(f"    nu(0.10) >= 0.30: n={hi.sum()}  pior vies MEDIANO {np.median(w[hi]):.3f}  "
              f"min {w[hi].min():.4f}  |  com dano ZERO: {100*np.mean(w[hi]==0):.1f}%")
