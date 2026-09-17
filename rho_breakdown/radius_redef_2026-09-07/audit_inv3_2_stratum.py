#!/usr/bin/env python3
"""inverse-3 audit 2: (a) is the support of nu exactly {r_val finite}?
(b) inside the r_val = 1 stratum -- the whole reason this definition exists --
what resolution does nu actually have, and does it track DAMAGE?
(c) GAMING: pad K with claims that are Meek-IMPLIED by the rest.  Those change
neither G0 nor Z nor the analyst's epistemic state.  Does nu move?
"""
import sys, itertools, math
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from scipy.stats import spearmanr
from harness import make_problem, ball, r_val
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (is_valid_adjustment_set_mpdag, random_sem,
                                    adjusted_estimand, is_valid_adjustment_set_dag,
                                    optimal_adjustment_set_dag)
from audit_inv3_1_core import corruption_table, nu_profile, nu, nu_cond

QS = [0.01, 0.05, 0.10, 0.20]


def nu_all(K, cpdag, x, y, z):
    pr = dict(K=K, cpdag=cpdag, x=x, y=y, z=z)
    m, rows = corruption_table(pr)
    nu_j, a_j = nu_profile(m, rows)
    return m, nu_j, a_j, rows


if __name__ == "__main__":
    rng = np.random.default_rng(20260907)
    got, tried = 0, 0
    R = []
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
        rv = r_val(rb)
        # ---- damage: worst |relative bias| over corruptions that break Z ----
        sem = random_sem(pr["dag"], rng)
        true_eff = adjusted_estimand(sem, pr["x"], pr["y"],
                                     optimal_adjustment_set_dag(pr["dag"], pr["x"], pr["y"]))
        worst = 0.0
        seen_break = False
        for j, adm, brk, J in rows:
            if not (adm and brk):
                continue
            seen_break = True
            KJ = [(b, a) if i in J else (a, b) for i, (a, b) in enumerate(pr["K"])]
            g = apply_orientations(pr["cpdag"], KJ)
            for d in enumerate_dag_extensions(g, limit=40):
                # a world compatible with the corrupted knowledge
                s2 = random_sem(d, rng)
                zt = optimal_adjustment_set_dag(d, pr["x"], pr["y"])
                truth = adjusted_estimand(s2, pr["x"], pr["y"], zt)
                got_est = adjusted_estimand(s2, pr["x"], pr["y"], pr["z"])
                if abs(truth) > 1e-9:
                    worst = max(worst, abs(got_est - truth) / abs(truth))
        # ---- gaming: pad K with Meek-implied (no-op) claims ------------------
        g0 = pr["g0"]
        Kset = {tuple(e) for e in pr["K"]}
        implied = [e for e in g0.directed_edges
                   if tuple(sorted(e)) in set(pr["undirected"]) and tuple(e) not in Kset]
        Kpad = list(pr["K"]) + [tuple(e) for e in implied]
        gpad = apply_orientations(pr["cpdag"], Kpad)
        pad_ok = gpad is not None and set(gpad.directed_edges) == set(g0.directed_edges)
        nu_pad = None
        if pad_ok and implied:
            mp, nujp, ajp, rowsp = nu_all(Kpad, pr["cpdag"], pr["x"], pr["y"], pr["z"])
            nu_pad = ({q: nu(q, mp, nujp) for q in QS}, mp,
                      {q: nu_cond(q, mp, rowsp) for q in QS})
        R.append(dict(m=m, nu_j=nu_j, rval=rv, worst=worst, brk=seen_break,
                      B1=sum(1 for j, a, b, _ in rows if j == 1 and a and b),
                      nu={q: nu(q, m, nu_j) for q in QS},
                      nuc={q: nu_cond(q, m, rows) for q in QS},
                      n_implied=len(implied), nu_pad=nu_pad))

    print(f"problemas: {got}")

    # ---- (a) support -------------------------------------------------------
    print("\n[a] SUPORTE: nu(q) > 0  vs  r_val finito")
    z0 = np.array([r["nu"][0.10] == 0 for r in R])
    inf = np.array([r["rval"] is None for r in R])
    print(f"   nu=0 em {z0.mean()*100:.1f}%;  r_val=UNREACHED em {inf.mean()*100:.1f}%;  "
          f"concordancia {100*np.mean(z0 == inf):.2f}%")
    print(f"   nu=0 e r_val finito: {int(np.sum(z0 & ~inf))} | nu>0 e r_val=UNREACHED: {int(np.sum(~z0 & inf))}")
    print("   -> onde discordam, quem 'nao e degenerado'?")

    # ---- (b) inside the r_val = 1 stratum ---------------------------------
    S = [r for r in R if r["rval"] == 1]
    print(f"\n[b] ESTRATO r_val = 1  (n = {len(S)})")
    for q in QS:
        v = np.array([r["nu"][q] for r in S])
        print(f"   q={q:<5} distintos {len(set(np.round(v,12))):3d}  min {v.min():.4f} "
              f"med {np.median(v):.4f} max {v.max():.4f}  | =max em {100*np.mean(v==v.max()):.1f}%")
    key = Counter((r["m"], r["B1"]) for r in S)
    print(f"   celulas (m,|B1|) distintas: {len(key)}  -> {dict(key)}")
    w = np.array([r["worst"] for r in S]); v = np.array([r["nu"][0.10] for r in S])
    ok = w > 0
    print(f"   spearman(nu(0.10), pior vies relativo) no estrato: "
          f"{spearmanr(v[ok], w[ok]).statistic:+.4f}  (n={ok.sum()})")
    wa = np.array([r["worst"] for r in R]); va = np.array([r["nu"][0.10] for r in R])
    oka = wa > 0
    print(f"   spearman global (todos os problemas):            "
          f"{spearmanr(va[oka], wa[oka]).statistic:+.4f}  (n={oka.sum()})")
    print(f"   pior vies relativo | Z quebra: mediana {np.median(wa[oka]):.3f} "
          f"p90 {np.percentile(wa[oka],90):.3f} max {wa[oka].max():.1f}")

    # ---- (c) gaming with no-op implied claims ------------------------------
    G = [r for r in R if r["nu_pad"] is not None]
    print(f"\n[c] GAMING: adicionar claims MEEK-IMPLICADOS (G0 e Z inalterados). n = {len(G)}")
    for q in QS:
        a = np.array([r["nu"][q] for r in G]); b = np.array([r["nu_pad"][0][q] for r in G])
        ac = np.array([r["nuc"][q] for r in G]); bc = np.array([r["nu_pad"][2][q] for r in G])
        mv = a != b
        ratio = b[a > 0] / a[a > 0]
        print(f"   q={q:<5} nu muda em {100*mv.mean():5.1f}% | razao nu_pad/nu: "
              f"mediana {np.median(ratio):.4f} min {ratio.min():.4f} max {ratio.max():.4f}"
              f" || nu_declarada muda em {100*np.mean(np.abs(ac-bc)>1e-12):.1f}%")
    print(f"   claims implicados adicionados: mediana {np.median([r['n_implied'] for r in G]):.0f} "
          f"max {max(r['n_implied'] for r in G)}")
