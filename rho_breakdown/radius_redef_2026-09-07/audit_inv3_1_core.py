#!/usr/bin/env python3
"""inverse-3 (fragility polynomial nu(q)) audit 1: is it non-degenerate, and is
it anything other than a rescaled count of singly-load-bearing claims?

nu_j := #{J admissible, |J|=j, Z not valid in G_J} / #{J admissible, |J|=j}
nu(q) := sum_j C(m,j) q^j (1-q)^{m-j} nu_j
Also computed:
  a_j       : admissible fraction of layer j  (the reweighting they admit to)
  nu_cond(q): the probability the DEFINITION SAYS IT IS -- each claim independently
              wrong w.p. q, conditioned on admissibility GLOBALLY, not per layer.
  |B1|      : the number of single claims whose reversal alone breaks Z
"""
import sys, itertools, math
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import make_problem, ball, r_val
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.evaluate import is_valid_adjustment_set_mpdag

QS = [0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50]


def corruption_table(pr):
    """For every J subset of K: admissible?  breaks Z?"""
    K, cp, x, y, z = pr["K"], pr["cpdag"], pr["x"], pr["y"], pr["z"]
    m = len(K)
    rows = []                                   # (|J|, admissible, breaks, J)
    for j in range(m + 1):
        for J in itertools.combinations(range(m), j):
            KJ = [(b, a) if i in J else (a, b) for i, (a, b) in enumerate(K)]
            g = apply_orientations(cp, KJ)
            if g is None:
                rows.append((j, False, False, J))
            else:
                rows.append((j, True, not is_valid_adjustment_set_mpdag(g, x, y, z), J))
    return m, rows


def nu_profile(m, rows):
    nu_j, a_j = [], []
    for j in range(m + 1):
        lay = [r for r in rows if r[0] == j]
        adm = [r for r in lay if r[1]]
        a_j.append(len(adm) / len(lay))
        nu_j.append((sum(r[2] for r in adm) / len(adm)) if adm else 0.0)
    return nu_j, a_j


def nu(q, m, nu_j):
    return sum(math.comb(m, j) * q**j * (1 - q)**(m - j) * nu_j[j] for j in range(m + 1))


def nu_cond(q, m, rows):
    """What the prose model actually specifies: P(break | admissible)."""
    num = den = 0.0
    for j, adm, brk, _ in rows:
        if not adm:
            continue
        w = q**j * (1 - q)**(m - j)
        den += w
        num += w * brk
    return num / den if den else 0.0


if __name__ == "__main__":
    rng = np.random.default_rng(20260907)
    got, tried = 0, 0
    R = []
    while got < 400 and tried < 40000:
        tried += 1
        pr = make_problem(rng)
        if pr is None:
            continue
        rows_ball = ball(pr)
        if not rows_ball:
            continue
        m, rows = corruption_table(pr)
        if m == 0:
            continue
        got += 1
        nu_j, a_j = nu_profile(m, rows)
        B1 = sum(1 for j, adm, brk, _ in rows if j == 1 and adm and brk)
        R.append(dict(m=m, nu_j=nu_j, a_j=a_j, B1=B1,
                      rval=r_val(rows_ball),
                      nu={q: nu(q, m, nu_j) for q in QS},
                      nuc={q: nu_cond(q, m, rows) for q in QS}))

    print(f"problemas: {got} de {tried} sorteios;  |K| = {Counter(r['m'] for r in R)}")

    # ---- 1. resolution: how many distinct values does nu really take? --------
    print("\n[1] RESOLUCAO de nu(q) sobre os problemas  (r_val toma 2-3 valores)")
    for q in QS:
        v = np.array([r["nu"][q] for r in R])
        print(f"   q={q:<5} distintos {len(set(np.round(v,12))):3d}   "
              f"min {v.min():.4f}  mediana {np.median(v):.4f}  max {v.max():.4f}   "
              f"=0 em {100*np.mean(v==0):.1f}%")
    c = Counter("UNREACHED" if r["rval"] is None else r["rval"] for r in R)
    print(f"   r_val: {dict(c)}")

    # ---- 2. is nu just a rescaled |B1|? -------------------------------------
    print("\n[2] nu(q) e so uma reescala de |B1| (n. de claims isoladamente load-bearing)?")
    B1 = np.array([r["B1"] for r in R], float)
    M  = np.array([r["m"] for r in R], float)
    from scipy.stats import spearmanr
    for q in QS:
        v = np.array([r["nu"][q] for r in R])
        rho = spearmanr(v, B1).statistic
        # exact first-order predictor: q*(1-q)^(m-1) * m * nu_1  == q*(1-q)^(m-1)*|B1|/a1frac
        approx = np.array([q * (1-q)**(r["m"]-1) * r["m"] * r["nu_j"][1] for r in R])
        rel = np.abs(approx - v) / np.maximum(v, 1e-12)
        print(f"   q={q:<5} spearman(nu,|B1|) = {rho:6.4f} | "
              f"erro relativo do termo j=1 sozinho: mediana {np.median(rel):.4f} p90 {np.percentile(rel,90):.4f}")
    # rank invariance in q
    from scipy.stats import kendalltau
    lo = np.array([r["nu"][0.01] for r in R]); hi = np.array([r["nu"][0.50] for r in R])
    print(f"   kendall tau entre a ordenacao em q=0.01 e q=0.50: {kendalltau(lo,hi).statistic:.4f}")
    print("   -> se ~1, o 'dial continuo' q nao separa nada: reordena zero problemas")

    # ---- 3. new degeneracy: nu_j == 1 for all j>=1 --------------------------
    print("\n[3] DEGENERACAO NOVA: nu_j = 1 para todo j>=1  =>  nu(q) = 1-(1-q)^m, funcao so de m")
    sat = [r for r in R if all(r["nu_j"][j] == 1.0 for j in range(1, r["m"]+1))]
    print(f"   saturados: {len(sat)}/{len(R)} = {100*len(sat)/len(R):.1f}%")
    rv1 = [r for r in R if r["rval"] == 1]
    sat1 = [r for r in rv1 if all(r["nu_j"][j] == 1.0 for j in range(1, r["m"]+1))]
    print(f"   dentro do estrato r_val=1 ({len(rv1)}): saturados {len(sat1)} = "
          f"{100*len(sat1)/max(len(rv1),1):.1f}%")
    # variance of nu(0.1) explained by (m,|B1|) alone
    v = np.array([r["nu"][0.10] for r in R])
    key = [(r["m"], r["B1"]) for r in R]
    grp = {}
    for k, val in zip(key, v):
        grp.setdefault(k, []).append(val)
    within = sum(np.var(g)*len(g) for g in grp.values()) / len(v)
    print(f"   variancia de nu(0.10) explicada so por (m,|B1|): "
          f"{100*(1 - within/np.var(v)):.2f}%   ({len(grp)} celulas)")

    # ---- 4. the Bernstein form is NOT the declared probability --------------
    print("\n[4] A FORMULA NAO E A PROBABILIDADE DECLARADA (renormalizacao por camada)")
    aj_drop = [r["a_j"] for r in R]
    maxm = max(r["m"] for r in R)
    for j in range(maxm + 1):
        vals = [a[j] for a in aj_drop if len(a) > j]
        if vals:
            print(f"   a_{j} (fracao admissivel na camada {j}): media {np.mean(vals):.4f}")
    for q in [0.05, 0.10, 0.20]:
        a = np.array([r["nu"][q] for r in R]); b = np.array([r["nuc"][q] for r in R])
        d = a - b
        print(f"   q={q}: nu_Bernstein - nu_declarada  media {d.mean():+.5f}  "
              f"max {d.max():+.5f}  >0 em {100*np.mean(d>1e-12):.1f}%  "
              f"| erro relativo mediano {np.median(np.abs(d)/np.maximum(b,1e-12)):.4f}")
