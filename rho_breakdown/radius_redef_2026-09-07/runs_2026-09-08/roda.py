#!/usr/bin/env python3
"""Roda os cinco experimentos e grava tudo em resultados.json.
Nenhum numero desta apostila e digitado a mao: todos saem daqui."""
import json, io, contextlib, sys
from collections import Counter
import numpy as np
import harness, confirm, signradius, stratum, local_hedge

out = {}

# --- Exp 1: a distribuicao de r_val -----------------------------------------
rng = np.random.default_rng(20260907)
radii, sizes, got, tried = [], [], 0, 0
while got < 400 and tried < 20000:
    tried += 1
    pr = harness.make_problem(rng)
    if pr is None: continue
    rows = harness.ball(pr)
    if not rows: continue
    got += 1
    radii.append(harness.r_val(rows)); sizes.append(len(rows))
c = Counter("UNREACHED" if r is None else r for r in radii)
out["exp1"] = dict(n=got, tried=tried, space_mean=float(np.mean(sizes)),
                   space_max=int(max(sizes)),
                   dist={str(k): int(v) for k, v in c.items()},
                   frac_1=c.get(1, 0)/got, frac_unreached=c.get("UNREACHED", 0)/got,
                   frac_binary=(c.get(1,0)+c.get("UNREACHED",0))/got)

# --- Exp 2: a fronteira cobertura-largura, tres bracos ----------------------
out["exp2"] = {}
for nf in (0, 1, 2):
    rows, tr = confirm.run(N=300, n_false=nf)
    inf = [r for r in rows if r["wb"] > 1e-9]
    arm = dict(n=len(rows), tried=tr, n_inf=len(inf), frac_inf=len(inf)/len(rows), budgets={})
    for j in range(4):
        hs = [r for r in rows if j in r["hulls"] and r["hulls"][j]]
        if not hs: continue
        cov = float(np.mean([h["hulls"][j][0]-1e-12 <= h["tau"] <= h["hulls"][j][1]+1e-12 for h in hs]))
        ii = [h for h in hs if h["wb"] > 1e-9]
        arm["budgets"][str(j)] = dict(
            cover=cov,
            ratio=float(np.mean([(h["hulls"][j][1]-h["hulls"][j][0])/h["wb"] for h in ii])) if ii else None,
            atom=float(np.mean([abs((h["hulls"][j][1]-h["hulls"][j][0])-h["wb"])<1e-9 for h in ii])) if ii else None)
    arm["blanket_cover"] = float(np.mean([r["blanket"][0]-1e-12 <= r["tau"] <= r["blanket"][1]+1e-12 for r in rows]))
    dis = [r for r in inf if not r["agree"]]
    share = [r["wdir"]/r["wb"] for r in inf]
    arm["direction"] = dict(disputed=len(dis), of=len(inf), frac=len(dis)/len(inf),
                            width_share=float(np.mean(share)),
                            frac_zero=float(np.mean([s < 1e-9 for s in share])))
    out["exp2"][str(nf)] = arm


# --- Exp 3: o raio de nulidade (a traducao de Rosenbaum) ---------------------
out["exp3"] = {}
for nf in (0, 1):
    rows, tried = signradius.run(N=400, n_false=nf)
    eff = [r for r in rows if abs(r["tau"]) > 1e-9]
    arm = dict(n=len(rows), n_efeito=len(eff))
    for name, key in (("nulidade", "rnull"), ("sinal", "rsign")):
        c = Counter("nunca" if r[key] is None else r[key] for r in eff)
        arm[name] = {str(k): c[k] / len(eff) for k in c}
        arm[name + "_interior"] = sum(v for k, v in c.items() if k not in ("nunca", 0)) / len(eff)
    out["exp3"][str(nf)] = arm

# --- Exp 4: o estrato informativo, e a troca de lugar ------------------------
import math
def _resumo(rs):
    n = len(rs)
    d = dict(n=n)
    for name, key in (("r_val", "rval"), ("nulidade", "rnull")):
        c = Counter("nunca" if r[key] is None else r[key] for r in rs)
        ps = [v / n for v in c.values()]
        H = -sum(p * math.log2(p) for p in ps if p > 0)
        d[name] = dict(dist={str(k): c[k] / n for k in c}, entropia=H,
                       valores_efetivos=2 ** H,
                       niveis_5pc=sum(1 for p in ps if p >= 0.05))
    return d

srows, stried = stratum.run(N=1800, n_false=1)
sinf = [r for r in srows if r["wb"] > 1e-9]
snon = [r for r in srows if r["wb"] <= 1e-9]
out["exp4"] = dict(n=len(srows), tried=stried, frac_inf=len(sinf) / len(srows),
                   todos=_resumo(srows), informativos=_resumo(sinf),
                   nao_informativos=_resumo(snon))

# --- Exp 5: a lei de localidade, medida em cobertura -------------------------
out["exp5"] = {}
POL = ("PLAIN", "RH1_far", "RH1_near", "RH_1", "LOCAL", "RH_2", "BLANKET")
for nf in (1, 2):
    rows, tried = local_hedge.run(N=600, n_false=nf)
    arm = dict(n=len(rows), tried=tried,
               n_near=float(np.mean([r["n_near"] for r in rows])),
               nK=float(np.mean([r["nK"] for r in rows])), politicas={})
    for name in POL:
        cov = float(np.mean([r["pol"][name][0] - 1e-12 <= r["tau"] <= r["pol"][name][1] + 1e-12
                             for r in rows]))
        w = float(np.mean([(r["pol"][name][1] - r["pol"][name][0]) / r["wb"] for r in rows]))
        if name == "RH_1":       sz = np.mean([1 + r["nK"] for r in rows])
        elif name == "RH1_near": sz = np.mean([1 + r["n_near"] for r in rows])
        elif name == "RH1_far":  sz = np.mean([1 + r["nK"] - r["n_near"] for r in rows])
        elif name == "RH_2":     sz = np.mean([1 + r["nK"] + r["nK"] * (r["nK"] - 1) / 2 for r in rows])
        elif name == "LOCAL":    sz = np.mean([2 ** r["n_near"] for r in rows])
        elif name == "FAR":      sz = np.mean([2 ** (r["nK"] - r["n_near"]) for r in rows])
        else:                    sz = 1.0
        arm["politicas"][name] = dict(cobre=cov, largura=w, enumeracoes=float(sz))
    out["exp5"][str(nf)] = arm

json.dump(out, open("resultados.json", "w"), indent=1)
print("exp1  r_val quase binario :", round(out["exp1"]["frac_binary"], 4))
print("exp2  fronteira, 1 erro    :", {j: round(b["ratio"], 3) for j, b in out["exp2"]["1"]["budgets"].items()})
print("exp3  nulidade 'nunca'     :", round(out["exp3"]["1"]["nulidade"].get("nunca", 0), 3))
print("exp4  informativos         :", round(out["exp4"]["frac_inf"], 4),
      "| r_val val.efetivos", round(out["exp4"]["informativos"]["r_val"]["valores_efetivos"], 2),
      "| nulidade", round(out["exp4"]["informativos"]["nulidade"]["valores_efetivos"], 2))
print("exp5  localidade, 1 erro   :",
      {k: (round(v["cobre"], 3), round(v["largura"], 3)) for k, v in out["exp5"]["1"]["politicas"].items()})
