#!/usr/bin/env python3
"""A direcao do erro: o raio calculado de Sigma-hat erra para o lado seguro?

Concordancia e' a metrica errada para a banda. A banda e' CONSERVADORA por
construcao (o casco alargado atinge zero mais cedo), entao ela deve discordar do
oraculo -- a pergunta e' se ela alguma vez discorda para o lado PERIGOSO, isto e',
declarando um raio MAIOR do que o verdadeiro, o que afirmaria mais robustez do
que existe.
"""
import json, sys
from collections import Counter
import numpy as np

rows = json.load(open("sweep.json"))
INF = float("inf")
f = lambda v: INF if v is None else v

for cov in (1.0, 0.5):
    rs = [r for r in rows if r["coverage"] == cov]
    if not rs: continue
    print(f"\n{'='*88}\n  COBERTURA {cov}  ({len(rs)} instancias)\n{'='*88}")
    print(f"{'estimador':>10} {'n':>7} {'igual':>8} {'seguro (<)':>11} "
          f"{'PERIGOSO (>)':>13}")
    for kind, lab in (("pt", "pontual"), ("bd", "banda")):
        for n in (200, 1000, 5000, 20000):
            eq = lo = hi = 0
            for r in rs:
                t = f(r["r0"])
                for v in r[f"{kind}_n{n}"]:
                    v = f(v)
                    eq += v == t; lo += v < t; hi += v > t
            tot = eq + lo + hi
            print(f"{lab:>10} {n:>7} {100*eq/tot:7.1f}% {100*lo/tot:10.1f}% "
                  f"{100*hi/tot:12.1f}%")
        print()
    c = Counter("inf" if r["r0"] is None else r["r0"] for r in rs)
    print("  distribuicao de r_0 (oraculo): "
          + "  ".join(f"{k}: {v} ({100*v/len(rs):.0f}%)"
                      for k, v in sorted(c.items(), key=lambda kv: (kv[0] == "inf", kv[0]))))
    fin = [r for r in rs if r["r0"] is not None]
    if fin:
        st = sum(1 for r in fin if r["structural"])
        print(f"  entre os r_0 finitos, o zero e' ESTRUTURAL em {st}/{len(fin)} "
              f"({100*st/len(fin):.0f}%) -- alguma extensao da efeito exatamente nulo,")
        print("  isto e', nela X nem sequer e' causa de Y.")
