#!/usr/bin/env python3
"""O raio de nulidade: em que orcamento de retracao o casco passa a conter zero?

Motivacao. r_val pergunta "o conjunto de ajuste ainda e valido?", um evento
discreto que medimos ser quase binario. Rosenbaum nunca perguntou isso: ele
pergunta quanto vies oculto seria preciso para VIRAR A CONCLUSAO. A traducao
direta para este cenario e: quantas afirmacoes do especialista teriam de ser
retiradas antes de eu deixar de poder afirmar que existe um efeito?
"""
import sys, itertools, json
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
import numpy as np
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.evaluate import random_sem, adjusted_estimand
import confirm
from confirm import problem, theta, hull, parents


def radii(pr, sem, cache, jmax=4):
    """Para cada orcamento: o casco, se contem zero, e se contem os dois sinais."""
    cpdag, x, y, K = pr["cpdag"], pr["x"], pr["y"], pr["K"]
    out = {}
    for j in range(jmax + 1):
        vals = []
        for r in range(j + 1):
            for S in itertools.combinations(range(len(K)), r):
                Kp = [K[i] for i in range(len(K)) if i not in S]
                g = apply_orientations(cpdag, Kp) if Kp else cpdag
                if g is None: continue
                vals += theta(g, sem, x, y, cache)
        h = hull(vals)
        if h is None: continue
        out[j] = dict(lo=h[0], hi=h[1],
                      has_zero=bool(h[0] <= 0.0 <= h[1]),
                      both_signs=bool(h[0] < -1e-12 and h[1] > 1e-12))
    return out


def first(d, key):
    for j in sorted(d):
        if d[j][key]:
            return j
    return None


def run(N=400, seed=20260908, n_false=0):
    rng = np.random.default_rng(seed)
    rows, tried = [], 0
    while len(rows) < N and tried < 40000:
        tried += 1
        pr = problem(rng, n_false=n_false)
        if pr is None: continue
        dag, x, y = pr["dag"], pr["x"], pr["y"]
        sem = random_sem(dag, rng); cache = {}
        try:
            tau = adjusted_estimand(sem, x, y, parents(dag, x))
        except Exception:
            continue
        d = radii(pr, sem, cache)
        if not d: continue
        rows.append(dict(tau=tau, d=d, rnull=first(d, "has_zero"),
                         rsign=first(d, "both_signs")))
    return rows, tried


if __name__ == "__main__":
    for nf in (0, 1):
        rows, tried = run(n_false=nf)
        eff = [r for r in rows if abs(r["tau"]) > 1e-9]     # so onde ha efeito a perder
        print(f"\n{'='*62}\n  {nf} afirmacao(oes) invertida(s) | {len(rows)} problemas, "
              f"{len(eff)} com efeito nao nulo\n{'='*62}")
        for name, key in (("RAIO DE NULIDADE  (casco passa a conter 0)", "rnull"),
                          ("RAIO DE SINAL     (casco contem os dois sinais)", "rsign")):
            c = Counter("nunca" if r[key] is None else r[key] for r in eff)
            tot = len(eff)
            print(f"\n{name}")
            for k in sorted(c, key=lambda v: (v == "nunca", v)):
                print(f"    = {str(k):<7} {c[k]:4d}  {100*c[k]/tot:5.1f}%")
            interior = sum(v for k, v in c.items() if k not in ("nunca", 0))
            print(f"    -> massa no INTERIOR (nem 0 nem nunca): {100*interior/tot:.1f}%")
