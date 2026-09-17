#!/usr/bin/env python3
"""Carrega uma rede gaussiana REAL do corpus do Chris como um LinearSEM.

ecoli70, arth150, magic-niab e magic-irri sao redes bayesianas GAUSSIANAS do
bnlearn: os coeficientes nao sao simulados, sao ajustados a dados reais
(ecoli70 = microarranjos de E. coli, Schaefer & Strimmer 2005). Isso da a unica
coisa que faltava do nosso lado: estrutura real COM parametros lineares reais.
"""
import json, sys
from pathlib import Path
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
import numpy as np
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.example import dag_to_cpdag, knowledge_to_recover
from bkrobust.demo.evaluate import LinearSEM

CORPUS = Path("/Users/josecosta/bkrobust/results/axisa3/networks/example_models")


def load(name):
    d = json.load(open(CORPUS / f"{name}.json"))
    nodes = sorted(d["nodes"])
    arcs = [(a, b) for a, b in d["arcs"]]
    dag = MPDAG(nodes=nodes, directed=arcs, undirected=[])
    W = {}       # (parent, child) -> coeficiente ajustado
    var = {}     # no -> variancia residual ajustada
    for v, cpd in d["cpds"].items():
        var[v] = float(cpd["variance"][0])
        for p, c in cpd["coefficients"].items():
            if p == "(Intercept)":
                continue
            W[(p, v)] = float(c[0])
    return dag, W, var


def sem_of(dag, W, var):
    """LinearSEM com os coeficientes REAIS (nao um random_sem)."""
    return LinearSEM(dag=dag,
                     weights={e: W[e] for e in W if e in set(dag.directed_edges)},
                     noise_var=dict(var))


if __name__ == "__main__":
    for name in ("ecoli70", "arth150", "magic-niab", "magic-irri"):
        try:
            dag, W, var = load(name)
        except Exception as e:
            print(f"{name:12s} FALHOU: {e}"); continue
        c = dag_to_cpdag(dag)
        und = len(c.undirected_edges)
        print(f"{name:12s} n={len(dag.nodes):4d} arcos={len(dag.directed_edges):4d} "
              f"nao-orientadas={und:4d} ({100*und/max(1,len(dag.directed_edges)):.1f}%) "
              f"coefs={len(W)}")
