#!/usr/bin/env python3
"""Audit of weighted-0: r_K, the contradiction-budget breakdown number."""
import sys, itertools, math, time
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
import numpy as np

from bkrobust.demo.graph import MPDAG
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.space import (enumerate_space, represented_dags, covering_pairs,
                                 neighbour_graph, bfs_distances)
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag,
                                    is_valid_adjustment_set_mpdag,
                                    is_valid_adjustment_set_dag,
                                    optimal_adjustment_set_dag,
                                    adjusted_estimand, asymptotic_variance, random_sem)

def random_dag(rng, n, p):
    order = list(rng.permutation(n))
    edges = set()
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < p:
                edges.add((f"v{order[i]}", f"v{order[j]}"))
    return MPDAG(nodes=[f"v{i}" for i in range(n)], directed=edges, undirected=[])

def m_of(d, K):
    """number of elicited claims the DAG d contradicts"""
    return sum(1 for (a, b) in K if d.is_directed_edge(b, a))

def r_K(fibre, K, x, y, z):
    best = math.inf
    for d in fibre:
        if not is_valid_adjustment_set_dag(d, x, y, z):
            best = min(best, m_of(d, K))
    return best

def r_val(cpdag, g0, x, y, z):
    space = enumerate_space(cpdag)
    reps = represented_dags(space)
    covers = covering_pairs(space, reps)
    nbrs = neighbour_graph(space, covers)
    g0i = next((g for g in space if g == g0), None)
    if g0i is None: return None, None
    dist = bfs_distances(nbrs, g0i)
    bad = [dist[g] for g in space if g in dist and dist[g] > 0
           and not is_valid_adjustment_set_mpdag(g, x, y, z)]
    return (min(bad) if bad else math.inf), len(space)

def make_problem(rng, n=6, p=0.35, k_claims=3, n_false=0):
    """n_false = how many of the elicited claims are DELIBERATELY reversed."""
    dag = random_dag(rng, n, p)
    if not dag.directed_edges: return None
    cpdag = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cpdag.undirected_edges)
    if not (1 <= len(und) <= 6): return None
    nodes = sorted(dag.nodes)
    cand = [(x, y) for x in nodes for y in nodes if x != y]
    rng.shuffle(cand)
    for x, y in cand:
        idx = rng.permutation(len(und))[:min(k_claims, len(und))]
        if len(idx) <= n_false: continue
        K, truth = [], []
        for j, t in enumerate(idx):
            a, b = und[t]
            true_dir = (a, b) if (a, b) in dag.directed_edges else (b, a)
            if j < n_false:
                K.append((true_dir[1], true_dir[0]))   # LIE
            else:
                K.append(true_dir)
            truth.append(true_dir)
        g0 = apply_orientations(cpdag, K)
        if g0 is None: continue
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None: continue
        # PRECONDITION (P)
        if not is_valid_adjustment_set_mpdag(g0, x, y, z): continue
        fibre = enumerate_dag_extensions(cpdag)
        return dict(dag=dag, cpdag=cpdag, x=x, y=y, K=[tuple(e) for e in K],
                    g0=g0, z=frozenset(z), und=und, fibre=fibre, n_false=n_false)
    return None

def gen(rng, want, **kw):
    out, tried = [], 0
    while len(out) < want and tried < 60000:
        tried += 1
        pr = make_problem(rng, **kw)
        if pr is not None: out.append(pr)
    return out, tried
