#!/usr/bin/env python3
"""Random-problem generator for the radius study.

One problem = (true DAG, its CPDAG, a query (X,Y), an elicited claim set K, the
Meek closure G0 = Meek(Chat, K), the adjustment set Z = O*(G0), and the enumerated
space around G0). Everything is exact: no sampling of data, population covariance.
"""
import sys, itertools
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
import numpy as np

from bkrobust.demo.graph import MPDAG
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.space import (enumerate_space, represented_dags, covering_pairs,
                                 neighbour_graph, bfs_distances)
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag,
                                    is_valid_adjustment_set_mpdag,
                                    is_valid_adjustment_set_dag)


def random_dag(rng, n, p):
    """DAG on v0..v{n-1} in a random topological order."""
    order = list(rng.permutation(n))
    edges = set()
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < p:
                edges.add((f"v{order[i]}", f"v{order[j]}"))
    return MPDAG(nodes=[f"v{i}" for i in range(n)], directed=edges, undirected=[])


def make_problem(rng, n=6, p=0.35, k_claims=3):
    """Returns a dict, or None if the draw is not usable."""
    dag = random_dag(rng, n, p)
    if not dag.directed_edges:
        return None
    cpdag = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cpdag.undirected_edges)
    if not (1 <= len(und) <= 6):          # keep 3^|und| enumerable
        return None

    # a query with a non-trivial adjustment problem
    nodes = sorted(dag.nodes)
    cand = [(x, y) for x in nodes for y in nodes if x != y]
    rng.shuffle(cand)
    for x, y in cand:
        # K: orient a random subset of the undirected edges, TRUTHFULLY
        idx = rng.permutation(len(und))[:min(k_claims, len(und))]
        K = []
        for t in idx:
            a, b = und[t]
            K.append((a, b) if (a, b) in dag.directed_edges else (b, a))
        g0 = apply_orientations(cpdag, K)
        if g0 is None:
            continue
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None:
            continue                       # not amenable: no adjustment set
        return dict(dag=dag, cpdag=cpdag, x=x, y=y, K=[tuple(e) for e in K],
                    g0=g0, z=frozenset(z), undirected=und)
    return None


def ball(problem):
    """The space, the distances from G0, and validity of Z at each element."""
    space = enumerate_space(problem["cpdag"])
    reps = represented_dags(space)
    covers = covering_pairs(space, reps)
    nbrs = neighbour_graph(space, covers)
    g0 = next((g for g in space if g == problem["g0"]), None)
    if g0 is None:
        return None
    dist = bfs_distances(nbrs, g0)
    x, y, z = problem["x"], problem["y"], problem["z"]
    rows = []
    for g in space:
        d = dist.get(g)
        if d is None:
            continue
        rows.append((d, bool(is_valid_adjustment_set_mpdag(g, x, y, z))))
    return rows


def r_val(rows):
    """The current definition: first shell containing an invalid element."""
    bad = [d for d, ok in rows if not ok and d > 0]
    return min(bad) if bad else None       # None = UNREACHED


if __name__ == "__main__":
    rng = np.random.default_rng(20260907)
    got, tried, radii, sizes = 0, 0, [], []
    while got < 400 and tried < 20000:
        tried += 1
        pr = make_problem(rng)
        if pr is None:
            continue
        rows = ball(pr)
        if not rows:
            continue
        got += 1
        radii.append(r_val(rows))
        sizes.append(len(rows))
    c = Counter("UNREACHED" if r is None else r for r in radii)
    print(f"problemas usaveis: {got} de {tried} sorteios")
    print(f"tamanho medio do espaco: {np.mean(sizes):.1f}  (max {max(sizes)})")
    print("distribuicao de r_val:")
    for k in sorted(c, key=lambda v: (v == 'UNREACHED', v)):
        print(f"   r_val = {k:<10} {c[k]:3d}   {100*c[k]/got:5.1f}%")
    deg = 100 * c.get(1, 0) / got
    print(f"\n=> r_val = 1 em {deg:.1f}% dos problemas  (o artefato que queremos remover)")
