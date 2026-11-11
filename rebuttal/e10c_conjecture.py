"""E10c: exhaustive test of "valid at G0 and amenable at G => valid at G".

For random small CPDAGs (n = 4..7, ER, submitted generator), every knowledge
state G0 of the CPDAG, every ordered pair (X, Y) amenable at G0 and every set
Z valid at G0 (GAC, all subsets of the other nodes), check every state G in
the retraction up-set of G0 (K_G a closed subset of K_G0): whenever (X, Y) is
amenable at G, is Z valid at G? A single violation shows that the radius can
be set-specific under the GAC; none supports r_val = amenability radius.
"""

from __future__ import annotations

import itertools
import json
import multiprocessing as mp
from collections import Counter

import numpy as np

from common import OUT, apply_orientations, is_gac_valid_mpdag, k_of

from bkrobust.demo.example import dag_to_cpdag
from bkrobust.mpdag_criterion import is_amenable
from bkrobust.search.space_fixed import enumerate_space_correct
from bkrobust.synth.generators import erdos_renyi_dag


def work(seed):
    rng = np.random.default_rng(seed)
    n = 4 + seed % 4
    dag = erdos_renyi_dag(n, rng, float(rng.uniform(0.4, 0.9)))
    cp = dag_to_cpdag(dag)
    c = Counter()
    ex = []
    if not cp.undirected_edges or len(cp.undirected_edges) > 7:
        return c, ex
    space = enumerate_space_correct(cp)
    for g0 in space:
        kg0 = k_of(cp, g0)
        ups = [g for g in space if k_of(cp, g) <= kg0]
        for x, y in itertools.permutations(g0.nodes, 2):
            if not is_amenable(g0, x, y):
                continue
            others = [v for v in g0.nodes if v not in (x, y)]
            for m in range(len(others) + 1):
                for zz in itertools.combinations(others, m):
                    z = frozenset(zz)
                    if not is_gac_valid_mpdag(g0, x, y, z):
                        continue
                    for g in ups:
                        if g is g0 or not is_amenable(g, x, y):
                            continue
                        c["checks"] += 1
                        if not is_gac_valid_mpdag(g, x, y, z):
                            c["violations"] += 1
                            if len(ex) < 3:
                                ex.append({"cpdag": cp.edge_string(), "G0": g0.edge_string(),
                                           "G": g.edge_string(), "x": x, "y": y, "z": sorted(z)})
    return c, ex


if __name__ == "__main__":
    tot = Counter()
    exs = []
    with mp.Pool(2) as P:
        for c, ex in P.imap_unordered(work, range(400)):
            tot.update(c)
            exs += ex
    out = {"checks": tot["checks"], "violations": tot["violations"], "examples": exs[:5]}
    p = OUT / "e10" / "conjecture_check.json"
    p.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
