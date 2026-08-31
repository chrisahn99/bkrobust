"""
Robustness sweep against Limitation 5.1: does the O*-locality mechanism survive
denser graphs and larger p? (Herman et al. arXiv:2503.17037 dense-graph caveat.)

Usage: python sweep_density.py [n_per_cell] [out.json]
"""
import json, sys
from multiprocessing import Pool
import numpy as np
from run_linear import analyse_one

CELLS = [(p, d) for p in (6, 8, 10) for d in (1.5, 2.5, 4.0, 6.0)]


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    out = sys.argv[2] if len(sys.argv) > 2 else "../results/density_sweep.json"
    jobs, tags = [], []
    base = 3_000_000
    for ci, (p, d) in enumerate(CELLS):
        for k in range(n * 4):
            jobs.append((base + ci * 100_000 + k, p, d))
            tags.append((p, d))
    with Pool(16) as pool:
        res = list(pool.imap(analyse_one, jobs, chunksize=4))
    cells = {}
    for (p, d), r in zip(tags, res):
        if r is None:
            continue
        key = f"p{p}_deg{d}"
        cells.setdefault(key, [])
        if len(cells[key]) < n:
            cells[key].append(r)
    json.dump(cells, open(out, "w"))
    for k, v in cells.items():
        print(k, len(v))


if __name__ == "__main__":
    main()
