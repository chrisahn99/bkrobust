"""The frozen core must be bit-reproducible across runs and PYTHONHASHSEED.

Two bugs of exactly this shape were found in the inherited code: RNG consumed in
frozenset-iteration order, and floating-point summation in set order. Both made
results stable within a process and silently different between runs. This test
extends that regression guard to the new core, as required.
"""

from __future__ import annotations

import subprocess
import sys

HASH_SEEDS = ("0", "1", "12345")


def _run(code: str) -> str:
    outs = set()
    for hs in HASH_SEEDS:
        env = {"PYTHONPATH": "src", "PYTHONHASHSEED": hs, "PATH": "/usr/bin:/bin"}
        res = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, env=env, check=True
        )
        outs.add(res.stdout.strip())
    assert len(outs) == 1, f"non-deterministic across PYTHONHASHSEED: {len(outs)} distinct outputs"
    return next(iter(outs))


def test_space_construction_is_hash_invariant():
    """Element order, cover set and neighbour degrees must not move."""
    code = (
        "from bkrobust.core.spacelib import build_space;"
        "from bkrobust.demo.example import dag_to_cpdag;"
        "from bkrobust.demo.scenario import true_dag;"
        "s = build_space(dag_to_cpdag(true_dag()));"
        "print([g.edge_string() for g in s.elements]);"
        "print(len(s.covers));"
        "print(sorted((g.edge_string(), len(n)) for g, n in s.neighbours.items()))"
    )
    _run(code)


def test_radii_and_distances_are_hash_invariant():
    code = (
        "from bkrobust.core.spacelib import build_space, distances_from, radius;"
        "from bkrobust.core.oracle import is_valid;"
        "from bkrobust.demo.example import dag_to_cpdag;"
        "from bkrobust.demo.meek import apply_orientations;"
        "from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag;"
        "from bkrobust.demo.scenario import true_dag, scenarios, TREATMENT, OUTCOME;"
        "c = dag_to_cpdag(true_dag()); s = build_space(c);"
        "out = [];"
        "\nfor lab, sc in sorted(scenarios().items()):\n"
        "    g0 = apply_orientations(c, sc['knowledge']);\n"
        "    z = frozenset(optimal_adjustment_set_mpdag(g0, TREATMENT, OUTCOME) or ());\n"
        "    d = distances_from(s, g0);\n"
        "    r, w = radius(s, d, lambda g: not is_valid(z, g, TREATMENT, OUTCOME));\n"
        "    out.append((lab, sorted(z), r, w.edge_string(), sorted(d.values())))\n"
        "print(out)"
    )
    _run(code)


def test_bias_stats_are_hash_invariant():
    """Bias is where the two inherited RNG-order bugs actually bit."""
    code = (
        "import numpy as np;"
        "from bkrobust.core.oracle import bias_stats;"
        "from bkrobust.demo.example import dag_to_cpdag;"
        "from bkrobust.demo.scenario import true_dag, TREATMENT, OUTCOME;"
        "c = dag_to_cpdag(true_dag());"
        "b = bias_stats(frozenset({'Age','Smoke'}), c, TREATMENT, OUTCOME,"
        "               np.random.default_rng(0), 25);"
        "print(repr(b))"
    )
    _run(code)
