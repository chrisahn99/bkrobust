"""Frozen shared core for the synthetic-ensemble and search-heuristic session.

This package is the single normative source for:

* the **radius convention** (:mod:`bkrobust.core.conventions`),
* the **validity / optimality / bias oracle** (:mod:`bkrobust.core.oracle`),
* **space construction and distance** (:mod:`bkrobust.core.spacelib`),
* the **instance schema** (:mod:`bkrobust.core.instance`),
* the **results contract** (:mod:`bkrobust.core.resultsio`).

It is frozen: downstream modules consume these signatures and may not change
them. Changes are made by the orchestrator, who then re-runs affected work and
records why.

It deliberately reuses the primitives in :mod:`bkrobust.demo` -- MPDAG, Meek's
rules, DAG-extension enumeration, adjustment-set validity, the linear-Gaussian
SEM -- rather than reimplementing them. Those were cross-checked against
independent implementations (networkx d-separation and chordality, brute-force
Markov equivalence classes over 300 random DAGs) and reimplementation would
throw that validation away.
"""

from __future__ import annotations
