"""Rebuild the published-network states of the paper's release records.

Rows of ``release/per_query_records.csv`` whose knowledge is in
``release/knowledge_states.json`` (the ten elicited conditions) are rebuilt:
CPDAG and data-generating DAG from the submitted loaders, K from the release,
checked against the row's ``k_sha256``, G0 = Meek(CPDAG, K).
"""

from __future__ import annotations

import csv
import hashlib
import json
from functools import lru_cache

from common import RELEASE, apply_orientations

from bkrobust.robustness import real_survival as rs


@lru_cache(maxsize=None)
def knowledge():
    return json.loads((RELEASE / "knowledge_states.json").read_text())


@lru_cache(maxsize=None)
def network(name):
    return rs.load_network(name)


def k_sha(k):
    return hashlib.sha256(json.dumps(sorted(map(list, k))).encode()).hexdigest()


def rows(arms=None):
    with open(RELEASE / "per_query_records.csv") as f:
        out = list(csv.DictReader(f))
    if arms is not None:
        out = [r for r in out if r["arm"] in arms]
    return out


_G0 = {}


def state(arm, net):
    """(dag, cpdag, k, g0, sha16) or None if the arm/network is not released."""
    key = (arm, net)
    if key in _G0:
        return _G0[key]
    ks = knowledge()
    if arm not in ks or net not in ks[arm]["networks"]:
        _G0[key] = None
        return None
    k = [tuple(e) for e in ks[arm]["networks"][net]["k"]]
    dag, cpdag = network(net)
    g0 = apply_orientations(cpdag, k) if k else cpdag
    _G0[key] = (dag, cpdag, k, g0, k_sha(k)[:16])
    return _G0[key]
