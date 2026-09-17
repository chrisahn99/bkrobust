"""M1: the knowledge interval, the null radius and the validity radius, with bands.

For every frozen frame row of the swept networks (and the author-declared queries
of the applied DAGs), under each knowledge condition, this computes the analyst
graph, the committed set, the retraction ball over the claims inside the
treatment's chain component, the effects each graph in the ball leaves possible
(semi-local IDA, ``bkrobust.estimation.knowledge_interval``), and from them the
Table 1 procedures at the population and at two sample sizes over twenty seeds:
PLAIN, FAR1, NEAR1, I1, I2, BLANKET and the oracle confounding bound, and (M1b) the
policy rows: the radius stop rule STOP1 and STOP2, and the leave-one-out screens SCREEN
and SCREEN2, whose certificates are computed from observables only. The truth
enters to draw parameters and samples and, after the intervals are computed, to
score them.

Panel B uses the fitted coefficients of the four real-coefficient networks; Panel A
draws semi-synthetic linear-Gaussian parameters. Every choice below is fixed in
``results/interval/PREREGISTRATION.md``.

Writes ``results/interval/rows.csv``, ``replicates.csv.gz``, ``profiles.jsonl`` and
``panel_summary.json``; ``experiments/interval_tables.py`` reads them.

    python experiments/interval_panel.py --workers 7
    python experiments/interval_panel.py --networks asia --out /tmp/smoke
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import itertools
import json
import math
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT / "src", ROOT / "experiments"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import ledger_sweep  # noqa: E402
from bkrobust.benchmarks.audit import is_valid_adjustment_set_gac_dag  # noqa: E402
from bkrobust.benchmarks.describe import parse_file  # noqa: E402
from bkrobust.core.conventions import NO_RETRACTABLE_EDGES  # noqa: E402
from bkrobust.demo.evaluate import optimal_adjustment_set_dag  # noqa: E402
from bkrobust.demo.example import knowledge_to_recover  # noqa: E402
from bkrobust.demo.graph import MPDAG, undirected_components  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.estimation.knowledge_interval import (  # noqa: E402
    Z975,
    Ball,
    BallElement,
    ParentSets,
    chain_component,
    closure_by_components,
    component_graph,
    contains_zero,
    local_claims,
    partial_r2,
    possible_parent_sets,
    regression_beta_se,
    regression_moments,
    retraction_ball,
    y_possible_descendant,
    zero_tolerance,
)
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag  # noqa: E402
from bkrobust.mpdag_criterion.optimal import CommittedSet, committed_adjustment_set  # noqa: E402
from bkrobust.mpdag_criterion.paths import possible_descendants  # noqa: E402

MODELS = ROOT / "results" / "axisa3" / "networks" / "example_models"
FRAME = ROOT / "results" / "frame" / "frame.jsonl"
KNOWLEDGE = ROOT / "results" / "elicit" / "knowledge.json"
LEDGER = ROOT / "results" / "ledger" / "rows.csv"
DECLARED = ROOT / "results" / "frame" / "declared_rows.csv"
OUT = ROOT / "results" / "interval"

PANEL_B = ("ecoli70", "arth150", "magic-niab", "magic-irri")
COEF_RANGE = (0.5, 1.5)
NOISE_VAR = 1.0
N_OBS = (1000, 20000)
SEEDS = 20
MAX_DEPTH = 3
SUBSET_BUDGET = 300
PARENT_CAP = 4096
CONDITIONS = ("CTRL_b0", "CTRL_b1", "CTRL_b2", "D_LLM", "D_LLM_72B", "A_ORACLE")
TABLE1_CONDITIONS = ("CTRL_b0", "CTRL_b1", "CTRL_b2", "D_LLM", "D_LLM_72B")
PROFILE_CONDITIONS = ("A_ORACLE", "CTRL_b0")
# every other condition keeps a profile only where Figure 1 may draw from
PROFILE_TIERS = ("T2", "T3")
PROCS = ("PLAIN", "FAR1", "NEAR1", "I0", "I1", "I2", "I3", "BLANKET")
REPLICATE_PROCS = ("PLAIN", "FAR1", "NEAR1", "I1", "I2", "BLANKET", "CH")
# M1b policy rows, appended so every M1 column keeps its position
POLICY_PROCS = ("SCREEN", "SCREEN2", "STOP1", "STOP2")
REPLICATE_PROCS = (*REPLICATE_PROCS, *POLICY_PROCS)
CENSORED = 4  # r0 code for never, >3 and censored


def stable_seed(*parts: object) -> int:
    """A seed from the first 16 hex digits of the sha256 of the ``:``-joined parts."""
    return int(hashlib.sha256(":".join(str(p) for p in parts).encode()).hexdigest()[:16], 16)


# ----------------------------------------------------------------------------- parameters


def sem_parameters(
    net: str, dag: MPDAG
) -> tuple[dict[tuple[str, str], float], dict[str, float], str]:
    """Coefficients and noise variances: fitted for Panel B, drawn for Panel A."""
    if net in PANEL_B:
        d = json.loads((MODELS / f"{net}.json").read_text())
        w: dict[tuple[str, str], float] = {}
        var: dict[str, float] = {}
        for v, cpd in d["cpds"].items():
            var[v] = float(cpd["variance"][0])
            for p, c in cpd["coefficients"].items():
                if p != "(Intercept)":
                    w[(p, v)] = float(c[0])
        if set(w) != set(dag.directed_edges):
            raise ValueError(f"{net}: fitted coefficients do not match the parsed edges")
        return w, var, "fitted"
    rng = np.random.default_rng(stable_seed("interval-sem-v1", net))
    w = {}
    for e in sorted(dag.directed_edges):
        sign = float(rng.choice(np.array([-1.0, 1.0])))
        w[e] = sign * float(rng.uniform(*COEF_RANGE))
    return w, {v: NOISE_VAR for v in dag.nodes}, "semi-synthetic"


def topological_order(dag: MPDAG) -> list[str]:
    """Kahn's order with sorted ties, so sampling never depends on set iteration."""
    indeg = {v: len(dag.parents(v)) for v in dag.nodes}
    children: dict[str, list[str]] = collections.defaultdict(list)
    for a, b in dag.directed_edges:
        children[a].append(b)
    ready = sorted(v for v, d in indeg.items() if d == 0)
    order: list[str] = []
    while ready:
        v = ready.pop(0)
        order.append(v)
        for c in sorted(children[v]):
            indeg[c] -= 1
            if indeg[c] == 0:
                ready.append(c)
        ready.sort()
    return order


def sample_covariances(
    net: str,
    dag: MPDAG,
    w: dict[tuple[str, str], float],
    var: dict[str, float],
    idx: dict[str, int],
    n: int,
) -> np.ndarray:
    """Sample covariances of ancestral draws of every node, one per seed, ``(SEEDS, p, p)``."""
    order = topological_order(dag)
    parents = {v: sorted(dag.parents(v)) for v in dag.nodes}
    out = np.empty((SEEDS, len(idx), len(idx)))
    for seed in range(SEEDS):
        rng = np.random.default_rng(stable_seed("interval-data-v1", net, n, seed))
        data = np.empty((n, len(idx)))
        for v in order:
            col = rng.normal(0.0, math.sqrt(var[v]), size=n)
            if parents[v]:
                pi = [idx[p] for p in parents[v]]
                col += data[:, pi] @ np.array([w[(p, v)] for p in parents[v]])
            data[:, idx[v]] = col
        out[seed] = np.cov(data, rowvar=False)
    return out


# ----------------------------------------------------------------------------- per row


class Effects:
    """Coefficient and standard error per parent set, at the population and each ``n``."""

    def __init__(
        self, sig_pop: np.ndarray, sig: dict[int, np.ndarray], idx: dict[str, int], x: str, y: str
    ) -> None:
        self.sig_pop, self.sig, self.idx, self.x, self.y = sig_pop, sig, idx, x, y
        self.cache: dict[frozenset[str], dict] = {}

    def of(self, z: frozenset[str]) -> dict:
        """``{"pop": beta, n: (beta_hat, se)}`` for adjustment on ``z``."""
        if z not in self.cache:
            beta, _ = regression_beta_se(self.sig_pop[None], self.idx, self.x, self.y, z, None)
            rec: dict = {"pop": float(beta[0])}
            for n, s in self.sig.items():
                rec[n] = regression_beta_se(s, self.idx, self.x, self.y, z, n)
            self.cache[z] = rec
        return self.cache[z]


def interval_of(eff: Effects, sets: list[frozenset[str]]) -> dict:
    """Population hull, and per ``n`` the band and plug-in hull over seeds."""
    if not sets:
        return {}
    recs = [eff.of(s) for s in sets]
    pops = [r["pop"] for r in recs]
    out: dict = {"pop": (min(pops), max(pops))}
    for n in eff.sig:
        b = np.stack([r[n][0] for r in recs])
        s = np.stack([r[n][1] for r in recs])
        out[n] = {
            "band": ((b - Z975 * s).min(axis=0), (b + Z975 * s).max(axis=0)),
            "plug": (b.min(axis=0), b.max(axis=0)),
        }
    return out


def plain_interval(eff: Effects, z: frozenset[str]) -> dict:
    """The committed set's interval: a point at the population, the CI per seed."""
    r = eff.of(z)
    out: dict = {"pop": (r["pop"], r["pop"])}
    for n in eff.sig:
        b, s = r[n]
        out[n] = {"band": (b - Z975 * s, b + Z975 * s), "plug": (b, b)}
    return out


def first_zero(
    intervals: list[dict | None], key: str, n: int | None, tol: float, n_loc: int
) -> tuple[np.ndarray | int, str]:
    """The null radius over budgets 0..3 from per-budget intervals.

    ``intervals[r]`` is ``None`` when budget ``r`` was not enumerated. Returns the
    radius (an int for the population, an array over seeds otherwise) with codes
    0..3 or :data:`CENSORED`, and a status for the population case.
    """
    if n is None:
        for r, iv in enumerate(intervals):
            if iv is None:
                return CENSORED, f"censored_at_depth_{r}"
            lo, hi = iv["pop"]
            if contains_zero(lo, hi, tol):
                return r, "exact"
        return CENSORED, ("never" if n_loc <= MAX_DEPTH else ">3")
    res = np.full(SEEDS, CENSORED)
    done = np.zeros(SEEDS, dtype=bool)
    for r, iv in enumerate(intervals):
        if iv is None:
            break
        lo, hi = iv[n][key]
        hit = (lo <= tol) & (hi >= -tol) & ~done
        res = np.where(hit, r, res)
        done |= hit
    return res, ""


def name_claims(claims: list[tuple[str, str]], names: dict[str, str]) -> str:
    """``a->b and c->d`` in the network's variable names."""
    return " and ".join(f"{names.get(a, a)}->{names.get(b, b)}" for a, b in claims)


def stop_certified(rec: dict, r: int) -> tuple[bool, str]:
    """Whether the radius stop rule vouches for the committed set at budget ``r``.

    Reads the row's validity radius exactly as M1 stores it. Returns the decision and
    how it was reached, so censored certificates can be counted.
    """
    if rec.get("committed_verdict") not in ("closed_form", "canonical"):
        return False, "no_committed_set"
    v, st = rec.get("r_val"), str(rec.get("r_val_status") or "")
    if v not in ("", None):
        return int(v) > r, "exact"
    if st in ("unreached", NO_RETRACTABLE_EDGES):
        return True, st
    if st.startswith("gt_"):
        return int(st.split("_")[1]) >= r, st
    if st.startswith("censored_at_depth_"):
        depth = int(st.rsplit("_", 1)[1])
        return depth > r, st if depth > r else f"refused_{st}"
    return False, f"refused_{st}"


def result_key(cs: CommittedSet) -> tuple:
    """What the screens compare: the committed set, or the verdict class when there is none."""
    if cs.z is not None:
        return ("set", frozenset(cs.z))
    return (
        "none",
        "y_not_possible_descendant" if cs.verdict == "y_not_possible_descendant" else "no_set",
    )


def process_network(net: str) -> dict:
    """Every row and condition of one network; returns rows, replicates, profiles, summary."""
    t_net = time.perf_counter()
    dag, cpdag, names = ledger_sweep.load(net)
    nodes = list(dag.nodes)
    idx = {v: i for i, v in enumerate(nodes)}
    w, var, params = sem_parameters(net, dag)
    bmat = np.zeros((len(nodes), len(nodes)))
    for (a, b), c in w.items():
        bmat[idx[b], idx[a]] = c
    total = np.linalg.inv(np.eye(len(nodes)) - bmat)
    sig_pop = (total * np.array([var[v] for v in nodes])) @ total.T
    sig = {n: sample_covariances(net, dag, w, var, idx, n) for n in N_OBS}
    t_sample = time.perf_counter() - t_net

    know = json.loads(KNOWLEDGE.read_text())
    k_arm = {
        a: [tuple(e) for e in know[a]["networks"].get(net, {}).get("k", [])]
        for a in ("D_LLM", "D_LLM_72B")
    }
    k_true = [(a, b) if dag.is_directed_edge(a, b) else (b, a) for (a, b) in k_arm["D_LLM"]]
    k_oracle = sorted(knowledge_to_recover(dag, cpdag))
    g0_cache: dict[tuple, MPDAG | None] = {}
    components = undirected_components(cpdag)
    local_closures: dict = {}

    def closure(k: list) -> MPDAG | None:
        key = tuple(sorted(k))
        if key not in g0_cache:
            g0_cache[key] = apply_orientations(cpdag, k) if k else cpdag
        return g0_cache[key]

    ledger = {
        (r["x"], r["y"], r["arm"]): r
        for r in csv.DictReader(LEDGER.open())
        if r["network"] == net and r["rep"] == "0"
    }
    queries = [
        (r["x"], r["y"], "frame", i)
        for i, r in enumerate(json.loads(line) for line in FRAME.open())
        if r["network"] == net
    ]
    declared = {r["network"]: r for r in csv.DictReader(DECLARED.open()) if r["status"] == "ok"}
    if net in declared:
        path = next(p for p in sorted(MODELS.iterdir()) if p.name.startswith(net))
        parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
        ex = [v for v, r in parsed.roles.items() if "exposure" in r]
        oc = [v for v, r in parsed.roles.items() if "outcome" in r]
        queries.append((ex[0], oc[0], "declared", -1))

    counts: collections.Counter = collections.Counter()
    rows: list[dict] = []
    reps: list[dict] = []
    profiles: list[dict] = []
    ps_cache: dict[tuple, ParentSets] = {}

    for x, y, source, frame_index in queries:
        comp = chain_component(cpdag, x)
        pa_ext = frozenset(cpdag.parents(x))
        base = component_graph(cpdag, comp) if comp is not None else None
        tau = float(total[idx[y], idx[x]])
        tol = zero_tolerance(sig_pop, idx, x, y)
        eff = Effects(sig_pop, sig, idx, x, y)

        def parent_sets(local: MPDAG | None, x: str = x) -> ParentSets:
            if local is None:
                return ParentSets((frozenset(),), (), 1, False)
            key = (x, local.key())
            if key not in ps_cache:
                ps_cache[key] = possible_parent_sets(local, x, PARENT_CAP)
            return ps_cache[key]

        # ---- the class
        ps_class = parent_sets(base)
        class_sets = [pa_ext | p for p in ps_class.local]
        if comp is None:
            toward = [y in possible_descendants(cpdag, x)]
        else:
            toward = [y_possible_descendant(cpdag, comp, g, x, y) for g in ps_class.closed]
        blanket = interval_of(eff, class_sets) if not ps_class.censored else {}
        common = {
            "panel": "B" if net in PANEL_B else "A",
            "network": net,
            "tier": ledger_sweep.TIER.get(net, "T1"),
            "params": params,
            "source": source,
            "frame_index": frame_index,
            "x": x,
            "y": y,
            "x_name": names.get(x, x),
            "y_name": names.get(y, y),
            "tau": tau,
            "tol": tol,
            "component_size": len(comp) if comp else 0,
            "n_siblings": len(cpdag.neighbors(x)),
            "class_parent_sets": len(class_sets),
            "class_censored": int(ps_class.censored),
        }
        if ps_class.censored:
            counts["class_parent_sets_censored"] += 1
            for cond in CONDITIONS:
                rows.append({**common, "condition": cond, "status": "class_parent_sets_censored"})
            continue
        c_lo, c_hi = blanket["pop"]
        restricted = [s for s, t in zip(class_sets, toward, strict=True) if t]
        if restricted:
            r_lo, r_hi = interval_of(eff, restricted)["pop"]
            r_width = r_hi - r_lo
        else:
            r_lo = r_hi = float("nan")
            r_width = 0.0
        width = c_hi - c_lo
        informative = width > tol
        common.update(
            {
                "class_lo": c_lo,
                "class_hi": c_hi,
                "class_width": width,
                "informative": int(informative),
                "restricted_lo": r_lo,
                "restricted_hi": r_hi,
                "restricted_share": (r_width / width) if informative else float("nan"),
                "dispute_share": (1 - r_width / width) if informative else float("nan"),
                "class_covers_tau": int(c_lo - tol <= tau <= c_hi + tol),
            }
        )

        for cond in CONDITIONS:
            rec: dict = {**common, "condition": cond}
            if source == "declared" and cond != "A_ORACLE":
                continue
            # ---- the knowledge for this row
            b_target = int(cond[-1]) if cond.startswith("CTRL") else None
            reversed_claims: list = []
            if cond.startswith("CTRL"):
                k_loc_true = local_claims(cpdag, k_true, comp)
                k_loc = list(k_loc_true)
                if b_target:
                    if len(k_loc_true) < b_target:
                        rec["status"] = "b_unreachable_too_few_local_claims"
                    else:
                        ok = []
                        for sub in itertools.combinations(range(len(k_loc_true)), b_target):
                            flipped = [
                                (e[1], e[0]) if i in sub else e for i, e in enumerate(k_loc_true)
                            ]
                            if base is not None and apply_orientations(base, flipped) is not None:
                                ok.append(sub)
                        if not ok:
                            rec["status"] = "b_unreachable_no_consistent_reversal"
                        else:
                            rng = np.random.default_rng(
                                stable_seed("interval-reverse-v1", net, x, y, b_target)
                            )
                            pick = ok[int(rng.integers(len(ok)))]
                            reversed_claims = [(k_loc_true[i][1], k_loc_true[i][0]) for i in pick]
                            k_loc = sorted(
                                (e[1], e[0]) if i in pick else e for i, e in enumerate(k_loc_true)
                            )
                    if "status" in rec:
                        counts[f"{cond}:{rec['status']}"] += 1
                        rows.append(rec)
                        continue
                loc_set = set(k_loc_true)
                k_full = [e for e in k_true if e not in loc_set] + k_loc
            elif cond == "A_ORACLE":
                k_full = list(k_oracle)
                k_loc = local_claims(cpdag, k_full, comp)
            else:
                k_full = list(k_arm[cond])
                k_loc = local_claims(cpdag, k_full, comp)
            g0 = closure(k_full)
            rec.update(
                {
                    "b_target": b_target,
                    "b_realised": len(reversed_claims) if b_target is not None else None,
                    "reversed": name_claims(reversed_claims, names),
                    "n_k": len(k_full),
                    "n_k_loc": len(k_loc),
                    "claims_loc": name_claims(k_loc, names),
                    "n_false_local": sum(1 for a, b in k_loc if not dag.is_directed_edge(a, b)),
                    "n_false_total": sum(1 for a, b in k_full if not dag.is_directed_edge(a, b)),
                }
            )
            if g0 is None:
                rec["status"] = "knowledge_inconsistent"
                counts[f"{cond}:knowledge_inconsistent"] += 1
                rows.append(rec)
                continue

            # ---- the committed set and the validity radius, observables only
            cs = committed_adjustment_set(g0, x, y)
            rec["committed_verdict"] = cs.verdict
            rec["z"] = " ".join(sorted(names.get(v, v) for v in cs.z)) if cs.z else ""
            rec["z_size"] = len(cs.z) if cs.z is not None else None
            arm = {"CTRL_b0": "A_TRUE"}.get(cond, cond)
            lr = ledger.get((x, y, arm)) if source == "frame" else None
            from_ledger = cond in ("CTRL_b0", "D_LLM", "D_LLM_72B", "A_ORACLE")
            if from_ledger and lr is not None:
                rec["r_val"] = lr.get("r_claim", "")
                rec["r_val_status"] = lr.get("r_claim_status") or lr.get("status", "")
                rec["r_val_source"] = "ledger"
                rec["ledger_verdict_match"] = int(lr.get("committed_verdict") == cs.verdict)
            elif cs.z is not None:
                r_claim, r_status, _, _ = ledger_sweep.breaking_depths(cpdag, k_full, x, y, cs.z)
                rec["r_val"] = "" if r_claim is None else r_claim
                rec["r_val_status"] = r_status
                rec["r_val_source"] = "computed"
            else:
                rec["r_val"] = ""
                rec["r_val_status"] = cs.verdict
                rec["r_val_source"] = "computed"

            # ---- the ball
            if base is None:
                ball = Ball((BallElement((), None),), MAX_DEPTH, "no_component")
            else:
                ball = retraction_ball(base, k_loc, MAX_DEPTH, SUBSET_BUDGET)
            rec["ball_status"] = ball.status
            rec["ball_depth_done"] = ball.depth_done
            rec["ball_subsets"] = len(ball.elements) - 1
            if ball.status != "exhaustive":
                counts[f"{cond}:ball_{ball.status}"] += 1
            el_sets: list[list[frozenset[str]]] = []
            censored = False
            for el in ball.elements:
                ps = parent_sets(el.graph) if el.graph is not None or base is None else None
                if ps is None:
                    raise RuntimeError(f"{net} {x} {y} {cond}: retraction left an inconsistency")
                censored |= ps.censored
                el_sets.append([pa_ext | p for p in ps.local])
            if censored:
                rec["status"] = "parent_sets_censored"
                counts[f"{cond}:parent_sets_censored"] += 1
                rows.append(rec)
                continue
            rec["status"] = "ok"
            el_iv = [interval_of(eff, s) for s in el_sets]
            n_loc = len(k_loc)

            def union(select: list[int], el_sets: list = el_sets, eff: Effects = eff) -> dict:
                sets = sorted({s for i in select for s in el_sets[i]}, key=sorted)
                return interval_of(eff, sets)

            budgets: list[dict | None] = []
            for r in range(MAX_DEPTH + 1):
                r_eff = min(r, n_loc)
                if r_eff > ball.depth_done:
                    budgets.append(None)
                    continue
                sel = [i for i, el in enumerate(ball.elements) if len(el.retracted) <= r_eff]
                budgets.append(union(sel))
            near = [i for i, (a, b) in enumerate(k_loc) if x in (a, b) or y in (a, b)]
            far = [i for i in range(n_loc) if i not in near]
            single = {
                el.retracted[0]: i for i, el in enumerate(ball.elements) if len(el.retracted) == 1
            }
            procs: dict[str, dict | None] = {
                "I0": budgets[0],
                "I1": budgets[1],
                "I2": budgets[2],
                "I3": budgets[3],
                "BLANKET": blanket,
                "FAR1": union([0] + [single[i] for i in far])
                if ball.depth_done >= 1 or not far
                else None,
                "NEAR1": union([0] + [single[i] for i in near])
                if ball.depth_done >= 1 or not near
                else None,
            }
            if cs.z is not None:
                procs["PLAIN"] = plain_interval(eff, cs.z)
                rec["plain_kind"] = "committed"
            else:
                procs["PLAIN"] = budgets[0]
                rec["plain_kind"] = "fallback_budget0"

            for p in PROCS:
                iv = procs.get(p)
                rec[f"pop_{p}_lo"] = iv["pop"][0] if iv else None
                rec[f"pop_{p}_hi"] = iv["pop"][1] if iv else None

            # ---- radii
            r0_pop, r0_status = first_zero(budgets, "", None, tol, n_loc)
            rec["r0_pop"] = r0_pop
            rec["r0_pop_status"] = r0_status
            if r0_status == "exact":
                zsets = {
                    s
                    for i, el in enumerate(ball.elements)
                    if len(el.retracted) <= min(r0_pop, n_loc)
                    for s in el_sets[i]
                }
                rec["r0_structural_zero"] = int(any(abs(eff.of(s)["pop"]) <= tol for s in zsets))
            for n in N_OBS:
                for key in ("band", "plug"):
                    arr, _ = first_zero(budgets, key, n, tol, n_loc)
                    rec[f"r0_{key}_n{n}"] = ";".join(str(int(v)) for v in arr)

            # ---- the oracle confounding bound (feasibility M2), scored after the freeze
            ch: dict | None = None
            if cs.z is not None:
                z = frozenset(cs.z)
                if y in dag.descendants(x):
                    v_true: set | None = optimal_adjustment_set_dag(dag, x, y)
                elif y in dag.parents(x):
                    v_true = None
                else:
                    v_true = set(dag.parents(x))
                rec["z_valid_at_truth"] = int(is_valid_adjustment_set_gac_dag(dag, x, y, z))
                rec["z_has_true_descendant"] = int(bool(z & dag.descendants(x)))
                if v_true is None:
                    rec["ch_defined"] = 0
                    rec["ch_reason"] = "no_valid_set_at_truth"
                else:
                    u = frozenset(v_true) - z - {x, y}
                    defined = is_valid_adjustment_set_gac_dag(dag, x, y, z | u)
                    rec["ch_defined"] = int(defined)
                    rec["ch_reason"] = "omitted_variable" if defined else "not_omitted_variable"
                    r2x = partial_r2(sig_pop, idx, x, u, z)
                    r2y = partial_r2(sig_pop, idx, y, u, z | {x})
                    rec["ch_r2_x"], rec["ch_r2_y"] = r2x, r2y
                    if r2x < 1:
                        k = math.sqrt(r2y * r2x / (1 - r2x))
                        adj = math.sqrt((1 - r2y) / (1 - r2x))
                        beta, s_y, s_x = regression_moments(sig_pop[None], idx, x, y, z)
                        bound = k * math.sqrt(max(s_y[0], 0.0) / s_x[0])
                        ch = {"pop": (beta[0] - bound, beta[0] + bound)}
                        for n in N_OBS:
                            bh, sy, sx = regression_moments(sig[n], idx, x, y, z)
                            _, se = regression_beta_se(sig[n], idx, x, y, z, n)
                            half = k * np.sqrt(np.maximum(sy, 0.0) / sx) + Z975 * se * adj
                            ch[n] = {"band": (bh - half, bh + half)}
                        rec["pop_CH_lo"], rec["pop_CH_hi"] = ch["pop"]

            # ---- M1b: the radius stop rule and the leave-one-out screens, observables only
            if cond in TABLE1_CONDITIONS and source == "frame":
                key0 = result_key(cs)
                z0 = frozenset(cs.z) if cs.z is not None else None
                n_k = len(k_full)
                single_change = single_change_z_valid = False
                for i in range(n_k):
                    g1 = closure_by_components(
                        cpdag,
                        [e for j, e in enumerate(k_full) if j != i],
                        components,
                        local_closures,
                    )
                    if g1 is None:
                        raise RuntimeError(f"{net} {x} {y} {cond}: retraction inconsistent")
                    if result_key(committed_adjustment_set(g1, x, y)) != key0:
                        single_change = True
                        if z0 is not None and is_gac_valid_mpdag(g1, x, y, z0):
                            single_change_z_valid = True
                pair_change = pair_change_z_valid = False
                pair_status = "not_needed" if single_change else "exhaustive"
                if not single_change and n_k >= 2:
                    if n_k + math.comb(n_k, 2) > SUBSET_BUDGET:
                        pair_status = "censored"
                    else:
                        for i, j in itertools.combinations(range(n_k), 2):
                            g2 = closure_by_components(
                                cpdag,
                                [e for t, e in enumerate(k_full) if t not in (i, j)],
                                components,
                                local_closures,
                            )
                            if g2 is None:
                                raise RuntimeError(f"{net} {x} {y} {cond}: pair inconsistent")
                            if result_key(committed_adjustment_set(g2, x, y)) != key0:
                                pair_change = True
                                if z0 is not None and is_gac_valid_mpdag(g2, x, y, z0):
                                    pair_change_z_valid = True
                                break
                cert1, how1 = stop_certified(rec, 1)
                cert2, how2 = stop_certified(rec, 2)
                trig2 = single_change or pair_change or pair_status == "censored"
                sets_b1 = {
                    s
                    for i, el in enumerate(ball.elements)
                    if len(el.retracted) <= 1
                    for s in el_sets[i]
                }
                sets_b2 = {
                    s
                    for i, el in enumerate(ball.elements)
                    if len(el.retracted) <= 2
                    for s in el_sets[i]
                }
                rval_exact = rec.get("r_val") not in ("", None)
                rec.update(
                    {
                        "stop1_certified": int(cert1),
                        "stop1_how": how1,
                        "stop2_certified": int(cert2),
                        "stop2_how": how2,
                        "screen_triggered": int(single_change),
                        "screen_change_z_valid": int(single_change_z_valid),
                        "screen2_triggered": int(trig2),
                        "screen2_pair_status": pair_status,
                        "screen2_pair_change_z_valid": int(pair_change_z_valid),
                        "interaction": int(
                            z0 is not None
                            and (
                                (rval_exact and int(rec["r_val"]) == 2)
                                or (not single_change and pair_change)
                            )
                        ),
                        "sets_b1_equal_b2": int(sets_b1 == sets_b2),
                    }
                )
                procs["STOP1"] = procs["PLAIN"] if cert1 else budgets[1]
                procs["STOP2"] = procs["PLAIN"] if cert2 else budgets[2]
                procs["SCREEN"] = budgets[1] if single_change else procs["PLAIN"]
                procs["SCREEN2"] = budgets[2] if trig2 else procs["PLAIN"]
                for p in POLICY_PROCS:
                    iv = procs.get(p)
                    rec[f"pop_{p}_lo"] = iv["pop"][0] if iv else None
                    rec[f"pop_{p}_hi"] = iv["pop"][1] if iv else None

            rows.append(rec)

            # ---- replicates for Table 1 (informative rows only)
            if informative and cond in TABLE1_CONDITIONS and source == "frame":
                key = {"network": net, "x": x, "y": y, "condition": cond}
                pop_rec = {**key, "n": 0, "seed": 0, "tau": tau}
                for p in REPLICATE_PROCS:
                    iv = ch if p == "CH" else procs.get(p)
                    pop_rec[f"{p}_lo"] = iv["pop"][0] if iv else ""
                    pop_rec[f"{p}_hi"] = iv["pop"][1] if iv else ""
                reps.append(pop_rec)
                for n in N_OBS:
                    for seed in range(SEEDS):
                        rr = {**key, "n": n, "seed": seed, "tau": tau}
                        for p in REPLICATE_PROCS:
                            iv = ch if p == "CH" else procs.get(p)
                            if iv:
                                lo, hi = iv[n]["band"]
                                rr[f"{p}_lo"] = float(lo[seed])
                                rr[f"{p}_hi"] = float(hi[seed])
                            else:
                                rr[f"{p}_lo"] = rr[f"{p}_hi"] = ""
                        reps.append(rr)

            # ---- per-budget profile for Table 2 and Figure 1
            wanted = cond in PROFILE_CONDITIONS or common["tier"] in PROFILE_TIERS
            if wanted and (cs.z is not None or source == "declared"):
                n_top = N_OBS[-1]
                prof: dict = {
                    "network": net,
                    "panel": common["panel"],
                    "params": params,
                    "source": source,
                    "frame_index": frame_index,
                    "condition": cond,
                    "x": names.get(x, x),
                    "y": names.get(y, y),
                    "x_id": x,
                    "y_id": y,
                    "tau": tau,
                    "n_k": len(k_full),
                    "n_k_loc": n_loc,
                    "n_false_local": rec["n_false_local"],
                    "n_false_total": rec["n_false_total"],
                    "claims": [[names.get(a, a), names.get(b, b)] for a, b in k_full],
                    "claims_local": [[names.get(a, a), names.get(b, b)] for a, b in k_loc],
                    "committed_verdict": cs.verdict,
                    "z": sorted(names.get(v, v) for v in cs.z) if cs.z is not None else None,
                    "r_val": rec.get("r_val"),
                    "r_val_status": rec.get("r_val_status"),
                    "r0_pop": r0_pop,
                    "r0_pop_status": r0_status,
                    "r0_band_seeds": [int(v) for v in rec[f"r0_band_n{n_top}"].split(";")],
                    "informative": int(informative),
                    "class_pop": [c_lo, c_hi],
                }
                if cs.z is not None:
                    b_hat, s_hat = eff.of(frozenset(cs.z))[n_top]
                    prof["estimate"] = float(b_hat[0])
                    prof["estimate_ci"] = [
                        float(b_hat[0] - Z975 * s_hat[0]),
                        float(b_hat[0] + Z975 * s_hat[0]),
                    ]
                    prof["estimate_pop"] = eff.of(frozenset(cs.z))["pop"]
                budget_rows = []
                prev: tuple[float, float] | None = None
                breaking = None
                r0_seed0 = int(rec[f"r0_band_n{n_top}"].split(";")[0])
                for r in range(MAX_DEPTH + 1):
                    iv = budgets[r]
                    if iv is None:
                        budget_rows.append({"r": r, "enumerated": False})
                        continue
                    r_eff = min(r, n_loc)
                    sel = [i for i, el in enumerate(ball.elements) if len(el.retracted) <= r_eff]
                    los = [(float(el_iv[i][n_top]["band"][0][0]), i) for i in sel if el_iv[i]]
                    his = [(float(el_iv[i][n_top]["band"][1][0]), i) for i in sel if el_iv[i]]
                    lo_w = min(los)[1]
                    hi_w = max(his)[1]
                    band = (float(iv[n_top]["band"][0][0]), float(iv[n_top]["band"][1][0]))
                    budget_rows.append(
                        {
                            "r": r,
                            "enumerated": True,
                            "pop": [float(iv["pop"][0]), float(iv["pop"][1])],
                            "band_seed0": list(band),
                            "witness_lo": [
                                [names.get(a, a), names.get(b, b)]
                                for a, b in (k_loc[j] for j in ball.elements[lo_w].retracted)
                            ],
                            "witness_hi": [
                                [names.get(a, a), names.get(b, b)]
                                for a, b in (k_loc[j] for j in ball.elements[hi_w].retracted)
                            ],
                        }
                    )
                    if r == r0_seed0 and r > 0 and prev is not None:
                        side = 0 if prev[0] > 0 else 1
                        for i in sel:
                            lo_i = float(el_iv[i][n_top]["band"][0][0])
                            hi_i = float(el_iv[i][n_top]["band"][1][0])
                            if (side == 0 and lo_i <= tol) or (side == 1 and hi_i >= -tol):
                                breaking = [
                                    [names.get(a, a), names.get(b, b)]
                                    for a, b in (k_loc[j] for j in ball.elements[i].retracted)
                                ]
                                break
                    prev = band
                prof["budgets"] = budget_rows
                prof["breaking_claims_seed0"] = breaking
                profiles.append(prof)
        counts["rows"] += 1

    summary = {
        "network": net,
        "panel": "B" if net in PANEL_B else "A",
        "params": params,
        "nodes": len(nodes),
        "edges": len(dag.directed_edges),
        "undirected": len(cpdag.undirected_edges),
        "sem_seed": None if net in PANEL_B else stable_seed("interval-sem-v1", net),
        "abs_coef_min": min(abs(c) for c in w.values()) if w else None,
        "abs_coef_max": max(abs(c) for c in w.values()) if w else None,
        "max_variance": float(sig_pop.diagonal().max()),
        "seconds_sampling": round(t_sample, 1),
        "seconds": round(time.perf_counter() - t_net, 1),
        "counts": dict(counts),
        "k_oracle": len(k_oracle),
        "k_true": len(k_true),
    }
    return {"rows": rows, "reps": reps, "profiles": profiles, "summary": summary}


def main(argv: list[str] | None = None) -> None:
    """Run every network, in parallel, and write the outputs."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--networks", default="")
    ap.add_argument("--out", default="")
    args = ap.parse_args(argv)
    out = Path(args.out) if args.out else OUT
    out.mkdir(parents=True, exist_ok=True)
    nets = sorted({r["network"] for r in csv.DictReader(LEDGER.open())})
    if args.networks:
        nets = [n for n in nets if n in set(args.networks.split(","))]
    # the heaviest first so the pool balances
    nets.sort(key=lambda n: (n != "arth150", n))
    t0 = time.perf_counter()
    results = []
    with mp.get_context("spawn").Pool(max(1, args.workers)) as pool:
        for res in pool.imap_unordered(process_network, nets):
            s = res["summary"]
            print(f"  {s['network']:16s} {s['seconds']:7.1f} s  {s['counts']}", flush=True)
            results.append(res)
    results.sort(key=lambda r: r["summary"]["network"])

    rows = [row for r in results for row in r["rows"]]
    fields: list[str] = []
    for row in rows:
        for k in row:
            if k not in fields:
                fields.append(k)
    with (out / "rows.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=fields)
        wr.writeheader()
        wr.writerows(rows)
    reps = [row for r in results for row in r["reps"]]
    with gzip.open(out / "replicates.csv.gz", "wt", newline="") as fh:
        rep_fields = list(reps[0]) if reps else []
        wr = csv.DictWriter(fh, fieldnames=rep_fields)
        wr.writeheader()
        wr.writerows(reps)
    with (out / "profiles.jsonl").open("w") as fh:
        for r in results:
            for p in r["profiles"]:
                fh.write(json.dumps(p) + "\n")
    summary = {
        "seconds": round(time.perf_counter() - t0, 1),
        "networks": [r["summary"] for r in results],
        "settings": {
            "panel_b": list(PANEL_B),
            "coef_range_abs": list(COEF_RANGE),
            "coef_sign": "uniform on {-1, +1}, drawn before the magnitude, sorted edge order",
            "noise_var": NOISE_VAR,
            "sem_seed_rule": "first 16 hex of sha256('interval-sem-v1:<network>')",
            "data_seed_rule": "first 16 hex of sha256('interval-data-v1:<network>:<n>:<seed>')",
            "reverse_seed_rule": "first 16 hex of sha256('interval-reverse-v1:<net>:<x>:<y>:<b>')",
            "n_obs": list(N_OBS),
            "seeds": SEEDS,
            "max_depth": MAX_DEPTH,
            "subset_budget": SUBSET_BUDGET,
            "parent_candidate_cap": PARENT_CAP,
            "conditions": list(CONDITIONS),
        },
        "rows": len(rows),
        "replicates": len(reps),
    }
    (out / "panel_summary.json").write_text(json.dumps(summary, indent=1))
    print(f"{len(rows)} rows, {len(reps)} replicates, {summary['seconds']} s")


if __name__ == "__main__":
    main()
