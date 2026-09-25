"""Split ``optimal_set_undefined`` into "O* is None" vs "O* is the empty set".

Context (see scratchpad BRIEF.md and the finding this script chases): in
``results/axis_robustness_llm/analysis_units.csv`` (5,400 units), 1,542 units
carry ``status == "optimal_set_undefined"``. That status is produced by
:func:`bkrobust.robustness.real_survival.commit_z_star`, whose test is
``if not z:`` -- true both when
:func:`bkrobust.demo.evaluate.optimal_adjustment_set_mpdag` returns ``None``
(no DAG extension of G0 agrees on a single O*, i.e. genuinely unidentified)
and when it returns ``set()`` (O* is identified as empty). Those two cases
were conflated, and ``commit_z_star`` never even reaches its later GAC
validity check for the ``set()`` case, because ``if not z`` fires on ``z ==
set()`` too. This script splits the 1,542 units back apart, and for every
unit where O* agrees to the empty set across DAG extensions, separately
checks whether that empty set is actually GAC-valid at G0 -- **it is not
always**, see "An important nuance" below.

Method
------
For each unit's ``(condition, network)`` we rebuild G0 = Meek(cpdag, K) once
(the elicited K read from ``results/elicit/knowledge.json``, exactly as
:func:`bkrobust.robustness.llm_survival.build_state` does) and reuse it, and
its DAG-extension enumeration, across every query of that group.

Why enumeration is not the problem here: ``build_g0`` already screens every
G0 to at most ``MAX_G0_UNDIRECTED_FOR_EXTENSIONS = 12`` undirected edges
before a unit can reach ``commit_z_star`` at all (see
``bkrobust.benchmarks.measure``), so ``enumerate_dag_extensions(G0)`` is
bounded (<= 4096 extensions) and empirically finishes in well under a second.
We still cache it once per ``(condition, network)`` group (not once per unit)
and apply a per-unit wall-clock budget (``PER_UNIT_TIME_LIMIT_S``) around the
per-extension ``optimal_adjustment_set_dag`` comparison as a safety net -- a
unit that exceeds it is recorded ``"undetermined"`` rather than blocking the
run.

An important nuance found while building this: O*-agreement does not imply
validity
---------------------------------------------------------------------------
The initial plan was to use the polynomial certificate
``bkrobust.mpdag_criterion.is_valid_mpdag(g0, x, y, frozenset())`` to *skip*
enumeration for units where it is False, on the reasoning "O*_D is always a
valid adjustment set of D, so if all extensions agree O*_D = empty, empty
must be valid everywhere". That reasoning is wrong, and a concrete
counterexample from this very dataset shows it: ``D_LLM_32B`` /
``Kampen_2014`` / ``(AIS, AFF)`` has G0 with 2 DAG extensions; both agree
``optimal_adjustment_set_dag(d, x, y) == set()``, yet the enumeration oracle
(:func:`bkrobust.core.oracle.is_valid`) says the empty set is **not** a valid
adjustment set in either extension, and
``bkrobust.mpdag_criterion.is_valid_mpdag`` agrees (False). The
``optimal_adjustment_set_dag`` formula (``pa(cn) \\ (cn u {x})``) is only
*proved* valid under conditions ``commit_z_star``'s pipeline never checks
before declaring O* "identified" by cross-extension agreement -- agreement of
the structural formula is not the same as the formula being valid. So this
script does **not** use that shortcut for the None/empty split (it enumerates
directly, which is cheap here regardless); it uses
``mpdag_criterion.is_valid_mpdag`` for its intended, sound purpose instead --
as a polynomial **validity check** on the empty set, once O* has already been
found to agree to empty by enumeration. Units land in three buckets instead
of two:

* ``"none"`` -- O* is genuinely unidentified (extensions disagree, or G0 has
  no extension at all).
* ``"empty_valid"`` -- O* agrees to the empty set across extensions, and that
  empty set is GAC-valid at G0. This is the case the paper can safely treat
  as "no adjustment needed"; r_val is computed for these.
* ``"empty_invalid"`` -- O* agrees to the empty set across extensions, but it
  is **not** valid at G0. Had ``commit_z_star`` not short-circuited on
  ``not z``, this unit would have gone on to fail its own validity check and
  land in ``"z_invalid_at_g0"``, not become a usable committed Z. These units
  are not "restorable" empty-valid units and are reported separately.

Ordering and budget
--------------------
Groups whose condition is one of the eight real-naming conditions (the
paper's pooled panel) are processed first; the scrambled-naming control
groups follow. A wall-clock budget (``TOTAL_BUDGET_S``, ~55 minutes) stops
the run before it would overrun; if hit, whatever has been processed so far
is written out and the run stops cleanly.

Nothing here calls an LLM or any external API -- it is pure local graph
computation over already-cached network files and the already-written
elicitation bundle, so it is expected to be fast (minutes, not the full
budget).

Usage::

    PYTHONPATH=src nohup .venv/bin/python experiments/llm_status_split_v2.py \
        > results/axis_robustness_llm_v2/status_split.stdout.log 2>&1 &
"""

from __future__ import annotations

import csv
import json
import time
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
UNITS_CSV = ROOT / "results/axis_robustness_llm/analysis_units.csv"
OUT_DIR = ROOT / "results/axis_robustness_llm_v2"
OUT_JSON = OUT_DIR / "status_split.json"
LOG_PATH = OUT_DIR / "status_split.progress.log"

PER_UNIT_TIME_LIMIT_S = 20.0
R_VAL_TIME_LIMIT_S = 60.0
TOTAL_BUDGET_S = 55 * 60.0

_START = time.time()


def _log(msg: str) -> None:
    line = f"[{time.time() - _START:8.1f}s] {msg}"
    print(line, flush=True)
    with LOG_PATH.open("a") as fh:
        fh.write(line + "\n")


def classify_group(g0, group_units, is_valid_mpdag, enumerate_dag_extensions, optimal_adjustment_set_dag):
    """Classify every (x, y) unit of one (condition, network) group.

    Returns a list of per-unit dicts (without condition/network/naming, the
    caller adds those) and summary counters for the log line.
    """
    t_enum0 = time.time()
    try:
        extensions = enumerate_dag_extensions(g0)
        enum_error = None
    except Exception as exc:  # noqa: BLE001
        extensions, enum_error = None, repr(exc)
    enum_elapsed = time.time() - t_enum0

    out = []
    n_timeout = 0
    for u in group_units:
        x, y = u["x"], u["y"]
        t0 = time.time()

        if enum_error is not None:
            classification, method = "undetermined", f"enumeration_error:{enum_error}"
        elif not extensions:
            classification, method = "none", "no_extensions"
        else:
            sets = []
            timed_out = False
            for d in extensions:
                if time.time() - t0 > PER_UNIT_TIME_LIMIT_S:
                    timed_out = True
                    break
                sets.append(frozenset(optimal_adjustment_set_dag(d, x, y)))
            if timed_out:
                classification, method = "undetermined", "enumeration_timeout"
                n_timeout += 1
            else:
                first = sets[0]
                if not all(s == first for s in sets[1:]):
                    classification, method = "none", "enumeration_disagree"
                elif len(first) == 0:
                    try:
                        valid = is_valid_mpdag(g0, x, y, frozenset())
                    except Exception as exc:  # noqa: BLE001
                        classification, method = "undetermined", f"validity_check_error:{exc!r}"
                    else:
                        classification = "empty_valid" if valid else "empty_invalid"
                        method = "enumeration+validity"
                else:
                    # Unexpected: a unit with status optimal_set_undefined
                    # implies commit_z_star saw a falsy z, so an agreed
                    # nonempty O* here is a genuine anomaly worth flagging.
                    classification, method = "anomaly_nonempty_agree", "enumeration"

        elapsed = time.time() - t0
        out.append({
            "x": x, "y": y, "naming": u.get("naming"),
            "classification": classification, "method": method,
            "elapsed_s": round(elapsed, 4),
        })

    return out, {
        "n_extensions": len(extensions) if extensions is not None else None,
        "enum_elapsed_s": round(enum_elapsed, 3),
        "n_timeout": n_timeout,
        "enum_error": enum_error,
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_PATH.write_text("")
    _log("starting llm_status_split_v2")

    from bkrobust.core.conventions import UNREACHED
    from bkrobust.demo.evaluate import optimal_adjustment_set_dag
    from bkrobust.demo.meek import enumerate_dag_extensions
    from bkrobust.hybrid import breakdown_radius
    from bkrobust.mpdag_criterion import is_valid_mpdag
    from bkrobust.robustness import llm_survival as ls
    from bkrobust.robustness import real_survival as rs

    bundle, bundle_sha = ls.load_bundle()
    real_conditions = set(ls.panel_conditions(bundle, naming=ls.REAL_NAMING))
    _log(f"elicitation bundle sha256={bundle_sha[:16]} real_conditions={sorted(real_conditions)}")

    rows = list(csv.DictReader(UNITS_CSV.open()))
    units = [r for r in rows if r["status"] == "optimal_set_undefined"]
    _log(f"loaded {len(rows)} analysis rows; {len(units)} carry optimal_set_undefined")

    groups: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for u in units:
        groups[(u["condition"], u["network"])].append(u)

    def group_sort_key(k: tuple[str, str]) -> tuple[int, str, str]:
        cond, net = k
        return (0 if cond in real_conditions else 1, cond, net)

    ordered_groups = sorted(groups, key=group_sort_key)
    _log(f"{len(ordered_groups)} (condition, network) groups; real-naming groups first")

    results: list[dict] = []
    method_counts: Counter = Counter()
    class_counts: Counter = Counter()
    stopped_early = False

    for gi, (condition, network) in enumerate(ordered_groups):
        if time.time() - _START > TOTAL_BUDGET_S:
            _log(f"TOTAL_BUDGET_S={TOTAL_BUDGET_S:.0f}s exceeded before group {gi}/{len(ordered_groups)}"
                 f" ({condition}, {network}); stopping early")
            stopped_early = True
            break

        group_units = groups[(condition, network)]
        try:
            k, _elicit_cols = ls.elicited_k(bundle, condition, network)
            dag, cpdag = rs.load_network(network)
            g0, g0_reason = rs.build_g0(cpdag, k)
        except Exception as exc:  # noqa: BLE001
            _log(f"group ({condition}, {network}): FAILED to build g0: {exc!r}")
            for u in group_units:
                results.append({
                    "condition": condition, "network": network, "x": u["x"], "y": u["y"],
                    "naming": u.get("naming"), "classification": "undetermined",
                    "method": "g0_build_error", "elapsed_s": 0.0,
                })
            continue

        if g0 is None:
            _log(f"group ({condition}, {network}): g0 is None (reason={g0_reason}) "
                 f"for {len(group_units)} units carrying optimal_set_undefined -- unexpected")
            for u in group_units:
                results.append({
                    "condition": condition, "network": network, "x": u["x"], "y": u["y"],
                    "naming": u.get("naming"), "classification": "undetermined",
                    "method": "g0_none", "elapsed_s": 0.0,
                })
            continue

        unit_results, summary = classify_group(
            g0, group_units, is_valid_mpdag, enumerate_dag_extensions, optimal_adjustment_set_dag,
        )
        for r in unit_results:
            r["condition"], r["network"] = condition, network
            results.append(r)
            class_counts[r["classification"]] += 1
            method_counts[r["method"].split(":")[0]] += 1

        _log(f"group {gi + 1}/{len(ordered_groups)} ({condition}, {network}): "
             f"{len(group_units)} units, |K|={len(k)}, g0_undirected={len(g0.undirected_edges)}, "
             f"n_extensions={summary['n_extensions']} (enum {summary['enum_elapsed_s']}s), "
             f"timeouts={summary['n_timeout']}")

        if (gi + 1) % 10 == 0:
            _flush(results, bundle_sha, real_conditions, stopped_early=False, partial=True)

    _log(f"classification done: {dict(class_counts)}")
    _log(f"methods used: {dict(method_counts)}")

    # --- r_val for empty_valid units -----------------------------------------
    empty_units = [r for r in results if r["classification"] == "empty_valid"]
    _log(f"computing r_val for {len(empty_units)} empty_valid units")

    by_group: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for r in empty_units:
        by_group[(r["condition"], r["network"])].append(r)

    n_rval_done = 0
    for (condition, network), group_rows in by_group.items():
        if time.time() - _START > TOTAL_BUDGET_S:
            _log("TOTAL_BUDGET_S exceeded during r_val phase; stopping r_val computation early")
            break
        try:
            k, _ = ls.elicited_k(bundle, condition, network)
            dag, cpdag = rs.load_network(network)
            g0, _reason = rs.build_g0(cpdag, k)
        except Exception as exc:  # noqa: BLE001
            for r in group_rows:
                r["r_val"], r["r_status"] = None, f"g0_rebuild_error:{exc!r}"
            n_rval_done += len(group_rows)
            continue
        for r in group_rows:
            x, y = r["x"], r["y"]
            try:
                out = breakdown_radius(
                    cpdag, None, x, y, frozenset(), g0=g0, time_limit_s=R_VAL_TIME_LIMIT_S,
                )
                if not out.exact:
                    r["r_val"], r["r_status"] = None, "censored"
                elif out.radius == UNREACHED:
                    r["r_val"], r["r_status"] = float("inf"), out.method
                else:
                    r["r_val"], r["r_status"] = out.radius, out.method
            except Exception as exc:  # noqa: BLE001
                r["r_val"], r["r_status"] = None, f"error:{exc!r}"
            n_rval_done += 1
        _log(f"r_val: {n_rval_done}/{len(empty_units)} done (group {condition}/{network}, n={len(group_rows)})")

    _flush(results, bundle_sha, real_conditions, stopped_early=stopped_early, partial=False)
    _log("done")


def _flush(results, bundle_sha, real_conditions, *, stopped_early: bool, partial: bool) -> None:
    class_counts = Counter(r["classification"] for r in results)
    method_counts = Counter(r["method"].split(":")[0] for r in results)
    real_units = [r for r in results if r["condition"] in real_conditions]
    real_class_counts = Counter(r["classification"] for r in real_units)

    empty_valid_scored = [r for r in results if r["classification"] == "empty_valid" and "r_val" in r]
    r_vals = [r["r_val"] for r in empty_valid_scored]
    finite = [v for v in r_vals if v is not None and v != float("inf")]
    n_inf = sum(1 for v in r_vals if v == float("inf"))

    payload = {
        "meta": {
            "bundle_sha256": bundle_sha,
            "real_conditions": sorted(real_conditions),
            "per_unit_time_limit_s": PER_UNIT_TIME_LIMIT_S,
            "r_val_time_limit_s": R_VAL_TIME_LIMIT_S,
            "total_budget_s": TOTAL_BUDGET_S,
            "stopped_early": stopped_early,
            "partial": partial,
            "n_units_total": len(results),
            "wall_s": round(time.time() - _START, 1),
        },
        "classification_counts_overall": dict(class_counts),
        "classification_counts_real_conditions": dict(real_class_counts),
        "method_counts": dict(method_counts),
        "r_val_summary": {
            "n_empty_valid_scored": len(empty_valid_scored),
            "n_inf": n_inf,
            "n_finite": len(finite),
            "share_r_eq_1": (sum(1 for v in finite if v == 1) / len(finite)) if finite else None,
            "max_finite": max(finite) if finite else None,
        },
        "units": results,
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2, default=str))


if __name__ == "__main__":
    main()
