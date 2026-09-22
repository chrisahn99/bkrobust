"""Reviewer Point 7 audit: what the 33 zero-radius rows in the final table actually are.

The reviewer's objection: 33 rows in ``results/final_table/instances.jsonl`` have
``r_val = 0``, and the objection reads this as "the optimal adjustment set is
empty and already invalid at G0", which would be a precondition violation of
the method.

Two facts, already established before this script was written and re-verified
by it (see :func:`check_real_corpus_baseline`):

* The 33 rows live only in ``results/final_table/`` (the LLM-elicited-knowledge
  table). The 831-row real corpus at ``results/axisa3/instances.jsonl``
  (``admissible == True`` rows) has zero instances with ``radius == 0`` and zero
  with an empty optimal set; its minimum finite radius is 1.
* All 33 are already labelled ``method == "degenerate"`` by
  :mod:`bkrobust.hybrid` (see ``hybrid.py`` around line 224) and are already
  excluded from every median in ``experiments/final_table.py``
  (:func:`summarise_network` splits ``degenerate`` out of ``ok`` before any
  ``_cell``/median is computed). So the 33 are not silently averaged in; the
  question this script answers is finer-grained: *what kind of "empty and
  invalid" are they, and is any of it a bug?*

This script does not trust the committed label. For every non-informative
instance in the final table -- the 33 ``r_val == 0`` "degenerate" rows and the
27 ``status == "error"`` ("no adjustment set is identified at G0") rows -- it
rebuilds the exact (network, K, G0) triple ``experiments/final_table.py`` built,
and recomputes, from the graph alone, using the same functions the pipeline
itself calls (:func:`bkrobust.demo.evaluate.optimal_adjustment_set_mpdag`,
:func:`bkrobust.mpdag_criterion.is_valid_mpdag`, :func:`bkrobust.mpdag_criterion.
is_amenable`, :func:`bkrobust.mpdag_criterion.why_invalid`):

* Is an optimal adjustment set identified at G0 at all (agrees across every DAG
  extension of G0)?
* If so, is *that* set GAC-valid at G0 -- checked independently, not inferred
  from ``r_val``?
* Is G0 amenable relative to (x, y)?

From those three graph-level facts alone (never from the committed ``r_val`` or
``status``) each instance is placed into exactly one of three classes:

    VALID_EMPTY    -- Z is identified, |Z| == 0, and Z passes GAC at G0.
                      A legitimate, well-posed query; a radius is meaningful.
    NO_SET_EXISTS  -- no optimal adjustment set is identified at G0 at all
                      (``optimal_adjustment_set_mpdag`` returns None: either G0
                      has no consistent DAG extension, or the extensions
                      disagree on the optimal set). Unidentified, not
                      "radius 0".
    INVALID_RETURNED -- an adjustment set (possibly empty) was identified but
                      fails GAC at G0. This is the only class that would be a
                      genuine implementation defect, and is reported as one if
                      it is ever non-empty.

Also computed, as a cross-check that is not part of the three-way split: for
every *informative* row (``status == "ok"`` and ``r_val != 0``), whether its
committed Z is empty. If some are, that independently demonstrates that
"|Z| == 0" and "the query is degenerate" are different properties -- a network
can have a legitimately empty, valid optimal adjustment set with a nonzero
radius.

Usage::

    PYTHONPATH=src /usr/bin/python3 experiments/zero_radius_audit.py

Writes ``results/zero_radius/instances.jsonl``, ``results/zero_radius/
instances.csv`` and ``results/zero_radius/summary.json``. Does not touch
``experiments/final_table.py``, ``src/bkrobust/hybrid.py``, or anything under
``results/final_table/`` -- it only reads them.
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from bkrobust.benchmarks.describe import parse_file, to_mpdag  # noqa: E402
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag  # noqa: E402
from bkrobust.demo.example import dag_to_cpdag  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions  # noqa: E402
from bkrobust.mpdag_criterion import is_amenable, is_valid_mpdag, why_invalid  # noqa: E402
from bkrobust.mpdag_criterion.paths import possibly_causal_paths  # noqa: E402

NETWORKS_DIR = ROOT / "results" / "axisa3" / "networks" / "example_models"
KNOWLEDGE_PATH = ROOT / "results" / "elicit" / "knowledge.json"
FINAL_TABLE_DIR = ROOT / "results" / "final_table"
FINAL_TABLE_EPS_GT1_DIR = ROOT / "results" / "final_table_eps_gt1"
AXISA3_INSTANCES = ROOT / "results" / "axisa3" / "instances.jsonl"
OUT_DIR = ROOT / "results" / "zero_radius"

COMMITTED_COUNTS = {"degenerate": 33, "no_identified_z": 27, "timeout": 34, "informative": 66, "total": 160}


def log(msg: str) -> None:
    print(f"[zero_radius_audit] {msg}", file=sys.stderr, flush=True)


# --------------------------------------------------------------------------------
# Loading (mirrors experiments/final_table.py exactly, read-only)
# --------------------------------------------------------------------------------


def load_instances(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text().splitlines() if line]


def load_knowledge(condition: str) -> dict[str, list[tuple[str, str]]]:
    data = json.loads(KNOWLEDGE_PATH.read_text())
    nets = data[condition]["networks"]
    return {net: [tuple(e) for e in info["k"]] for net, info in nets.items()}


def build_network(name: str) -> dict[str, MPDAG] | None:
    """Parse one real network and build its (dag, cpdag), exactly as final_table.py does."""
    for path in sorted(NETWORKS_DIR.iterdir()):
        if not path.is_file() or path.name.startswith("."):
            continue
        try:
            parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
        except Exception:
            continue
        if parsed.bidirected or parsed.name != name:
            continue
        dag = to_mpdag(parsed)
        cpdag = dag_to_cpdag(dag)
        return {"dag": dag, "cpdag": cpdag}
    return None


# --------------------------------------------------------------------------------
# The three-way classification
# --------------------------------------------------------------------------------


def classify_instance(
    cpdag: MPDAG, k: list[tuple[str, str]], x: str, y: str
) -> dict[str, Any]:
    """Rebuild G0 and recompute the three-way class from the graph alone.

    Returns a dict of every intermediate fact, not just the final class label,
    so the record is auditable on its own.
    """
    g0 = apply_orientations(cpdag, k)
    if g0 is None:
        # Knowledge inconsistent with the CPDAG -- did not happen in the
        # committed sweep (it would have surfaced as a different error message
        # than the one all 27 "error" rows carry), but handled defensively.
        return {
            "g0_built": False,
            "amenable": None,
            "z_identified": None,
            "z": None,
            "z_size": None,
            "z_valid_at_g0": None,
            "why_invalid": None,
            "n_extensions": None,
            "recomputed_class": "NO_SET_EXISTS",
            "recompute_note": "apply_orientations(cpdag, k) returned None: knowledge inconsistent with cpdag",
        }

    amenable = is_amenable(g0, x, y)
    extensions = enumerate_dag_extensions(g0)
    z_opt = optimal_adjustment_set_mpdag(g0, x, y)

    if z_opt is None:
        reason = (
            "no_consistent_extension"
            if not extensions
            else "extensions_disagree_on_optimal_set"
        )
        return {
            "g0_built": True,
            "amenable": amenable,
            "z_identified": False,
            "z": None,
            "z_size": None,
            "z_valid_at_g0": None,
            "why_invalid": None,
            "n_extensions": len(extensions),
            "recomputed_class": "NO_SET_EXISTS",
            "recompute_note": reason,
        }

    z = frozenset(z_opt)
    valid = is_valid_mpdag(g0, x, y, z)
    reason = "" if valid else why_invalid(g0, x, y, z)

    if valid and len(z) == 0:
        cls = "VALID_EMPTY"
    elif valid:
        # Identified, non-empty, and GAC-valid at G0 -- not a "non-informative"
        # verdict at all. If this is reached for a row the committed pipeline
        # called degenerate/error, that is itself a finding worth flagging.
        cls = "VALID_NONEMPTY_UNEXPECTED"
    else:
        cls = "INVALID_RETURNED"

    # Diagnostic only, not part of the classification: when the reason is
    # "open_noncausal_path", is that because X has literally no possibly
    # causal path to Y at all in G0 (so *every* X-Y path is non-causal by
    # definition, and the Henckel O-set formula's cn(x, y) = empty degenerates
    # to O = empty without ever checking whether the empty set blocks the
    # confounding paths that remain)? Cheap to check only in this branch.
    has_causal_path = None
    if reason == "open_noncausal_path":
        has_causal_path = len(possibly_causal_paths(g0, x, y)) > 0

    return {
        "g0_built": True,
        "amenable": amenable,
        "z_identified": True,
        "z": sorted(z),
        "z_size": len(z),
        "z_valid_at_g0": valid,
        "why_invalid": reason,
        "n_extensions": len(extensions),
        "has_possibly_causal_path_x_to_y": has_causal_path,
        "recomputed_class": cls,
        "recompute_note": "",
    }


# --------------------------------------------------------------------------------
# The real-corpus baseline (re-verified, not assumed)
# --------------------------------------------------------------------------------


def check_real_corpus_baseline() -> dict[str, Any]:
    """Re-verify the orchestrator's claim about results/axisa3/instances.jsonl.

    Claim: among the 831 admissible rows, zero have radius == 0, zero have
    z_size == 0, and the minimum finite radius is 1.
    """
    log(f"re-verifying real-corpus baseline from {AXISA3_INSTANCES} ...")
    n_admissible = 0
    zero_radius = 0
    empty_z = 0
    min_radius = None
    with AXISA3_INSTANCES.open() as f:
        for line in f:
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("_meta"):
                continue
            if not row.get("admissible"):
                continue
            n_admissible += 1
            radius = row["radius"]
            if radius != -1:
                if min_radius is None or radius < min_radius:
                    min_radius = radius
                if radius == 0:
                    zero_radius += 1
            if row.get("z_size") == 0:
                empty_z += 1
    result = {
        "n_admissible": n_admissible,
        "n_zero_radius": zero_radius,
        "n_empty_z": empty_z,
        "min_finite_radius": min_radius,
    }
    log(f"real-corpus baseline: {result}")
    ok = n_admissible == 831 and zero_radius == 0 and empty_z == 0 and min_radius == 1
    log(f"matches orchestrator's claim (831 admissible, min radius 1, 0 zero-radius, 0 empty-Z)? {ok}")
    result["matches_orchestrator_claim"] = ok
    return result


# --------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------


def main() -> None:
    t_start = time.perf_counter()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    real_corpus_check = check_real_corpus_baseline()

    log(f"loading {FINAL_TABLE_DIR / 'instances.jsonl'} ...")
    ft_rows = load_instances(FINAL_TABLE_DIR / "instances.jsonl")
    log(f"  {len(ft_rows)} rows")

    eps_gt1_path = FINAL_TABLE_EPS_GT1_DIR / "instances.jsonl"
    eps_gt1_rows: list[dict[str, Any]] = []
    eps_gt1_consistency: dict[str, Any] = {"checked": False}
    if eps_gt1_path.exists():
        log(f"loading {eps_gt1_path} for a schema/consistency cross-check ...")
        eps_gt1_rows = load_instances(eps_gt1_path)
        same_schema = bool(eps_gt1_rows) and set(eps_gt1_rows[0]) >= {
            "network", "x", "y", "condition", "r_val", "status",
        }
        ft_by_key = {(r["network"], r["x"], r["y"], r["condition"]): r for r in ft_rows}
        eg_by_key = {(r["network"], r["x"], r["y"], r["condition"]): r for r in eps_gt1_rows}
        same_keys = set(ft_by_key) == set(eg_by_key)
        diffs = []
        if same_keys:
            for key in ft_by_key:
                a, b = ft_by_key[key], eg_by_key[key]
                if (a["status"], a["r_val"]) != (b["status"], b["r_val"]):
                    diffs.append(
                        {
                            "network": key[0], "x": key[1], "y": key[2],
                            "final_table": {"status": a["status"], "r_val": a["r_val"]},
                            "final_table_eps_gt1": {"status": b["status"], "r_val": b["r_val"]},
                        }
                    )
        eps_gt1_consistency = {
            "checked": True,
            "same_schema": same_schema,
            "same_instance_keys": same_keys,
            "n_rows": len(eps_gt1_rows),
            "n_status_or_rval_diffs": len(diffs),
            "diffs": diffs,
        }
        log(
            f"  final_table_eps_gt1: same instance set as final_table? {same_keys}; "
            f"{len(diffs)} status/r_val disagreement(s) (see summary.json for detail)"
        )
    else:
        log(f"  {eps_gt1_path} not found; skipping cross-check")

    # The audit scope: every non-informative instance in the canonical
    # final_table/instances.jsonl -- the 33 r_val==0 "degenerate" rows and the
    # 27 status=="error" ("no adjustment set identified") rows. Informative
    # rows (status=="ok" and r_val != 0) are out of scope for the three-way
    # class, but scanned separately below for the VALID_EMPTY cross-check.
    degenerate_rows = [r for r in ft_rows if r["status"] == "ok" and r["r_val"] == 0]
    error_rows = [r for r in ft_rows if r["status"] == "error"]
    timeout_rows = [r for r in ft_rows if r["status"] == "timeout"]
    informative_rows = [r for r in ft_rows if r["status"] == "ok" and r["r_val"] != 0]

    # UNREACHED (r_val == -1) is a sentinel, never a radius (core.conventions.
    # UNREACHED; see PAPER_INTEGRATION_FINAL_TABLE.md Sec 5: "UNREACHED is
    # never averaged"). It lives *inside* the "informative" bucket above
    # (status == "ok" and r_val != 0 is true for both finite radii and -1), and
    # is entirely orthogonal to the three-way VALID_EMPTY / NO_SET_EXISTS /
    # INVALID_RETURNED split this script computes: that split is scoped to the
    # 33 degenerate (r_val == 0) and 27 error ("no adjustment set identified")
    # rows only, none of which can be -1 by construction. Reported here purely
    # so "66 informative" is never misread as 66 homogeneous finite radii.
    informative_unreached_rows = [r for r in informative_rows if r["r_val"] == -1]
    informative_finite_rows = [r for r in informative_rows if r["r_val"] != -1]

    log(
        f"committed buckets in final_table/instances.jsonl: "
        f"degenerate={len(degenerate_rows)} error={len(error_rows)} "
        f"timeout={len(timeout_rows)} informative={len(informative_rows)} "
        f"(of which UNREACHED sentinel={len(informative_unreached_rows)}, "
        f"finite radius={len(informative_finite_rows)}) "
        f"total={len(ft_rows)}"
    )
    non_informative = degenerate_rows + error_rows
    log(f"audit scope: {len(non_informative)} non-informative instances")

    unexpected_error_messages = sorted(
        {r["error"] for r in error_rows if r["error"] != "no adjustment set is identified at G0; nothing to certify"}
    )
    if unexpected_error_messages:
        log(f"WARNING: unexpected error messages among status=='error' rows: {unexpected_error_messages}")

    # Recompute, network by network (cache the (dag, cpdag) build per network).
    condition = "D_LLM"
    log(f"loading knowledge for condition={condition} ...")
    knowledge = load_knowledge(condition)

    networks_needed = sorted({r["network"] for r in non_informative} | {r["network"] for r in informative_rows})
    net_cache: dict[str, dict[str, MPDAG] | None] = {}
    for net in networks_needed:
        cond_mismatch = any(
            r["network"] == net and r["condition"] != condition for r in (non_informative + informative_rows)
        )
        if cond_mismatch:
            log(f"WARNING: {net} has a row with condition != {condition!r}; audit assumes condition={condition!r}")
        t0 = time.perf_counter()
        net_cache[net] = build_network(net)
        log(f"  built {net}: {'ok' if net_cache[net] else 'FAILED'} ({time.perf_counter() - t0:.2f}s)")

    records: list[dict[str, Any]] = []
    for i, row in enumerate(non_informative, 1):
        net, x, y = row["network"], row["x"], row["y"]
        built = net_cache.get(net)
        committed_label = "degenerate" if row["status"] == "ok" else "no_identified_z"
        base = {
            "network": net,
            "x": x,
            "y": y,
            "condition": row["condition"],
            "committed_status": row["status"],
            "committed_r_val": row["r_val"],
            "committed_label": committed_label,
        }
        if built is None:
            base.update(
                {
                    "recomputed_class": "ERROR_COULD_NOT_BUILD_NETWORK",
                    "recompute_note": f"failed to rebuild network {net!r}",
                }
            )
            records.append(base)
            log(f"[{i}/{len(non_informative)}] {net} {x}->{y}: FAILED to rebuild network")
            continue
        t0 = time.perf_counter()
        try:
            result = classify_instance(built["cpdag"], knowledge[net], x, y)
        except Exception as exc:  # keep going; record the failure as data
            result = {
                "g0_built": None,
                "amenable": None,
                "z_identified": None,
                "z": None,
                "z_size": None,
                "z_valid_at_g0": None,
                "why_invalid": None,
                "n_extensions": None,
                "recomputed_class": "ERROR_RECOMPUTE_RAISED",
                "recompute_note": f"{type(exc).__name__}: {exc}",
            }
        dt = time.perf_counter() - t0
        base.update(result)
        base["seconds"] = round(dt, 4)
        records.append(base)
        log(
            f"[{i}/{len(non_informative)}] {net} {x}->{y}: committed={committed_label} "
            f"recomputed={result['recomputed_class']} amenable={result.get('amenable')} "
            f"z_size={result.get('z_size')} valid={result.get('z_valid_at_g0')} ({dt:.3f}s)"
        )

    # Cross-check: among informative rows, how many have a committed Z that is
    # empty? (Not part of the three-way split -- these are already r_val != 0,
    # i.e. already informative. This just demonstrates that |Z|==0 does not by
    # itself imply degeneracy.) We recompute Z the same way for these too.
    log("cross-check: scanning informative rows for committed adjustment sets with |Z| == 0 ...")
    informative_empty_z: list[dict[str, Any]] = []
    for row in informative_rows:
        net, x, y = row["network"], row["x"], row["y"]
        built = net_cache.get(net)
        if built is None:
            continue
        g0 = apply_orientations(built["cpdag"], knowledge[net])
        if g0 is None:
            continue
        z_opt = optimal_adjustment_set_mpdag(g0, x, y)
        if z_opt is not None and len(z_opt) == 0:
            informative_empty_z.append(
                {"network": net, "x": x, "y": y, "r_val": row["r_val"]}
            )
    log(
        f"  {len(informative_empty_z)}/{len(informative_rows)} informative rows have a "
        f"legitimately empty, valid, nonzero-radius adjustment set"
    )

    # Tally the three-way split.
    counts = {"VALID_EMPTY": 0, "NO_SET_EXISTS": 0, "INVALID_RETURNED": 0}
    other = {}
    for r in records:
        c = r["recomputed_class"]
        if c in counts:
            counts[c] += 1
        else:
            other[c] = other.get(c, 0) + 1

    # Cross-tab against the committed label, to see whether recomputation
    # agrees with what the pipeline already believed.
    crosstab: dict[str, dict[str, int]] = {}
    for r in records:
        lbl = r["committed_label"]
        cls = r["recomputed_class"]
        crosstab.setdefault(lbl, {}).setdefault(cls, 0)
        crosstab[lbl][cls] += 1

    invalid_returned_detail = [
        {
            "network": r["network"], "x": r["x"], "y": r["y"],
            "z": r["z"], "z_size": r["z_size"], "why_invalid": r["why_invalid"],
            "has_possibly_causal_path_x_to_y": r.get("has_possibly_causal_path_x_to_y"),
            "committed_label": r["committed_label"], "committed_r_val": r["committed_r_val"],
        }
        for r in records
        if r["recomputed_class"] == "INVALID_RETURNED"
    ]
    invalid_returned_why = {}
    for r in records:
        if r["recomputed_class"] == "INVALID_RETURNED":
            invalid_returned_why[r["why_invalid"]] = invalid_returned_why.get(r["why_invalid"], 0) + 1
    invalid_returned_no_causal_path = sum(
        1 for r in records
        if r["recomputed_class"] == "INVALID_RETURNED"
        and r.get("has_possibly_causal_path_x_to_y") is False
    )
    unexpected_detail = [
        {
            "network": r["network"], "x": r["x"], "y": r["y"],
            "z": r["z"], "z_size": r["z_size"],
            "committed_label": r["committed_label"], "committed_r_val": r["committed_r_val"],
        }
        for r in records
        if r["recomputed_class"] == "VALID_NONEMPTY_UNEXPECTED"
    ]

    reconciliation = {
        "committed_counts": COMMITTED_COUNTS,
        "recomputed_degenerate_rows_in_final_table": len(degenerate_rows),
        "recomputed_error_rows_in_final_table": len(error_rows),
        "recomputed_timeout_rows_in_final_table": len(timeout_rows),
        "recomputed_informative_rows_in_final_table": len(informative_rows),
        "recomputed_informative_unreached_sentinel_rows": len(informative_unreached_rows),
        "recomputed_informative_finite_radius_rows": len(informative_finite_rows),
        "recomputed_total_rows": len(ft_rows),
        "matches_committed_counts": (
            len(degenerate_rows) == COMMITTED_COUNTS["degenerate"]
            and len(error_rows) == COMMITTED_COUNTS["no_identified_z"]
            and len(timeout_rows) == COMMITTED_COUNTS["timeout"]
            and len(informative_rows) == COMMITTED_COUNTS["informative"]
            and len(ft_rows) == COMMITTED_COUNTS["total"]
        ),
    }

    unreached_sentinel = {
        "note": (
            "r_val == -1 is the UNREACHED sentinel (bkrobust.core.conventions.UNREACHED), "
            "never a radius and never averaged. It occurs only inside the 'informative' "
            "bucket (status=='ok' and r_val != 0), alongside finite-radius rows, and is "
            "orthogonal to this script's VALID_EMPTY/NO_SET_EXISTS/INVALID_RETURNED "
            "three-way split, whose scope is exactly the 33 degenerate (r_val==0) and 27 "
            "error ('no adjustment set identified') rows -- neither can be -1 by "
            "construction, so no UNREACHED row is or could be folded into any of the "
            "three classes."
        ),
        "n_informative_total": len(informative_rows),
        "n_informative_unreached_sentinel": len(informative_unreached_rows),
        "n_informative_finite_radius": len(informative_finite_rows),
        "unreached_rows": [
            {"network": r["network"], "x": r["x"], "y": r["y"]} for r in informative_unreached_rows
        ],
    }

    summary = {
        "generated_by": "experiments/zero_radius_audit.py",
        "condition": condition,
        "interpreter_note": (
            "This script imports only bkrobust.demo.* and bkrobust.mpdag_criterion.* "
            "(graph/GAC recomputation at G0) -- never bkrobust.epsilon.certify, which uses "
            "zip(..., strict=True) (Python 3.10+ only) and cannot import under the "
            "PYTHONPATH=src /usr/bin/python3 (3.9.6) interpreter this script is run with. "
            "results/axisa3/instances.jsonl (the 831-row real-corpus baseline re-verified "
            "below) was also produced under 3.9.6-compatible code. results/final_table/ and "
            "results/final_table_eps_gt1/ (read-only inputs here, not recomputed) were "
            "produced by experiments/final_table.py under .venv/bin/python (3.10+), since "
            "that script does call certify(). This script never calls certify() and never "
            "needs the venv."
        ),
        "real_corpus_baseline_check": real_corpus_check,
        "unreached_sentinel": unreached_sentinel,
        "final_table_eps_gt1_consistency_check": eps_gt1_consistency,
        "reconciliation_against_committed_counts": reconciliation,
        "three_way_class_counts": counts,
        "other_recomputed_classes": other,
        "crosstab_committed_label_vs_recomputed_class": crosstab,
        "n_invalid_returned": counts["INVALID_RETURNED"],
        "any_genuine_precondition_violation": counts["INVALID_RETURNED"] > 0,
        "invalid_returned_why_breakdown": invalid_returned_why,
        "invalid_returned_no_causal_path_x_to_y_count": invalid_returned_no_causal_path,
        "invalid_returned_detail": invalid_returned_detail,
        "unexpected_valid_nonempty_detail": unexpected_detail,
        "informative_rows_with_empty_committed_z": {
            "n": len(informative_empty_z),
            "of_informative": len(informative_rows),
            "rows": informative_empty_z,
        },
        "n_audited": len(records),
        "seconds_total": round(time.perf_counter() - t_start, 2),
    }

    instances_jsonl = OUT_DIR / "instances.jsonl"
    with instances_jsonl.open("w") as f:
        for r in records:
            f.write(json.dumps(r) + "\n")
    log(f"wrote {instances_jsonl} ({len(records)} rows)")

    csv_path = OUT_DIR / "instances.csv"
    fieldnames = [
        "network", "x", "y", "condition", "committed_status", "committed_r_val",
        "committed_label", "g0_built", "amenable", "z_identified", "z", "z_size",
        "z_valid_at_g0", "why_invalid", "n_extensions", "has_possibly_causal_path_x_to_y",
        "recomputed_class", "recompute_note", "seconds",
    ]
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for r in records:
            row = dict(r)
            if row.get("z") is not None:
                row["z"] = ";".join(row["z"])
            writer.writerow(row)
    log(f"wrote {csv_path}")

    summary_path = OUT_DIR / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    log(f"wrote {summary_path}")

    log("=" * 70)
    log(f"THREE-WAY SPLIT of the {len(records)} non-informative instances:")
    log(f"  VALID_EMPTY      = {counts['VALID_EMPTY']}")
    log(f"  NO_SET_EXISTS    = {counts['NO_SET_EXISTS']}")
    log(f"  INVALID_RETURNED = {counts['INVALID_RETURNED']}")
    if other:
        log(f"  other/unexpected = {other}")
    log(f"genuine precondition violation found? {counts['INVALID_RETURNED'] > 0}")
    log("=" * 70)


if __name__ == "__main__":
    main()
