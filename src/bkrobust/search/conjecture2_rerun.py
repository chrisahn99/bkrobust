r"""Re-run of the Conjecture 2 exhaustive sweep on the CORRECTED space.

Session 2 / Axis B. Task 0 of this session found that the inherited
``bkrobust.demo.space.enumerate_space`` drops legitimate knowledge states: it
filters through ``is_valid_mpdag``, which demands chordality of the
*undirected subgraph* -- a CPDAG property, not an MPDAG property. Background
knowledge can place a directed edge inside what would otherwise be one
undirected component; the chord that would make the component chordal is then
present but directed, so the undirected subgraph looks chordless and a real
state is dropped. See :mod:`bkrobust.search.space_fixed` for the full account
and the fixpoint predicate that replaces the chordality filter.

The previous session's ~2M-comparison Conjecture 2 sweep
(``results/search/conjectures/n{3,4,5}.json``) ran on the defective space, so
its "no counterexample" claim needs re-deriving on the corrected one. This
module does that:

* :func:`run_study_corrected` mirrors
  :func:`bkrobust.search.conjecture_study.run_study` structurally (same
  guards, same ``c1_every`` sampling, same adjustment-set cap) but builds
  spaces with :func:`bkrobust.search.space_fixed.build_space_correct`, and
  additionally reports how many space elements and (G0, X, Y, Z) combinations
  exist only because of the correction.
* :func:`identify_density_gap` names the CPDAGs the corrected sweep still had
  to skip for size, exactly as the old sweep's skip counts are reported
  rather than assumed.
* :func:`run_dense_closure` runs those CPDAGs to completion (or records
  precisely why one could not be), closing the gap the old sweep left at
  n=5's densest instances -- exactly where a counterexample would be most
  likely to live.

Determinism and crash-safety (see the module's tests):

* No global RNG is used anywhere in this module.
* Every ordering that reaches output is either an already-sorted tuple
  (``space.elements``, ``all_cpdags(n)``) or explicitly ``sorted(...)`` here;
  nothing depends on ``set``/``dict`` iteration order, which is what makes the
  output byte-identical across ``PYTHONHASHSEED``.
* Results are written **one line per CPDAG**, flushed and fsynced
  immediately, so a kill loses at most the CPDAG in flight. A second call
  with ``resume=True`` (the default) picks up exactly where the first left
  off, replaying already-written deltas to reconstruct the running
  ``c1_every`` sample counter before continuing.
"""

from __future__ import annotations

import itertools
import json
import os
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from bkrobust.core.oracle import is_valid
from bkrobust.core.spacelib import Space, distances_from
from bkrobust.demo.evaluate import (
    all_valid_adjustment_sets_mpdag,
    optimal_adjustment_set_mpdag,
)
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.space import enumerate_space  # the OLD (defective) space
from bkrobust.search.conjecture_study import all_cpdags, all_dags
from bkrobust.search.conjectures import (
    Conjecture2Result,
    check_conjecture1,
    check_conjecture2,
    check_no_empty_extensions,
)
from bkrobust.search.space_fixed import build_space_correct, enumerate_space_correct

#: The additive fields tracked per CPDAG and summed into the running totals.
#: A fixed list literal, never derived from set/dict iteration, so the sum
#: order is identical regardless of PYTHONHASHSEED.
DELTA_FIELDS: tuple[str, ...] = (
    "n_spaces",
    "n_combos",
    "n_radius_comparisons",
    "c1_violations",
    "c2_violations",
    "empty_extension_spaces",
    "skipped_no_valid_z",
    "skipped_z_fails_at_g0",
    "skipped_cpdag_too_many_undirected",
    "skipped_cpdag_space_too_big",
    "skipped_cpdag_no_undirected",
    "c1_checks_run",
    "c1_checks_skipped",
    "n_states_added_vs_old",
    "n_combos_new",
)


@dataclass
class CorrectedStudyTotals:
    """Running totals for the corrected-space sweep.

    Mirrors :class:`bkrobust.search.conjecture_study.StudyTotals` field for
    field (so old and new are directly comparable), plus two fields that only
    make sense once there is an "old" space to compare against:
    ``n_states_added_vs_old`` (space elements that exist only because the
    chordality filter no longer drops them) and ``n_combos_new``
    (``(G0, X, Y, Z)`` combinations touching at least one such element, either
    as ``G0`` or as a witness).
    """

    n_dags: int = 0
    n_cpdags: int = 0
    n_spaces: int = 0
    n_combos: int = 0
    n_radius_comparisons: int = 0
    c1_violations: int = 0
    c2_violations: int = 0
    empty_extension_spaces: int = 0
    skipped_no_valid_z: int = 0
    skipped_z_fails_at_g0: int = 0
    skipped_cpdag_too_many_undirected: int = 0
    skipped_cpdag_space_too_big: int = 0
    skipped_cpdag_no_undirected: int = 0
    c1_checks_run: int = 0
    c1_checks_skipped: int = 0
    n_states_added_vs_old: int = 0
    n_combos_new: int = 0


def _new_deltas() -> dict[str, int]:
    return {f: 0 for f in DELTA_FIELDS}


def _add_deltas(totals: CorrectedStudyTotals, deltas: dict[str, int]) -> None:
    for f in DELTA_FIELDS:
        setattr(totals, f, getattr(totals, f) + deltas.get(f, 0))


# --- the per-space sweep, shared by the main run and the dense closure -----


def _sweep_one_space(
    cpdag: MPDAG,
    space: Space,
    old_strings: frozenset[str],
    *,
    n: int,
    combo_offset: int,
    c1_every: int,
    all_z: bool,
    max_z_per_g0: int,
) -> tuple[dict[str, int], list[dict[str, Any]]]:
    """Run the full (G0, X, Y, Z) sweep over one already-built space.

    This is the inner loop of :func:`bkrobust.search.conjecture_study.run_study`,
    restructured to (a) build on the CORRECTED space and (b) additionally flag
    combinations that exist only because of the correction.

    Args:
        cpdag: The CPDAG the space refines.
        space: The corrected space, already built.
        old_strings: ``edge_string()`` of every element of the OLD (defective)
            space for this CPDAG -- used to flag which combinations are new.
        n: Number of nodes (for record-keeping only).
        combo_offset: The GLOBAL running combo count immediately before this
            CPDAG, so the ``c1_every`` sample lands on the same combinations
            regardless of how the sweep is chunked across processes/resumes.
        c1_every: Sample rate for the O(|space|^2) Conjecture-1 guard check.
        all_z: Test every valid adjustment set rather than only ``O*``.
        max_z_per_g0: Cap on adjustment sets tested per ``G0``.

    Returns:
        ``(deltas, counterexamples)``.
    """
    deltas = _new_deltas()
    deltas["n_spaces"] = 1
    counterexamples: list[dict[str, Any]] = []
    if not check_no_empty_extensions(space):
        deltas["empty_extension_spaces"] = 1

    nodes = list(cpdag.nodes)
    combo_count = combo_offset

    for g0 in space.elements:
        dists = distances_from(space, g0)
        g0_is_new = g0.edge_string() not in old_strings
        for x, y in itertools.permutations(nodes, 2):
            if all_z:
                zs = all_valid_adjustment_sets_mpdag(g0, x, y)[:max_z_per_g0]
            else:
                o = optimal_adjustment_set_mpdag(g0, x, y)
                zs = [] if o is None else [frozenset(o)]
            if not zs:
                deltas["skipped_no_valid_z"] += 1
                continue
            for zset in zs:
                z = frozenset(zset)
                if not is_valid(z, g0, x, y):
                    deltas["skipped_z_fails_at_g0"] += 1
                    continue
                deltas["n_combos"] += 1
                combo_count += 1

                def fails(
                    g: MPDAG,
                    _z: frozenset[str] = z,
                    _x: str = x,
                    _y: str = y,
                ) -> bool:
                    return not is_valid(_z, g, _x, _y)

                c2: Conjecture2Result = check_conjecture2(space, g0, fails, dists)
                deltas["n_radius_comparisons"] += 1

                wit_full_new = c2.witness_full is not None and c2.witness_full not in old_strings
                wit_up_new = c2.witness_up is not None and c2.witness_up not in old_strings
                if g0_is_new or wit_full_new or wit_up_new:
                    deltas["n_combos_new"] += 1

                if not c2.holds:
                    deltas["c2_violations"] += 1
                    counterexamples.append(
                        {
                            "conjecture": 2,
                            "n": n,
                            "cpdag": cpdag.edge_string(),
                            "g0": g0.edge_string(),
                            "g0_is_new_state": g0_is_new,
                            "x": x,
                            "y": y,
                            "z": sorted(z),
                            "r_full": c2.r_full,
                            "r_up": c2.r_up,
                            "witness_full": c2.witness_full,
                            "witness_up": c2.witness_up,
                            "space_size": c2.n_space,
                        }
                    )

                # Conjecture 1 is a theorem; sampled the same way, and at the
                # same rate, as the old sweep -- see conjecture_study.run_study.
                deltas["c1_checks_skipped"] += 1
                if combo_count % c1_every != 0:
                    continue
                deltas["c1_checks_skipped"] -= 1
                deltas["c1_checks_run"] += 1
                c1 = check_conjecture1(space, fails)
                if not c1.holds:
                    deltas["c1_violations"] += 1
                    counterexamples.append(
                        {
                            "conjecture": 1,
                            "n": n,
                            "cpdag": cpdag.edge_string(),
                            "g0": g0.edge_string(),
                            "x": x,
                            "y": y,
                            "z": sorted(z),
                            "violations": list(c1.violations[:3]),
                        }
                    )
    return deltas, counterexamples


def _record_skeleton(cpdag: MPDAG, n: int) -> dict[str, Any]:
    return {
        "cpdag": cpdag.edge_string(),
        "n": n,
        "k_undirected": len(cpdag.undirected_edges),
        "skip_reason": None,
        "space_size_old": None,
        "space_size_new": None,
        "n_states_added": 0,
        "counterexamples": [],
        "deltas": _new_deltas(),
    }


# --- deliverable 1: the main corrected-space sweep --------------------------


def _process_one_cpdag_main(
    cpdag: MPDAG,
    *,
    n: int,
    combo_offset: int,
    max_space: int,
    max_undirected: int,
    c1_every: int,
    all_z: bool,
    max_z_per_g0: int,
) -> dict[str, Any]:
    """One CPDAG's contribution to the main corrected sweep, with the old guards.

    Mirrors ``run_study``'s per-CPDAG guard order exactly: the undirected-edge
    count guard fires before any space is enumerated (both old and new spaces
    share the same undirected-edge count, since that is a property of the
    CPDAG itself), and the size guard fires on the CORRECTED enumeration
    (cheap; ``3**k``) before the cubic covering-relation build -- the same
    two-stage discipline the old sweep used, now applied to the space actually
    being built.
    """
    rec = _record_skeleton(cpdag, n)

    if not cpdag.undirected_edges:
        rec["deltas"]["skipped_cpdag_no_undirected"] = 1
        rec["skip_reason"] = "no_undirected_edges"
        return rec
    if len(cpdag.undirected_edges) > max_undirected:
        rec["deltas"]["skipped_cpdag_too_many_undirected"] = 1
        rec["skip_reason"] = (
            f"k={len(cpdag.undirected_edges)} undirected edges > max_undirected={max_undirected}"
        )
        return rec

    new_elements = enumerate_space_correct(cpdag)
    rec["space_size_new"] = len(new_elements)
    if len(new_elements) > max_space:
        rec["deltas"]["skipped_cpdag_space_too_big"] = 1
        rec["skip_reason"] = (
            f"corrected space has {len(new_elements)} elements > max_space={max_space}"
        )
        return rec

    old_elements = enumerate_space(cpdag)
    old_strings = frozenset(g.edge_string() for g in old_elements)
    new_strings = frozenset(g.edge_string() for g in new_elements)
    rec["space_size_old"] = len(old_elements)
    n_added = len(new_strings - old_strings)
    rec["n_states_added"] = n_added

    space = build_space_correct(cpdag)
    deltas, counterexamples = _sweep_one_space(
        cpdag,
        space,
        old_strings,
        n=n,
        combo_offset=combo_offset,
        c1_every=c1_every,
        all_z=all_z,
        max_z_per_g0=max_z_per_g0,
    )
    deltas["n_states_added_vs_old"] = n_added
    rec["deltas"] = deltas
    rec["counterexamples"] = counterexamples
    return rec


def _read_checkpoint(jsonl_path: Path) -> dict[str, dict[str, Any]]:
    done: dict[str, dict[str, Any]] = {}
    if not jsonl_path.exists():
        return done
    with jsonl_path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            done[rec["cpdag"]] = rec
    return done


def _append_record(fh: Any, rec: dict[str, Any]) -> None:
    fh.write(json.dumps(rec, sort_keys=True))
    fh.write("\n")
    fh.flush()
    os.fsync(fh.fileno())


def run_study_corrected(
    n: int,
    out_dir: str | Path,
    *,
    max_space: int = 200,
    max_undirected: int = 6,
    c1_every: int = 50,
    all_z: bool = True,
    max_z_per_g0: int = 8,
    verbose: bool = True,
    resume: bool = True,
) -> tuple[CorrectedStudyTotals, list[dict[str, Any]]]:
    """Run the Conjecture 2 sweep for ``n`` nodes on the CORRECTED space.

    Structurally identical to
    :func:`bkrobust.search.conjecture_study.run_study` -- same guards, same
    ``c1_every`` sampling, same adjustment-set cap -- so the two are like for
    like. The one behavioural difference is the space itself: every space is
    built with :func:`bkrobust.search.space_fixed.build_space_correct` rather
    than :func:`bkrobust.demo.space.enumerate_space` /
    :func:`bkrobust.core.spacelib.build_space`.

    Results are checkpointed one line per CPDAG to
    ``<out_dir>/n<n>_corrected.jsonl``, flushed and fsynced after every line.
    With ``resume=True`` (default), a prior partial run's lines are read back
    first and their deltas replayed into the running totals -- including the
    running combo counter that drives ``c1_every`` sampling -- before the
    sweep continues from the first CPDAG not yet checkpointed. A kill loses at
    most the one CPDAG in flight.

    When every CPDAG in ``all_cpdags(n)`` has a checkpointed record (whether
    from this call or a prior one), the summary file
    ``<out_dir>/n<n>_corrected.json`` is written -- see
    :func:`write_corrected_summary`.

    Returns:
        ``(totals, counterexamples)`` for the full sweep (replayed lines
        included).
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    jsonl_path = out / f"n{n}_corrected.jsonl"

    cpdags = all_cpdags(n)
    done = _read_checkpoint(jsonl_path) if resume else {}

    totals = CorrectedStudyTotals()
    totals.n_dags = len(all_dags(n))
    totals.n_cpdags = len(cpdags)
    counterexamples: list[dict[str, Any]] = []

    # Replay already-checkpointed records, IN CPDAG ORDER, to reconstruct the
    # running totals (in particular the combo counter that drives c1_every)
    # exactly as if the sweep had run start to finish in one process. `done`
    # is always a contiguous prefix of `cpdags` because the sweep below only
    # ever writes records in that order, so stopping at the first gap and
    # counting replayed records are the same thing.
    n_replayed = 0
    for cpdag in cpdags:
        rec = done.get(cpdag.edge_string())
        if rec is None:
            break
        _add_deltas(totals, rec["deltas"])
        counterexamples.extend(rec["counterexamples"])
        n_replayed += 1

    start = time.time()
    mode = "a" if jsonl_path.exists() and resume else "w"
    if mode == "w" and jsonl_path.exists():
        jsonl_path.unlink()  # resume=False: start clean, never silently append
    with jsonl_path.open(mode) as fh:
        for ci, cpdag in enumerate(cpdags):
            if ci < n_replayed:
                continue
            rec = _process_one_cpdag_main(
                cpdag,
                n=n,
                combo_offset=totals.n_combos,
                max_space=max_space,
                max_undirected=max_undirected,
                c1_every=c1_every,
                all_z=all_z,
                max_z_per_g0=max_z_per_g0,
            )
            _add_deltas(totals, rec["deltas"])
            counterexamples.extend(rec["counterexamples"])
            _append_record(fh, rec)
            if verbose and (ci + 1) % 25 == 0:
                print(
                    f"  n={n} corrected: cpdag {ci + 1}/{len(cpdags)} "
                    f"combos={totals.n_combos} c1v={totals.c1_violations} "
                    f"c2v={totals.c2_violations} [{time.time() - start:.0f}s]",
                    flush=True,
                )

    if len(_read_checkpoint(jsonl_path)) == len(cpdags):
        write_corrected_summary(n, out_dir, runtime_seconds=time.time() - start)

    return totals, counterexamples


def write_corrected_summary(
    n: int,
    out_dir: str | Path,
    *,
    runtime_seconds: float | None = None,
) -> Path:
    """Aggregate ``n<n>_corrected.jsonl`` into the final ``n<n>_corrected.json``.

    Field-compatible with ``results/search/conjectures/n{3,4,5}.json``, plus
    ``n_states_added_vs_old`` and ``n_combos_new``.
    """
    out = Path(out_dir)
    jsonl_path = out / f"n{n}_corrected.jsonl"
    done = _read_checkpoint(jsonl_path)
    cpdags = all_cpdags(n)

    totals = CorrectedStudyTotals()
    totals.n_dags = len(all_dags(n))
    totals.n_cpdags = len(cpdags)
    counterexamples: list[dict[str, Any]] = []
    for cpdag in cpdags:
        rec = done.get(cpdag.edge_string())
        if rec is None:
            continue
        _add_deltas(totals, rec["deltas"])
        counterexamples.extend(rec["counterexamples"])

    payload: dict[str, Any] = {
        "n": n,
        "scope": (
            f"all CPDAGs on {n} nodes with <=6 undirected edges and corrected space <=200 elements"
        ),
        "totals": asdict(totals),
        "counterexamples": counterexamples,
    }
    if runtime_seconds is not None:
        payload["runtime_seconds"] = runtime_seconds
    path = out / f"n{n}_corrected.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True))
    return path


# --- deliverable 2: the n=5 density gap -------------------------------------


def identify_density_gap(
    n: int,
    out_dir: str | Path,
) -> list[dict[str, Any]]:
    """The CPDAGs the main corrected sweep had to skip for size, from its checkpoint.

    Reads ``<out_dir>/n<n>_corrected.jsonl`` (which must already exist -- run
    :func:`run_study_corrected` first) and returns every record whose
    ``skip_reason`` is a density skip (too many undirected edges, or corrected
    space too big) -- i.e. everything except the trivial "no undirected
    edges" skip, which carries no knowledge states to examine at all.
    """
    out = Path(out_dir)
    jsonl_path = out / f"n{n}_corrected.jsonl"
    done = _read_checkpoint(jsonl_path)
    cpdags = all_cpdags(n)
    gap: list[dict[str, Any]] = []
    for cpdag in cpdags:
        rec = done.get(cpdag.edge_string())
        if rec is None:
            continue
        reason = rec.get("skip_reason")
        if reason and reason != "no_undirected_edges":
            gap.append(
                {
                    "cpdag": rec["cpdag"],
                    "k_undirected": rec["k_undirected"],
                    "space_size_new": rec["space_size_new"],
                    "reason": reason,
                }
            )
    return gap


def _process_one_cpdag_dense(
    cpdag_str: str,
    cpdag: MPDAG,
    *,
    n: int,
    combo_offset: int,
    max_space_dense: int,
    c1_every: int,
    all_z: bool,
    max_z_per_g0: int,
) -> dict[str, Any]:
    """One CPDAG's contribution to the dense closure -- no undirected/size skip.

    The only remaining reason to not fully examine a CPDAG here is genuine
    infeasibility: the corrected space is large enough that building the
    O(N^3) covering relation is impractical within this session, which is
    recorded explicitly (never silently dropped) via ``max_space_dense``.
    """
    rec = _record_skeleton(cpdag, n)
    rec["cpdag"] = cpdag_str

    new_elements = enumerate_space_correct(cpdag)
    rec["space_size_new"] = len(new_elements)
    if len(new_elements) > max_space_dense:
        rec["skip_reason"] = (
            f"INFEASIBLE: corrected space has {len(new_elements)} elements; "
            f"the covering relation is built by brute-force enumeration of "
            f"represented-DAG subset comparisons, O(N^3) = "
            f"{len(new_elements) ** 3:.3e} comparisons, which exceeds this "
            f"session's practical time budget (max_space_dense={max_space_dense}). "
            "Left unexamined; see 'infeasible' in the n5_dense.json summary."
        )
        return rec

    old_elements = enumerate_space(cpdag)
    old_strings = frozenset(g.edge_string() for g in old_elements)
    new_strings = frozenset(g.edge_string() for g in new_elements)
    rec["space_size_old"] = len(old_elements)
    n_added = len(new_strings - old_strings)
    rec["n_states_added"] = n_added

    space = build_space_correct(cpdag)
    deltas, counterexamples = _sweep_one_space(
        cpdag,
        space,
        old_strings,
        n=n,
        combo_offset=combo_offset,
        c1_every=c1_every,
        all_z=all_z,
        max_z_per_g0=max_z_per_g0,
    )
    deltas["n_states_added_vs_old"] = n_added
    rec["deltas"] = deltas
    rec["counterexamples"] = counterexamples
    return rec


def run_dense_closure(
    n: int,
    gap: list[dict[str, Any]],
    out_dir: str | Path,
    *,
    max_space_dense: int = 2000,
    c1_every: int = 50,
    all_z: bool = True,
    max_z_per_g0: int = 8,
    verbose: bool = True,
    resume: bool = True,
) -> tuple[CorrectedStudyTotals, list[dict[str, Any]], list[dict[str, Any]]]:
    """Run every CPDAG in ``gap`` to completion, however expensive, or record why not.

    Args:
        n: Number of nodes.
        gap: The CPDAGs to close, from :func:`identify_density_gap` -- each
            item must carry a ``cpdag`` edge-string.
        out_dir: Directory for incremental output
            (``<out_dir>/n<n>_dense.jsonl``).
        max_space_dense: A CPDAG whose corrected space exceeds this many
            elements is recorded as infeasible rather than attempted -- the
            O(N^3) covering-relation build makes anything past a few thousand
            elements impractical in pure Python within a session. Every
            infeasible CPDAG is named, with its size and the reason, never
            silently dropped.
        c1_every: Sample rate for the O(|space|^2) Conjecture-1 guard check.
        all_z: Test every valid adjustment set rather than only ``O*``.
        max_z_per_g0: Cap on adjustment sets tested per ``G0``.
        verbose: Print progress after every CPDAG (there are few of these,
            and each can be slow, unlike the main sweep's every-25 cadence).
        resume: As in :func:`run_study_corrected` -- checkpoint after every
            CPDAG, replay on restart.

    Returns:
        ``(totals, counterexamples, infeasible)`` where ``infeasible`` lists
        the CPDAGs that could not be completed, each with its name, ``k``,
        space size and reason.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    jsonl_path = out / f"n{n}_dense.jsonl"

    # Look CPDAGs up by name against all_cpdags(n) so we get the actual MPDAG
    # object (the checkpoint only stores edge_string), in a stable order.
    by_string = {c.edge_string(): c for c in all_cpdags(n)}
    target_strings = sorted(g["cpdag"] for g in gap)

    done = _read_checkpoint(jsonl_path) if resume else {}
    totals = CorrectedStudyTotals()
    totals.n_dags = len(all_dags(n))
    totals.n_cpdags = len(target_strings)
    counterexamples: list[dict[str, Any]] = []

    n_replayed = 0
    for s in target_strings:
        rec = done.get(s)
        if rec is None:
            break
        _add_deltas(totals, rec["deltas"])
        counterexamples.extend(rec["counterexamples"])
        n_replayed += 1

    start = time.time()
    mode = "a" if jsonl_path.exists() and resume else "w"
    if mode == "w" and jsonl_path.exists():
        jsonl_path.unlink()
    with jsonl_path.open(mode) as fh:
        for ci, s in enumerate(target_strings):
            if ci < n_replayed:
                continue
            cpdag = by_string[s]
            rec = _process_one_cpdag_dense(
                s,
                cpdag,
                n=n,
                combo_offset=totals.n_combos,
                max_space_dense=max_space_dense,
                c1_every=c1_every,
                all_z=all_z,
                max_z_per_g0=max_z_per_g0,
            )
            _add_deltas(totals, rec["deltas"])
            counterexamples.extend(rec["counterexamples"])
            _append_record(fh, rec)
            if verbose:
                status = "INFEASIBLE" if rec["skip_reason"] else "done"
                print(
                    f"  n={n} dense: cpdag {ci + 1}/{len(target_strings)} "
                    f"[{status}] size={rec['space_size_new']} "
                    f"combos={totals.n_combos} c2v={totals.c2_violations} "
                    f"[{time.time() - start:.0f}s]",
                    flush=True,
                )

    all_done = _read_checkpoint(jsonl_path)
    infeasible = [
        {
            "cpdag": r["cpdag"],
            "k_undirected": r["k_undirected"],
            "space_size_new": r["space_size_new"],
            "reason": r["skip_reason"],
        }
        for r in all_done.values()
        if r.get("skip_reason")
    ]
    infeasible.sort(key=lambda r: r["cpdag"])

    if len(all_done) == len(target_strings):
        write_dense_summary(n, out_dir, runtime_seconds=time.time() - start)

    return totals, counterexamples, infeasible


def write_dense_summary(
    n: int,
    out_dir: str | Path,
    *,
    runtime_seconds: float | None = None,
) -> Path:
    """Aggregate ``n<n>_dense.jsonl`` into ``n<n>_dense.json``."""
    out = Path(out_dir)
    jsonl_path = out / f"n{n}_dense.jsonl"
    done = _read_checkpoint(jsonl_path)

    totals = CorrectedStudyTotals()
    totals.n_dags = len(all_dags(n))
    totals.n_cpdags = len(done)
    counterexamples: list[dict[str, Any]] = []
    infeasible: list[dict[str, Any]] = []
    for key in sorted(done):
        rec = done[key]
        if rec.get("skip_reason"):
            infeasible.append(
                {
                    "cpdag": rec["cpdag"],
                    "k_undirected": rec["k_undirected"],
                    "space_size_new": rec["space_size_new"],
                    "reason": rec["skip_reason"],
                }
            )
            continue
        _add_deltas(totals, rec["deltas"])
        counterexamples.extend(rec["counterexamples"])

    payload: dict[str, Any] = {
        "n": n,
        "scope": (
            "the density gap: every CPDAG the main corrected sweep skipped "
            "for size (>6 undirected edges, or corrected space >200 "
            "elements), run to completion with a generous per-CPDAG size "
            "budget; CPDAGs still too large to complete are listed under "
            "'infeasible' rather than silently dropped"
        ),
        "totals": asdict(totals),
        "counterexamples": counterexamples,
        "infeasible": infeasible,
        "n_cpdags_attempted": len(done),
        "n_cpdags_completed": len(done) - len(infeasible),
        "n_cpdags_infeasible": len(infeasible),
    }
    if runtime_seconds is not None:
        payload["runtime_seconds"] = runtime_seconds
    path = out / f"n{n}_dense.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True))
    return path


# --- deliverable 3: old vs new side by side ----------------------------------


def _load_json(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text())


def build_old_vs_new(
    ns: tuple[int, ...] = (3, 4, 5),
    *,
    old_dir: str | Path = "results/search/conjectures",
    new_dir: str | Path = "results/axisb2/conjecture2",
) -> dict[str, Any]:
    """Side-by-side comparison of the OLD (defective-space) and NEW (corrected) sweeps.

    Reads the old numbers from ``old_dir/n<n>.json`` verbatim -- they are
    never recomputed here -- and the new numbers from
    ``new_dir/n<n>_corrected.json`` (and, for n=5, ``new_dir/n5_dense.json``
    if present, folded in so the comparison reflects the CLOSED density gap,
    not just the guarded sweep).

    Returns:
        A dict with one entry per ``n`` plus an overall verdict on whether
        Conjecture 2's status changed.
    """
    old_dir = Path(old_dir)
    new_dir = Path(new_dir)
    per_n: dict[str, Any] = {}
    any_counterexample = False

    for n in ns:
        old_path = old_dir / f"n{n}.json"
        new_path = new_dir / f"n{n}_corrected.json"
        entry: dict[str, Any] = {"n": n}
        if old_path.exists():
            old = _load_json(old_path)
            entry["old"] = {
                "scope": old.get("scope"),
                "totals": old.get("totals"),
                "n_counterexamples": len(old.get("counterexamples", [])),
                "runtime_seconds": old.get("runtime_seconds"),
            }
        else:
            entry["old"] = None

        if new_path.exists():
            new = _load_json(new_path)
            entry["new"] = {
                "scope": new.get("scope"),
                "totals": new.get("totals"),
                "n_counterexamples": len(new.get("counterexamples", [])),
                "runtime_seconds": new.get("runtime_seconds"),
            }
            any_counterexample = any_counterexample or bool(new.get("counterexamples"))
        else:
            entry["new"] = None

        dense_path = new_dir / f"n{n}_dense.json"
        if dense_path.exists():
            dense = _load_json(dense_path)
            entry["dense_closure"] = {
                "scope": dense.get("scope"),
                "totals": dense.get("totals"),
                "n_counterexamples": len(dense.get("counterexamples", [])),
                "n_cpdags_attempted": dense.get("n_cpdags_attempted"),
                "n_cpdags_completed": dense.get("n_cpdags_completed"),
                "n_cpdags_infeasible": dense.get("n_cpdags_infeasible"),
                "infeasible": dense.get("infeasible", []),
            }
            any_counterexample = any_counterexample or bool(dense.get("counterexamples"))
        else:
            entry["dense_closure"] = None

        per_n[str(n)] = entry

    verdict = (
        "COUNTEREXAMPLE FOUND -- Conjecture 2 is FALSE on the corrected space; "
        "see the counterexamples fields above for the witnesses."
        if any_counterexample
        else (
            "No counterexample found on the corrected space either, across every "
            "n=3,4,5 combination examined (including the closed n=5 density gap "
            "where completed). Conjecture 2's status is UNCHANGED: still an "
            "empirically-supported conjecture, not a theorem, and the scope of "
            "that support is stated exactly (see any 'infeasible' entries above "
            "for what remains unexamined)."
        )
    )

    return {"per_n": per_n, "verdict": verdict}


def write_old_vs_new(
    ns: tuple[int, ...] = (3, 4, 5),
    *,
    old_dir: str | Path = "results/search/conjectures",
    new_dir: str | Path = "results/axisb2/conjecture2",
) -> Path:
    """Compute and write ``<new_dir>/old_vs_new.json``."""
    comparison = build_old_vs_new(ns, old_dir=old_dir, new_dir=new_dir)
    out = Path(new_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "old_vs_new.json"
    path.write_text(json.dumps(comparison, indent=2, sort_keys=True))
    return path


# --- reproducing a recorded counterexample from its fields alone ------------


def reproduce_conjecture2(rec: dict[str, Any]) -> Conjecture2Result:
    """Recompute a Conjecture 2 record's result from only its recorded fields.

    Rebuilds the corrected space from ``cpdag`` and re-derives ``r_full``,
    ``r_up`` and both witnesses independently from ``g0``, ``x``, ``y`` and
    ``z`` alone -- the check that a saved counterexample (or any other saved
    combination) depends on nothing beyond what was recorded.

    Args:
        rec: A dict carrying ``cpdag`` and ``g0`` as :class:`MPDAG` objects
            (not edge strings -- an edge string alone cannot be parsed back
            into a graph without its node set, so callers look the objects up
            via :func:`bkrobust.search.conjecture_study.all_cpdags` first),
            plus ``x``, ``y`` (node names) and ``z`` (an iterable of node
            names).

    Returns:
        A freshly computed :class:`Conjecture2Result`.
    """
    cpdag: MPDAG = rec["cpdag"]
    g0: MPDAG = rec["g0"]
    x: str = rec["x"]
    y: str = rec["y"]
    z = frozenset(rec["z"])
    space = build_space_correct(cpdag)
    dists = distances_from(space, g0)

    def fails(g: MPDAG, _z: frozenset[str] = z, _x: str = x, _y: str = y) -> bool:
        return not is_valid(_z, g, _x, _y)

    return check_conjecture2(space, g0, fails, dists)
