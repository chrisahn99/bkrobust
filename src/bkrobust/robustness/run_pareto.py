r"""Driver for the Phase 2 Pareto-frontier run.

Usage (always ``PYTHONPATH=src``):

    python3 -m bkrobust.robustness.run_pareto --mode smoke
    python3 -m bkrobust.robustness.run_pareto --mode full
    python3 -m bkrobust.robustness.run_pareto --mode single --base-index 0 \
        --out results/axis_robustness/pareto_proposals_check.csv

See ``results/axis_robustness/PREREGISTRATION.md`` Section 9 and this
project's PARETO_RUN_NOTES.md for the design and the acceptance-criteria
evidence gathered by running this script.
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

from bkrobust.core.resultsio import ResultWriter, write_manifest
from bkrobust.robustness import pareto as P

RESULTS_DIR = Path("results/axis_robustness")
ROOT_SEED = 0


def build_grid(limit_bases: int | None = None) -> list[P.BaseCpdag]:
    return P.select_base_cpdags(root_seed=ROOT_SEED, limit=limit_bases)


def run_one_base(base: P.BaseCpdag) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """Process a single base CPDAG end to end. Returns (proposal_rows, sem_rows, info)."""
    t0 = time.perf_counter()
    sems = P.draw_sems(base, ROOT_SEED, n=P.N_SEM_DRAWS)

    b2_biases = P.b2_sanity_check(base, sems)
    b2_max_abs = max(abs(b) for b in b2_biases)
    b2_ok = b2_max_abs < 1e-8

    menu = P.truthful_orientation_menu(base)
    desc_x = base.dag_obj.descendants(base.treatment) | {base.treatment}
    anc_y = base.dag_obj.ancestors(base.outcome) | {base.outcome}
    proposals = P.sample_proposals(base, menu, desc_x, anc_y, root_seed=ROOT_SEED, cap=P.PROPOSAL_CAP)

    proposal_rows: list[dict[str, Any]] = []
    sem_rows: list[dict[str, Any]] = []
    for prop in proposals:
        result = P.measure_proposal(base, prop, sems, root_seed=ROOT_SEED)
        proposal_rows.append(result.row)
        sem_rows.extend(result.sem_rows)

    P.compute_fronts(proposal_rows)

    elapsed = time.perf_counter() - t0
    info = {
        "base_id": base.base_id,
        "n_undirected": len(menu),
        "n_proposals": len(proposals),
        "b2_max_abs_bias": b2_max_abs,
        "b2_ok": b2_ok,
        "elapsed_s": elapsed,
    }
    return proposal_rows, sem_rows, info


def write_csvs(
    out_dir: Path,
    proposals_name: str,
    sems_name: str,
    all_proposal_rows: list[dict[str, Any]],
    all_sem_rows: list[dict[str, Any]],
) -> tuple[Path, Path]:
    prop_path = out_dir / proposals_name
    sem_path = out_dir / sems_name
    if prop_path.exists():
        prop_path.unlink()
    if sem_path.exists():
        sem_path.unlink()
    with ResultWriter(prop_path) as w:
        for row in all_proposal_rows:
            w.write(row)
    with ResultWriter(sem_path) as w:
        for row in all_sem_rows:
            w.write(row)
    return prop_path, sem_path


def sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cmd_smoke(args: argparse.Namespace) -> int:
    bases = build_grid(limit_bases=1)
    if not bases:
        print("SMOKE: no base CPDAG admitted at all -- cannot even smoke-test.")
        return 1
    base = bases[0]
    proposal_rows, sem_rows, info = run_one_base(base)
    print(f"SMOKE base_id={info['base_id']}")
    print(f"  undirected edges in Cpdag: {info['n_undirected']}")
    print(f"  proposals sampled: {info['n_proposals']}")
    print(f"  B2 sanity: max |bias(true optimal set)| over {P.N_SEM_DRAWS} SEM draws = {info['b2_max_abs_bias']:.3e} ({'OK' if info['b2_ok'] else 'FAIL'})")
    print(f"  wall time: {info['elapsed_s']:.2f} s")
    full_grid = build_grid(limit_bases=None)
    n_full = len(full_grid)
    projected = info["elapsed_s"] * n_full
    print(f"  full grid would be {n_full} base CPDAGs; projected time = {projected:.1f} s ({projected/60:.1f} min)")
    scratch = Path("/private/tmp/claude-502/-Users-ahn-Documents-Research-iclr27-bkrobust/616b70e8-9605-41df-a870-61f608f1764b/scratchpad")
    scratch.mkdir(parents=True, exist_ok=True)
    write_csvs(scratch, "pareto_proposals_smoke.csv", "pareto_sems_smoke.csv", proposal_rows, sem_rows)
    print(f"  smoke CSVs written under {scratch} for inspection (not final outputs)")
    return 0


def cmd_single(args: argparse.Namespace) -> int:
    bases = build_grid(limit_bases=args.base_index + 1)
    if len(bases) <= args.base_index:
        print(f"base index {args.base_index} not available (only {len(bases)} bases found)")
        return 1
    base = bases[args.base_index]
    proposal_rows, sem_rows, info = run_one_base(base)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists():
        out_path.unlink()
    with ResultWriter(out_path) as w:
        for row in proposal_rows:
            w.write(row)
    print(f"single-base run: base_id={info['base_id']} -> {out_path}")
    print(f"sha256={sha256_of(out_path)}")
    return 0


def cmd_full(args: argparse.Namespace) -> int:
    t0 = time.perf_counter()
    n_bases = args.limit_bases
    bases = build_grid(limit_bases=n_bases)
    print(f"FULL run: {len(bases)} base CPDAGs")

    all_proposal_rows: list[dict[str, Any]] = []
    all_sem_rows: list[dict[str, Any]] = []
    per_base_info: list[dict[str, Any]] = []
    b2_failures: list[dict[str, Any]] = []

    for base in bases:
        proposal_rows, sem_rows, info = run_one_base(base)
        per_base_info.append(info)
        if not info["b2_ok"]:
            b2_failures.append(info)
        all_proposal_rows.extend(proposal_rows)
        all_sem_rows.extend(sem_rows)
        print(
            f"  {info['base_id']}: {info['n_proposals']} proposals, "
            f"B2 max|bias|={info['b2_max_abs_bias']:.2e}, {info['elapsed_s']:.2f}s"
        )

    if b2_failures:
        print("B2 SANITY CHECK FAILED for the following bases -- STOPPING, not writing final outputs:")
        for f in b2_failures:
            print(f"  {f['base_id']}: max|bias|={f['b2_max_abs_bias']:.3e}")
        return 2

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    prop_path, sem_path = write_csvs(
        RESULTS_DIR, "pareto_proposals.csv", "pareto_sems.csv", all_proposal_rows, all_sem_rows
    )

    total_elapsed = time.perf_counter() - t0

    # manifest: write_manifest() hardcodes "manifest.json" in the target
    # directory, which is the OTHER worker's file. Write it into a scratch
    # directory and move the result to pareto_manifest.json instead.
    scratch = Path("/private/tmp/claude-502/-Users-ahn-Documents-Research-iclr27-bkrobust/616b70e8-9605-41df-a870-61f608f1764b/scratchpad/manifest_tmp")
    scratch.mkdir(parents=True, exist_ok=True)
    grid = {
        "component_sizes": [6, 8, 10],
        "seps_per_size": 4,
        "accepted_per_cell": 2,
        "n_bases": len(bases),
        "proposal_cap": P.PROPOSAL_CAP,
        "n_sem_draws": P.N_SEM_DRAWS,
        "n_bias_reps": P.N_BIAS_REPS,
        "bias_flip_rate": P.BIAS_FLIP_RATE,
        "radius_time_limit_s": P.RADIUS_TIME_LIMIT_S,
    }
    extra = {
        "gate": P.GATE_NAME,
        "n_proposal_rows": len(all_proposal_rows),
        "n_sem_rows": len(all_sem_rows),
        "total_elapsed_s": total_elapsed,
        "per_base": per_base_info,
    }
    manifest_tmp_path = write_manifest(scratch, seed=ROOT_SEED, grid=grid, extra=extra)
    manifest_final_path = RESULTS_DIR / "pareto_manifest.json"
    shutil.copyfile(manifest_tmp_path, manifest_final_path)

    print(f"wrote {prop_path} ({len(all_proposal_rows)} rows)")
    print(f"wrote {sem_path} ({len(all_sem_rows)} rows)")
    print(f"wrote {manifest_final_path}")
    print(f"total wall time: {total_elapsed:.1f} s")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["smoke", "full", "single"], required=True)
    ap.add_argument("--limit-bases", type=int, default=None)
    ap.add_argument("--base-index", type=int, default=0)
    ap.add_argument("--out", type=str, default="results/axis_robustness/pareto_proposals_single.csv")
    args = ap.parse_args(argv)

    if args.mode == "smoke":
        return cmd_smoke(args)
    if args.mode == "single":
        return cmd_single(args)
    return cmd_full(args)


if __name__ == "__main__":
    sys.exit(main())
