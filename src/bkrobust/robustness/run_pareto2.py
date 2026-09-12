r"""Driver for the Phase 2 continuation run: P6' (achievable efficiency) and
P7' (aggressive background knowledge on real networks).

Usage (always ``PYTHONPATH=src``):

    python3 -m bkrobust.robustness.run_pareto2 --mode smoke
    python3 -m bkrobust.robustness.run_pareto2 --mode census
    python3 -m bkrobust.robustness.run_pareto2 --mode full
    python3 -m bkrobust.robustness.run_pareto2 --mode single --base-index 0 \
        --out /tmp/check.csv

See ``results/axis_robustness/PREREGISTRATION.md`` Appendix A.3 and this
project's ``PARETO2_RUN_NOTES.md`` for the design and the acceptance-criteria
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

from bkrobust.core.resultsio import ResultWriter, write_manifest
from bkrobust.robustness import pareto as P
from bkrobust.robustness import pareto2 as P2

RESULTS_DIR = Path("results/axis_robustness")
ROOT_SEED = 0
SCRATCH = Path(
    "/private/tmp/claude-502/-Users-ahn-Documents-Research-iclr27-bkrobust/"
    "616b70e8-9605-41df-a870-61f608f1764b/scratchpad"
)


def build_grid(limit_bases: int | None = None) -> list[P.BaseCpdag]:
    """The identical base-CPDAG grid round 1 used (same root seed, same
    selection function), so Task A strata here line up one-to-one with
    round 1's 24 strata."""
    return P.select_base_cpdags(root_seed=ROOT_SEED, limit=limit_bases)


def sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_rows(path: Path, rows: list[dict[str, Any]]) -> Path:
    if path.exists():
        path.unlink()
    with ResultWriter(path) as w:
        for row in rows:
            w.write(row)
    return path


# ---------------------------------------------------------------------------
# smoke
# ---------------------------------------------------------------------------


def cmd_smoke(args: argparse.Namespace) -> int:
    bases = build_grid(limit_bases=1)
    if not bases:
        print("SMOKE: no base CPDAG admitted at all -- cannot even smoke-test.")
        return 1
    base = bases[0]
    prows, mrows, info = P2.run_one_base_task_a(base, ROOT_SEED)
    print(f"SMOKE (Task A) base_id={info['base_id']}")
    print(f"  undirected edges in Cpdag: {info['n_undirected_cpdag']}")
    print(f"  menu members: {info['n_menu_members']}")
    print(f"  proposals sampled: {info['n_proposals']}")
    print(f"  status counts: ok={info['n_ok']} no_menu_member_valid={info['n_no_menu_member_valid']} "
          f"avar_undefined={info['n_avar_undefined']} proposal_contradictory={info['n_proposal_contradictory']}")
    print(f"  distinct achievable_utility_median: {info['n_distinct_achievable_utility']} "
          f"(varies={info['utility_varies']})")
    print(f"  distinct winning_member_index: {info['n_distinct_winning_member']}")
    print(f"  C2 max|bias| over GAC-valid menu members x SEM draws: {info['max_abs_bias_check']:.3e}")
    print(f"  wall time: {info['elapsed_s']:.2f} s")

    full_grid = build_grid(limit_bases=None)
    n_full = len(full_grid)
    projected_a = info["elapsed_s"] * n_full
    print(f"  full Task A grid would be {n_full} base CPDAGs; projected time = "
          f"{projected_a:.1f} s ({projected_a / 60:.1f} min)")

    SCRATCH.mkdir(parents=True, exist_ok=True)
    write_rows(SCRATCH / "pareto2_proposals_smoke.csv", prows)
    write_rows(SCRATCH / "pareto2_menus_smoke.csv", mrows)
    print(f"  smoke CSVs written under {SCRATCH} for inspection (not final outputs)")

    # Cheap Task B smoke: load networks + census one small network, timed
    # separately since it must never call fast_gate.
    t0 = time.perf_counter()
    networks, load_exclusions = P2.load_networks()
    t_load = time.perf_counter() - t0
    print(f"  Task B: loaded {len(networks)} networks ({len(load_exclusions)} excluded) in {t_load:.2f}s")
    if networks:
        one = sorted(networks, key=lambda n: n.name)[0]
        t1 = time.perf_counter()
        row, on_path = P2.census_one_network(one, ROOT_SEED)
        t_census = time.perf_counter() - t1
        print(f"  Task B: census of one network ({one.name}) took {t_census:.3f}s -> {row}")
    return 0


# ---------------------------------------------------------------------------
# single (determinism check, C3)
# ---------------------------------------------------------------------------


def cmd_single(args: argparse.Namespace) -> int:
    bases = build_grid(limit_bases=args.base_index + 1)
    if len(bases) <= args.base_index:
        print(f"base index {args.base_index} not available (only {len(bases)} bases found)")
        return 1
    base = bases[args.base_index]
    prows, mrows, info = P2.run_one_base_task_a(base, ROOT_SEED)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    write_rows(out_path, prows)
    print(f"single-base run: base_id={info['base_id']} -> {out_path}")
    print(f"sha256={sha256_of(out_path)}")
    return 0


# ---------------------------------------------------------------------------
# census only
# ---------------------------------------------------------------------------


def cmd_census(args: argparse.Namespace) -> int:
    t0 = time.perf_counter()
    networks, load_exclusions = P2.load_networks()
    rows, on_path_by_network = P2.run_census(networks, ROOT_SEED)
    elapsed = time.perf_counter() - t0
    print(f"CENSUS: {len(networks)} networks loaded, {len(load_exclusions)} excluded, {elapsed:.2f}s")
    for r in load_exclusions:
        print(f"  excluded: {r}")
    for r in sorted(rows, key=lambda r: r["network"]):
        print(f"  {r['network']}: n_nodes={r['n_nodes']} n_undirected={r['n_undirected_edges']} "
              f"n_sampled={r['n_pairs_sampled']} n_on_path={r['n_pairs_with_on_path_undirected_edge']}")
    total_on_path = sum(r["n_pairs_with_on_path_undirected_edge"] for r in rows)
    print(f"  total on-path pairs across corpus: {total_on_path}")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    census_path = write_rows(RESULTS_DIR / "pareto2_aggressive_census.csv", rows)
    print(f"wrote {census_path} ({len(rows)} rows)")
    return 0


# ---------------------------------------------------------------------------
# full
# ---------------------------------------------------------------------------


def cmd_full(args: argparse.Namespace) -> int:
    t0 = time.perf_counter()
    n_bases = args.limit_bases
    bases = build_grid(limit_bases=n_bases)
    print(f"FULL run: Task A over {len(bases)} base CPDAGs")

    all_proposal_rows: list[dict[str, Any]] = []
    all_menu_rows: list[dict[str, Any]] = []
    per_base_info: list[dict[str, Any]] = []

    for base in bases:
        prows, mrows, info = P2.run_one_base_task_a(base, ROOT_SEED)
        per_base_info.append(info)
        all_proposal_rows.extend(prows)
        all_menu_rows.extend(mrows)
        print(
            f"  {info['base_id']}: {info['n_proposals']} proposals, "
            f"ok={info['n_ok']}, distinct_utility={info['n_distinct_achievable_utility']}, "
            f"C2 max|bias|={info['max_abs_bias_check']:.2e}, {info['elapsed_s']:.2f}s"
        )

    c2_failures = [i for i in per_base_info if i["max_abs_bias_check"] >= 1e-8]
    if c2_failures:
        print("C2 MENU SANITY CHECK FAILED for the following bases -- STOPPING, not writing final outputs:")
        for f in c2_failures:
            print(f"  {f['base_id']}: max|bias|={f['max_abs_bias_check']:.3e}")
        return 2

    n_strata_varying = sum(1 for i in per_base_info if i["utility_varies"])
    print(f"DECISIVE CHECK (Task A / P6'): strata with >1 distinct achievable_utility_median: "
          f"{n_strata_varying}/{len(per_base_info)}")

    print("Task B: loading real networks and running the structural census")
    tb0 = time.perf_counter()
    networks, load_exclusions = P2.load_networks()
    census_rows, on_path_by_network = P2.run_census(networks, ROOT_SEED)
    tb_census = time.perf_counter() - tb0
    total_on_path = sum(r["n_pairs_with_on_path_undirected_edge"] for r in census_rows)
    print(f"  {len(networks)} networks loaded ({len(load_exclusions)} excluded), "
          f"census in {tb_census:.2f}s; total on-path pairs = {total_on_path}")

    followup_prows: list[dict[str, Any]] = []
    followup_mrows: list[dict[str, Any]] = []
    gate_exclusions: list[dict[str, Any]] = []
    followup_info: dict[str, Any] = {}
    followup_targets: list[tuple[str, list[tuple[str, str]]]] = []
    if total_on_path > 0:
        followup_targets = P2.select_followup_targets(census_rows, on_path_by_network)
        print(f"  follow-up: gating pairs from {len(followup_targets)} smallest on-path-bearing networks: "
              f"{[t[0] for t in followup_targets]}")
        networks_by_name = {n.name: n for n in networks}
        followup_prows, followup_mrows, gate_exclusions, followup_info = P2.run_followup(
            networks_by_name, followup_targets, ROOT_SEED
        )
        print(f"  follow-up: {followup_info['n_bases_measured']} (network,x,y) bases measured, "
              f"{followup_info['n_gate_excluded']} gate-excluded, {followup_info['elapsed_s']:.2f}s")
    else:
        print("  no network in the corpus has any on-path undirected edge: P7' is structurally "
              "unavailable across the whole corpus. Not running the follow-up stage.")

    all_proposal_rows.extend(followup_prows)
    all_menu_rows.extend(followup_mrows)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    prop_path = write_rows(RESULTS_DIR / "pareto2_proposals.csv", all_proposal_rows)
    menu_path = write_rows(RESULTS_DIR / "pareto2_menus.csv", all_menu_rows)
    census_path = write_rows(RESULTS_DIR / "pareto2_aggressive_census.csv", census_rows)

    total_elapsed = time.perf_counter() - t0

    scratch_manifest_dir = SCRATCH / "manifest2_tmp"
    scratch_manifest_dir.mkdir(parents=True, exist_ok=True)
    grid = {
        "task_a": {
            "n_bases": len(bases),
            "proposal_cap": P.PROPOSAL_CAP,
            "n_sem_draws": P.N_SEM_DRAWS,
            "radius_time_limit_s": P.RADIUS_TIME_LIMIT_S,
            "menu_cap": P2.MENU_CAP,
            "menu_perturb_cap": P2.MENU_PERTURB_CAP,
        },
        "task_b": {
            "networks_dir": str(P2.NETWORKS_DIR),
            "excluded_files": sorted(P2.EXCLUDED_NETWORK_FILES),
            "census_pair_cap": P2.CENSUS_PAIR_CAP,
            "followup_max_networks": P2.FOLLOWUP_MAX_NETWORKS,
            "followup_pairs_per_network": P2.FOLLOWUP_PAIRS_PER_NETWORK,
        },
    }
    extra = {
        "gate": P2.GATE_NAME,
        "n_proposal_rows": len(all_proposal_rows),
        "n_menu_rows": len(all_menu_rows),
        "n_census_rows": len(census_rows),
        "total_elapsed_s": total_elapsed,
        "task_a_per_base": per_base_info,
        "task_a_n_strata_varying": n_strata_varying,
        "task_a_n_strata_total": len(per_base_info),
        "task_b_load_exclusions": load_exclusions,
        "task_b_total_on_path_pairs": total_on_path,
        "task_b_followup_targets": [t[0] for t in followup_targets],
        "task_b_followup_gate_exclusions": gate_exclusions,
        "task_b_followup_info": followup_info,
    }
    manifest_tmp_path = write_manifest(scratch_manifest_dir, seed=ROOT_SEED, grid=grid, extra=extra)
    manifest_final_path = RESULTS_DIR / "pareto2_manifest.json"
    shutil.copyfile(manifest_tmp_path, manifest_final_path)

    print(f"wrote {prop_path} ({len(all_proposal_rows)} rows)")
    print(f"wrote {menu_path} ({len(all_menu_rows)} rows)")
    print(f"wrote {census_path} ({len(census_rows)} rows)")
    print(f"wrote {manifest_final_path}")
    print(f"total wall time: {total_elapsed:.1f} s")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["smoke", "census", "full", "single"], required=True)
    ap.add_argument("--limit-bases", type=int, default=None)
    ap.add_argument("--base-index", type=int, default=0)
    ap.add_argument("--out", type=str, default="/tmp/pareto2_single.csv")
    args = ap.parse_args(argv)

    if args.mode == "smoke":
        return cmd_smoke(args)
    if args.mode == "single":
        return cmd_single(args)
    if args.mode == "census":
        return cmd_census(args)
    return cmd_full(args)


if __name__ == "__main__":
    sys.exit(main())
