"""Driver for the hops-axis (P5) check: cost census, matched subsample, dual-axis
sampling, and the P5 sign-agreement test between `d_claims` and `d_hops`.

See :mod:`bkrobust.robustness.hops` for the corruption/scoring logic and
``results/axis_robustness/PREREGISTRATION.md`` §2, §6 (P5) for the design this
must not deviate from. This module is orchestration and I/O only, and owns
only files prefixed ``hops_`` under ``results/axis_robustness/``.

Usage::

    PYTHONPATH=src python3 -m bkrobust.robustness.run_hops census \
        --sizes 4,5,6,7,8,9,10 --out-dir results/axis_robustness

    PYTHONPATH=src python3 -m bkrobust.robustness.run_hops sweep \
        --component-sizes 3,4,5,6 --coverages 0.5,1.0 \
        --base-wrongness 0.0,0.10,0.25 --instances-per-cell 15 \
        --seed-budget 150 --n-reps 100 --out-dir results/axis_robustness

No global RNG is touched anywhere in this module -- every draw happens inside
:mod:`bkrobust.robustness.hops`, which threads explicit
``np.random.Generator`` objects throughout (via ``survival.derived_seed``).
"""

from __future__ import annotations

import argparse
import json
import shutil
import tempfile
import time
from pathlib import Path
from typing import Any

from scipy.stats import kendalltau

from bkrobust.core.conventions import UNREACHED
from bkrobust.core.resultsio import ResultWriter, write_manifest
from bkrobust.core.spacelib import distances_from
from bkrobust.robustness import hops as hp
from bkrobust.robustness import survival as sv

CENSUS_FIELDS = [
    "component_size",
    "separation",
    "coverage",
    "triangle_prob",
    "seeds_tried",
    "n_undirected_cpdag",
    "n_space_elements",
    "build_seconds",
    "status",
]

SAMPLES_FIELDS = [
    "instance_id",
    "component_size",
    "separation",
    "coverage",
    "base_wrongness",
    "gate",
    "d_claims",
    "rep",
    "seed",
    "status",
    "gac_status",
    "in_space",
    "d_hops",
    "symdiff_proxy_not_distance",
    "survived",
    "g_equals_g0",
]

INSTANCES_FIELDS = [
    "instance_id",
    "component_size",
    "separation",
    "coverage",
    "base_wrongness",
    "seed",
    "x",
    "y",
    "n_k",
    "z_size",
    "z_star",
    "s0_ok",
    "r_val",
    "r_status",
    "r_method",
    "r_oracle",
    "r_exact",
    "r_assumes",
    "r_witness",
    "r_search_seconds",
    "r_ladder_seconds",
    "r_total_seconds",
    "wall_until_timeout_s",
    "space_n_elements",
    "space_build_seconds",
    "n_samples",
    "n_ok",
    "n_unreachable_in_space",
    "n_not_in_space_bug",
    "n_contradictory",
    "h_max",
    "AUC_claims",
    "AUC_abs_claims",
    "n_d_claims_evaluated",
    "AUC_hops",
    "n_d_hops_evaluated",
    "n_symdiff_compared",
    "n_symdiff_eq_dhops",
    "n_symdiff_gt_dhops",
    "n_symdiff_lt_dhops",
    "frac_symdiff_eq_dhops",
    "frac_symdiff_gt_dhops",
    "frac_symdiff_lt_dhops",
]

P5_FIELDS = [
    "stratum",
    "coverage",
    "base_wrongness",
    "n_instances",
    "n_tau_claims_defined",
    "tau_claims",
    "tau_claims_pvalue",
    "n_tau_hops_defined",
    "tau_hops",
    "tau_hops_pvalue",
    "sign_agree",
    "verdict",
]


def _fill(row: dict[str, Any], fields: list[str]) -> dict[str, Any]:
    return {k: row.get(k, "") for k in fields}


def _parse_float_list(s: str) -> list[float]:
    return [float(x) for x in s.split(",") if x.strip() != ""]


def _parse_int_list(s: str) -> list[int]:
    return [int(x) for x in s.split(",") if x.strip() != ""]


def write_hops_manifest(out_dir: Path, *, seed: int, grid: dict[str, Any], extra: dict[str, Any]) -> Path:
    """``write_manifest`` hardcodes ``manifest.json``, owned by another worker.

    Writes to a temp dir, then renames/copies to ``hops_manifest.json`` in
    ``out_dir`` -- never touches ``out_dir/manifest.json``.
    """
    tmp = Path(tempfile.mkdtemp(prefix="hops_manifest_"))
    write_manifest(tmp, seed=seed, grid=grid, extra=extra)
    dest = out_dir / "hops_manifest.json"
    shutil.move(str(tmp / "manifest.json"), str(dest))
    shutil.rmtree(tmp, ignore_errors=True)
    return dest


# --- census -------------------------------------------------------------------


def cmd_census(args: argparse.Namespace) -> int:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sizes = _parse_int_list(args.sizes)
    path = out_dir / "hops_cost_census.csv"

    rows = []
    with ResultWriter(path, resume=args.resume) as w:
        for size in sizes:
            row = hp.census_one_size(
                size,
                seed_budget=args.seed_budget,
                triangle_prob=args.triangle_prob,
                timeout_s=args.timeout_s,
            )
            rows.append(row)
            w.write(_fill(row, CENSUS_FIELDS))
            print(json.dumps(row, default=str))

    print(f"wrote {path}")
    return 0


# --- sweep (steps 2-4) ---------------------------------------------------------


def run_cell(
    component_size: int,
    separation: int,
    coverage: float,
    base_wrongness: float,
    *,
    instances_per_cell: int,
    seed_budget: int,
    n_reps: int,
    r_time_limit_s: float,
    samples_w: ResultWriter,
    instances_w: ResultWriter,
    exclusions: dict[str, int],
    bug_rows: list[dict[str, Any]],
    instance_records: list[dict[str, Any]],
    dclaims_dhops_pool: list[dict[str, Any]],
) -> int:
    """Fill one flip-arm cell up to quota, sampling both axes for each instance."""
    n_accepted = 0
    for seed in range(seed_budget):
        if n_accepted >= instances_per_cell:
            break
        inst, reason = sv.build_flip_instance(component_size, separation, coverage, base_wrongness, seed)
        if inst is None:
            exclusions[reason] = exclusions.get(reason, 0) + 1
            continue
        n_accepted += 1

        s0_ok = sv.is_gac_valid_mpdag(inst["g0"], inst["x"], inst["y"], inst["z_star"])

        space, space_build_seconds = hp.build_space_for_instance(inst["cpdag"])
        space_elements = frozenset(space.elements)
        g0 = inst["g0"]
        g0_in_space = g0 in space_elements
        if not g0_in_space:
            bug_rows.append({"instance_id": inst["instance_id"], "what": "g0_not_in_own_space"})
        dist = distances_from(space, g0)
        d_hops_of_g0 = dist.get(g0)

        n_k = inst["n_k"]
        rows: list[dict[str, Any]] = []
        for d in range(1, n_k + 1):
            for rep in range(n_reps):
                row = hp.sample_dual_axis_state(
                    inst, d, rep, space_elements=space_elements, dist=dist
                )
                rows.append(row)
                samples_w.write(_fill(row, SAMPLES_FIELDS))
                if row["status"] == "not_in_space_BUG":
                    bug_rows.append(
                        {"instance_id": inst["instance_id"], "what": "corrupted_state_not_in_space",
                         "d_claims": d, "rep": rep}
                    )
                dclaims_dhops_pool.append(
                    {"d_claims": row["d_claims"], "status": row["status"], "d_hops": row["d_hops"]}
                )

        n_ok = sum(1 for r in rows if r["status"] == "ok")
        n_unreach = sum(1 for r in rows if r["status"] == "unreachable_in_space")
        n_bug = sum(1 for r in rows if r["status"] == "not_in_space_BUG")
        n_contra = sum(1 for r in rows if r["status"] == "corrupted_k_contradictory")

        curve_c = hp.claims_curve(rows)
        auc_c = sv.auc_frac(curve_c, n_k)
        auc_c_abs = sv.auc_abs(curve_c)
        curve_h = hp.hops_curve(rows)
        auc_h, h_max = hp.auc_hops(curve_h)
        gap = hp.symdiff_gap_stats(rows)

        rv = sv.compute_r_val(inst["cpdag"], g0, inst["x"], inst["y"], inst["z_star"], time_limit_s=r_time_limit_s)

        irow: dict[str, Any] = {
            "instance_id": inst["instance_id"],
            "component_size": component_size,
            "separation": separation,
            "coverage": coverage,
            "base_wrongness": base_wrongness,
            "seed": seed,
            "x": inst["x"],
            "y": inst["y"],
            "n_k": n_k,
            "z_size": inst["z_size"],
            "z_star": ",".join(sorted(inst["z_star"])),
            "s0_ok": s0_ok,
            **rv,
            "space_n_elements": len(space),
            "space_build_seconds": space_build_seconds,
            "n_samples": len(rows),
            "n_ok": n_ok,
            "n_unreachable_in_space": n_unreach,
            "n_not_in_space_bug": n_bug,
            "n_contradictory": n_contra,
            "h_max": h_max,
            "AUC_claims": auc_c,
            "AUC_abs_claims": auc_c_abs,
            "n_d_claims_evaluated": len(curve_c),
            "AUC_hops": auc_h,
            "n_d_hops_evaluated": len(curve_h),
            **gap,
        }
        instances_w.write(_fill(irow, INSTANCES_FIELDS))

        instance_records.append(
            {
                "instance_id": inst["instance_id"],
                "component_size": component_size,
                "separation": separation,
                "coverage": coverage,
                "base_wrongness": base_wrongness,
                "s0_ok": s0_ok,
                "d_hops_of_g0": d_hops_of_g0,
                "g0_in_space": g0_in_space,
                "r_val": rv.get("r_val"),
                "r_status": rv.get("r_status"),
                "AUC_claims": auc_c,
                "AUC_hops": auc_h,
                "h_max": h_max,
            }
        )
    return n_accepted


def cmd_sweep(args: argparse.Namespace) -> int:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    component_sizes = _parse_int_list(args.component_sizes)
    coverages = _parse_float_list(args.coverages)
    base_wrongness = _parse_float_list(args.base_wrongness)

    samples_path = out_dir / "hops_samples.csv"
    instances_path = out_dir / "hops_instances.csv"

    exclusions: dict[str, int] = {}
    bug_rows: list[dict[str, Any]] = []
    instance_records: list[dict[str, Any]] = []
    dclaims_dhops_pool: list[dict[str, Any]] = []
    cell_counts: dict[str, int] = {}

    t_start = time.perf_counter()

    with ResultWriter(samples_path, resume=args.resume) as samples_w, ResultWriter(
        instances_path, resume=args.resume
    ) as instances_w:
        for c in sorted(component_sizes):
            seps = sorted(set(sv.separations_for(c)))
            for s in seps:
                for cov in sorted(coverages):
                    for b in sorted(base_wrongness):
                        n = run_cell(
                            c,
                            s,
                            cov,
                            b,
                            instances_per_cell=args.instances_per_cell,
                            seed_budget=args.seed_budget,
                            n_reps=args.n_reps,
                            r_time_limit_s=args.r_time_limit_s,
                            samples_w=samples_w,
                            instances_w=instances_w,
                            exclusions=exclusions,
                            bug_rows=bug_rows,
                            instance_records=instance_records,
                            dclaims_dhops_pool=dclaims_dhops_pool,
                        )
                        cell_counts[f"c{c}_s{s}_cov{cov}_b{b}"] = n
                        print(
                            f"cell c={c} s={s} cov={cov} b={b}: accepted={n} "
                            f"elapsed={time.perf_counter() - t_start:.1f}s"
                        )

    elapsed = time.perf_counter() - t_start

    # --- P5: per-stratum Kendall tau-b -----------------------------------------
    p5_path = out_dir / "hops_p5.csv"
    strata: dict[tuple[float, float], list[dict[str, Any]]] = {}
    for rec in instance_records:
        key = (rec["coverage"], rec["base_wrongness"])
        strata.setdefault(key, []).append(rec)

    def _tau(pairs: list[tuple[float, float]]) -> tuple[Any, Any]:
        if len(pairs) < 4:
            return "", ""
        xs = [p[0] for p in pairs]
        ys = [p[1] for p in pairs]
        if len(set(xs)) < 2 or len(set(ys)) < 2:
            return "", ""
        res = kendalltau(xs, ys, variant="b")
        return res.correlation, res.pvalue

    p5_rows: list[dict[str, Any]] = []
    all_pairs_claims: list[tuple[float, float]] = []
    all_pairs_hops: list[tuple[float, float]] = []

    for key in sorted(strata):
        cov, bw = key
        recs = strata[key]
        pairs_claims = [
            (r["AUC_claims"], r["r_val"])
            for r in recs
            if r["AUC_claims"] != "" and r["r_val"] != "" and r["r_val"] != UNREACHED
        ]
        pairs_hops = [
            (r["AUC_hops"], r["r_val"])
            for r in recs
            if r["AUC_hops"] != "" and r["r_val"] != "" and r["r_val"] != UNREACHED
        ]
        all_pairs_claims.extend(pairs_claims)
        all_pairs_hops.extend(pairs_hops)

        tau_c, p_c = _tau(pairs_claims)
        tau_h, p_h = _tau(pairs_hops)
        if tau_c == "" or tau_h == "":
            verdict = "underpowered"
            sign_agree = ""
        else:
            sign_agree = (tau_c >= 0) == (tau_h >= 0)
            verdict = "P5_pass" if sign_agree else "P5_fail"

        p5_rows.append(
            {
                "stratum": f"cov{cov}_bw{bw}",
                "coverage": cov,
                "base_wrongness": bw,
                "n_instances": len(recs),
                "n_tau_claims_defined": len(pairs_claims),
                "tau_claims": tau_c,
                "tau_claims_pvalue": p_c,
                "n_tau_hops_defined": len(pairs_hops),
                "tau_hops": tau_h,
                "tau_hops_pvalue": p_h,
                "sign_agree": sign_agree,
                "verdict": verdict,
            }
        )

    tau_c_all, p_c_all = _tau(all_pairs_claims)
    tau_h_all, p_h_all = _tau(all_pairs_hops)
    if tau_c_all == "" or tau_h_all == "":
        verdict_all = "underpowered"
        sign_agree_all = ""
    else:
        sign_agree_all = (tau_c_all >= 0) == (tau_h_all >= 0)
        verdict_all = "P5_pass" if sign_agree_all else "P5_fail"
    p5_rows.append(
        {
            "stratum": "ALL_POOLED",
            "coverage": "",
            "base_wrongness": "",
            "n_instances": len(instance_records),
            "n_tau_claims_defined": len(all_pairs_claims),
            "tau_claims": tau_c_all,
            "tau_claims_pvalue": p_c_all,
            "n_tau_hops_defined": len(all_pairs_hops),
            "tau_hops": tau_h_all,
            "tau_hops_pvalue": p_h_all,
            "sign_agree": sign_agree_all,
            "verdict": verdict_all,
        }
    )

    with ResultWriter(p5_path, resume=args.resume) as p5_w:
        for row in p5_rows:
            p5_w.write(_fill(row, P5_FIELDS))

    # --- D2 checks ---------------------------------------------------------
    s0_failures = [r["instance_id"] for r in instance_records if not r["s0_ok"]]
    g0_dhops_nonzero = [
        r["instance_id"] for r in instance_records if r["d_hops_of_g0"] not in (0, None)
    ]

    relationship = hp.dclaims_dhops_relationship(dclaims_dhops_pool)
    n_status_pool: dict[str, int] = {}
    for r in dclaims_dhops_pool:
        n_status_pool[r["status"]] = n_status_pool.get(r["status"], 0) + 1
    relationship_path = out_dir / "hops_dclaims_dhops.csv"
    with ResultWriter(relationship_path, resume=args.resume) as rel_w:
        for d in sorted(relationship):
            rel_w.write(
                _fill(
                    {"d_claims": d, **relationship[d]},
                    ["d_claims", "n", "mean_d_hops", "sd_d_hops"],
                )
            )

    grid = {
        "component_sizes": component_sizes,
        "coverages": coverages,
        "base_wrongness": base_wrongness,
        "n_reps": args.n_reps,
        "instances_per_cell": args.instances_per_cell,
        "seed_budget": args.seed_budget,
        "r_time_limit_s": args.r_time_limit_s,
    }
    manifest_path = write_hops_manifest(
        out_dir,
        seed=0,
        grid=grid,
        extra={
            "gate": hp.GATE_NAME,
            "elapsed_seconds": elapsed,
            "cell_counts": cell_counts,
            "exclusions": exclusions,
            "n_bug_rows": len(bug_rows),
            "bug_rows_sample": bug_rows[:50],
            "n_s0_failures": len(s0_failures),
            "s0_failures": s0_failures,
            "n_g0_dhops_nonzero_bug": len(g0_dhops_nonzero),
            "g0_dhops_nonzero_bug": g0_dhops_nonzero,
            "sample_status_counts": n_status_pool,
        },
    )

    summary = {
        "elapsed_seconds": elapsed,
        "cell_counts": cell_counts,
        "exclusions": exclusions,
        "n_instances": len(instance_records),
        "n_bug_rows": len(bug_rows),
        "n_s0_failures": len(s0_failures),
        "n_g0_dhops_nonzero_bug": len(g0_dhops_nonzero),
    }
    (out_dir / "_run_summary_hops.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))
    print(f"wrote {samples_path}")
    print(f"wrote {instances_path}")
    print(f"wrote {p5_path}")
    print(f"wrote {manifest_path}")
    print(f"elapsed_seconds={elapsed:.2f}")
    return 0


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    pc = sub.add_parser("census", help="Step 1: time build_space per component_size.")
    pc.add_argument("--sizes", default="4,5,6,7,8,9,10")
    pc.add_argument("--out-dir", default="results/axis_robustness")
    pc.add_argument("--seed-budget", type=int, default=200)
    pc.add_argument("--triangle-prob", type=float, default=0.3)
    pc.add_argument("--timeout-s", type=float, default=30.0)
    pc.add_argument("--resume", action="store_true")
    pc.set_defaults(func=cmd_census)

    ps = sub.add_parser("sweep", help="Steps 2-4: matched subsample, dual-axis sampling, P5.")
    ps.add_argument("--out-dir", default="results/axis_robustness")
    ps.add_argument("--component-sizes", default="4,5,6")
    ps.add_argument("--coverages", default="0.5,1.0")
    ps.add_argument("--base-wrongness", default="0.0,0.10,0.25")
    ps.add_argument("--instances-per-cell", type=int, default=15)
    ps.add_argument("--seed-budget", type=int, default=150)
    ps.add_argument("--n-reps", type=int, default=100)
    ps.add_argument("--r-time-limit-s", type=float, default=60.0)
    ps.add_argument("--resume", action="store_true")
    ps.set_defaults(func=cmd_sweep)

    return p


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
