"""Paired network-of-instances resampling on the N=1000 synthetic rerun.

Reviewer point 6 asks for differences between predictors to be evaluated with a
*paired* resampling design rather than inferred by eye from separate confidence
intervals. ``src/bkrobust/robustness/paired_resample.py`` implements that test
and was applied to the 831 real rows (``results/axis_robustness_real_p6/``).
This script applies the identical machinery to the synthetic rerun produced by
``run_survival_p6.py`` (``results/axis_robustness_p6/survival_instances.csv``),
so the two corpora are analysed the same way.

Clustering. The real corpus clusters on ``network`` because its 831 rows sit on
25 networks and within-network radius variation is thin. The synthetic corpus
has no networks: instances are independent draws from the component generator,
so each instance is its own cluster and the test reduces to an ordinary paired
bootstrap over instances. That is the correct analogue, not a weaker one.

Degeneracies are **measured per stratum** by :func:`degeneracy_note`, not
asserted, because they do not hold uniformly and mislabelling a real comparison
as an identity would understate it. What the data shows:

  * ``k_g0 == n_k`` on every row at coverage 1.0, where the analyst's claims are
    Meek-closed -- there it is a relabelling and its delta is byte-identical to
    ``n_k``'s. At coverage 0.5 the two differ and the comparison is real.
    ``k_g0`` is definitionally ``shd_cpdag`` in both cases.
  * ``r_val == separation`` exactly at coverage 1.0 with zero base wrongness,
    because the generator draws instances at target separations. A delta of
    exactly 0 there is an identity, not a null result. Elsewhere it is a genuine
    (and small) comparison.
  * ``shd_truth == 0`` identically at coverage 1.0 with zero base wrongness
    (``G0`` is then the ground truth), so it is predictor-constant there and the
    pair is reported ``undefined``.
  * ``r_claim`` is right-censored on many rows (240 of 480 in one stratum);
    censored values are dropped from the matched sample and counted in the note,
    never imputed. It coincides with ``r_val`` on the uncensored rows of some
    strata and not of others.

Usage::

    PYTHONPATH=src /usr/bin/python3 experiments/paired_synth_p6.py
"""

from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from bkrobust.robustness.paired_resample import paired_comparison  # noqa: E402

IN_CSV = ROOT / "results" / "axis_robustness_p6" / "survival_instances.csv"
OUT_DIR = ROOT / "results" / "axis_robustness_p6_paired"

#: Endpoints, matched to the arm, exactly as the committed analysis pairs them.
ENDPOINT = {"flip": "AUC_frac", "tiered": "AUC_frac"}

BASELINES = ["shd_truth", "n_k", "k_g0", "separation", "phi_1", "r_claim"]

N_BOOT = 10000
SEED = 0

COLUMNS = [
    "stratum", "arm", "predictor_a", "predictor_b", "endpoint", "status",
    "stratum_n", "n", "n_excluded", "n_clusters",
    "tau_a_point", "tau_b_point", "delta_point", "ci_lo_2p5", "ci_hi_97p5",
    "frac_delta_gt_0", "n_boot", "seed", "degeneracy_note",
]

def degeneracy_note(srows: list[dict[str, str]], baseline: str) -> str:
    """Report a degeneracy only where the stratum's own rows exhibit it.

    Identities are *measured* here rather than asserted, because they do not
    hold uniformly: ``k_g0 == n_k`` holds at coverage 1.0 (where the analyst's
    claims are Meek-closed) but not at coverage 0.5, and labelling the latter an
    identity would understate a real comparison.

    Args:
        srows: The stratum's instance rows.
        baseline: The baseline column being compared against ``r_val``.

    Returns:
        A note describing the degeneracy, or the empty string if none holds.
    """
    def col(name: str) -> list[str]:
        return [r[name] for r in srows if r[name] not in ("", "None")]

    if baseline == "k_g0":
        if col("k_g0") and col("k_g0") == col("n_k"):
            return "k_g0 == n_k on every row of this stratum; a relabelling, not an independent baseline."
        return "k_g0 is definitionally shd_cpdag, but differs from n_k here."
    if baseline == "separation":
        pairs = [(r["r_val"], r["separation"]) for r in srows
                 if r["r_val"] not in ("", "None") and r["separation"] not in ("", "None")]
        if pairs and all(a == b for a, b in pairs):
            return ("r_val == separation exactly on every row: the generator draws instances at "
                    "target separations, so a zero delta is an identity, not a null result.")
        return ""
    if baseline == "r_claim":
        pairs = [(r["r_val"], r["r_claim"]) for r in srows
                 if r["r_val"] not in ("", "None") and r["r_claim"] not in ("", "None")]
        n_cens = sum(1 for r in srows if r.get("r_claim_censored") == "1")
        note = f"{n_cens} of {len(srows)} rows right-censored and dropped, not imputed."
        if pairs and all(a == b for a, b in pairs):
            note = "r_val == r_claim on every uncensored row; a zero delta is an identity. " + note
        return note
    return ""


def stratum_of(row: dict[str, str]) -> str:
    """Name the pre-registered stratum a row belongs to."""
    if row["arm"] == "flip":
        return f"flip_cov{row['coverage']}_bw{row['base_wrongness']}"
    return f"tiered_nt{row['n_tiers']}"


def main() -> None:
    """Run every r_val-vs-baseline paired comparison, per stratum, and write them."""
    t0 = time.perf_counter()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    with IN_CSV.open() as fh:
        rows = list(csv.DictReader(fh))
    print(f"[paired_synth] loaded {len(rows)} instances from {IN_CSV}", flush=True)

    strata: dict[str, list[dict[str, str]]] = {}
    for r in rows:
        strata.setdefault(stratum_of(r), []).append(r)
    print(f"[paired_synth] {len(strata)} strata", flush=True)

    out_rows = []
    for stratum in sorted(strata):
        srows = strata[stratum]
        arm = srows[0]["arm"]
        endpoint = ENDPOINT[arm]
        for b in BASELINES:
            res = paired_comparison(
                srows, "r_val", b, endpoint, "instance_id",
                n_boot=N_BOOT, seed=SEED,
            )
            out_rows.append({
                "stratum": stratum,
                "arm": arm,
                "predictor_a": "r_val",
                "predictor_b": b,
                "endpoint": endpoint,
                "status": res.get("status", "ok"),
                "stratum_n": len(srows),
                "n": res.get("n"),
                "n_excluded": res.get("n_excluded"),
                "n_clusters": res.get("n_clusters"),
                "tau_a_point": res.get("tau_a_point"),
                "tau_b_point": res.get("tau_b_point"),
                "delta_point": res.get("delta_point"),
                "ci_lo_2p5": res.get("ci_lo_2p5"),
                "ci_hi_97p5": res.get("ci_hi_97p5"),
                "frac_delta_gt_0": res.get("frac_delta_gt_0"),
                "n_boot": N_BOOT,
                "seed": SEED,
                "degeneracy_note": degeneracy_note(srows, b),
            })
            print(
                f"[paired_synth] {stratum:24s} r_val vs {b:11s} "
                f"status={out_rows[-1]['status']} delta={out_rows[-1]['delta_point']}",
                flush=True,
            )

    out_csv = OUT_DIR / "paired_diff_synth.csv"
    with out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(out_rows)

    manifest = {
        "source": str(IN_CSV.relative_to(ROOT)),
        "n_instances": len(rows),
        "n_strata": len(strata),
        "cluster_col": "instance_id",
        "n_boot": N_BOOT,
        "seed": SEED,
        "baselines": BASELINES,
        "endpoint_by_arm": ENDPOINT,
        "elapsed_seconds": time.perf_counter() - t0,
        "degeneracies_note": "measured per stratum from the rows themselves; see degeneracy_note column",
    }
    (OUT_DIR / "manifest.json").write_text(json.dumps(manifest, indent=1))
    print(f"[paired_synth] wrote {out_csv} ({len(out_rows)} rows) "
          f"in {manifest['elapsed_seconds']:.1f}s", flush=True)


if __name__ == "__main__":
    main()
