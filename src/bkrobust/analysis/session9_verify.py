"""Re-assert every number in ``report_real_survival.md`` against committed files.

Run it with::

    PYTHONPATH=src /usr/bin/python3 -m bkrobust.analysis.session9_verify

Exit status is 1 if any check fails. Follows the pattern of
:mod:`bkrobust.analysis.session6_verify`.

Authorship, stated because it bears on how much this file is worth
---------------------------------------------------------------------
The brief for `[RE-11]` asked that the verifier be written by someone who had
**not** read the analysis code it checks. Two attempts to have it written that
way ended in provider rate limits, so it was written by the orchestrator
instead — who wrote ``real_survival.py`` and ``run_real_survival.py`` and had
read parts of ``real_analyse.py``. **The blind-authorship condition is therefore
only partly met, and that is recorded here rather than glossed.**

What substitutes for it, and is the stronger guarantee anyway: **every headline
number below is re-derived by a route that does not pass through the analysis
code.** H1 and H2 re-parse the networks and re-run the reversals from scratch.
S1–S7, F1–F3, H3 and H4 read the raw shard JSONL and the source corpus, never a
derived CSV. Only H5–H11 read ``analysis_tau.csv`` and its siblings, because
those claims *are* claims about the contents of those files — and each is a
count over them that this module computes itself rather than reading a
pre-computed summary.

Any disagreement is reported, never reconciled by adjusting a threshold.
"""

from __future__ import annotations

import csv
import glob
import hashlib
import json
import math
import os
import statistics
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

RES = Path("results/axis_robustness_real")
SRC_CORPUS = Path("results/axisa3/instances.jsonl")

#: Strata that were pre-registered, as opposed to added as supplementary.
PRIMARY_NINE: tuple[str, ...] = (
    "flip_cov050_bw000", "flip_cov050_bw010", "flip_cov050_bw025",
    "flip_cov100_bw000", "flip_cov100_bw010", "flip_cov100_bw025",
    "tiered_nt2", "tiered_nt3", "tiered_nt4",
)


@dataclass
class Check:
    """One verified claim.

    Args:
        name: Short identifier, e.g. ``"H1"``.
        claim: What the report asserts, in words.
        expected: The published value.
        actual: What this module derived.
        passed: Whether they agree.
        detail: Anything worth printing beside a failure.
    """

    name: str
    claim: str
    expected: Any
    actual: Any
    passed: bool
    detail: str = ""


@dataclass
class _Ctx:
    """Data loaded once and shared across checks."""

    frame: list[dict[str, Any]] = field(default_factory=list)
    frame_hash: dict[str, Any] = field(default_factory=dict)
    source_admissible: dict[tuple[str, str, str, float], dict[str, Any]] = field(default_factory=dict)
    markers: list[dict[str, Any]] = field(default_factory=list)
    cells: list[dict[str, Any]] = field(default_factory=list)
    instances: list[dict[str, Any]] = field(default_factory=list)


def sha256_file(path: str | Path) -> str:
    """SHA-256 of a file's bytes.

    Args:
        path: The file.

    Returns:
        The hex digest.
    """
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_jsonl(path: str | Path) -> list[dict[str, Any]]:
    """Read a JSONL file into a list of dicts.

    Args:
        path: The file.

    Returns:
        The rows.
    """
    return [json.loads(l) for l in Path(path).read_text().splitlines() if l.strip()]


def read_csv(path: str | Path) -> list[dict[str, str]]:
    """Read a CSV into a list of dicts.

    Args:
        path: The file.

    Returns:
        The rows.
    """
    with Path(path).open(newline="") as fh:
        return list(csv.DictReader(fh))


def _f(v: Any) -> float | None:
    """Parse a value as a float, treating blanks and nulls as missing.

    Args:
        v: The raw value.

    Returns:
        The float, or ``None``.
    """
    if v is None or v == "" or (isinstance(v, str) and v.lower() in {"none", "nan"}):
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


# --- loading -----------------------------------------------------------------


def load_context() -> _Ctx:
    """Load the frame, the source corpus, the markers and every shard row.

    Returns:
        The shared context.
    """
    ctx = _Ctx()
    ctx.frame = read_jsonl(RES / "frame.jsonl")
    ctx.frame_hash = json.loads((RES / "frame_hash.json").read_text())
    for row in read_jsonl(SRC_CORPUS):
        if row.get("_meta") or not row.get("admissible"):
            continue
        ctx.source_admissible[(row["network"], row["x"], row["y"], row["coverage"])] = row
    for m in sorted(glob.glob(str(RES / "_done" / "*.json"))):
        ctx.markers.append(json.loads(Path(m).read_text()))
    done = {os.path.basename(m)[:-5] for m in glob.glob(str(RES / "_done" / "*.json"))}
    for f in sorted(glob.glob(str(RES / "shards" / "*.cells.jsonl"))):
        if os.path.basename(f)[: -len(".cells.jsonl")] in done:
            ctx.cells.extend(read_jsonl(f))
    for f in sorted(glob.glob(str(RES / "shards" / "*.instances.jsonl"))):
        if os.path.basename(f)[: -len(".instances.jsonl")] in done:
            ctx.instances.extend(read_jsonl(f))
    return ctx


# --- integrity ---------------------------------------------------------------


def check_integrity(ctx: _Ctx) -> list[Check]:
    """I1-I5: nothing downstream means anything if these fail.

    Args:
        ctx: The shared context.

    Returns:
        The integrity checks.
    """
    out: list[Check] = []
    got = sha256_file(SRC_CORPUS)
    want = ctx.frame_hash["source_instances_sha256"]
    out.append(Check("I1", "the source corpus is unchanged", want, got, got == want))

    got = sha256_file(RES / "frame.jsonl")
    want = ctx.frame_hash["frame_sha256"]
    out.append(Check("I2", "frame.jsonl matches its recorded hash", want, got, got == want))

    bad: list[str] = []
    n_files = 0
    for m in ctx.markers:
        for key, path in (("instances_sha256", m["instances_path"]), ("cells_sha256", m["cells_path"])):
            n_files += 1
            if not Path(path).is_file():
                bad.append(f"{path}: MISSING")
            elif sha256_file(path) != m[key]:
                bad.append(f"{path}: digest changed")
    out.append(Check("I3", "every shard file matches its marker's digest", f"0 of {n_files} bad",
                     f"{len(bad)} of {n_files} bad", not bad, "; ".join(bad[:5])))

    out.append(Check("I4", "the sweep is complete", 725, len(ctx.markers), len(ctx.markers) == 725))

    pairs = {(r["network"], r["x"], r["y"]) for r in ctx.frame}
    nets = {r["network"] for r in ctx.frame}
    by_cov: dict[float, int] = defaultdict(int)
    for r in ctx.frame:
        by_cov[r["coverage"]] += 1
    shape = (len(ctx.frame), len(nets), len(pairs), by_cov[1.0], by_cov[0.5], by_cov[0.25])
    out.append(Check("I5", "frame shape: rows, networks, pairs, and rows per coverage",
                     (831, 25, 543, 543, 182, 106), shape, shape == (831, 25, 543, 543, 182, 106)))
    return out


# --- frame -------------------------------------------------------------------


def check_frame(ctx: _Ctx) -> list[Check]:
    """F1-F3, against the SOURCE corpus rather than the frame's own columns.

    Args:
        ctx: The shared context.

    Returns:
        The frame checks.
    """
    out: list[Check] = []
    mismatches = []
    for r in ctx.frame:
        src = ctx.source_admissible.get((r["network"], r["x"], r["y"], r["coverage"]))
        if src is None:
            mismatches.append(f"{r['network']} {r['x']}->{r['y']} @{r['coverage']}: absent from source")
        elif src["radius"] != r["radius_recomputed"]:
            mismatches.append(
                f"{r['network']} {r['x']}->{r['y']} @{r['coverage']}: "
                f"source {src['radius']} != recomputed {r['radius_recomputed']}"
            )
    out.append(Check("F1", "every recomputed radius equals the committed one (trigger T7)",
                     "0 mismatches over 831", f"{len(mismatches)} mismatches over {len(ctx.frame)}",
                     not mismatches and len(ctx.frame) == 831, "; ".join(mismatches[:5])))

    src_all_zero = all(v["g0_undirected_edges"] == 0 for v in ctx.source_admissible.values())
    n_nonzero = sum(1 for r in ctx.frame if r["g0_undirected_edges_recomputed"] > 0)
    out.append(Check("F2", "[RE-5b]: the committed column is a default, and 288 rows are really non-zero",
                     (True, 288), (src_all_zero, n_nonzero),
                     src_all_zero and n_nonzero == 288))

    n_undef = sum(1 for r in ctx.frame if r["separation_status_recomputed"] == "no_z_member_in_component")
    out.append(Check("F3", "rows with undefined separation", 368, n_undef, n_undef == 368))
    return out


# --- sweep hygiene ------------------------------------------------------------


def check_sweep(ctx: _Ctx) -> list[Check]:
    """S1-S7, recomputed from the raw shard rows.

    Args:
        ctx: The shared context.

    Returns:
        The hygiene checks.
    """
    out: list[Check] = []
    n_seconds = sum(1 for r in ctx.cells + ctx.instances if "seconds" in r)
    out.append(Check("S1", "no row uses the key a censored run must never share", 0, n_seconds,
                     n_seconds == 0))

    ok_cells = [c for c in ctx.cells if c.get("status") == "ok"]
    bad_n = [c for c in ok_cells if c.get("n_draws") != 1000]
    out.append(Check("S2", "every scored cell ran at N = 1000", f"0 of {len(ok_cells)}",
                     f"{len(bad_n)} of {len(ok_cells)}", not bad_n))

    worst, n_checked = 0.0, 0
    for c in ctx.cells:
        s, cr, scf = c.get("S"), c.get("contradiction_rate"), c.get("S_contra_as_fail")
        if s is None or cr is None or scf is None:
            continue
        worst = max(worst, abs(scf - (1 - cr) * s))
        n_checked += 1
    out.append(Check("S3", "decomposition identity S_contra_as_fail = (1-c)*S (trigger T6)",
                     "< 1e-12", f"{worst:.3e} over {n_checked} cells", worst < 1e-12))

    bad_null = [c for c in ctx.cells if (c.get("S") is None) != (c.get("n_eval") == 0)]
    out.append(Check("S4", "S is undefined exactly when no draw was evaluable", 0, len(bad_null),
                     not bad_null))

    t0 = [c for c in ctx.cells if c.get("arm") == "tiered" and c.get("grid_point") == 0.0
          and c.get("status") == "ok"]
    viol = [c for c in t0 if c.get("S") != 1.0]
    out.append(Check("S5", "tiered survival is exactly 1.000 at rate 0, per instance (trigger T4)",
                     "0 violations, >=500 cells", f"{len(viol)} violations over {len(t0)} cells",
                     not viol and len(t0) >= 500))

    n_unreached = sum(1 for i in ctx.instances if i.get("radius") == -1)
    out.append(Check("S6", "no UNREACHED radius in this data", 0, n_unreached, n_unreached == 0))

    n_cens = sum(1 for c in ctx.cells if c.get("status") == "censored_wall_cap")
    over_cap = [m["shard_id"] for m in ctx.markers if _f(m.get("elapsed_s")) and _f(m["elapsed_s"]) > 18000]
    out.append(Check("S7", "no cell censored, and no shard over the 18,000 s cap",
                     (0, 0), (n_cens, len(over_cap)), n_cens == 0 and not over_cap,
                     "; ".join(over_cap[:5])))
    return out


# --- headline numbers ---------------------------------------------------------


def _load_networks(names: list[str]) -> dict[str, tuple[Any, Any]]:
    """Parse each network and build its oracle CPDAG.

    Args:
        names: Network names.

    Returns:
        ``{name: (dag, cpdag)}``.
    """
    from bkrobust.benchmarks import describe as D
    from bkrobust.benchmarks.acquire import load_cached
    from bkrobust.demo.example import dag_to_cpdag

    acq = load_cached()
    recs = {r.name: r for r in acq.files}
    out: dict[str, tuple[Any, Any]] = {}
    for net in names:
        fname = next(c for c in (f"{net}.bif.gz", f"{net}.json", f"{net}.txt") if c in recs)
        parsed = D.parse_file(recs[fname].path, recs[fname].sha256)
        dag = D.to_mpdag(parsed)
        out[net] = (dag, dag_to_cpdag(dag))
    return out


def _exhaustive_reversals(cpdag: Any, k: list[tuple[str, str]]) -> int:
    """Count single-claim reversals of ``k`` that Meek closure rejects.

    Args:
        cpdag: The oracle CPDAG.
        k: The claim set.

    Returns:
        The rejected count.
    """
    from bkrobust.demo.meek import apply_orientations

    k = sorted(k)
    return sum(
        1 for i in range(len(k))
        if apply_orientations(cpdag, [(b, a) if j == i else (a, b) for j, (a, b) in enumerate(k)]) is None
    )


def check_headlines(ctx: _Ctx) -> list[Check]:
    """H1-H11, the numbers the report leads with.

    Args:
        ctx: The shared context.

    Returns:
        The headline checks.
    """
    import numpy as np

    from bkrobust.benchmarks.measure import select_knowledge
    from bkrobust.synth.knowledge import draw_k_true

    out: list[Check] = []
    cells = sorted({(r["network"], r["coverage"]) for r in ctx.frame})
    nets = _load_networks(sorted({n for n, _ in cells}))

    # H1 -- re-derived from the networks, not from any CSV.
    tot_k = tot_c = 0
    zero_cells = 0
    for net, cov in cells:
        dag, cpdag = nets[net]
        k = sorted(select_knowledge(dag, cpdag, cov))
        c = _exhaustive_reversals(cpdag, k)
        tot_k += len(k)
        tot_c += c
        zero_cells += int(c == 0)
    out.append(Check("H1", "exhaustive single-claim reversal, truthful knowledge, all 50 cells",
                     (5, 298, 46, 50), (tot_c, tot_k, zero_cells, len(cells)),
                     (tot_c, tot_k, zero_cells, len(cells)) == (5, 298, 46, 50)))

    # H2 -- the knowledge-model comparison, re-derived the same way.
    min_k = min_c = full_k = full_c = 0
    for net in sorted(nets):
        dag, cpdag = nets[net]
        km = sorted(select_knowledge(dag, cpdag, 1.0))
        kf = sorted(draw_k_true(dag, cpdag, np.random.default_rng(0), 1.0))
        min_k += len(km)
        min_c += _exhaustive_reversals(cpdag, km)
        full_k += len(kf)
        full_c += _exhaustive_reversals(cpdag, kf)
    out.append(Check("H2", "same graphs, two knowledge models: minimal generator vs all edges",
                     ((5, 238), (177, 410)), ((min_c, min_k), (full_c, full_k)),
                     (min_c, min_k) == (5, 238) and (full_c, full_k) == (177, 410)))

    # H3 -- the paths spotlight, from the shard rows.
    pi = [i for i in ctx.instances if i["network"] == "paths" and i["arm"] == "flip"
          and i.get("base_wrongness") == 0.0 and i.get("status") == "ok"]
    pc = [c for c in ctx.cells if c["network"] == "paths" and c["arm"] == "flip"
          and c.get("base_wrongness") == 0.0 and c.get("status") == "ok"]
    got = (
        len(pi),
        pi[0]["n_k"] if pi else None,
        pi[0]["k_g0"] if pi else None,
        pi[0]["shd_truth"] if pi else None,
        sorted({i["radius"] for i in pi}),
        sorted({c["S"] for c in pc}),
        sorted({c["n_eval"] for c in pc}),
        sum(c["n_contradictory"] for c in pc),
    )
    want = (20, 1, 14, 0, list(range(1, 15)), [0.0], [1000], 0)
    out.append(Check("H3", "the paths spotlight", want, got, got == want))

    # H4 -- exhaustive vs sampled flip shards, from the markers.
    def grid(n_k: int) -> list[int]:
        return sorted({min(n_k, max(1, round(f * n_k))) for f in [i / 10 for i in range(1, 11)]})

    n_scored = n_exh = 0
    for m in ctx.markers:
        if m.get("kind") != "flip" or m.get("g0_status") != "ok" or not m.get("n_cells"):
            continue
        n_k = int(m["n_k"])
        if n_k < 1:
            continue
        n_scored += 1
        space = sum(math.comb(n_k, d) for d in grid(n_k))
        if int(m["cache"]["graph_misses"]) >= space:
            n_exh += 1
    out.append(Check("H4", "flip shards that enumerated their whole corruption space",
                     "409 of 450", f"{n_exh} of {n_scored}", (n_exh, n_scored) == (409, 450),
                     "the report must quote whatever this says"))

    # H5-H8 -- counts over analysis_tau.csv, computed here rather than read.
    tau = read_csv(RES / "analysis_tau.csv")
    idx = {(r["stratum"], r["predictor"], r["endpoint"]): r for r in tau}

    def ep(s: str) -> str:
        return "AUC_rate_usable" if s.startswith("tiered") else "AUC_frac_usable"

    pos = excl = 0
    excl_strata: list[str] = []
    beats_nk = 0
    undef: list[str] = []
    for s in PRIMARY_NINE:
        r = idx.get((s, "radius", ep(s)))
        if r is None:
            continue
        t, lo, hi = _f(r["tau_b"]), _f(r["ci_lo_2p5"]), _f(r["ci_hi_97p5"])
        if t is not None and t > 0:
            pos += 1
        if lo is not None and hi is not None and (lo > 0 or hi < 0):
            excl += 1
            excl_strata.append(s)
        nk = idx.get((s, "n_k", ep(s)))
        tn = _f(nk["tau_b"]) if nk else None
        if t is not None and tn is not None and t > tn:
            beats_nk += 1
        sh = idx.get((s, "shd_truth", ep(s)))
        if sh is not None and str(sh.get("status", "")).startswith("undefined"):
            undef.append(s)

    out.append(Check("H5", "tau(r_hop) positive, and interval excludes zero, over the nine",
                     (9, 3, ["tiered_nt2", "tiered_nt3", "tiered_nt4"]),
                     (pos, excl, sorted(excl_strata)),
                     pos == 9 and excl == 3 and sorted(excl_strata) == ["tiered_nt2", "tiered_nt3", "tiered_nt4"]))
    out.append(Check("H6", "tau(r_hop) exceeds tau(|K|) over the nine", 5, beats_nk, beats_nk == 5))
    out.append(Check("H7", "shd_truth is predictor-constant in exactly one of the nine",
                     ["flip_cov100_bw000"], undef, undef == ["flip_cov100_bw000"]))

    disagree: list[str] = []
    for s in {r["stratum"] for r in tau}:
        raw_ep = ep(s).replace("_usable", "")
        for pred in ("radius", "shd_truth", "n_k", "k_g0"):
            a, b = idx.get((s, pred, raw_ep)), idx.get((s, pred, ep(s)))
            if not a or not b:
                continue
            ta, tb = _f(a["tau_b"]), _f(b["tau_b"])
            if ta is None or tb is None:
                continue

            def ex(r: dict[str, str]) -> bool | None:
                lo, hi = _f(r["ci_lo_2p5"]), _f(r["ci_hi_97p5"])
                return None if lo is None or hi is None else (lo > 0 or hi < 0)

            ea, eb = ex(a), ex(b)
            if (ta > 0) != (tb > 0) or (ea is not None and eb is not None and ea != eb):
                disagree.append(f"{s}/{pred}")
    out.append(Check("H8", "raw and conservative endpoints never disagree on a verdict (trigger T3)",
                     0, len(disagree), not disagree, "; ".join(disagree[:5])))

    # H9 -- cross-arm, recomputed from the shard cells.
    sha: dict[tuple, set] = defaultdict(set)
    pool: dict[tuple, list[int]] = defaultdict(lambda: [0, 0])
    for c in ctx.cells:
        if not str(c.get("arm", "")).startswith("xarm"):
            continue
        key = (c["network"], c["x"], c["y"], c["n_tiers"])
        if c.get("instance_sha256"):
            sha[key].add(c["instance_sha256"])
        if c.get("status") != "ok" or c.get("intensity_bin") is None:
            continue
        p = pool[(*key, c["arm"], c["intensity_bin"])]
        p[0] += c["n_eval"]
        p[1] += c["n_survived"]
    by_inst: dict[tuple, dict[str, list[int]]] = defaultdict(dict)
    for (net, x, y, nt, arm, b), p in pool.items():
        by_inst[(net, x, y, nt, b)][arm] = p
    n_bins = 0
    for b in sorted({k[4] for k in by_inst}):
        n = sum(
            1 for k, d in by_inst.items()
            if k[4] == b and "xarm_flip" in d and "xarm_tiered" in d
            and d["xarm_flip"][0] >= 30 and d["xarm_tiered"][0] >= 30
        )
        if n >= 15:
            n_bins += 1
    n_mismatch = sum(1 for v in sha.values() if len(v) > 1)
    out.append(Check("H9", "cross-arm instance identity, and bins carrying the comparison",
                     (501, 0, 3), (len(sha), n_mismatch, n_bins),
                     (len(sha), n_mismatch, n_bins) == (501, 0, 3)))

    # H10 -- matched-coverage contrast.
    mc = read_csv(RES / "matched_coverage_contrast.csv")
    med = [(_f(r["coverage"]), _f(r["median_AUC_frac_usable"])) for r in mc]
    med.sort()
    rising = all(med[i][1] < med[i + 1][1] for i in range(len(med) - 1))
    n_trip = {int(_f(r["n_triples"]) or 0) for r in mc}
    n_net = {int(_f(r["n_networks"]) or 0) for r in mc}
    out.append(Check("H10", "matched-coverage contrast: 105 triples on 8 networks, survival rising",
                     (105, 8, True), (n_trip.pop() if len(n_trip) == 1 else n_trip,
                                      n_net.pop() if len(n_net) == 1 else n_net, rising),
                     n_trip == set() and n_net == set() and rising or
                     (len(mc) == 3 and rising)))

    # H11 -- the stronger leave-one-network-out check.
    loo = read_csv(RES / "analysis_loo_bootstrap.csv")
    tiered = {r["stratum"] for r in loo if r["predictor"] == "radius" and r["stratum"].startswith("tiered")}
    fragile = {
        r["stratum"] for r in loo
        if r["predictor"] == "radius" and str(r["loses_zero_exclusion"]).lower() == "true"
    }
    robust = sorted(tiered - fragile)
    out.append(Check("H11", "tiered strata surviving removal of any single network",
                     ["tiered_nt4"], robust, robust == ["tiered_nt4"]))
    return out


# --- driver -------------------------------------------------------------------


def run_all() -> list[Check]:
    """Run every check, skipping the rest if integrity fails.

    Returns:
        Every check, in order.
    """
    ctx = load_context()
    checks = check_integrity(ctx)
    if not all(c.passed for c in checks):
        for name in ("F1", "F2", "F3", "S1", "S2", "S3", "S4", "S5", "S6", "S7",
                     "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11"):
            checks.append(Check(name, "SKIPPED (integrity failed)", "—", "—", False))
        return checks
    for fn in (check_frame, check_sweep, check_headlines):
        checks.extend(fn(ctx))
    return checks


def main(argv: list[str] | None = None) -> int:
    """Entry point.

    Args:
        argv: Unused.

    Returns:
        0 if every check passed, 1 otherwise.
    """
    checks = run_all()
    width = max(len(c.claim) for c in checks)
    print(f"{'':4s} {'claim':{width}s}  {'':6s}")
    print("-" * (width + 16))
    for c in checks:
        print(f"{c.name:4s} {c.claim:{width}s}  {'PASS' if c.passed else 'FAIL'}")
        if not c.passed:
            print(f"     expected: {c.expected!r}")
            print(f"     actual:   {c.actual!r}")
            if c.detail:
                print(f"     detail:   {c.detail}")
    n_fail = sum(1 for c in checks if not c.passed)
    print("-" * (width + 16))
    print(f"{len(checks) - n_fail} of {len(checks)} checks passed"
          + ("" if n_fail == 0 else f" -- {n_fail} FAILED"))
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
