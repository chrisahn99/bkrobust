"""Assert every quoted number in ``report_axisb_oracle.md`` against the data.

Run from the repository root::

    python -m bkrobust.analysis.session4_verify

Session 2's lesson was that a rebuild is not a verification. Each anchor is
checked twice over: that the claimed string is actually PRESENT in the rendered
markdown, so a silently no-op edit fails loudly, and that it matches what the
committed files say. Exits non-zero on any disagreement, printing all of them.
"""

# ruff: noqa: RUF001
# The report quotes multiplication signs; the anchors must match it byte for byte.

from __future__ import annotations

import json
import pathlib
import re
import statistics as st
import sys

TXT = pathlib.Path("report_axisb_oracle.md").read_text()
RES = pathlib.Path("results/axisb4")


def _jsonl(name: str) -> list[dict]:
    return [json.loads(line) for line in (RES / name).read_text().splitlines() if line.strip()]


def _json(name: str) -> dict:
    return json.loads((RES / name).read_text())


fails: list[str] = []


def check(label: str, computed: object, anchor: str) -> None:
    """Assert the anchor is present in the report AND agrees with the data."""
    if anchor not in TXT:
        fails.append(f"MISSING from report: {label!r} -> {anchor!r}")
    elif computed is not None and str(computed).replace(",", "") not in anchor.replace(",", ""):
        fails.append(f"MISMATCH {label!r}: data says {computed}, report says {anchor!r}")


# --- the profile ------------------------------------------------------------
prof = [r for r in _jsonl("profile_local_up.jsonl") if "total_s" in r and "timed_out_at_s" not in r]
sat = [r for r in prof if r["side"] == "SAT"]
unsat = [r for r in prof if r["side"] == "UNSAT"]


def share(rs: list[dict], key: str) -> float:
    """Percentage of total wall-clock spent in one component."""
    tot = sum(r["total_s"] for r in rs)
    return 100 * sum(r[key] for r in rs) / tot if tot else 0.0


check("profiled instances", len(prof), "115 instances profiled")
check("oracle share all", f"{share(prof, 'oracle_s'):.1f}", "| 25.6% |")
check("cover share all", f"{share(prof, 'cover_ext_s'):.1f}", "| **70.4%** |")
check("oracle share unsat", f"{share(unsat, 'oracle_s'):.1f}", "25.2%")
check("cover share unsat", f"{share(unsat, 'cover_ext_s'):.1f}", "**71.3%**")
check("enum share all", f"{share(prof, 'oracle_s') + share(prof, 'cover_ext_s'):.1f}", "**96.1%**")
check(
    "enum share unsat", f"{share(unsat, 'oracle_s') + share(unsat, 'cover_ext_s'):.1f}", "**96.5%**"
)
o = share(unsat, "oracle_s") / 100
c = share(unsat, "cover_ext_s") / 100
check("amdahl criterion", f"{1 / (1 - o):.1f}", "**1.3×**")
check("amdahl lemma O", f"{1 / (1 - c):.1f}", "**3.5×**")
check("amdahl both", f"{1 / (1 - o - c):.1f}", "**28.7×**")
check("sat total", f"{sum(r['total_s'] for r in sat):.2f}", "0.69 s")
check("unsat total", f"{sum(r['total_s'] for r in unsat):.2f}", "22.74 s")

# --- Lemma O ----------------------------------------------------------------
order = _json("order_identification.json")
cov = _json("fast_covers_differential.json")
check("order pairs", order["pairs"], "80,480 ordered pairs")
assert order["disagree"] == 0
check("cover sets", cov["cover_sets_compared"], "1,588 cover sets")
check("radii compared", cov["radii_compared"], "1,044 radii")
assert cov["cover_sets_identical"] == cov["cover_sets_compared"]
assert cov["radii_agree_with_bfs"] == cov["radii_compared"]

speed = _jsonl("fast_vs_frozen.jsonl")
both = [
    r
    for r in speed
    if "total_s" in r["slow"]
    and "total_s" in r["fast"]
    and "timed_out_at_s" not in r["slow"]
    and "timed_out_at_s" not in r["fast"]
]
resc = [
    r
    for r in speed
    if "timed_out_at_s" in r["slow"]
    and "total_s" in r["fast"]
    and "timed_out_at_s" not in r["fast"]
]
s_side = [r for r in both if r["fast"]["side"] == "SAT"]
u_side = [r for r in both if r["fast"]["side"] == "UNSAT"]


def agg(rs: list[dict]) -> float:
    """Time-weighted speedup: total slow time over total fast time."""
    return sum(r["slow"]["total_s"] for r in rs) / max(sum(r["fast"]["total_s"] for r in rs), 1e-9)


def med(rs: list[dict]) -> float:
    """Median per-instance speedup."""
    return st.median(r["slow"]["total_s"] / max(r["fast"]["total_s"], 1e-9) for r in rs)


check("speed pairs", len(speed), "144 instance pairs")
check("speed both", len(both), "| all | 136 |")
check("speed agg all", f"{agg(both):.2f}", "**3.77×**")
check("speed med all", f"{med(both):.2f}", "1.51×")
check("speed sat", len(s_side), "| SAT side | 67 |")
check("speed unsat", len(u_side), "| UNSAT side | 69 |")
check("speed agg unsat", f"{agg(u_side):.2f}", "**3.79×**")
check("rescued", len(resc), "**0 radius disagreements.** 2 instances")
check("frozen timeouts", sum(1 for r in speed if "timed_out_at_s" in r["slow"]), "| 8 |")
check("fast timeouts", sum(1 for r in speed if "timed_out_at_s" in r["fast"]), "| 6 |")

# --- criterion --------------------------------------------------------------
crit = _json("criterion_agreement.json")
s2 = _json("stage2_rerun.json")
strata = s2["strata"]
stress = _json("criterion_stress_walks.json")
check("crit cases", crit["n_cases"], "74,568")
assert crit["n_disagree"] == 0
check("crit nonam", crit["non_amenable_cases"], "28,464")
check("stage2 cases", s2["n_cases"], "| E1 vs enumeration oracle | 74,568 |")
assert s2["e1_vs_enumeration_agree"] == s2["n_cases"]
assert s2["e1_vs_criterion_agree"] == s2["n_cases"]
check("stage2 nonam", strata["NON-amenable/total"], "| **28,464** |")
assert strata["NON-amenable/agree"] == strata["NON-amenable/total"]
check("stress cases", f"{stress['n_cases']:,}", "**1,131,460 cases, 0 disagreements**")
assert stress["n_disagree"] == 0

pathenum = _jsonl("oracle_envelope_pathenum.jsonl")
pe_to = [r for r in pathenum if "timed_out_at_s" in r["criterion"]]
pe_ok = [r for r in pathenum if "largest_component" in r["criterion"]]
check(
    "pathenum worst",
    f"{max(r['criterion']['seconds'] for r in pe_ok):.2f}",
    "**2.35 seconds per call**",
)
check(
    "pathenum timeouts",
    f"{len(pe_to)} of {len(pathenum)} instances",
    "**60 s cap on 12 of 81 instances**",
)

# --- hybrid envelope --------------------------------------------------------
hd = _json("hybrid_differential.json")
check("hybrid diff", hd["n_instances"], "1,044 instances at n = 4")
assert hd["hybrid_criterion_agrees_with_bfs"] == hd["n_instances"]
assert hd["hybrid_enumeration_agrees_with_bfs"] == hd["n_instances"]

env = _jsonl("hybrid_envelope.jsonl")


def ok(e: dict) -> bool:
    """Whether an entry is a real measurement rather than a censored run."""
    return "total_s" in e and "timed_out_at_s" not in e


hyb = [r["hybrid"] for r in env if ok(r["hybrid"])]
froz_to = [r for r in env if "timed_out_at_s" in r["frozen"]]
pair = [r for r in env if ok(r["hybrid"]) and ok(r["frozen"])]
check("env rows", len(env), "**256 / 256**")
check("env total", len(env), "256 instances, n = 8…24")
assert len(hyb) == len(env), (len(hyb), len(env))
check("frozen env completed", len(env) - len(froz_to), "| 247 / 256 |")
check("frozen env timeouts", len(froz_to), "| **9** |")
check(
    "env agg",
    f"{agg([{'slow': r['frozen'], 'fast': r['hybrid']} for r in pair]):.2f}",
    "aggregate **3.29×**",
)
check(
    "env med",
    f"{med([{'slow': r['frozen'], 'fast': r['hybrid']} for r in pair]):.2f}",
    "median **2.03×**",
)
check(
    "max component",
    max(e["largest_component"] for e in hyb),
    "**largest undirected component solved: 11 vertices**",
)
check("max k", max(e["k_g0"] for e in hyb), "**largest \\|K_{G₀}\\| solved: 33 orientations**")
check("max n", max(e["n"] for e in hyb), "**largest n solved: 24**")
check("slowest", f"{max(e['total_s'] for e in hyb):.2f}", "**slowest single instance: 21.62 s**")
disp = {}
for e in hyb:
    disp[e["method"]] = disp.get(e["method"], 0) + 1
check("dispatch search", disp.get("local_up_fast", 0), "148 instances answered by the search")
check("dispatch ladder", disp.get("e1_ladder", 0), "108 by the")
for r in froz_to:
    assert r["hybrid"]["side"] == "UNSAT", r
    assert r["hybrid"]["method"] == "e1_ladder", r
check("all rescued are unsat ladder", None, "all nine are **UNSAT**")

# per-component envelope table
by: dict[int, list[dict]] = {}
for e in hyb:
    by.setdefault(e["largest_component"], []).append(e)
for comp, rows_n, kmed, kmax, m, w in (
    (2, 49, 1, 3, "0.0005", "0.33"),
    (3, 58, 3, 5, "0.0011", "6.29"),
    (4, 52, 5, 10, "0.0559", "9.55"),
    (5, 50, 6, 13, "0.0618", "12.32"),
    (6, 33, 11, 14, "0.0289", "21.62"),
    (7, 7, 17, 21, "0.0713", "9.86"),
    (8, 4, 20, 24, "1.1023", "4.47"),
    (9, 2, 32, 33, "0.0046", "0.0055"),
    (11, 1, 14, 14, "0.3885", "0.3885"),
):
    g = by.get(comp, [])
    if len(g) != rows_n:
        fails.append(f"component {comp}: {len(g)} rows, report says {rows_n}")
        continue
    got_k = st.median(e["k_g0"] for e in g)
    if f"{got_k:.0f}" != str(kmed) or max(e["k_g0"] for e in g) != kmax:
        fails.append(f"component {comp}: k {got_k}/{max(e['k_g0'] for e in g)} vs {kmed}/{kmax}")
    if f"{st.median(e['total_s'] for e in g):.4f}" != m:
        fails.append(f"component {comp}: median {st.median(e['total_s'] for e in g)} vs {m}")
    got_w = max(e["total_s"] for e in g)
    if f"{got_w:.2f}" != w and f"{got_w:.4f}" != w:
        fails.append(f"component {comp}: worst {got_w} vs {w}")

# --- figures referenced exist ----------------------------------------------
for fig in re.findall(r"!\[[^\]]*\]\(([^)]+)\)", TXT):
    if not pathlib.Path(fig).exists():
        fails.append(f"figure missing on disk: {fig}")

print(f"anchors checked; {len(fails)} problem(s)")
for f in fails:
    print("  ", f)
sys.exit(1 if fails else 0)
