"""Assert every quoted number in ``report_real_graphs.md`` against the data.

Run from the repository root::

    python -m bkrobust.analysis.session6_verify

Each anchor is checked twice over: that the claimed string is actually PRESENT in
the rendered markdown, so a silently no-op edit fails loudly, and that it agrees
with the committed files. Sessions 3, 4 and 5 each caught a stale or
double-counted figure this way.
"""

from __future__ import annotations

import csv
import json
import pathlib
import re
import statistics as st
import sys

TXT = pathlib.Path("report_real_graphs.md").read_text()
RES = pathlib.Path("results/axisa3")
fails: list[str] = []


def check(label: str, computed: object, anchor: str) -> None:
    """Assert the anchor appears in the report AND agrees with the data."""
    if anchor not in TXT:
        fails.append(f"MISSING from report: {label!r} -> {anchor!r}")
    elif computed is not None and str(computed).replace(",", "") not in anchor.replace(",", ""):
        fails.append(f"MISMATCH {label!r}: data says {computed}, report says {anchor!r}")


S = json.loads((RES / "analysis_summary.json").read_text())
desc = list(csv.DictReader(open(RES / "descriptive_structure.csv")))
per = list(csv.DictReader(open(RES / "per_network.csv")))
manifest = json.loads((RES / "networks" / "acquisition_manifest.json").read_text())

# --- acquisition ------------------------------------------------------------
sha = json.dumps(manifest)
check("sdist sha", None, "aef361e0858bbb1de839c54b940b203170609e1822aff37fc6853e715478255a")
if "aef361e0858bbb1de839c54b940b203170609e1822aff37fc6853e715478255a" not in sha:
    fails.append("sdist sha256 in the report is not the one in the manifest")
check(
    "bif count",
    sum(1 for d in desc if d["source_format"] == "bif"),
    "| BIF benchmark networks | 24 |",
)
check(
    "dagitty count",
    sum(1 for d in desc if d["source_format"] == "dagitty"),
    "| dagitty DAGs from applied papers | 12 (11 usable) |",
)
check(
    "bnjson count",
    sum(1 for d in desc if d["source_format"] == "bnjson"),
    "| bnjson gene/plant networks | 4 |",
)

# --- descriptive ------------------------------------------------------------
dags = [d for d in desc if d["is_dag"] == "True"]
fr = [100 * float(d["undirected_fraction"]) for d in dags]
check("n dags", len(dags), "across 39 parsed DAGs")
check("median undirected", f"{st.median(fr):.1f}", "| median undirected fraction | **10.7%** |")
zero = sorted(d["network"] for d in dags if float(d["undirected_fraction"]) == 0)
check("fully compelled", len(zero), "**6 networks**")
for z in zero:
    if z not in TXT:
        fails.append(f"fully-compelled network {z} not named in the report")
mx = [int(d["max_component_size"]) for d in dags]
check("median max comp", f"{st.median(mx):.0f}", "| median max chain component | **4** |")
check("largest comp", max(mx), "| largest chain component | **85** (pathfinder) |")
check("gt6", sum(1 for v in mx if v > 6), "| networks with max component > 6 | **11 of 39** |")

# --- the measurement --------------------------------------------------------
sep, rad, law, tr = S["separation"], S["radius"], S["law"], S["tractability"]
check("pairs screened", f"{S['pairs_screened']:,}", "**115,974 ordered pairs screened**")
check("admissible", S["admissible_instances"], "**831 admissible instances**")
check("nets with adm", S["networks_with_admissible"], "**25 networks**")
check("censored", S["censored_network_runs"], "3 network-coverage runs")
check("intractable", tr["o_g0_extensions_intractable"], "**280**")
check("sep ge2 pair", sep["ge2_pair_weighted_pct"], "**45.4% of measurable pairs**")
check("sep ge2 net", sep["ge2_network_weighted_pct"], "**35.6%**")
check("sep max", sep["max"], "maximum of\n**7**")
check("sep undefined", sep["undefined"], "**368 of 831**")
check("sep undefined pct", sep["undefined_pct"], "(44.3%)")
check("r1 pair", rad["r1_pair_weighted_pct"], "**61.1%**")
check("r1 net", rad["r1_network_weighted_pct"], "**71.3%**")
check("r max", rad["max"], "maximum **14**")
if rad["unreached"] != 0:
    fails.append("report says no UNREACHED instances but the data has some")
check("no unreached", None, "**no UNREACHED instances at all**")
for cov, key in (("1.00", "1.0"), ("0.50", "0.5"), ("0.25", "0.25")):
    c = S["by_coverage"][key]
    check(f"cov {cov}", f"{c['r1_pct']}", f"| {cov} | {c['admissible']} | {c['r1_pct']}% |")

# --- the law ----------------------------------------------------------------
check("law agree", f"{law['agree']} of {law['tested']}", f"**{law['agree']} of {law['tested']}")
check("law pct", law["agree_pct"], "(72.8%)**")
check("law exceptions", law["exceptions"], "The 126 exceptions")
check(
    "law greater",
    law["exception_direction"]["radius_greater_than_law"],
    "110 have a radius *larger*",
)
check("law less", law["exception_direction"]["radius_less_than_law"], "only\n16 smaller")
check(
    "shape 2/1", law["exception_shapes"].get("r=2,min=1"), "`r = 2` where the law says 1 (52 cases)"
)
check("shape 3/2", law["exception_shapes"].get("r=3,min=2"), "`r = 3` where it says 2 (40)")

# --- tractability -----------------------------------------------------------
check("max comp radius", tr["max_component_with_completed_radius"], "**85**")
check("slowest", tr["slowest_seconds"], "**0.16 s**")
check("inexact", tr["inexact_instances"], "**0 inexact instances**")

# --- per-network table rows -------------------------------------------------
for p in per:
    if int(p["admissible"]) == 0:
        continue
    row = (
        f"| {p['network']} | {p['admissible']} | {p['separation_measurable']} | "
        f"{p['sep_ge2_pct'] if p['sep_ge2_pct'] != '' else '—'}%"
    )
    if row not in TXT:
        fails.append(f"per-network row missing or altered: {p['network']}")

for fig in re.findall(r"!\[[^\]]*\]\(([^)]+)\)", TXT):
    if not pathlib.Path(fig).exists():
        fails.append(f"figure missing on disk: {fig}")

print(f"anchors checked; {len(fails)} problem(s)")
for f in fails:
    print("  ", f)
sys.exit(1 if fails else 0)
