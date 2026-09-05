"""Assert every quoted number in ``report_axisb_sat.md`` against the committed data.

Run from the repository root::

    python -m bkrobust.analysis.session3_verify

Exits non-zero on the first disagreement, and prints all of them.

Session 2's lesson: a rebuild is not a verification. Each anchor is asserted to
be PRESENT in the rendered text as well as correct, so a silently no-op edit or
a stale citation fails loudly instead of passing.
"""

# ruff: noqa: RUF001
# The report quotes multiplication signs; the anchors must match it byte for byte.

import json
import pathlib
import re
import statistics as st
import sys

TXT = pathlib.Path("report_axisb_sat.md").read_text()


def _jsonl(name: str) -> list[dict]:
    path = pathlib.Path("results/axisb3") / name
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


S = _jsonl("scaling.jsonl")
H = _jsonl("highk.jsonl")
X = json.loads(pathlib.Path("results/axisb3/cross_encoding_n4.json").read_text())
fails = []


def check(label: str, computed: object, anchor: str) -> None:
    """Assert the anchor string appears in the report AND matches the data."""
    if anchor not in TXT:
        fails.append(f"MISSING from report: {label!r} -> {anchor!r}")
    elif computed is not None and str(computed).replace(",", "") not in anchor.replace(",", ""):
        fails.append(f"MISMATCH {label!r}: data says {computed}, report says {anchor!r}")


# --- cross-encoding at n = 4 ------------------------------------------------
check("n=4 instances", X["n_instances"], "1,044")
for k in ("e1_bfs", "e1_up", "e2_bfs", "e1_e2", "e3_bfs"):
    assert X["agreements"][k] == 1044, (k, X["agreements"][k])
check("n=4 walks", X["walks_certified"], "792")
assert X["walks_certified"] == X["walks_monotone"] == 792
assert X["c2_counterexamples"] == []
assert X["e2_optima_with_down_leg"] == 0

# --- scaling ---------------------------------------------------------------
check("scaling rows", len(S), "360 instances")
bfs = [r for r in S if "bfs" in r]
check("bfs scope", len(bfs), "346 of 360")
assert all(r["bfs"]["radius"] == r["e1"]["radius"] for r in bfs)
assert all(r["local_up"]["radius"] == r["e1"]["radius"] for r in S if "local_up" in r)
e3s = [r for r in S if "e3" in r]
check("scaling e3", len(e3s), "129")
assert all(r["e3"]["certified"] and r["e3"]["monotone"] for r in e3s)
assert all(r["e3"]["radius"] == r["e1"]["radius"] for r in e3s)
deg = sum(1 for r in S if r["e1"]["radius"] == -1)
check("degeneracy", f"{100 * deg / len(S):.1f}", "64.2%")
check("bfs not built", len(S) - len(bfs), "14 of 360 rows")

for n, build, unsat, sat, med, worst in (
    (8, "0.004", "0.001", "0.000", "0.006", "0.015"),
    (10, "0.007", "0.002", "0.002", "0.012", "0.074"),
    (12, "0.014", "0.007", "0.000", "0.025", "0.134"),
    (15, "0.031", "0.019", "0.000", "0.054", "0.374"),
    (20, "0.087", "0.053", "0.000", "0.136", "1.497"),
):
    g = [r for r in S if r["n"] == n]
    for key, want in (("build_s", build), ("unsat_s", unsat), ("sat_s", sat)):
        got = f"{st.median(r['e1'][key] for r in g):.3f}"
        if got != want:
            fails.append(f"n={n} {key}: data {got} != report {want}")
    if f"{st.median(r['e1']['total_s'] for r in g):.3f}" != med:
        fails.append(f"n={n} median")
    if f"{max(r['e1']['total_s'] for r in g):.3f}" != worst:
        fails.append(f"n={n} worst")
    check(f"n={n} worst row", None, f"| {worst} s |")

# --- high-k hunt -----------------------------------------------------------
check("hunt rows", len(H), "144 instances")
check("max k", max(r["k_g0"] for r in H), "|K_{G₀}| = 51")
check("max component", max(r["largest_component"] for r in H), "**11 vertices**")
ok = [r for r in H if "timed_out_at_s" not in r["local_up"]]
check("hunt agreement", len(ok), "**117 of 117**")
assert all(r["local_up"]["radius"] == r["e1"]["radius"] for r in ok)
near = [r for r in H if r["e1"]["radius"] not in (-1, -2)]
far = [r for r in H if r["e1"]["radius"] == -1]
check("SAT-side rows", len(near), "| 95 |")
check("UNSAT-side rows", len(far), "| 49 |")
_near_ok = [r for r in near if "timed_out_at_s" not in r["local_up"]]
check(
    "SAT-side up median",
    f"{st.median(r['local_up']['total_s'] for r in _near_ok):.4f}",
    "0.0027 s",
)
check("SAT-side E1 median", f"{st.median(r['e1']['total_s'] for r in near):.3f}", "0.099 s")
check("SAT-side E1 worst", f"{max(r['e1']['total_s'] for r in near):.3f}", "4.117 s")
check("UNSAT-side E1 median", f"{st.median(r['e1']['total_s'] for r in far):.3f}", "1.865 s")
check("UNSAT-side E1 worst", f"{max(r['e1']['total_s'] for r in far):.3f}", "8.745 s")
farok = [r for r in far if "timed_out_at_s" not in r["local_up"]]
check("UNSAT up median", f"{st.median(r['local_up']['total_s'] for r in farok):.3f}", "10.459 s")
check("UNSAT completed", len(farok), "23 that finished")
check("UNSAT timeouts", sum(1 for r in far if "timed_out_at_s" in r["local_up"]), "**26 / 49**")
check("SAT timeouts", sum(1 for r in near if "timed_out_at_s" in r["local_up"]), "1 / 95")
ratios = sorted((r["local_up"]["total_s"] / r["e1"]["total_s"], r) for r in farok)
check("median speedup", f"{st.median(x for x, _ in ratios):.1f}", "8.8×")
check("max speedup", f"{ratios[-1][0]:.1f}", "117.9×")
check("max speedup k", ratios[-1][1]["k_g0"], "k = 12: E1 0.129 s against 15.212 s")
w = max((r for r in far if "timed_out_at_s" in r["local_up"]), key=lambda r: r["k_g0"])
check("worst timeout", f"n = {w['n']} and k = {w['k_g0']}", f"n = {w['n']} and k = {w['k_g0']}")
check("worst timeout E1", f"{w['e1']['total_s']:.3f}", "8.745 s")
check("hunt degeneracy", f"{100 * len(far) / len(H):.1f}", "34.0%")
he3 = [r for r in H if "e3" in r]
check("hunt e3", len(he3), "| 86 |")
for _r in he3:
    assert _r["e3"]["certified"] and _r["e3"]["monotone"]
    assert _r["e3"]["radius"] == _r["e1"]["radius"]
check("total walks", 792 + len(e3s) + len(he3), "1,007")
big = max(he3, key=lambda r: r["k_g0"])
check("largest e3", f"{big['n']}", "**n = 18, `|K_{G₀}| = 51`")
check("largest e3 vars", f"{big['e3']['n_variables']:,}", "2,125 variables")
check("largest e3 time", f"{big['e3']['total_s']:.3f}", "1.355 s")
check("largest e3 density", f"{big['largest_component_density']:.2f}", "density\n0.93")
radii = sorted({r["e1"]["radius"] for r in H})
assert 5 in radii and 6 in radii and 8 in radii, radii
check("uncovered radii", None, "radii 5, 6 and 8 exist")
check(
    "k=33 conflicts",
    f"{st.median(r['e1']['conflicts'] for r in H if r['k_g0'] == 33):,.0f}",
    "118,839",
)

# --- crossover vs BFS table ------------------------------------------------
for m, rows_n, size, b, e, ratio in (
    (3, 40, "18", "0.002", "0.017", "0.10"),
    (4, 62, "30", "0.009", "0.026", "0.32"),
    (5, 58, "60", "0.043", "0.040", "1.08"),
    (6, 30, "98", "0.147", "0.021", "6.88"),
    (7, 22, "294", "0.830", "0.064", "13.02"),
    (8, 22, "210", "1.696", "0.033", "51.00"),
    (9, 4, "515", "13.172", "0.245", "53.84"),
):
    g = [r for r in bfs if r["m_undirected"] == m]
    if len(g) != rows_n:
        fails.append(f"m={m} rows: {len(g)} != {rows_n}")
    got_b = st.median(r["bfs"]["build_s"] for r in g)
    got_e = st.median(r["e1"]["total_s"] for r in g)
    for got, want, lab in (
        (f"{got_b:.3f}", b, "bfs"),
        (f"{got_e:.3f}", e, "e1"),
        (f"{got_b / got_e:.2f}", ratio, "ratio"),
        (f"{st.median(r['bfs']['space_size'] for r in g):.0f}", size, "size"),
    ):
        if got != want:
            fails.append(f"m={m} {lab}: data {got} != report {want}")

# --- figures referenced exist ----------------------------------------------
for fig in re.findall(r"!\[[^\]]*\]\(([^)]+)\)", TXT):
    if not pathlib.Path(fig).exists():
        fails.append(f"figure missing on disk: {fig}")

print(f"anchors checked; {len(fails)} problem(s)")
for f in fails:
    print("  ", f)
sys.exit(1 if fails else 0)
