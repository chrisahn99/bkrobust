"""Assert every quoted number in ``report_saturation_gac.md`` against the data.

Run from the repository root::

    python -m bkrobust.analysis.session5_verify

Each anchor is checked twice over: that the claimed string is actually PRESENT
in the rendered markdown, so a silently no-op edit fails loudly, and that it
agrees with the committed files. Exits non-zero on any disagreement, printing
all of them. Sessions 3 and 4 both caught a stale figure this way.
"""

from __future__ import annotations

import json
import pathlib
import re
import statistics as st
import sys

TXT = pathlib.Path("report_saturation_gac.md").read_text()
RES = pathlib.Path("results/axisa2")
UNREACHED = -1
fails: list[str] = []


def _jsonl(name: str) -> list[dict]:
    return [json.loads(line) for line in (RES / name).read_text().splitlines() if line.strip()]


def _json(name: str) -> dict:
    return json.loads((RES / name).read_text())


def check(label: str, computed: object, anchor: str) -> None:
    """Assert the anchor appears in the report AND agrees with the data."""
    if anchor not in TXT:
        fails.append(f"MISSING from report: {label!r} -> {anchor!r}")
    elif computed is not None and str(computed).replace(",", "") not in anchor.replace(",", ""):
        fails.append(f"MISMATCH {label!r}: data says {computed}, report says {anchor!r}")


# --- Phase 1: the GAC -------------------------------------------------------
gac = _json("gac_agreement.json")
c1, c2, c3, c4 = (
    gac["check1_dag_gac_vs_backdoor"],
    gac["check2_mpdag_vs_enumeration"],
    gac["check3_optimal_set"],
    gac["check4_performance"],
)
check("dag graphs", c1["n_graphs"], "all 29,824 labelled DAGs")
check("dag cases", c1["n_cases"], "4,711,024")
check("dag disagree", c1["n_disagree"], "227,224")
assert c1["n_unexplainable_disagreements"] == 0
assert c1["implication_failures"] == 0
check("mpdag cases", c2["n_cases"], "74,568")
assert c2["n_disagree"] == 0
check("mpdag nonamenable", c2["non_amenable_cases"], "28,464")
check("perf", f"{c4['worst_ms_at_n20']:.2f}", "1.07 ms per call at n = 20")
assert c4["bar_met"]
check("oset cases", c3["n_cases"], "8,154")
assert c3["n_disagree"] == 0
check("oset bothvalid", c3["n_both_valid"], "3,360 both-valid")
check("oset bothinvalid", c3["n_both_invalid"], "4,794 both-invalid")

h9 = _json("h9_definitional.json")
assert h9["invariant_violations"] == 0
o, a = h9["optimal_set"], h9["arbitrary_z"]
check("h9 opt n", o["n_instances"], "| `Z = O(G₀)` | 5,304 |")
check("h9 opt pct", f"{100 * o['radii_identical'] / o['n_instances']:.2f}", "**100.00%**")
check("h9 opt r1", o["backdoor_r1"], "**0** of **2,790**")
assert o["of_those_gac_strictly_larger"] == 0
check("h9 arb n", a["n_instances"], "| `Z` arbitrary | 16,926 |")
check("h9 arb pct", f"{100 * a['radii_identical'] / a['n_instances']:.2f}", "70.86%")
check("h9 arb r1", a["backdoor_r1"], "1,308 of 6,588")
check(
    "h9 arb frac", f"{100 * a['of_those_gac_strictly_larger'] / a['backdoor_r1']:.2f}", "**19.85%**"
)
check("h9 total", o["n_instances"] + a["n_instances"], "22,230 instances at n = 3, 4")

# --- Phase 2: the generator -------------------------------------------------
pilot = _json("generator_pilot.json")
tot = pilot["totals"]
check("pilot accepted", f"{tot['accepted']:,} of {tot['attempts']:,}", "1,584 of 1,584")
check(
    "max separation",
    pilot["separation_range"]["max_realised_separation_overall"],
    "max realised separation is 11",
)
fid = pilot["realisation_fidelity"]
assert fid["intended_c_equals_realised_c"] == fid["n_attempts"]
assert fid["intended_s_equals_realised_s"] == fid["n_attempts"]

# --- Phase 3: the census ----------------------------------------------------
cen = _jsonl("census.jsonl")
acc = [r for r in cen if r.get("accepted") and "backdoor" in r]
check("census rows", len(cen), "1,584 instances over the full")
check("census measured", len(acc), "1,572 measured")
assert not [r for r in cen if r.get("censored")], "censored rows present"
assert not [r for r in cen if r.get("error")], "error rows present"
assert not [r for r in acc if not r["invariant_gac_ge_backdoor"]]
assert not [r for r in acc if r["backdoor"]["radius"] != r["gac"]["radius"]]
check("census gac diff", None, "**0 of 1,572 instances\ndiffer**")

full = [r for r in acc if r["coverage"] == 1.0]
eq = sum(1 for r in full if r["backdoor"]["radius"] == r["realised"]["realised_separation"])
check("r equals s", f"{eq} of {len(full)}", "**`r_val` equals `s` in 792 of 792 instances**")
assert eq == len(full)
for s in sorted({r["realised"]["realised_separation"] for r in full}):
    g = [r for r in full if r["realised"]["realised_separation"] == s]
    frac = 100 * sum(1 for r in g if r["backdoor"]["radius"] == 1) / len(g)
    want = 100.0 if s == 1 else 0.0
    if abs(frac - want) > 1e-9:
        fails.append(f"census s={s}: r=1 in {frac:.1f}%, expected {want}")
check("r1 at s>=2", None, "**0% at every separation ≥ 2**")
# flat in c at fixed s, coverage 1.0
for s in sorted({r["realised"]["realised_separation"] for r in full}):
    meds = {
        c: st.median(
            r["backdoor"]["radius"]
            for r in full
            if r["realised"]["realised_separation"] == s
            and r["realised"]["realised_component_size"] == c
        )
        for c in sorted(
            {
                r["realised"]["realised_component_size"]
                for r in full
                if r["realised"]["realised_separation"] == s
            }
        )
    }
    if len(set(meds.values())) > 1:
        fails.append(f"census s={s}: median radius NOT flat in c: {meds}")
check("flat in c", None, "median radius is identical across every")
mn = sum(
    1
    for r in acc
    if r["backdoor"]["radius"]
    == min(r["realised"]["realised_separation"], r["realised"]["realised_k_g0"])
)
check("min law", f"{mn} of {len(acc)}", "in 1,572 of 1,572 instances, 100.00%")
assert mn == len(acc)
cov05 = [r for r in acc if r["coverage"] == 0.5]
eq05 = sum(1 for r in cov05 if r["backdoor"]["radius"] == r["realised"]["realised_separation"])
check("cov05 r=s", f"{100 * eq05 / len(cov05):.1f}", "55.5%")
check(
    "cov05 binding",
    sum(1 for r in cov05 if r["realised"]["realised_k_g0"] < r["realised"]["realised_separation"]),
    "binding term in 347 of 780",
)

# --- The control ------------------------------------------------------------
rc = _jsonl("random_control.jsonl")
racc = [r for r in rc if r.get("accepted") and "backdoor" in r]
rnd = [r for r in racc if r["generator"] != "backdoor"]
des = [r for r in racc if r["generator"] == "backdoor"]
check("control attempted", len(rc), "**2,880 instances attempted")
check("control usable", len(racc), "2,334 usable**")
check("control errors", sum(1 for r in rc if r.get("error")), "60 errored")
rej = sum(1 for r in rc if not r.get("accepted") and not r.get("error"))
check(
    "control rejected", f"{rej} rejected ({100 * rej / len(rc):.1f}%)", "**486 rejected (16.9%)**"
)


def r1(g: list[dict]) -> float:
    """Percentage of a group at radius exactly 1."""
    return 100 * sum(1 for r in g if r["backdoor"]["radius"] == 1) / max(len(g), 1)


check("random pooled n", len(rnd), "| **1,914** |")
check("random pooled r1", f"{r1(rnd):.1f}", "**79.6%**")
for name, n_want, pct_want in (
    ("er_dense", 432, "85.0%"),
    ("er_medium", 389, "79.2%"),
    ("scale_free", 403, "77.9%"),
    ("block", 377, "78.0%"),
    ("er_sparse", 313, "77.0%"),
):
    g = [r for r in rnd if r["generator"] == name]
    if len(g) != n_want or f"{r1(g):.1f}%" != pct_want:
        fails.append(f"{name}: n={len(g)} r1={r1(g):.1f}%, report says {n_want} / {pct_want}")
check("designed r1", f"{r1(des):.0f}", "**0% at `r = 1`.**")
meas = [r for r in rnd if r["realised"]["realised_separation"] is not None]
sepd: dict[int, int] = {}
for r in meas:
    sepd[r["realised"]["realised_separation"]] = (
        sepd.get(r["realised"]["realised_separation"], 0) + 1
    )
check("sep 1", sepd.get(1), "| **1,362** | 192 | 29 | 4 |")
check("max nat sep", max(sepd), "Maximum natural separation is 4")
check("sep1 share", f"{100 * sepd[1] / len(meas):.1f}", "85.8% of measurable separations")
nz = len(rnd) - len(meas)
check(
    "no z in comp",
    f"**{100 * nz / len(rnd):.1f}%** of instances ({nz})",
    "**17.1%** of instances (327)",
)
small = sum(1 for r in rnd if r["realised"]["x_component_size"] <= 6)
check("comp<=6", f"{small:,} of the {len(rnd):,}", "**1,776 of the 1,914**")
check("comp<=6 pct", f"{100 * small / len(rnd):.1f}", "**92.8%**")
chk = [r for r in racc if r.get("gac_checked") and "gac" in r]
exact = [r for r in chk if not r["gac"].get("budget_exhausted")]
check("gac spotcheck", f"{len(chk)} of the {len(racc):,}", "469 of the 2,334")
check("gac exact", len(exact), "Of the 457 with an exact GAC radius")
assert not [r for r in exact if r["backdoor"]["radius"] != r["gac"]["radius"]]
assert not [r for r in chk if not r.get("invariant_gac_ge_backdoor")]

# --- H10 --------------------------------------------------------------------
fr = [r for r in _jsonl("frontier.jsonl") if r.get("accepted") and "r_opt_backdoor" in r]
tb = sum(r["n_backdoor_candidates"] for r in fr)
tg = sum(r["n_gac_candidates"] for r in fr)
check("frontier instances", len(fr), "0 / 420")
check("bd pool", f"{tb:,}", "| 27,552 |")
check("gac pool", f"{tg:,}", "| 40,740 |")
check("gac only pct", f"{100 * (tg - tb) / tg:.1f}", "**32.4% of the GAC pool")
check("gac only n", f"{tg - tb:,}", "13,188 candidate sets")
assert sum(r["backdoor_beats_opt"] for r in fr) == 0
assert sum(r["gac_beats_opt"] for r in fr) == 0
check(
    "pool larger",
    sum(1 for r in fr if r["n_gac_only_candidates"] > 0),
    "strictly larger on 256 of the 420",
)
assert all(r["r_opt_backdoor"] == r["r_opt_gac"] for r in fr)

# --- the gate ---------------------------------------------------------------
fg = _json("fast_gate_differential.json")
check("gate cases", f"{fg['n_cases']:,}", "46,800 cases")
check("gate disagree", fg["n_disagree"], "**16 disagreements**")
check("gate pct", f"{100 * fg['n_disagree'] / fg['n_cases']:.3f}", "(0.034%)")

ext = [r for r in _jsonl("frontier_ext.jsonl") if r.get("accepted") and "r_opt_backdoor" in r]
etb = sum(r["n_backdoor_candidates"] for r in ext)
etg = sum(r["n_gac_candidates"] for r in ext)
check("ext instances", len(ext), "Committed state: **242 instances**")
check("ext bd", f"{etb:,}", "**115,299** back-door")
check("ext gac", f"{etg:,}", "**194,794** GAC candidates")
check("ext gaconly", f"{etg - etb:,}", "(**79,495** of them GAC-only)")
assert sum(r["gac_beats_opt"] for r in ext) == 0
assert sum(r["backdoor_beats_opt"] for r in ext) == 0
check("combined pool", f"{tg + etg:,}", "**235,534 GAC-admissible candidate sets**")

re_rows = _jsonl("r_eps.jsonl")
re_acc = [r for r in re_rows if r.get("accepted") and "r_val" in r]
check("reps rows", len(re_acc), "**0.0%** | 560 instances" if False else "560 instances")
assert not [r for r in re_rows if r.get("censored") or r.get("error")]
for e in ("0.01", "0.02", "0.05", "0.1", "0.2", "0.5"):
    if any(r["middle_regime"][e] for r in re_acc):
        fails.append(f"r_eps: middle regime present at eps={e}, report says 0%")
check("reps zero", None, "| middle regime (`r_ε < r_val`) | **0.0%**")
if max(r["bias_at_g0"] for r in re_acc) > 1e-9:
    fails.append("r_eps: nonzero bias at G0, report says exactly 0")
check("reps bias0", None, "exactly 0 in all 560 instances**")
_lt = sum(
    1
    for r in re_acc
    if r["r_eps"]["0.05"] != -1 and r["r_val"] != -1 and r["r_eps"]["0.05"] < r["r_val"]
)
_eq = sum(1 for r in re_acc if r["r_eps"]["0.05"] == r["r_val"])
_un = sum(1 for r in re_acc if r["r_eps"]["0.05"] == -1)
assert _lt == 0
check("reps eq", _eq, "`r_ε = r_val` in 505 cases")
check("reps un", _un, "UNREACHED in 26")
check("reps max rval", max(r["r_val"] for r in re_acc), "radii out to 8")

for fig in re.findall(r"!\[[^\]]*\]\(([^)]+)\)", TXT):
    if not pathlib.Path(fig).exists():
        fails.append(f"figure missing on disk: {fig}")

print(f"anchors checked; {len(fails)} problem(s)")
for f in fails:
    print("  ", f)
sys.exit(1 if fails else 0)
