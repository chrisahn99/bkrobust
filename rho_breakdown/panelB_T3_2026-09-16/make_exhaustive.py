#!/usr/bin/env python3
"""The exhaustive T3: every query, no real-coefficient rule, no cap of five.

Reads results/interval/profiles.jsonl under the recovering set (A_ORACLE) and emits
sections/_t3_exhaustive.tex, a landscape longtable of all 303 rows that commit to an
adjustment set, plus sections/_t3_noset.tex listing the 234 rows that do not.

The label helpers reproduce experiments/interval_tables.py exactly; the 26 rows T3
prints are cross-checked against table2.json cell by cell before anything is written.
"""
import csv
import json
import os
import re
from collections import Counter

SRC = os.path.expanduser("~/bkrobust/results/interval")
HERE = os.path.dirname(os.path.abspath(__file__))
CENSORED = 4

# ----------------------------------------------------------------- label helpers
def r0_seed_range(p):
    vals = p["r0_band_seeds"]
    lab = [">3" if v == CENSORED else str(v) for v in (min(vals), max(vals))]
    return lab[0] if lab[0] == lab[1] else f"{lab[0]}..{lab[1]}"


def r0_label(v, status="exact"):
    return str(v) if v != CENSORED else ("never" if status == "never" else ">3")


def rval_label(v, status):
    if v not in ("", None):
        return str(v)
    if status in ("unreached", "no_retractable_edges"):
        return "never"
    if status.startswith("gt_"):
        return ">" + status.split("_")[1]
    if status.startswith("censored_at_depth_"):
        return ">=" + status.rsplit("_", 1)[1]
    return status


# ----------------------------------------------------------------- latex helpers
def esc(s):
    return str(s).replace("_", "\\_").replace("&", "\\&").replace("%", "\\%")


def qesc(s):
    """Variable names in some networks are long, full of underscores, or CamelCase, none of
    which TeX will break. Allow a break after each underscore and at each case change."""
    out = re.sub(r"(?<=[a-z])(?=[A-Z])", lambda _: chr(92) + "allowbreak{}", str(s))
    return esc_keep(out).replace("\\_", "\\_\\allowbreak{}")


def esc_keep(s):
    """esc(), but without mangling the \\allowbreak{} we just inserted."""
    return s.replace("_", "\\_").replace("&", "\\&").replace("%", "\\%")


def num(v, d=3):
    return f"${v:+.{d}f}$".replace("+", "\\phantom{-}" if v >= 0 else "")


def iv(v):
    return f"[{num(v[0])}, {num(v[1])}]"


def sym(s):
    return esc(s).replace(">=", "$\\geq$").replace(">", "$>$")


def claims(p):
    v = p.get("breaking_claims_seed0")
    if not v:
        return "--"
    return " \\,/\\, ".join(f"{esc(a)}$\\to${esc(b)}" for a, b in v)


# ----------------------------------------------------------------- load
def load():
    P = [json.loads(l) for l in open(f"{SRC}/profiles.jsonl")]
    A = [p for p in P if p["condition"] == "A_ORACLE"]
    R = [r for r in csv.DictReader(open(f"{SRC}/rows.csv")) if r["condition"] == "A_ORACLE"]
    return A, R


def verify_against_published(A):
    """Every cell T3 prints must come out identical from this path."""
    T = json.load(open(f"{SRC}/table2.json"))
    idx = {(p["network"], p["x"], p["y"]): p for p in A}
    bad = 0
    for r in T["rows"]:
        x, y = [s.strip() for s in r["query"].split("->")]
        p = idx[(r["network"], x, y)]
        b = {d["r"]: d for d in p["budgets"]}
        checks = [
            (r["n_k"], p["n_k"]), (r["n_k_loc"], p["n_k_loc"]),
            (round(r["tau_true"], 9), round(p["tau"], 9)),
            (round(r["estimate"], 9), round(p["estimate"], 9)),
            (r["I1_band"], b[1]["band_seed0"]), (r["I2_band"], b[2]["band_seed0"]),
            (r["r0_pop"], r0_label(p["r0_pop"], p["r0_pop_status"])),
            (r["r0_band_seed0"], r0_label(p["r0_band_seeds"][0])),
            (r["r0_band_range"], r0_seed_range(p)),
            (r["r_val"], rval_label(p["r_val"], p["r_val_status"])),
            (r["breaking_claims"], p["breaking_claims_seed0"]),
        ]
        for a, c in checks:
            if a != c:
                print(f"  MISMATCH {r['network']} {r['query']}: {a!r} vs {c!r}")
                bad += 1
    print(f"verify: {len(T['rows'])} published rows, {bad} cell mismatches")
    return bad == 0


# ----------------------------------------------------------------- emit
COLSPEC = "@{}llp{2.15cm}ccrl cccc c l@{}"   # 13 columns
HEADER = (r"network & param. & query & $|K|$ & loc. & true & estimate [95\,\%] & "
          r"one revision & two revisions & no knowledge & $r_0$ & $r_{\mathrm{val}}$ & "
          r"first breaking claim \\")


def row_tex(p, printed):
    b = {d["r"]: d for d in p["budgets"]}
    mark = "$\\bullet$" if printed else ""
    par = "fitted" if p["params"] == "fitted" else "semi-s."
    return (f"{esc(p['network'])}{mark} & {par} & {qesc(p['x'])} $\\to$ {qesc(p['y'])} & "
            f"{p['n_k']} & {p['n_k_loc']} & {num(p['tau'])} & "
            f"{num(p['estimate'])} {iv(p['estimate_ci'])} & "
            f"{iv(b[1]['band_seed0'])} & {iv(b[2]['band_seed0'])} & {iv(p['class_pop'])} & "
            f"{sym(r0_label(p['r0_pop'], p['r0_pop_status']))} & "
            f"{sym(rval_label(p['r_val'], p['r_val_status']))} & {claims(p)} \\\\")


def key(p):
    """Two queries are drawn both as a declared pair and, independently, into the frame
    (Thoemmes_2013 x->y and mediator X->Y), so the source has to be part of the key."""
    return (p["network"], p["x"], p["y"], p["source"], p["frame_index"])


def printed_in_t3(A):
    PANEL_B = ("ecoli70", "arth150", "magic-niab", "magic-irri")
    out = {key(p) for p in A if p["source"] == "declared"}
    for net in PANEL_B:
        c = [q for q in A if q["network"] == net and q["source"] == "frame" and q["z"] is not None]
        c.sort(key=lambda q: (q["r0_pop"] == CENSORED, q["frame_index"]))
        out |= {key(q) for q in c[:5]}
    return out


def main():
    A, R = load()
    assert verify_against_published(A)
    mark = printed_in_t3(A)

    committed = [p for p in A if p["z"] is not None]
    declared = sorted([p for p in committed if p["source"] == "declared"],
                      key=lambda p: p["network"])
    frame = sorted([p for p in committed if p["source"] == "frame"],
                   key=lambda p: (p["params"] != "fitted", p["network"], p["frame_index"]))

    L = [r"\begin{landscape}", r"\begingroup", r"\fontsize{6.3}{7.3}\selectfont",
         r"\setlength{\tabcolsep}{1.5pt}", r"\renewcommand{\arraystretch}{1.04}",
         r"\begin{longtable}{" + COLSPEC + "}",
         r"\caption*{\textbf{T3-X, the exhaustive table.} Every query under the recovering set "
         r"whose analyst graph commits to an adjustment set: " + str(len(committed)) +
         r" rows over 27 networks. Neither rule of the printed T3 is applied here, so networks with "
         r"parameters we drew appear beside the four with coefficients fitted to real data "
         r"(\emph{param.}), and no per-network cap is used. A bullet after the network name marks "
         r"the 26 rows the printed T3 selects. \emph{loc.} is the number of asserted claims inside "
         r"the treatment's own chain component. The \emph{no knowledge} column is the interval over "
         r"the whole equivalence class, which is the ceiling every other interval is bounded by. "
         r"Estimates at $n=20\,000$, first seed; the two revision columns are bands at that seed, "
         r"the class column is the population.}\\",
         r"\toprule", HEADER, r"\midrule", r"\endfirsthead",
         r"\toprule", HEADER, r"\midrule", r"\endhead",
         r"\midrule \multicolumn{13}{r@{}}{\itshape continues on the next page}\\", r"\endfoot",
         r"\bottomrule", r"\endlastfoot",
         r"\multicolumn{13}{@{}l}{\itshape Queries the source papers declare.}\\[1pt]"]
    L += [row_tex(p, key(p) in mark) for p in declared]
    L += [r"\addlinespace",
          r"\multicolumn{13}{@{}l}{\itshape Sampled queries, networks with coefficients fitted to "
          r"real data.}\\[1pt]"]
    L += [row_tex(p, key(p) in mark)
          for p in frame if p["params"] == "fitted"]
    L += [r"\addlinespace",
          r"\multicolumn{13}{@{}l}{\itshape Sampled queries, networks whose parameters we drew "
          r"(structure from a published source).}\\[1pt]"]
    L += [row_tex(p, key(p) in mark)
          for p in frame if p["params"] != "fitted"]
    L += [r"\end{longtable}", r"\endgroup", r"\end{landscape}"]
    open(f"{HERE}/sections/_t3_exhaustive.tex", "w").write("\n".join(L) + "\n")

    # the rows with no adjustment set, listed so nothing is omitted
    noset = [r for r in R if r["committed_verdict"] != "closed_form"]
    per = Counter(r["network"] for r in noset)
    N = [r"\begingroup\footnotesize",
         r"\setlength{\tabcolsep}{5pt}",
         r"\begin{longtable}{@{}lrl@{}}",
         r"\toprule", r"network & queries & verdict on the analyst's own graph \\", r"\midrule",
         r"\endfirsthead", r"\toprule",
         r"network & queries & verdict on the analyst's own graph \\", r"\midrule", r"\endhead",
         r"\bottomrule", r"\endlastfoot"]
    for net in sorted(per):
        N.append(f"{esc(net)} & {per[net]} & the outcome cannot descend from the treatment \\\\")
    N += [f"\\midrule total & {len(noset)} & \\\\", r"\end{longtable}", r"\endgroup"]
    open(f"{HERE}/sections/_t3_noset.tex", "w").write("\n".join(N) + "\n")

    print(f"committed rows written: {len(committed)} "
          f"(declared {len(declared)}, fitted {sum(1 for p in frame if p['params']=='fitted')}, "
          f"semi-synthetic {sum(1 for p in frame if p['params']!='fitted')})")
    print(f"no-set rows listed: {len(noset)} over {len(per)} networks")
    print(f"marked as printed in T3: {sum(1 for p in committed if key(p) in mark)}")
    print(f"networks represented: {len({p['network'] for p in committed})}")


if __name__ == "__main__":
    main()
