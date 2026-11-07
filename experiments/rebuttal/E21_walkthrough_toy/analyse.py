"""E21 analysis: the registered radius per toy, its agreement across methods, and hand vs code.

The radii are exact, so the reported interval is the degenerate [r, r]; what can vary is
agreement, and that is what the table shows.

    python experiments/rebuttal/E21_walkthrough_toy/analyse.py --out results/rebuttal/E21
"""

from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path


def load(out: Path) -> dict[str, dict]:
    rows: dict[str, dict] = {}
    for f in sorted(glob.glob(str(out / "rows*.jsonl"))):
        for line in open(f):
            if line.strip():
                r = json.loads(line)
                rows[r["unit"]] = r  # a resumed run appends; the last row wins
    return rows


def hand_checks(r: dict) -> dict[str, object]:
    """Each hand prediction against what the code printed: True, False, or None (not predicted)."""
    h = r.get("hand", {})
    out: dict[str, object] = {}
    if "z" in h:
        out["z"] = sorted(h["z"]) == r.get("z")
    if "r_val" in h:
        out["r_val"] = h["r_val"] == r.get("r_val")
    if "r_claim" in h:
        out["r_claim"] = h["r_claim"] == r.get("r_claim")
    if h.get("n_space") is not None:
        out["n_space"] = h["n_space"] == r.get("n_space")
    if "n_upset" in h:
        out["n_upset"] = h["n_upset"] == r.get("n_upset")
    if "upset_valid" in h:
        got = {s["k"]: s["valid_gac"] for s in r.get("upset", [])}
        out["upset_valid"] = got == h["upset_valid"]
    if "witnesses" in h:
        out["witnesses"] = sorted(h["witnesses"]) == r.get("upset_witnesses")
    if "n_dags_below_g0" in h:
        below = r.get("covered_by_g0", [])
        out["dags_below_g0"] = (
            sum(b["is_dag"] for b in below) == h["n_dags_below_g0"]
            and sum(b["is_dag"] and b["valid_gac"] for b in below) == h["n_dags_below_g0_valid"]
        )
    if "non_states" in h:
        ghosts = {g["subset"] for g in r.get("non_states", [])}
        out["non_states"] = set(h["non_states"]) <= ghosts
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/rebuttal/E21")
    args = ap.parse_args(argv)
    out = Path(args.out)
    rows = load(out)
    if not rows:
        print(f"no rows in {out}")
        return 1

    res: dict[str, object] = {"experiment": "E21", "units": {}}
    print(f"{'unit':<20}{'Z':<16}{'r_val':>6}{'[lo,hi]':>9}{'agree':>7}{'r_claim':>8}"
          f"{'|space|':>8}{'legacy':>7}{'|up|':>5}{'r_bd':>6}  hand")
    any_fail = False
    for name, r in rows.items():
        if r.get("status") != "ok":
            print(f"{name:<20}ERROR: {r.get('consistency_failures')}")
            res["units"][name] = {"status": r.get("status"), "failures": r.get("consistency_failures")}
            any_fail = True
            continue
        radii = r.get("radii", {})
        agree = len(set(radii.values())) == 1
        rv = r.get("r_val")
        hc = hand_checks(r)
        wrong = [k for k, v in hc.items() if v is False]
        any_fail |= bool(r["consistency_failures"])
        rv_s = "inf" if rv == -1 else str(rv)
        print(f"{name:<20}{'{' + ','.join(r.get('z', [])) + '}':<16}{rv_s:>6}"
              f"{'[' + rv_s + ',' + rv_s + ']':>9}{str(agree):>7}{r.get('r_claim', ''):>8}"
              f"{r.get('n_space', ''):>8}{r.get('n_space_legacy_enumerate_space', ''):>7}"
              f"{r.get('n_upset', ''):>5}{r['hybrid']['backdoor']['radius']:>6}  "
              f"{'all match' if not wrong else 'MISMATCH ' + ','.join(wrong)}")
        res["units"][name] = {
            "status": "ok", "commit": r.get("commit"), "criterion": r.get("criterion"),
            "z": r.get("z"), "r_val": rv, "r_val_interval": [rv, rv], "radii": radii,
            "radii_agree": agree, "r_claim": r.get("r_claim"), "phi_1": r.get("phi_1"),
            "r_backdoor": r["hybrid"]["backdoor"]["radius"],
            "n_space": r.get("n_space"), "n_space_legacy": r.get("n_space_legacy_enumerate_space"),
            "n_upset": r.get("n_upset"), "upset": r.get("upset"),
            "upset_witnesses": r.get("upset_witnesses"),
            "space_witnesses": r.get("witness_knowledge"),
            "covered_by_g0": r.get("covered_by_g0"), "non_states": r.get("non_states"),
            "added_by_meek": r.get("added_by_meek"),
            "hand": r.get("hand"), "hand_checks": hc, "hand_mismatches": wrong,
            "consistency_failures": r["consistency_failures"],
        }
    for name, u in res["units"].items():
        if u.get("status") == "ok":
            print(f"\n{name}: up-set of G0 (k = orientations held; r = retractions from G0)")
            for s in u["upset"]:
                print(f"   r={s['retractions_from_g0']}  {'valid  ' if s['valid_gac'] else 'INVALID'}"
                      f"  {{{s['k']}}}")
            print(f"   nearest failures in the up-set: {u['upset_witnesses']}")
            if u["consistency_failures"]:
                print(f"   CONSISTENCY FAILURES: {u['consistency_failures']}")
    res["consistency_ok"] = not any_fail
    (out / "results.json").write_text(json.dumps(res, indent=1))
    print(f"\nwrote {out / 'results.json'}; consistency_ok={not any_fail}")
    return 0 if not any_fail else 2


if __name__ == "__main__":
    raise SystemExit(main())
