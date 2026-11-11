"""E21: compute the walkthrough toy and the Figure 2(a) toy exactly.

For each toy: enumerate every knowledge state of the CPDAG, the order distance
to G0 (Corollary 1(a)), validity of the committed Z = O(G0) under the GAC, the
hybrid radius, a brute-force radius over closed retraction subsets, all
nearest witnesses, the claim radius, and the join used in the walkthrough.
"""

from __future__ import annotations

from collections import Counter

from common import (MPDAG, OUT, apply_orientations, claim_radius, dist, dump,
                    hybrid, is_gac_valid_mpdag, k_of, rad_str,
                    radius_bruteforce_retractions)

from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag
from bkrobust.search.space_fixed import enumerate_space_correct

TOYS = {
    "walkthrough": dict(
        nodes="UXYW", und=[("U", "X"), ("X", "Y"), ("X", "W"), ("W", "Y")],
        k=[("U", "X"), ("X", "Y"), ("X", "W")], x="X", y="Y"),
    "figure2a": dict(
        nodes="AWXY", und=[("A", "W"), ("W", "X"), ("W", "Y"), ("X", "Y")],
        k=[("A", "W"), ("X", "Y")], x="X", y="Y"),
}


def run(name, spec):
    cpdag = MPDAG(list(spec["nodes"]), undirected=spec["und"])
    x, y = spec["x"], spec["y"]
    g0 = apply_orientations(cpdag, spec["k"])
    assert g0 is not None
    z = optimal_adjustment_set_mpdag(g0, x, y)
    out = {"toy": name, "cpdag": cpdag.edge_string(), "K": spec["k"],
           "G0": g0.edge_string(), "K_G0": sorted(k_of(cpdag, g0)),
           "Z": sorted(z) if z is not None else None}
    if z is None:
        out["note"] = "optimal set undefined at G0"
        return out
    z = frozenset(z)
    out["Z_valid_at_G0"] = is_gac_valid_mpdag(g0, x, y, z)
    out["Z_valid_at_Chat"] = is_gac_valid_mpdag(cpdag, x, y, z)
    space = enumerate_space_correct(cpdag)
    rows = []
    for g in space:
        d = dist(cpdag, g0, g)
        rows.append({"state": g.edge_string(), "K": sorted(k_of(cpdag, g)), "d": d,
                     "above_G0": k_of(cpdag, g) <= k_of(cpdag, g0),
                     "Z_valid": is_gac_valid_mpdag(g, x, y, z)})
    out["n_states"] = len(space)
    shells = Counter((r["d"], r["Z_valid"]) for r in rows)
    out["shells"] = {d: {"states": sum(v for (dd, _), v in shells.items() if dd == d),
                         "Z_invalid": shells.get((d, False), 0)}
                     for d in sorted({r["d"] for r in rows})}
    fail_d = [r["d"] for r in rows if not r["Z_valid"]]
    out["r_by_full_space"] = min(fail_d) if fail_d else "inf"
    out["nearest_failing_states"] = [r for r in rows if not r["Z_valid"]
                                     and r["d"] == out["r_by_full_space"]]
    h = hybrid(cpdag, g0, x, y, z)
    out["r_hybrid"] = rad_str(h.radius)
    out["hybrid_method"] = h.method
    out["hybrid_witness"] = h.witness
    rb, wits = radius_bruteforce_retractions(cpdag, g0, x, y, z)
    out["r_retraction_bruteforce"] = rb
    out["retraction_witnesses"] = wits
    rc, cw = claim_radius(cpdag, spec["k"], x, y, z)
    out["r_claim"] = rc
    out["claim_witnesses"] = cw
    out["up_set_states"] = sum(r["above_G0"] for r in rows)
    out["states_table"] = rows
    # Walkthrough join: of G0 and a nearest failing state.
    if out["nearest_failing_states"]:
        nf = apply_orientations(cpdag, out["nearest_failing_states"][0]["K"])
        j = apply_orientations(cpdag, sorted(k_of(cpdag, g0) & k_of(cpdag, nf)))
        out["join_G0_nearest"] = {"state": j.edge_string(),
                                  "d_to_G0": dist(cpdag, g0, j),
                                  "Z_valid": is_gac_valid_mpdag(j, x, y, z)}
    return out


if __name__ == "__main__":
    res = {n: run(n, s) for n, s in TOYS.items()}
    dump(res, OUT / "e21" / "toys.json")
    for n, r in res.items():
        print(f"== {n}: Z={r['Z']} |K_G0|={len(r['K_G0'])} states={r.get('n_states')} "
              f"r_full={r.get('r_by_full_space')} r_hybrid={r.get('r_hybrid')} "
              f"r_brute={r.get('r_retraction_bruteforce')} r_claim={r.get('r_claim')} "
              f"up_set={r.get('up_set_states')} valid_at_Chat={r.get('Z_valid_at_Chat')}")
        print("   G0:", r["G0"], " K_G0:", r["K_G0"])
        print("   shells:", r.get("shells"))
        print("   retraction witnesses:", r.get("retraction_witnesses"))
        print("   claim witnesses:", r.get("claim_witnesses"))
        print("   join:", r.get("join_G0_nearest"))
