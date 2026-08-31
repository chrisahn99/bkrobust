import json, sys

def f(d):
    if d is None: return "n/a"
    if not isinstance(d, dict): return str(d)
    if d.get("rate") is None: return "-- (n=0)"
    return "%.4f [%.4f,%.4f] k=%d n=%d" % (d["rate"], d["ci"][0], d["ci"][1], d["k"], d["n"])

for e in ("original", "licensed", "large"):
    a = json.load(open("../results/x1_analysis_%s.json" % e))
    h = a["header"]
    print("=" * 104)
    print("ENSEMBLE %s   n_analysed=%d  readable=%s" % (e, h["n_analysed"], h["verdict_readable"]))
    print("  stmt/arm rho=1 = %d | hop>=1 reachable = %d | unreachable = %d"
          % (h["n_stmt_rho1_per_arm"], h["n_hop_ge1_reachable_Suni_rho1"],
             h["n_hop_unreachable_Suni_rho1"]))
    print("  P(|O0|=0) = %s   mean|O0| = %.3f" % (f(h["P_O0_empty"]), h["mean_O0"]))
    p1 = a["P1_catchability"]
    print("  --- P1 catchability ---")
    print("   R meek-inconsistent rho=1   : %s" % f(p1["R_meek_inconsistent_rho1"]))
    print("   R meek-inconsistent rho<=4  : %s" % f(p1["R_meek_inconsistent_rho_le4"]))
    for arm in ("Suni", "Sloc"):
        d = p1[arm]
        print("   %-5s S1 caught (intrinsic) rho=1 : %s" % (arm, f(d["S1_caught_intrinsic_rho1"])))
        dd = d["S1_decomposition_rho1"]
        print("         decomp: conflict %s | cycle %s | nonext-acyclic %s"
              % (f(dd["conflict"]), f(dd["cycle"]), f(dd["non_extendable_acyclic"])))
        v = d["vstruct_vs_C_DIAGNOSTIC_rho1"]
        print("         vsC fires %s" % f(v["fires"]))
        print("             purely shielding artefact %s | any new unshielded %s"
              % (f(v["purely_collider_of_C_shielded_by_added_edge"]),
                 f(v["any_new_unshielded_collider"])))
        print("         mpdag_valid rho=1 %s" % f(d["mpdag_valid_rho1"]))
        print("         mpdag_valid by rho: " + "  ".join(
            "%s:%.4f(n=%d)" % (r, d["mpdag_valid_by_rho"][r]["rate"], d["mpdag_valid_by_rho"][r]["n"])
            for r in "1234"))
        print("         path_blowup %s" % f(d["path_blowup"]))
    g = a["gates"]
    print("   G3 skel-detect S-uni %s | S-loc(MEASURED) %s"
          % (f(g["G3_skeleton_detect_Suni"]), f(g["G3_skeleton_detect_Sloc_MEASURED"])))
    print("   G8 closure-order disagreement %s" % f(g["G8_closure_order_MEASUREMENT"]))
    print("   G12 unreachable: edge-attributable moves %d / %d ; raw moves %d  -> %s"
          % (g["G12_unreachable_structural_zero"]["edge_attributable_Ostar_moves"],
             g["G12_unreachable_structural_zero"]["n"],
             g["G12_unreachable_structural_zero"]["raw_moves_vs_G0"],
             g["G12_unreachable_structural_zero"]["verdict"]))
    print("   design effects (rho=1 silent): " + json.dumps(
        {k: round(v, 3) for k, v in a["design_effects_rho1_silent"].items()}))
