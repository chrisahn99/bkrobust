import json

def f(d):
    if d is None: return "n/a"
    if not isinstance(d, dict): return str(d)
    if d.get("rate") is None: return "-- (n=0)"
    return "%.4f [%.4f,%.4f] k=%d n=%d" % (d["rate"], d["ci"][0], d["ci"][1], d["k"], d["n"])

def rr(d):
    return "%.3f [%.3f,%.3f]" % (d["point"], d["ci"][0], d["ci"][1])

for e in ("original", "licensed", "large"):
    a = json.load(open("../results/x1_analysis_%s.json" % e))
    print("=" * 104)
    print("ENSEMBLE %s" % e)
    p2 = a["P2_damage"]
    print("  --- P2 damage (rho=1) ---")
    for k in ("PRIMARY_Suni_vs_R_reachable", "Suni_vs_R_all_incl_unreachable",
              "Suni_vs_R_hop0_only", "Suni_vs_R_mpdag_valid_only", "Sloc_vs_R_reachable"):
        b = p2[k]
        print("   [%s]" % k)
        print("      S(cond)=%s" % f(b["S_conditional"]))
        print("      R(cond)=%s   R(uncond)=%s   P(cons|R)=%s"
              % (f(b["R_conditional"]), f(b["R_unconditional"]), f(b["P_consistent_R"])))
        print("      RR=%s  -> %s" % (rr(b["RR_conditional"]), b["verdict_M1_riskratio"]))
        print("      diff(uncond)=%.4f [%.4f,%.4f]  free=%.4f class=%.4f"
              % (b["diff_unconditional"]["point"], b["diff_unconditional"]["ci"][0],
                 b["diff_unconditional"]["ci"][1],
                 b["diff_decomposition"]["free_component"],
                 b["diff_decomposition"]["class_component"]))
    d = p2["ADDED_Suni_vs_DROP_reachable"]
    print("   [ADDED Suni vs DROP]  S=%s  DROP=%s" % (f(d["S"]), f(d["DROP"])))
    print("      RR=%s -> %s" % (rr(d["RR"]), d["verdict"]))
    d = p2["ADDED_DROP_vs_R_reachable"]
    print("   [ADDED DROP vs R]     DROP=%s  R(cond)=%s  RR=%s"
          % (f(d["DROP"]), f(d["R_conditional"]), rr(d["RR"])))
    st = p2["M9_3_standardised_Suni_onto_R_hopdist"]
    print("   [M9.3 standardised S-uni onto R hop dist] = %.4f  (n_covered=%d)  weights=%s"
          % (st["value"], st["n_covered"], json.dumps({k: round(v, 4) for k, v in st["weights"].items()})))
    print("   hop strata S-uni rho=1 : " + " | ".join("%s:%s" % (k, f(v)) for k, v in p2["hop_strata_Suni_rho1"].items()))
    print("   hop strata R rho=1 cond: " + " | ".join("%s:%s" % (k, f(v)) for k, v in p2["hop_strata_R_rho1_conditional"].items()))
    pr = p2["PREREG_4_2_absolute_rule_DEMOTED"]
    print("   PREREG 4.2 absolute rule (DEMOTED): diff=%.4f disjoint=%s -> %s"
          % (pr["diff"], pr["disjoint"], pr["verdict"]))
    p3 = a["P3_A26"]
    print("  --- P3 A26 (S-loc) ---")
    c = p3["composition_Sloc_pool"]
    print("   composition: pure_spurious %s | reversal %s | query_pair %s | SCMs w/ >=1 query pair %s"
          % (f(c["pure_spurious"]), f(c["reversal_of_true_edge"]), f(c["query_pair"]),
             f(c["scms_with_ge1_query_pair"])))
    for k in ("pure_spurious_NON_QUERY_PRIMARY", "pure_spurious_QUERY_PAIR",
              "reversal_of_true_edge", "POOLED_all"):
        print("   %-34s rho=1 %s | rho<=4 %s"
              % (k, f(p3[k]["member_rho1"]), f(p3[k]["member_rho_le4"])))
    ib = p3["IN_RUN_reversal_baseline_member_rho_le4"]
    print("   IN-RUN reversal baseline rho<=4: uncond %s | cond %s"
          % (f(ib["unconditional"]), f(ib["conditional"])))
    av = p3["A26_verdict"]
    print("   A26 window [0.087,0.100]; primary %s ; overlaps=%s"
          % (f(av["primary_rate"]), av["overlaps_window"]))
    p4 = a["P4_locality"]
    print("  --- P4 locality (rho=1, S-uni) ---")
    print("   hop0 %s" % f(p4["Suni_hop0"]))
    print("   hop>=1 reachable %s" % f(p4["Suni_hop_ge1_reachable"]))
    print("   unreachable %s" % f(p4["Suni_hop_unreachable"]))
    print("   materiality ratio = %s ; n_hop>=1 = %d ; verdict %s (events=%d)"
          % (("%.4f" % p4["materiality_ratio"]) if p4["materiality_ratio"] else "n/a",
             p4["n_hop_ge1_denominator"], p4["verdict_M12"]["verdict"],
             p4["verdict_M12"]["n_events"]))
    print("   EDGE-ATTRIBUTABLE O* move at hop>=1 : %s" % f(p4["ADDED_edge_attributable_hop_ge1"]))
    print("   hopC x hopG (S-uni, silent):")
    for k, v in p4["Suni_hopC_x_hopG_silent_rho1"].items():
        if v["n"]: print("      %-28s %s" % (k, f(v)))
    print("   hopC x hopG (S-uni, EDGE-ATTRIBUTABLE O* move):")
    for k, v in p4["Suni_hopC_x_hopG_EDGE_ATTRIBUTABLE_Ostar_move_rho1"].items():
        if v["n"]: print("      %-28s %s" % (k, f(v)))
    p5 = a["P5_severity"]
    cd = p5["cliffs_delta_Suni_vs_R"]
    print("  --- P5 severity --- Cliff delta=%s ci=%s n_S=%d n_R=%d median_S=%s median_R=%s -> %s"
          % (("%.4f" % cd["delta"]) if cd["delta"] is not None else "n/a",
             ("[%.4f,%.4f]" % tuple(cd["ci"])) if cd["ci"][0] is not None else "n/a",
             cd["n_a"], cd["n_b"],
             ("%.4f" % cd.get("median_a", float("nan"))), ("%.4f" % cd.get("median_b", float("nan"))),
             p5["dominance"]))
    p6 = a["P6_movement"]["channels"]
    print("  --- P6 movement (censoring at n=20000, z=1.960) ---")
    key = "n=20000|z=1.960"
    for arm in ("R", "Suni", "Sloc", "Drop"):
        b = p6[arm][key]
        print("   %-5s rho*_se %s | rho*_ident %s | rho*_any %s"
              % (arm, f(b["censored_rho_star_se"]), f(b["censored_rho_star_ident"]),
                 f(b["censored_rho_star_any"])))
