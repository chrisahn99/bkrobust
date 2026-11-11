# E11 preregistration: survival ranking on queries chosen from the CPDAG alone (2026-11-06)

Written before any run. Not edited after the first full run; later notes go in a dated file beside it.

**Question.** Table 1's within-state ranking (Section 6.3, corrected 0.43 [0.07, 0.69]) was measured on queries admitted with five reads of the generating DAG (experiments/frame_build.py docstring). Does the radius still rank survival inside a fixed knowledge state when the queries are the frozen CPDAG-only frame (results/frame/frame.jsonl, 648 rows, 33 networks, sha256 6f05f044…43fc), which never reads the DAG? Serves reviewer point h8pC-W4b.

**Protocol.** The submitted elicited-knowledge survival pipeline, unchanged: K from results/elicit/knowledge.json (sha256 45119393…bfc7bb) for the eight real-name conditions; G0 = Meek closure (llm_survival.build_state); Z* = optimal set of G0, a valid empty set counted as committed (the corrected empty-set rule of real_survival.commit_z_star); survival = GAC validity of Z* under exactly d reversed claims, 1,000 draws per depth on the frac10 depth grid, draws seeded by the committed family llm|condition|network; endpoint AUC_frac (real_analyse._build_unit). Radius: breakdown_radius with criterion "gac" (registered) and, beside it, "backdoor".

**Registered statistic.** Within-state Kendall tau_b between the GAC radius and AUC_frac, computed as in Appendix G: states are (network, condition); units have status ok, a finite radius (not UNREACHED) and a defined AUC_frac; only states in which the radius varies enter; concordant, discordant and one-sided tie counts are summed over within-state pairs; tau_b = (C-D)/sqrt((C+D+Tx)(C+D+Ty)). 95% percentile interval from a network-cluster bootstrap, 10,000 resamples, seed 0, states selected once on the full sample.

**Prediction.** tau_b > 0 and its interval excludes zero; point estimate inside [0.07, 0.69].

**Outcome against the authors.** The interval includes zero (or lies below it). Then we say the within-state ranking holds on queries chosen with the truth and print both results in Section 6.3, with the CPDAG-only number first.

**Populations.** Primary: CPDAG-only frame x 8 real-name conditions. Every unit is written, whatever its status (n_k_zero, o_g0_extensions_intractable, optimal_set_undefined, z_invalid_at_g0, network_absent_from_bundle for pathfinder), and status counts are reported. Secondary, same code and seeds: (a) the CPDAG-only frame with the back-door radius; (b) the CPDAG-only frame without its 78 pairs that also lie in the Table 1 frame; (c) the Table 1 frame (results/axis_robustness_real/frame.jsonl, 543 distinct pairs) rescored under both criteria, the same-code reference; (d) per-condition values; (e) a within-state permutation p-value (2,000 shuffles, seed 0). The two scrambled-name conditions are not run by default.

**Exclusions.** Only those above, all fixed by the frozen eligibility rules. If fewer than 5 states show radius variation, the registered statistic is reported as undetermined and no direction is claimed.

**Consistency.** The run exits non-zero if a pair scored in both frames gets different radii or endpoints, if a GAC radius is 0 on a Z* valid at G0, or if the enumeration oracle gives a different GAC radius from the criterion on the units it checks. The analysis also reports how often the Table 1 frame rescored with the back-door radius reproduces the committed results/axis_robustness_llm units (plus the 94 corrected ones).

**Stopping rule.** One full run. A crashed run is resumed, never restarted with other settings. No change to frame, draws, grid, seeds or criterion after the first row is written.

## Deviations from the plan

1. The published 0.43 used radii from the back-door criterion (oracle "mpdag_criterion", commit 29e6c09, before the GAC switch a278b96); survival was GAC there too. The plan asks for "gac". So the GAC number is registered, the back-door number is reported beside it, and the Table 1 frame is rescored under both so the frame effect is not confounded with the criterion.
2. Separation is not computed: the CPDAG-only frame carries no separation column and the registered statistic does not need it.
3. The plan names run_llm_survival.py; its shard loop is reproduced in run.py with a frame adapter (row ids, component sizes) instead of being called, because it reads the Table 1 frame's row format only.
