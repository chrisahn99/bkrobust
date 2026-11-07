# E15 preregistration: adjustment identification distance as a survival baseline (2026-11-06)

Written before any run. Not edited after the first full run; later notes go in a dated file beside it.

**Question.** Does the adjustment identification distance (AID) of Henckel, Würtzen and Weichwald (2024), computed for each instance, rank survival as well as the breakdown radius? Serves h8pC-W3 (novelty against graph-space distances) and C2Bd-W3 (cheap proxies rank about as well).

**Populations.** The 480 untargeted instances (results/axis_robustness_unconfounded_v2/survival_instances.csv, sha256 4d911e7f…308a) and the 4,104 targeted instances (results/axis_robustness_p6/survival_instances.csv, sha256 1a5f44fe…de2; 2,730 flip, 1,374 tiered). Each graph is regenerated from its stored seed with the submitted builders (unconfounded_generator.build_unconfounded_tiered_instance, survival.build_flip_instance, survival.build_tiered_instance). Survival is the stored AUC_frac (1,000 draws per grid point); it is not recomputed. No seeds are drawn by this experiment.

**Consistency (the run exits non-zero if any fails).** Rebuild succeeds; rebuilt x equals stored x; rebuilt Z* equals stored Z*; the radius from hybrid.breakdown_radius with criterion "gac" is exact and equals the stored r_val; on the untargeted instances the rebuilt separation equals the stored one; the pair-level CPDAG AID of (x, Y) equals 1 exactly when the radius is finite (both say Z* fails somewhere in the class). On targeted instances separation is a generator input, so a separation measured from the tiered or flip Z* is recorded beside it and not required to match.

**AID.** For an ordered pair (t, v) the guess G0 answers: "no effect" if v is not a possible descendant of t in G0; "not identifiable" if G0 is not amenable for (t, v) or its optimal set differs across the extensions of G0; else the optimal set O(t, v, G0). Against a reference graph T it is a mistake when: "no effect" but v is a possible descendant of t in T; "not identifiable" but T is amenable for (t, v); a set that is not GAC-valid for (t, v) in T. Two references: the true DAG (oracle, not available to the analyst) and the CPDAG (available to the analyst). Three scopes: all ordered pairs (the published graph-level AID, count divided by n(n-1)), pairs with treatment x (row, divided by n-1), and the single pair (x, Y). gadjid's oset_aid, parent_aid and ancestor_aid are computed beside it as the published implementation; agreement of our oset count with gadjid's is reported.

**Registered statistic.** On the 480 untargeted instances pooled: Delta = tau_b(radius, AUC_frac) - tau_b(-AID_cpdag, AUC_frac), with AID_cpdag the graph-level normalised oset AID from G0 to the CPDAG (our implementation). The sign is fixed in advance: more mistakes should mean lower survival. Interval: 95% percentile, paired instance bootstrap, 10,000 resamples, seed 0 (paired_resample.paired_comparison, cluster = instance_id).

**Prediction.** Delta > 0 and its interval excludes zero.

**Outcome against the authors.** The interval of Delta includes zero or lies below it, or AID with the opposite sign has |tau_b| at least the radius's tau_b. Then we write that the AID available to the analyst ranks survival as well as the radius, claim no ranking advantage over it, and keep exactness, the witness and the certified region as the radius's contribution.

**Secondary, same code.** The same Delta for the oracle AID, the row and pair scopes, and the three gadjid variants; per n_tiers stratum; on the targeted instances per arm and per (arm, component size). Marginal tau_b of each predictor on pooled strata. Bootstrap resamples: 10,000 on the untargeted pool, 2,000 in every other stratum, to keep the analysis near an hour.

**Exclusions.** None by choice. A unit whose AID cannot be computed is written with its status and counted; a predictor with fewer than 30 defined units in a stratum is reported as undetermined.

**Stopping rule.** One full run. A crashed run is resumed, never restarted with other settings. No change to populations, definitions, criterion or bootstrap after the first row is written.

## Deviations from the plan

1. The plan says AID "per query". Henckel et al. define AID over all ordered pairs of a graph; each instance here carries one query, so the graph-level AID is per query. The row and pair scopes are added so that a query-local reading is also scored. The pair scope is nearly constant by construction (consistency check above), so it cannot rank.
2. gadjid's handling of a guess that is an MPDAG (not a CPDAG) could not be checked without running it; if it refuses G0, gadjid columns carry the error and the registered statistic, which does not depend on gadjid, is unaffected.
