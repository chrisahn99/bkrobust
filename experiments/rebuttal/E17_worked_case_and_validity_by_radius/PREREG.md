# E17 preregistration: worked applied case, and validity at truth by radius (2026-11-06)

Written before the first run; not edited after it. Serves reviewer point qQYf-S2 (practical motivation).

**Question.** In one published applied network, does one elicited knowledge state give radius 1 with a named claim to confirm for one query and a large radius (finite at least 2, or unreachable) for another? Across the corpus, is the committed set more often invalid at the true graph when the radius is 1, and is that gap larger for elicited knowledge (D_LLM) than for random orientations at matched |K| (D_RANDOM_MATCHED)?

**Units.** (a) worked case: every D_LLM query of Acid_1996 and Polzer_2012 listed in the released per_query_records.csv, recomputed with the submitted code. (b) table: the matched-control frame, 5 queries per network (seed 20260917, pathfinder skipped): D_LLM K from results/elicit/knowledge.json, and each of the 532 stored D_RANDOM_MATCHED draws (20 seeds per network, 2 on the six largest), exhausted draws excluded and counted. Per unit: G0 = Meek(CPDAG, K); Z = the optimal set at G0; radius = breakdown_radius with criterion "gac"; claim radius and witness = fewest asserted claims whose retraction breaks Z (depth at most 3, as the released records, 400 subsets); validity of Z at the data-generating DAG under the GAC (registered) and the back-door criterion (secondary).

**Registered statistic.** For each arm, share of informative units whose Z is GAC-invalid at truth, by radius bin {1, 2, >=3, unreachable}; and the gap G = share invalid at radius 1 minus share invalid at radius >= 2 or unreachable. Registered: Delta = G(D_LLM) - G(D_RANDOM_MATCHED), with a 95% percentile interval from 10,000 network-clustered paired bootstrap resamples (seed 20261106) over networks informative in both arms; units of one network are pooled across seeds. Informative = status ok and radius != 0. Radius 0 (degenerate), timeouts (120 s per unit), errors and unidentified Z are excluded from shares and reported as denominators. Resamples where a gap is undefined are dropped and counted. Secondary: the same table from the released records alone (r_hop bins, z_valid_at_truth), D_LLM rep 0 against D_RAND reps 0-2.

**Worked case rule.** First network in order Acid_1996, Polzer_2012 whose single D_LLM state has a query at radius 1 with a witness and a query with radius unreachable or >= 2. Report both queries: Z, radius, claim radius, witness, the witness claim's direction in the true graph, validity at truth. Labelled in the paper as read after the fact.

**Prediction.** The worked case exists in Acid_1996 (one claim, x4->x1, wrong at truth: radius-1 queries invalid at truth, unreachable queries valid). G > 0 in both arms. Delta's interval contains 0, as the matched control's share at radius 1 already did (+0.038 [-0.079, 0.134]).

**Against the authors.** (1) If G(D_RANDOM_MATCHED) is as large as G(D_LLM) (Delta interval containing 0), the table describes the real CPDAGs and the radius, not the knowledge: we will word it as "radius 1 marks the sets that fail when a claim is wrong, whatever the source of the claims" and make no claim about elicited knowledge. (2) If G(D_LLM) <= 0 with an interval below 0.05, radius 1 does not mark invalid sets in this corpus and the table is dropped; the worked case is then reported alone, as one case.

**Consistency checks (non-zero exit).** Criterion and enumeration oracles give the same radius where the CPDAG has at most 8 undirected edges; no unreachable unit is invalid at truth; a radius >= 1 unit invalid at truth has at least one closure orientation of G0 wrong at truth; every stored random draw is reproduced from its seed; on worked-case units the recomputed |Z|, radius and claim radius equal the released values.

**Stopping rule.** One full pass over all units. No re-run with other seeds, bins, time limits or criteria after seeing results; a crashed run is resumed, not restarted with new settings.

## Deviations from the plan

- The specification names validity for "the matched random state" from results/matched_control/arms/D_RANDOM_MATCHED, whose stored instances carry no Z and no validity, and whose r_val predates the switch of breakdown_radius to the GAC (2026-09-25). Both arms are therefore recomputed under the GAC on the stored draws; the stored r_val is kept in each row for reference only.
- The released records' own random arm (D_RAND, same edges as D_LLM) is reported as a secondary table; it is not the matched-|K| control.
