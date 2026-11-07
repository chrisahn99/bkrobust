# E13 preregistration: consistency checks on committed sets and the error budget
Dated 2026-11-06. Not to be edited after the first run (smoke runs excepted).
Serves reviewer qQYf, points M2 (which set Definition 3 commits) and Q1 (invalid set with a large radius).
## Question
(a) Are the 898 real-name units whose optimal set was recorded as undefined ("none" in
results/axis_robustness_llm_v2/status_split.json) unidentified at G0, i.e. not amenable?
(b) On every released row (release/per_query_records.csv), does the closed-form optimal set agree with
the optimal set read off by enumerating the DAG extensions of G0?
(c) On released rows whose committed set is invalid at the generating DAG, does that DAG reverse at
least r_val closure orientations of G0, as Corollary 1(c) requires?
## Registered statistics (criterion: GAC throughout)
(a) A = number of the 898 units with G0 amenable for (x, y); also how many of those get a GAC-valid
closed-form set. Share with a 95% Clopper-Pearson interval.
(b) B = number of closed_form rows whose enumeration is complete and where the extensions disagree on
O_D or agree on a set different from the committed Z; share with a Clopper-Pearson interval.
(c) C = number of ok rows invalid at the truth with c < bound, where c = |K_G0 \ dir(D_true)| (G0
orientations of CPDAG-undirected edges that the truth reverses) and bound = r_hop when exact, the
released hop_lower_bound when censored, infinity when r_hop is UNREACHED. Secondary check: where
Ext(CPDAG) can be enumerated (at most 256 extensions), the minimum of |K_G0 \ dir(D)| over extensions
D in which Z is invalid must equal r_hop (or exceed the released lower bound).
## Prediction
A = 0, B = 0, C = 0, and no mismatch in the secondary check.
## Against the authors
C > 0 or a secondary mismatch is a code error: no use of the error budget until it is found.
A > 0: App. H wording changes from "no optimal set" to "extensions disagree on O_D although the effect
is identified", and the excluded units are reported as identified. B > 0: the paper says the committed
set is the closed-form formula, not the enumerated optimal set, and reports B.
## Populations, seeds, exclusions
(a) the 898 real-condition "none" units; the 331 scrambled-name ones are reported apart. (b, c) all
8,308 released rows; D_RAND knowledge rebuilt with seed 20260913. Rows that fail the rebuild (hashes,
verdict, |Z|, counts) are excluded and counted. (b) excludes rows whose G0 has more than 4,096
extensions or whose enumeration exceeds 600 s, counted.
## Stopping rule
One full pass over every unit; no reruns except of units that raised an error.
## Deviations from the plan
None at registration. The committed Z is recomputed (the release stores only |Z|) and checked
against the released verdict and |Z|. r_val is the released r_hop (closure orientations), not
r_claim. Z is checked in every extension of G0 only when G0 has at most 512 extensions.
