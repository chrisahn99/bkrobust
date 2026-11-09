# E28 preregistration: joint radius over orientation retractions and edge additions near the query (2026-11-07)

Written 2026-11-07, frozen (sha256 in `_queue/PREREG_SHA256.txt`) before the first full run. Not edited after it;
changes go in `DEVIATIONS.md`. Serves C2Bd-W2, C2Bd-Q3, h8pC-W2, qQYf-Q1. Nothing from this experiment is printed in
the paper or the response; it writes a result folder and a summary.

**Disclosure before freezing.** The code was smoke-tested on 6 rows (3 invalid) before this file was frozen, and the
smoke output was read: all six rows had joint radius 1, and on the 3 invalid rows no true missing edge was a
distance-1 witness (0 of 3), with median 61 witness pairs against 2 orientations in K. The registered statistics,
bars and predictions below are the plan's, written before any code, and are not changed by that reading. The smoke
rows belong to the population and are scored again in the full run.

**Question.** If the radius search is extended from retracting the analyst's orientations to also adding edges near
the query, do the distance-1 additions that break the committed set name the edge whose absence actually broke it at
the truth, with a list short enough for an analyst to confirm?

**Population.** Every scored row of E22 (`results/rebuttal/E22/rows.shard0of1.jsonl`, sha256 f7929e06...bd34):
984 rows on learned PC and GES CPDAGs of ecoli70, arth150, magic-niab and magic-irri, arm A_ORACLE_R, 534 invalid at
the truth and 450 valid; the three rows E22 excluded (no DAG extension) stay out. The registered population is the
452 invalid rows whose G0 misses at least one true edge inside An_true({x, y} u Z) (E22 field `s3_missing_in_A`;
the other 82 invalid rows have only reversed edges there).

**Rebuild check.** For each row, G (the learned CPDAG), K (E5's `oracle_r_knowledge`) and Z
(`committed_adjustment_set(G0, x, y)`) are rebuilt with E5's functions; Z and |K| must equal the stored row, Z must
be GAC-valid in G0, and its validity in the true DAG must equal the stored one. Any mismatch is a consistency
failure (exit 2).

**Moves and state.** A move is (a) the retraction of one orientation of K, or (b) the addition of one directed edge
a -> b between two nodes non-adjacent in G and both within two hops of {x, y} u Z in the skeleton of G (both
directions of a pair are separate moves; a pair counts once in witness lists). The state after retractions R and
additions A is `meek_closure(G + A)` followed by `apply_orientations(., K minus R)`, the radius's own primitive; the
learned orientations of G are kept as learned. A state whose closure fails is inconsistent and is not a break. Z
breaks when it fails the GAC in the state (`is_gac_valid_mpdag` of the submitted code).

**Joint radius r_joint.** The fewest moves after which Z breaks. Depth 1 is exhaustive (every single move, time limit
1,800 s per row; an incomplete depth 1 is recorded) and lists every breaking move (the witnesses). Depth 2 runs only
when depth 1 has no break, in the order retraction pairs, retraction plus addition, addition pairs, stopping at the
first break or at 600 s (censored). Above 2 is recorded as "more than 2".

**Orientation-only radius r_orient.** The ledger's `breaking_depths(G, K, x, y, Z)` (retractions only, depth 3,
300-subset budget), computed on the same rows, for the per-row comparison.

**Registered statistics,** on the 452 rows, each with its share and a 95 % percentile interval from a bootstrap
that resamples learned graphs (2,000 draws, seed 20261107):
- S1: share with r_joint equal to 1 or 2 (censored and incomplete rows count as not). Bar 0.8.
- S2 (coverage, the plan's kill condition): share of rows where at least one true missing edge in An_true({x, y} u Z) is a
  distance-1 addition witness, in either direction (an incomplete depth 1 counts as not covered). Bar 0.8.
- S3 (burden): median number of distance-1 witness pairs against the median |K| on the same rows.

**Prediction (the plan's).** S1 at least 0.8; S2 at least 0.8; S3 median witness pairs at most median |K|. The plan
also predicts r_joint = 1 on most valid rows, which is not a death: the claim is coverage and burden, not ranking.

**Against the method (dies if).** S2 below 0.8, or S3 median witness pairs above median |K|. Either result is
analysed (which true edges were missed and why: outside the region, an inconsistent state, or a break that needs
several true edges together) and is not printed. S1 at or above 0.8 together with the same share on valid rows says
the joint radius has no ranking power; that is reported beside it.

**Reported with it, not registered.** S1 on all 534 invalid rows and on the 450 valid rows; share with r_joint = 1;
S2 with the true direction only; share of rows whose true missing edge lies in the search region; whether adding all
true missing edges in A (true direction) breaks Z, and the fewest of them (up to 3) that do; witness-pair quantiles
on invalid and valid rows; cross-tabulation of r_orient against r_joint by validity; AUC for invalidity of 1/r_orient,
1/r_joint (infinite = 0) and of the witness-pair count; the GAC reason of the first three addition witnesses; each by
network.

**Not covered.** Whether r_joint differs across sets valid in G0 (the plan's secondary statistic) is not computed. Reversed
edges are not moves; the 82 reversed-only rows enter only the unregistered lines.

**Consistency failures (exit 2).** A rebuild mismatch; Z not valid in G0; a breaking addition state where the added
edge is not directed as added or whose skeleton is not G's plus the added pair (first three per row); depth-1
retraction witnesses that disagree with r_orient = 1.

**Seeds, machine, stopping.** Deterministic except the bootstrap. The workstation, 6 shards on 6 cores (nice 10,
taskset, 2 GB address space each). One launch; a resume only reruns rows not written. No change of limits after the
first full run; any rerun is reported in DEVIATIONS.md.
