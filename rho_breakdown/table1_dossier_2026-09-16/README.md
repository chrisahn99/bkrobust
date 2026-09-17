# Table 1, every version of it — design dossier for the breakdown-radius paper

**2026-09-16.** Asked for by Zé: every possible Table 1 for the bkrobust ICLR 2027 paper, with the
claims stated, their significance and originality, and the comparison with the literature made
visible. Deliverable: `TABLE1-DOSSIER.pdf`, 15 pages, 13 candidates, every cell filled from a file
or computed here.

```
TABLE1-DOSSIER.tex          preamble + includes
sections/01_brief.tex       what this is, what it was built from, the nine requirements
sections/02_claims.tex      C1–C11: statement, evidence, significance, originality, falsifier
sections/03_comparison.tex  the named prior work and where each enters a float
sections/04_space.tex       the five design axes, and the combinations already dead
sections/05_candidates.tex  T1, T1b, T2, T3, T4, T5
sections/05b_candidates.tex T6–T13
sections/06_scorecard.tex   the two matrices, the recommendation, the blocking list
```

## Recommendation

**Table 1 is T11**: two panels, the certificate ladder over `over-certification` and the price
table with the screen row. Every cell is measured. It carries C3 (the unit) and C4 (the price),
and gives a number to four named comparisons inside the float: the model-oriented distance as the
revision-unit instrument, the claim count, the whole-class interval, and the leave-one-out screen.

Then T3 as Table 2, T8 as Table 3 (or its numbers in the results prose), T4 in the ranking section
with the screen column added, T5 as the first figure, T6 and T7 as setup and limitation, T9/T10/T13
in the appendix.

Two defensible alternatives, neither better: T2 alone if the merged paper keeps the diagnostic
framing; T4 with its fourth column if the ranking result is judged the strongest thing either side
has.

## Where the numbers come from

Nothing was written from memory. Sources, all read this session:

| candidate | source |
|---|---|
| T1, T1b | `~/bkrobust/results/interval/TABLE1_L95.md` |
| T2, T9, T13 | computed here from `~/bkrobust/results/ledger/rows.csv` (queries in the session transcript) |
| T3 | `~/bkrobust/results/interval/TABLE2_REPORTS.md` |
| T4 | the draft's `sections/06_simulation.tex`, table `tab:tau` |
| T5 | `~/bkrobust/results/ledger/TABLE8_DECISION.md`, `CHECKS.md` |
| T6 | `~/bkrobust/results/ledger/TABLE0_FUNNEL.md`, `funnel.json` |
| T7 | `~/bkrobust/results/pest/TABLE7.md` |
| T8, T11 | `~/bkrobust/results/interval/PREREGISTRATION.md` appendices A and B |
| T10 | `~/bkrobust/results/ledger/TABLE4.md` |

The aggregate certificate ladder of T2 and T11 Panel A was computed for this document:
3,014 certifiable rows over eleven deployment conditions and 27 networks, of which 902 committed
sets are invalid at the truth. Claim unit 0 over-certifications, graph-revision unit 39, free
separation rule 15, `|K|` 729, `k_g0` 791. This reproduces the ledger's published per-condition
convention (`dangerous` over certifiable rows; detection over the invalid ones).

## Design provenance

The dossier reuses rather than re-ideates, per the 2026-09-12 instruction. The paid design work it
reads: `paper-table-ideation` (wf_07538611-a07), `paper-table-selection` (wf_ed5b8f98-5d1),
`main-table-design` (wf_9198cc6a-2cf), `radius-redefinition-design` (wf_0cf50270-0b5),
`iclr-table-archetypes` (wf_a7abeacd-24f), and Part IV of the 2026-09-09 apostila —
210 agents, 26.2 M context tokens. Four candidates (T2, T3, T7, T8) could not have been written
before the 11–13 September measurements and are new here.

## Relation to this morning's draft review

`output/2026-09-16_chris-draft-review/` reviewed the same draft a few hours earlier. This dossier is
its companion and carries two of its findings: the gradedness point (Taeb et al. Prop. 9 proves the
general MPDAG poset is not graded for |V| ≥ 4; this space is graded because it freezes the skeleton,
and saying so turns the liability into positioning) into C2, and the licensed screen sentence
(the radius says *which depth* was needed) into C5 and T8.

One number needed reconciling. The review quotes the depth-one screen gap as **+0.305 [0.182, 0.410]**
and this dossier had **+0.256 [0.149, 0.346]**. Both are in `TABLE1_L95.md` lines 393–394 and both
are right: 0.305 is the population row, 0.256 the n = 1000 row (0.269 at n = 20000). Every quotation
of either now carries its sample size.

## Two defects found and fixed while writing

1. **A mixed population inside one row.** The first draft of T8 and T11 Panel B filled the
   two-wrong-claim cells of the screen and stop-rule rows from the M1b committed-rows population
   (59 rows) while the neighbouring cells came from the price-table population (131 rows). Both are
   now on one population per panel, with the refusal number moved out of the panel and its own
   denominator printed. This is constraint R5, broken in the act of writing the document that
   states it.
2. **The free rule conflated with the screen.** The 15 over-certifications belong to the
   treatment-to-set separation rule, not to the leave-one-out re-closure screen. The screen never
   over-certifies relative to the stop rule: on all 2,212 rows a refusal at depth *r* is matched by
   a screen triggering at depth at most *r*.

## Blocking, before the 18/09 abstract

- **The novelty search for C1 has never been run in the wording the claim uses.** The register
  entry `rho-breakdown-sensitivity-structural` carries `verified_NOT_RUN` deliberately: the three
  sweeps on record were pointed at latent-causal representation learning and at Taeb, not at the
  sensitivity-analysis lineage the claim excludes, and the framing rests on a statement about
  arXiv:2010.08611 that nobody here has read in full. Either run it, or narrow C1 to the scope
  actually searched and let the abstract say the narrower thing.
- **The author list is written nowhere**, and 18/09 locks it.

## Owed, not blocking

The screen column of T4 (half a day, the corruption states are stored); one validity oracle shared
by both sides, the generalised adjustment criterion rather than the back-door one; the
anti-exchange caveat replaced by the monotonicity audit's number; the sensitivity row of T1; the
corrupted-knowledge calibration on the frame, which is being rerun after the seed fix.

## Build

```bash
latexmk -pdf TABLE1-DOSSIER.tex      # 15 pages, zero overfull boxes
```

Prose gate: `paper_lint` reports 0 errors on every section but `02_claims.tex`, whose single error
is `I1` on the word `verified_NOT_RUN`. That is internal vocabulary, and it is deliberate here:
this is an internal document and constraint R9 says so. The remaining warnings are emphasis density
from bolded winning cells and the vault's own words (`register`, `falsifier`), which do not travel
into the paper.
