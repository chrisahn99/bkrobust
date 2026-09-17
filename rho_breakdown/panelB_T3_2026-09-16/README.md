# What Panel B measures, T3 in full, and the figure that should carry it

**2026-09-16.** Companion to `output/2026-09-16_table1-dossier/`. Questions from Zé: what exactly are the
numbers in T11 Panel B, the full T3, the best figure for T3, why there are only 26 rows, why the
one- and two-revision columns barely differ, and an exhaustive table with both selection rules off.

Deliverable: **`PANELB-T3.pdf`**, 25 pages, 11 figures, the exhaustive 303-row table.

```
PANELB-T3.tex          preamble + includes
sections/01_panelB.tex  the metric explained in three steps, with five figures
sections/02_t3.tex      the full 26-row table (landscape) + what it says
sections/03_figure.tex  the design constraint, five figure candidates, the recommendation
sections/04_why26.tex   why 26 rows: the funnel, the two selection rules, what they biased
sections/05_exhaustive.tex  why one and two revisions coincide; the exhaustive table
sections/_t3_exhaustive.tex GENERATED — all 303 committed rows, 13 columns, landscape longtable
sections/_t3_noset.tex      GENERATED — the 234 rows with no adjustment set, by network
make_exhaustive.py          GENERATED both; cross-checks all 26 published cells first
sections/_t3_table.tex  GENERATED — the whole T3 float, written from table2.json
make_figs.py            GENERATED all ten figures; reads the run files at plot time
figs/*.pdf              f1..f11
```

Rebuild: `python3 make_figs.py && python3 make_exhaustive.py && latexmk -pdf PANELB-T3.tex`

## The answer to question 1, in one paragraph

Each Panel B cell is `L95 / raw coverage`. Raw coverage is the fraction of problems on which the
interval the procedure actually reported contained the true effect. **L95 is not a width in effect
units.** For each problem and each draw, take λᵢ = the factor by which that interval would have to
be stretched, about its own midpoint, to reach the truth. Take the 95th percentile of λᵢ over the
whole cell, call it λ\*, multiply every interval by it, average the widths, and divide by the same
quantity for the procedure that discards the knowledge. So the cell answers: *forced to be honest,
how wide would this have to be, relative to using no background knowledge at all?*

Worked, for the two headline cells (one wrong claim, n = 1000, semi-synthetic panel, 7,180 replicates):

| procedure | λ\* | mean half-width | calibrated width | L95 |
|---|---|---|---|---|
| report the estimate as committed | 20.416 | 0.3855 | 15.7397 | **5.123** |
| interval over one revision | 1.000 | 1.3064 | 2.6129 | **0.851** |
| discard the knowledge, whole class | 1.000 | 1.5361 | 3.0721 | 1.000 |

A value **above 1** means the procedure would have to be wider than throwing the knowledge away, and
it is earned by being narrow and centred in the wrong place, because stretching happens about the
interval's own midpoint.

## Verification

`make_figs.py` reimplements the L95 pipeline from `replicates.csv.gz` + `rows.csv` and **reproduces
all 27 published Panel A cells of `TABLE1_L95.md` to three decimals** (`verify()`, asserted at every
run). Nothing in the report is transcribed from the published table; the numbers come out of the
replicate data.

The column filter that makes it reproduce: a problem is in a column iff it is `informative` **and**
`b_realised` is non-empty. That gives 391 / 359 / 131 rows, which are the counts printed under the
panel. My first attempt joined on `profiles.jsonl` (934 rows, a subset) and produced 0.180 instead
of 0.464 — worth knowing, because `profiles.jsonl` is the wrong join key for anything cell-level.

## The answer to question 2

The full T3 is the landscape table on p. 6: all 26 rows, 10 columns, generated from `table2.json`.
Three declared queries are excluded with their reasons printed.

What it says: **8 of the 9 certifiable declared queries are unbroken by any retraction to depth 3**,
while **16 of 17 informative sampled queries break at the first**. The mechanism is the sampling
gate — a query no retraction can touch is rejected before it is measured — so that contrast is an
applicability statement, not a result.

## The answer to question 3

🔑 **The design constraint nobody wrote down: the sensitivity curve is not a curve, it is one step.**
On 16 of the 17 informative rows, the interval after retracting **one** claim equals the whole-class
interval to four decimals. All the identifying power sits in a single claim. This is why the paper's
own pre-registered Figure 1 rule (flat, then a break) returned **no candidate** — that shape is not
in the data.

Five candidates drawn: F-A the one-claim cliff (`f6`), F-B the step law (`f7`), F-C the state strip
(`f8`), F-D the paradigm card (`f9`), F-E the cliff with its population (`f6`+`f10`).

**Recommendation: F-A with the population panel of F-E beside it**, F-D as the page-one teaser on a
different query, F-B as a sentence, F-C in the appendix. F-A is the only candidate carrying the
quantity, the consequence in the analyst's own units, and the comparison, at once. Its caption
sentence: *on these analyses, retracting a single claim returns the analyst to the interval they
would have had if they had asserted nothing.*

Two cautions before it becomes the main **claim** rather than the main **figure**, both in §3.5:
selection, which §4 measures (the width claim survives, the null-radius column does not), and that
the figure reads as showing the **null** radius while the paper's headline quantity is the
**validity** radius — they coincide on most of these rows and are different objects.

## The answer to question 4 (why only 26 rows)

T3 is a selection from a 303-row eligible population, by two rules.

| | stage, under the recovering set |
|---:|---|
| 39 | networks parsed |
| 33 | have a frame |
| 27 | are swept (6 exceed 500 nodes or the questionnaire cap) |
| 80,132 | candidate pairs from observables |
| 648 | sampled into the frame (20 per network cap) |
| 537 | rows: 528 sampled + 9 declared, all certifiable |
| 234 | say the outcome cannot descend from the treatment — no set to report on |
| **303** | **commit to a set: 9 declared + 294 sampled, over 27 networks** |
| 26 | printed in T3 |

**Rule one, defensible:** only 4 of the 27 networks have coefficients fitted to real data
(`ecoli70`, `arth150`, `magic-niab`, `magic-irri`). The other 23 carry parameters we drew, so a
"sensitivity report on a named analysis" would print an invented effect size. That excludes 257
eligible rows. It should be in the caption and currently is not.

**Rule two, NOT defensible as written:** on those 4 networks, 37 queries are eligible and the table
takes the first 5 per network (17; arth150 has only 2). The cap is arbitrary, but the *ordering* is
the problem — `cand.sort(key=lambda q: (q["r0_pop"] == CENSORED, q["frame_index"]))` puts queries
with a **finite null radius first**, which is the outcome the table exists to illustrate.

🔴 Measured: finite null radius on **76%** of the 17 printed, **15%** of the 20 dropped, **43%** of
the 37 eligible. The `r0` column of T3 is not representative. Fix: drop the `r0_pop == CENSORED`
term from the sort, or print all 37.

🟢 **What the selection did NOT bias** — the claim the recommended figure rests on:

| population | queries | informative | one retraction reaches the whole class |
|---|---|---|---|
| printed in T3 | 17 | 16 | 16 / 16 |
| dropped by the cap | 20 | 10 | 10 / 10 |
| full eligible pool | 37 | 26 | **26 / 26** |
| every committed row, 27 networks | 294 | 226 | **211 / 226 (93%)** |

So "one retraction returns you to the whole class" is the pattern of the corpus, not of the cap.
`f11_pool_dumbbell.pdf` redraws F-A on all 37, marking which are printed in T3; it also shows the
11 rows where nothing moves, which the table currently hides.

## The answer to question 5 (one vs two revisions)

The interval over *r* revisions is a union over knowledge states, so it grows with *r* and is
**bounded above by the whole-class interval** — retracting everything is what the class already
assumes. The two columns can only differ where the first revision has not yet hit that ceiling.

Across all 303 committed rows:

| rows | |
|---:|---|
| **281 (92.7%)** | identical at one and two revisions, band and population |
| 76 | class interval is a point — no asserted claim bears on the query at any budget |
| 205 | the first revision already returns the whole class |
| 22 | first revision strictly inside the class; on **20** of those the second widens it |

🔑 Where those 20 sit is the point: `insurance`, `Polzer_2012`, `Kampen_2014`, `sachs`,
`Schipf_2010`, `mediator`. All but `mediator` are **excluded from the printed T3 by the
real-coefficient rule**. The column looks empty in T3 because T3 drops the networks on which it is
full.

**Recommendation:** in T3, replace *two revisions* with *no knowledge at all* (the class interval).
On the saturated rows the reader then sees two identical intervals side by side, and that identity
is the finding rather than a redundancy.

## The answer to question 6 (the exhaustive table)

`T3-X`, §6: **all 303 rows** that commit to an adjustment set, over **all 27 networks**, with
neither rule applied — no real-coefficient restriction, no cap of five, no finite-r₀-first sort.
13 columns including the class interval; a bullet marks the 26 rows the printed T3 selects; three
blocks (declared / fitted-coefficient / semi-synthetic) so the dropped rule stays visible.

`make_exhaustive.py` re-derives every field from `profiles.jsonl` and **asserts all 26 published
rows agree cell by cell with `table2.json`** before writing anything.

Accounting: 303 committed + 234 with no set = 537 = the whole sweep under the recovering set. The
234 are listed by network in §6.1, so nothing is omitted.

⚠️ Two queries appear twice — `Thoemmes_2013` x→y and `mediator` X→Y are drawn both as a declared
pair and as an independently sampled frame row. Same query, two routes; both rows kept.

## Palette and chart rules

Data-viz reference instance: categorical slots 1–3 (`#2a78d6`, `#eb6834`, `#1baf7a`), a sequential
blue ramp for the ordered `b` variable, neutral gray for the reference procedure. Thin marks,
hairline solid grids, one axis per panel, legend on every ≥2-series figure, direct labels used
selectively. Every plotted number is computed at plot time.
