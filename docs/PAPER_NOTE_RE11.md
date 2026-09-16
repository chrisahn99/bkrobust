# What this discharges, and what it forces the paper to say differently

**Audience:** whoever drafts the ICLR 2027 paper next.
**Scope:** this note is about `[RE-11]` only. Every number it refers to lives in
`results/axis_robustness_real/` and is re-asserted by
`src/bkrobust/analysis/session9_verify.py`.

---

## 1. The TODO this discharges

`docs/PAPER_NARRATIVE.md` §6 marks **Table 4** — *"Worst-case radius against
average-case survival, same instances"* — as `TODO [RE-11]` for the real-network
version, and §8 answers the reviewer objection *"Synthetic only"* with *"True of
the survival phase today. Say so, and cite `[RE-11]` as the completion."*

The results section's TODO reads:

> *The ranking evaluation of §6.1 has not yet been run on these 831 real
> instances; survival and the paired cross-arm design exist only on synthetic
> structure. This is the experiment that turns the ranking claim into a claim
> about practice, and Table 1 must not be read as covering it.*

**That TODO can now be deleted.** It is replaced by `table_tau_real.md` and by
the sentences in §3 below. The experiment ran on the committed corpus at full
pre-registered precision, and the result is **weaker than the synthetic one**.

---

## 2. What the paper must now say differently

### 2.1 Result B needs a real-structure qualifier, not a real-structure victory

The draft's Result B says the worst-case radius ranks average-case survival where
the baselines do not do so consistently, and instructs the writer to report the
synthetic result honestly (8 of 9 intervals excluding zero, two unstable cells,
one confirmed null). **On real structure the same measurement is materially
weaker**, and the paper should say so in the same breath as the synthetic result
rather than in a later limitations paragraph.

This costs the paper very little, because §6.2 currently promises nothing about
real structure. It would cost the paper everything to be caught overstating it.

### 2.2 Result D's number does not transfer, and the reason is a finding

Result D currently says most orientation errors are self-revealing — Meek closure
fails about six times in ten. **On the committed real corpus that rate is near
zero**, and the reason is not the graphs. Holding the same 25 networks fixed and
changing only how the analyst's knowledge is written down moves the exhaustive
single-reversal contradiction rate by more than an order of magnitude.

So Result D is a property of **how the knowledge is phrased**, not of the
structure. The paper should state it that way. This is a strictly more
interesting claim than the current one, and it makes the diagnostic framing
*more* necessary rather than less: if errors do not announce themselves, pricing
them is the only thing left to do.

It also raises the stakes on `[RE-1]`, which asks exactly how real elicited
knowledge compares with a minimal generator. That item is no longer a tidying-up
exercise; it decides which of two very different numbers the paper is entitled
to print.

### 2.3 The two-radii result now has a real-structure, average-case demonstration

`docs/PAPER_NARRATIVE.md` §4 presents the geometric/claim distinction as a
measurement result and states the coincidence condition (`|K| == k_g0`). This
campaign supplies the missing half: **what the gap costs when it is there**,
measured against an average-case endpoint on real graphs.

`paths` is the figure. An analyst who made **one** statement, which Meek closure
turned into **fourteen** committed orientations, gets a hop radius of up to 14 —
and reversing that one statement destroys the adjustment set every time, on 1000
draws per row with zero contradictions. `SHD(G₀, truth)` is 0: the analyst is
perfectly correct. This is the over-promise §4 warns about, with a number
attached.

**Recommendation:** make this the paper's Figure 1, or a panel of it. §6 of the
narrative already reserves Figure 1 for *"one real instance — the first-failing
retraction named in domain terms, and the hop count that missed it"*. This is
that figure, and it now carries a survival curve as well as a hop count.

### 2.4 The certificate sentence must be restated in commitment units

`[RE-16]` already asks for this on `hybrid.py:225`. This campaign turns the
request into a measured consequence rather than a style note: on the corpus the
leverage runs to 17, and a certificate denominated in stated claims is wrong by
that factor where it is wrong.

---

## 3. Sentences a reviewer can check against committed files

Each of these is re-asserted by `session9_verify.py` against the committed CSVs.

All of these are checkable against `results/axis_robustness_real/`.

> On the 831 admissible instances of the real-network corpus, the breakdown
> radius ranks average-case survival of a committed adjustment set positively in
> all nine pre-registered strata, but its 95% interval — bootstrapped over the 25
> networks, which are the unit of analysis — excludes zero in only three of them,
> and only one of those three survives the removal of any single network.

> Decomposed, the pooled association is a mixture. **Within a fixed network and
> knowledge state — the comparison an analyst actually faces — the radius ranks
> survival with an interval excluding zero in six of the nine strata, at a mean
> Kendall τ_b between +0.30 and +0.84.** Across networks it is confounded with the
> size of the asserted knowledge set and with a resolution limit of the corruption
> model, and contributes little or nothing.

> Neither baseline is dominated. `SHD(G₀, truth)` outranks the radius in five of
> the eight strata where it is defined, and `|K|` in four of nine — though on this
> corpus both are constants within a network, so they rank networks rather than
> queries.

And, for Result D:

> On the real corpus, reversing a single truthful orientation claim produces
> knowledge inconsistent with the CPDAG in only 5 of 298 cases, measured
> exhaustively. That is a property of how the knowledge was written down rather
> than of the graphs: asserting every undirected edge instead of a minimal
> generating set, on the same 25 networks, raises the rate to 177 of 410.

---

## 4. What must NOT be said

- **Do not** claim the radius dominates the baselines on real structure. It does
  not, and the anchor table shows it.
- **Do not** quote the coverage-0.25 strata's significant *negative* τ as evidence
  that the radius ranks survival backwards. It is −0.775 between networks and
  +0.158 within them, and is produced by three `|K| = 1` networks whose endpoint
  is pinned at zero by the depth grid.
- **Do not** present the within-network figure as the pre-registered result. It is
  post-hoc, and `results/axis_robustness_real/PREREGISTRATION.md` Appendix J says
  so.
- **Do not** carry Result D's 0.627 into any sentence about real structure.
- **Do not** compare any real row with a synthetic row, or pool the two corpora.

---

## 5. The table that replaces the TODO

`table_tau_real.md`, the nine pre-registered strata, in the same order as the
synthetic anchor table so the two can be read side by side — **and never pooled**.
`docs/PAPER_NARRATIVE.md` §6's **Table 4** is now writable.

For **Figure 1**, which §6 reserves for *"one real instance — the first-failing
retraction named in domain terms, and the hop count that missed it"*: use `paths`.
It now carries a survival curve as well as a hop count, and the survival number is
an exhaustive enumeration rather than an estimate.

---

## 6. What this costs the paper, and what it buys

**Costs.** §8's reviewer defence *"Synthetic only — true of the survival phase
today"* can no longer be answered with a promise, and the honest answer is that
the real-structure result is weaker. Result D's headline number has to be
restated.

**Buys.** A real-structure result at all, with an explanation for its weakness
that is itself a finding; a quantitative demonstration of the two-radii
over-promise on a real network; a direct measurement bearing on `[RE-1]`; and a
pre-registration with ten appendices recording four falsified predictions and four
defects, which is the kind of record that makes a reviewer read the rest as
calibration rather than salesmanship.
