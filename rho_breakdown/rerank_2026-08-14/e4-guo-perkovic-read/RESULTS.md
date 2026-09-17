# E4 — THE READ THE FRAMING RESTS ON

**Date:** 2026-08-14 · **Verdict: SUPPORTED** (delta INTACT, three narrowings)
**Both papers fetched and read in full text. `ran = true`.**

---

## 0. Commands that produced everything below

```bash
cd /Users/josecosta/mugango/output/2026-08-14_latent-causal-iclr2027-rerank/e4-guo-perkovic-read
curl -sL -A "Mozilla/5.0" -o guo-perkovic-2010.08611.pdf   https://arxiv.org/pdf/2010.08611   # HTTP 200, 1,025,666 B
curl -sL -A "Mozilla/5.0" -o taeb-guo-henckel-2511.10625.pdf https://arxiv.org/pdf/2511.10625 # HTTP 200, 1,407,021 B
pdftotext -layout guo-perkovic-2010.08611.pdf   guo.txt    # 1,135 lines /  68,647 chars / 24 pp
pdftotext -layout taeb-guo-henckel-2511.10625.pdf taeb.txt  # 2,605 lines / 187,663 chars / 49 pp
```

**Provenance.** Guo & Perković, *Minimal enumeration of all possible total effects in a Markov
equivalence class*, arXiv:2010.08611**v3** (3 Mar 2021, math.ST), 24 pp, no journal-ref on the abs
page. Taeb, Guo & Henckel, *Model-oriented Graph Distances via Partially Ordered Sets*,
arXiv:2511.10625**v3** (30 Jun 2026), 49 pp — confirming the campaign's "v3 posted 30 Jun 2026,
17 days pre-campaign" and its "49 pages".

**Extraction integrity check** (a null here would invalidate every zero below):
`guo.txt` "the"×540, "graph"×80, "causal"×110 · `taeb.txt` "the"×1468, "graph"×397, "causal"×90.
Text is intact; the zeros reported below are real absences, not extraction failures.

**Evidence standard.** Every quote below is grep-reproducible from `guo.txt` / `taeb.txt` sitting
beside this file, with line numbers. No WebFetch summarizer was used for any load-bearing claim.

---

## 1. Pre-registered decision rule (written before either fetch)

See `PREREGISTRATION.md`, written and saved **before** the first `curl`.

> **REFUTED** iff arXiv:2010.08611 does ANY of:
> (a) defines a distance/radius/tolerance/ball/neighbourhood **over background-knowledge sets**;
> (b) considers K wrong/misspecified/possibly-false beyond a future-work line **and** computes a
>     downstream consequence;
> (c) reports a **bias magnitude** attributable to a wrong orientation constraint;
> (d) weights/ranks/filters enumerated members by relevance to the **optimal** adjustment set O*.
>
> **SUPPORTED** iff none of (a)–(d) fires.

**Result: none of (a)–(d) fires. → SUPPORTED.**

---

## 2. Q1–Q5 on Guo & Perković 2010.08611

### Q1. Do they perturb, relax, or question K — or is K assumed true throughout?

**K IS ASSUMED TRUE THROUGHOUT. It is never perturbed, doubted, or varied.**

Verbatim term counts over the full 24 pages:

| term | count | term | count |
|---|---|---|---|
| `perturb` | **0** | `sensitivity` | **0** |
| `misspecif` | **0** | `robust` | **0** |
| `incorrect` | **0** | `erroneous` | **0** |
| `wrong` | **0** | `inconsisten` | **0** |
| `false` | **0** | `conflict` | **0** |
| `possibly false` | **0** | `background knowledge` | 12 |

All 12 `background knowledge` hits are the **source** of the MPDAG refinement or bibliography.
The load-bearing one (guo.txt:63–78):

> "Often, we may have **additional background knowledge** on the underlying causal system. For
> example, we may know that A temporally precedes Y and therefore determine (or reveal) the edge
> orientation A → Y in CPDAG C. […] The background knowledge of pairwise causal relationships of
> this type **can be derived from field expertise**."

"determine (or reveal)" — knowledge *reveals* a true orientation. It is an input, taken as given.

The only failure mode they model is **Meek-inconsistency**, and it is inherited pseudocode, not a
contribution (guo.txt:705–717, Algorithm 2, credited to "Meek, 1995 and Algorithm 1 of Perković
et al., 2017"): `input: MPDAG G, set of background knowledge edge orientations R. output: MPDAG G'
or FAIL.` The 4 hits of `error` are the noise terms ε of the linear SEM (guo.txt:580, 621, 690, 831).

**Critical distinction the framing depends on, now verified.** IDGraphs *does* add orientations —
guo.txt:332–333 branches on `G1 = MPDAG(G, {A1 → V1})` and `G2 = MPDAG(G, {A1 ← V1})`. But those are
hypothesised orientations of the **problematic undirected edges** — the part of the graph that is
*unknown* — enumerated **inside** a fixed K. They are not reversals of anything the expert asserted.
The premise "they enumerate within a fixed assumed-true K and never perturb it" is **exactly right**.

### Q2. Any radius, tolerance, or distance over knowledge sets?

**NO. Not one.** `radius` 0 · `tolerance` 0 · `ball` 0 · `distance` 0 · `neighbourhood` 0.
The 2 hits of `neighborhood` (guo.txt:667, 669) are the graph-local neighbourhood of node A used by
local IDA — "the local versions of IDA […] only require the neighborhood information of A instead of
the whole MPDAG" — a statement about *graph locality*, not a metric on knowledge sets.

### Q3. Do they weight by relevance to O*, or is O* absent?

**O* IS ABSENT.** `optimal adjustment` = **0**. They never define O* = pa(cn(X,Y)) \ forb(X,Y).
Efficiency is explicitly outsourced (guo.txt:507–509): "Some recently developed efficient estimators
**can be employed**, including the semiparametric efficient estimator of Rotnitzky and Smucler (2020)
and the efficient least-squares estimator of Guo and Perković (2020)."

⚠️ **But Corollary 8 is closer to this territory than the campaign recorded, and this is a real
finding** (guo.txt:541–543):

> "**Corollary 8.** Then there are no two MPDAGs in L that share the same adjustment set relative to
> (A, Y). Further, if |A| = |Y| = 1, then there exists an adjustment set relative to (A, Y) for each
> MPDAG in L."

So *"enumerate the class and obtain one distinct adjustment set per member"* **is published**, for
the single-treatment/single-outcome case that is exactly the pilot's setting. This does **not** fire
trigger (d) — it is a *minimality/distinctness* guarantee, with no weighting, no ranking, no
relevance notion, and no O*. But it constrains how the ICLR paper may describe Chris's contribution:
**"O* per member" is not a deliverable, it is Corollary 8 plus Henckel's O*.** The contribution must
be stated as **O*-relevance** — *which asserted statements move O***  — not as *obtaining* an
adjustment set per member.

### Q4. What EXACTLY do they enumerate, and what do they report?

They enumerate a **minimal set L of MPDAGs** such that the total effect in each is identified as a
**distinct functional of the observed distribution** (Theorem 3, Corollary 7: "no two graphs in L
share the same formula Eq. (2)"). The output is a **de-duplicated set of identified functionals /
adjustment sets** — a set of *possible true effects*.

It is **not an interval**: `interval` = **0** occurrences.

The entire contribution is *removing duplicates* from IDA's multiset (abstract, guo.txt:26–29):

> "This resolves an issue with existing methods, which often report possible total effects with
> duplicates, namely those that are **numerically distinct due to sampling variability but are in
> fact causally identical**."

**Note this sentence carefully.** Sampling variability appears only as a **nuisance to be removed**,
so that the reported set is causally clean. It is never a *scale against which knowledge error is
measured*. That is the precise gap leg (3) occupies.

### Q5. Is the surviving delta intact, narrowed, or gone?

**INTACT on all three legs, with one framing narrowing.** Detail in §4.

---

## 3. Verification of the campaign's six claims about Taeb §6

§6 spans taeb.txt:1219–1262 (§7 Discussion begins 1263); page header "24" at 1233.

| # | Campaign claim | Verdict | Evidence |
|---|---|---|---|
| 1 | §6 is a one-page illustrative example | ✅ **CONFIRMED** | 44 lines = ~12 lines prose + Figure 7 + caption, all on p.24 of 49 |
| 2 | on a 4-node ADMG | ✅ **CONFIRMED** | nodes v1,v2,v3,v4; "the model-oriented poset of (𝔊, M_ADMGc)"; `ADMG`×55 |
| 3 | validity-survival, not bias-magnitude | ✅ **CONFIRMED** | "the valid adjustment sets are **color-coded**"; "{v3,v4} is the most robust: it **remains valid** whenever the other sets do" |
| 4 | "bias" appears 0 times in 49 pp | ✅ **CONFIRMED** | `grep -oi bias taeb.txt \| wc -l` → **0**. Substring match, so "biased"/"unbiased" are also absent |
| 5 | "optimal adjustment" appears 0 times | ✅ **CONFIRMED** | **0** |
| 6 | "ball-of-radius-d + tolerance-radius", called "this sensitivity analysis" verbatim | ⚠️ **HALF-CONFIRMED — the vocabulary is the campaign's paraphrase, not Taeb's** | see below |

### On claim 6 — the campaign's own wording is not verbatim

**`radius` = 0 occurrences. `tolerance radius` = 0 occurrences. `ball` = 1**, and that single hit is
in the *introduction* as a **future** application (taeb.txt:63–64): "Other **potential applications**
include constructing **ball-like confidence regions** for a graph and conducting sensitivity analyses
with respect to potential mistakes in the causal graph."

What **is** verbatim (taeb.txt:1231–1232):

> "…{v3, v4} is the most robust: it remains valid whenever the other sets do and is able to
> **tolerate misspecifications as large as 𝑑_L = 5**. **This sensitivity analysis** depends on the
> fact that 𝑑_L is a metric so that the set of graphs **{G ∈ 𝔊 : 𝑑_L(G, G0) ≤ 𝑑}** can be iteratively
> enlarged."

So the **construction** is substantively a ball of radius d and a tolerance level — but
🔴 **do NOT put "tolerance radius" or "ball of radius d" in quotation marks attributed to Taeb.**
Only *"this sensitivity analysis"* and *"tolerate misspecifications as large as d_L = 5"* are quotable.
This campaign was damaged once by a citation that said something its source did not.

### The material finding the campaign missed: Taeb §6 holds background knowledge FIXED

The ball is over **supergraphs** — *edge addition*, not orientation reversal (Figure 7 caption,
taeb.txt:1254–1256): "the least element is the specified graph and each graph above represents **a
supergraph containing extra directed or bidirected edges**"; §6 body: "misspecifications, **in
particular directed or bidirected edges missing from the specification**".

And the background knowledge is the **fixed restriction defining the class** (taeb.txt:1228–1229):
"𝔊 contains the specified graph and its supergraphs **subject to the restriction that v3, v4
temporally precede v1**."

That restriction is a *tiering* — background knowledge — and it is held **true and constant** while
the graph varies. All 7 `background knowledge` hits in Taeb are definitional ("an MPDAG […] subject
to certain background knowledge", taeb.txt:278) or bibliographic (Meek 1995). `perturb` = 0.

⇒ **Taeb §6 points sensitivity analysis at GRAPH MISSPECIFICATION, with background knowledge as a
fixed constraint. It does not point it at the background knowledge itself.**

**Consequence for the struck sentence.** The campaign struck *"sensitivity analysis has never been
pointed at graphical background knowledge"* as "refuted verbatim". That is **slightly overstated** —
Taeb perturbs the graph, not K. But the sentence must **stay struck**, for two independent reasons:
(i) a reviewer who knows Taeb will read a ball-over-graphs and a ball-over-K as the same idea and
will not grant the distinction as a *novelty*; (ii) resting the headline on a distinction that fine
is fragile rhetoric, which is precisely why the refinement phase already moved novelty onto the three
concrete deltas. Keep it struck; the reason is now *precision*, not *refutation*.

### Bonus finding: the whole distance literature Taeb surveys is COUNT-based, never magnitude

taeb.txt:134–139 — SID (Peters & Bühlmann 2015) "**counts the number** of causal effects that would
be wrongly inferred by adjusting for the parents"; "Henckel, Würtzen and Weichwald (2024)
**generalized this approach into a broader class of adjustment identification distances**" (= gadjid
2402.08616); Wahl & Runge (2025) separation distances are "**counts** of conditional independence
statements".

**This independently corroborates leg (2) from a second angle**: the field's own survey of
graph-distance-for-adjustment measures binary identification failures and counts them. Not one
measures the *magnitude* of the resulting error. This is a ready-made related-work paragraph.

---

## 4. Verdict on the three delta legs

| Leg | Status | Evidence |
|---|---|---|
| **(1) PERTURBED OBJECT** — ball over Meek-consistent **reversals** of a stated, **possibly-false** K | ✅ **INTACT, and sharper than recorded** | Guo: K assumed true, `perturb`/`misspecif`/`wrong`/`incorrect` all 0; enumeration is over the *unknown* undirected edges inside a fixed K. Taeb §6: perturbs the **graph** by **supergraph/edge-addition**, holding background knowledge fixed. Neither reverses an asserted orientation. |
| **(2) CURRENCY** — **bias magnitude** of the effect estimate | ✅ **INTACT, doubly corroborated** | `bias` = 0 in both papers (0/24 pp and 0/49 pp, substring match). Guo reports distinct *identified functionals*; Taeb §6 reports *binary validity*, colour-coded. Taeb's own survey confirms SID / gadjid / Wahl-Runge are all **counts**. |
| **(3) INDEXING by sample size n** | ✅ **INTACT, one framing narrowing** | Guo `sample size` = 0, `interval` = 0. Taeb's single `sample size` hit is a BIC table caption (n=11, sample size 1000), unrelated. **BUT** see narrowing N1. |

### Four narrowings that must be carried into the draft

**N1 — the identification-vs-estimation *separation* is Guo & Perković's framing; only the
*indexing* is ours.** guo.txt:651–653, Discussion:

> "our result can be viewed as **separating two sources of uncertainty — identification and
> estimation** — the two crucial steps in causal inference."

They separate the two conceptually and hand estimation to Henckel/Rotnitzky-Smucler. They never
compare them on a common scale (`sample size` 0, `interval` 0, `bias` 0). So ρ*(n) — *the smallest
number of expert errors whose induced bias exceeds the estimator's own sampling error* — remains
unclaimed. But the paper must **cite this sentence** when introducing the two-scale framing rather
than presenting the separation as its own idea.

**N2 — Corollary 8 already gives one distinct adjustment set per enumerated member** (single X,
single Y). Chris's contribution must be positioned as **O*-relevance** (which asserted statements
move O*) and as the *efficient/optimal* O* specifically (Henckel 1907.02435), **never** as
"we obtain an adjustment set per member" — which is Corollary 8, published 2021.

**N3 — Taeb's "tolerance radius" / "ball of radius d" is the campaign's paraphrase.** `radius` = 0.
Quote only *"this sensitivity analysis"* and *"tolerate misspecifications as large as d_L = 5"*, and
describe the construction as `{G ∈ 𝔊 : d_L(G, G0) ≤ d}` — which is what the paper actually writes.

**N4 — Taeb Eq. (8) is the nearest published neighbour to leg (3), and the campaign never saw it.**
taeb.txt:544–555, §4 on why the triangle inequality matters:

> "the triangle inequality is essential for proving properties like the **consistency of an
> estimator**. It allows one to **bound the total error** between an estimate and the truth by
> introducing an intermediate 'population-level' proxy […]
> **d_L(Ĝ_n, G) ≤ d_L(Ĝ_n, G_∞) + d_L(G_∞, G)   (8)**
> [under-braced in the PDF as] **variability** + **approx. error**
> where Ĝ_n is the estimate, G_∞ is the intermediate population graph, and G is the true graph."

This is an **n-indexed decomposition of estimation variability against structural error** — the same
*shape* of idea as ρ*(n). It does **not** fire any trigger and leg (3) survives, because Eq. (8) is
measured in **graph-distance units** (d_L) about a **graph estimate** Ĝ_n from structure learning,
whereas ρ*(n) is measured in **units of the effect estimator's own standard error** about a
**causal-effect estimate** — and `bias` is still 0 in the whole document, so no effect estimate
appears in it at all. But the paper must **cite Eq. (8)** when it introduces the two-scale
comparison. Writing "nobody has compared structural error against sampling error" without this
citation would be the third avoidable overstatement on this line.

---

## 5. What would still change this

This read closes decisive step #3 (open 28 days). It does **not** close:

- **Fang, Zhao, Liu & He 2207.05067** (73 pp, JMLR 26) — never read by anyone. If its
  decomposed-MPDAG result contains a *perturbation* rather than an *equivalence* statement, leg (1)
  partly falls. **This is now the highest-value unread item on the line.**
- **gadjid 2402.08616** — absent from the 1,336-note corpus, and Taeb cites it as "a broader class of
  adjustment identification distances". Its *class* may be broader than the campaign assumed.
- The non-arXiv statistics venues (Biometrika, JRSS-B, Annals, JASA) — swept by **no** agent. The
  machinery ρ perturbs was published in JRSS-B.
- Neither paper has a literature note in the vault. **2511.10625 still has none**, so
  `check_claims.py` remains structurally blind to the nearest competitor.

Per CLAUDE.md, do **not** record a `verified:` date on any claim scoped wider than
"arXiv + ML venues, through 2026-08-14".
