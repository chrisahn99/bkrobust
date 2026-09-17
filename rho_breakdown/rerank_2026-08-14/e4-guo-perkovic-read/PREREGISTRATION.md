# E4 PRE-REGISTRATION — written BEFORE fetching either paper

Written 2026-08-14, before any text of arXiv:2010.08611 or arXiv:2511.10625 was read
in this session. Timestamp of this file precedes the first curl call.

## The premise under test

The surviving framing of idea #1 (rho-breakdown radius) rests on one sentence:

> "Guo & Perkovic (2010.08611) enumerate within a FIXED, ASSUMED-TRUE background
> knowledge set K and never perturb it."

This premise has only ever been read at ABSTRACT level. It is decisive step #3 on
wiki/projects/rho-breakdown-adjustment.md and has been open for 28 days.

## The three delta legs the read can break

| Leg | Statement | Broken if 2010.08611 ... |
|-----|-----------|--------------------------|
| (1) PERTURBED OBJECT | a ball over Meek-consistent REVERSALS of statements in a stated, possibly-FALSE K | ... perturbs, relaxes, doubts, or defines a neighbourhood over K itself |
| (2) CURRENCY | BIAS MAGNITUDE of the effect estimate under wrong knowledge | ... reports a bias magnitude attributable to knowledge error (not merely a range of correct effects) |
| (3) INDEXING | rho* indexed by SAMPLE SIZE n — smallest number of expert errors whose induced bias exceeds sampling error | ... indexes anything by n against estimation noise |

## PRE-REGISTERED DECISION RULE

Let me define, before reading:

- **REFUTED** — declare the delta GONE — iff the paper does **ANY** of:
  - (a) defines a distance, radius, tolerance, ball, or neighbourhood over
    background-knowledge sets (not over graphs, not over DAGs in a MEC), **or**
  - (b) explicitly considers K being WRONG / MISSPECIFIED / possibly-false anywhere
    beyond a one-line "future work" mention, and computes a downstream consequence
    of that wrongness, **or**
  - (c) reports a bias magnitude of an effect estimate attributable to a wrong
    orientation constraint, **or**
  - (d) weights / ranks / filters enumerated members by relevance to the optimal
    adjustment set O*.

  Any ONE of (a)-(d) kills at least one of the three legs. If (a) OR (b) fires,
  leg (1) — the load-bearing one — is gone and the honest verdict is REFUTED.

- **SUPPORTED** — delta intact — iff **NONE** of (a)-(d) fires, i.e. K is assumed
  true throughout, the enumeration is over DAGs/effects consistent with a fixed K,
  no knowledge-space metric exists, O* is absent, and the output is a set/interval
  of possible TRUE effects rather than a bias under error.

- **INCONCLUSIVE** — iff the full text cannot be obtained, or if the paper's
  treatment is genuinely ambiguous on (a)-(d) after a full read.

## Secondary target: arXiv:2511.10625 section 6 (Taeb/Guo/Henckel)

Campaign's characterisation to VERIFY (not to trust):
1. Section 6 is a **one-page illustrative example**.
2. On a **4-node ADMG**.
3. Measures **validity-survival**, NOT bias-magnitude.
4. The word **"bias" appears 0 times** in the 49 pages.
5. The phrase **"optimal adjustment" appears 0 times**.
6. It contains a **ball-of-radius-d + tolerance-radius** construction and calls it
   "this sensitivity analysis" verbatim.

Each of the six is a checkable string/structure fact. I will report each as
CONFIRMED / REFUTED / UNVERIFIABLE with the verbatim evidence.

If (4) or (5) is REFUTED — i.e. "bias" or "optimal adjustment" DO appear — leg (2)
of the delta (bias-magnitude currency) is damaged and the overall verdict moves
toward REFUTED regardless of what 2010.08611 says.

## Evidence standard (campaign rule, enforced)

- pdftotext or the actual HTML text. **No WebFetch summarizer snippet may carry a
  load-bearing quote.** Idea #10 died on a hallucinated Melnychuk sentence; a
  summarizer separately INVERTED the Rotnitzky-Smucler theorem.
- Every quote below must be reproducible by grepping the extracted text file that
  is saved next to this one.
- If a fetch fails, say NOT_RUN. Do not reason from the abstract and present it as
  a full read — that is exactly the failure this experiment exists to correct.
