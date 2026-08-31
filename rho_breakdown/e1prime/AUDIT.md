# E1′ DESIGN AUDIT — written with `PREREG.md`, before implementation

Adversarial read of my own design. E1 was refuted 3/3 on *interpretation* while every number
reproduced; the cheapest way to repeat that is to skip this file.

---

## A1. The location bar can still be free — and the check is a number, not an argument

`n50 ≤ 20000` is a bar on *where* the transition happens. It is free if the transition has already
happened at the bottom of the grid.

> **Mitigation, binding:** report `censored(n = 50)` and `censored(n = 200)` explicitly. If
> `censored(200) ≤ 0.50`, the primary rule was **cleared without the n-axis doing any work**, and the
> honest finding is *"ρ\* was never degenerate; `sign` was the wrong criterion"* — a **different and
> stronger** paper that also makes E1's |K|=8 ensemble moot. That must be written that way, not
> smoothed into "SUPPORTED".

## A2. 🔴 "Analyst-facing" is scoped to *given the CPDAG*, and the CPDAG here is the oracle's

`C = dag_to_cpdag(D)` (`run_linear.py:80`). The analyst in the real world runs PC/GES on the same
`n` rows and gets `Ĉ ≠ C`. So the honest claim is:

> ρ\*_se is computable **from (C, K, Σ̂, n)** — no `tau`, no `D*`. It is *not* demonstrated to be
> stable under an **estimated** CPDAG.

**ARM C (Ĉ from discovery on the same n rows) is NOT RUN.** It is the first thing a reviewer asks
and the honest answer is "not tested", not "robust". Cost estimate: needs a PC implementation
(`causal-learn` absent from every local env per `RANKING.md`), ≈ 1 day + a package install. Written
here so the gap cannot read as coverage.

## A3. The yardstick is conservative — and it is conservative AGAINST us, which is the right sign

`d(ρ)` compares `est_m` to `est0` computed on the **same data**, so `Var(est_m − est0) < Var(est0)`.
Using `z·SE_n(est0)` as the threshold is therefore a **stricter** bar than the variance of the
difference would justify: it **under-detects** breakdowns and inflates censoring. Every censoring
number in ARM 1 is an upper bound. State it; do not "correct" it — a diagnostic that reports the
analyst's *own published standard error* as the yardstick is the one an analyst can actually act on.

## A4. Homoskedastic OLS SE is *exact* here, and that is a scope limit, not a strength

Everything is jointly Gaussian, so `Y | X, O` is linear with **constant** variance for *any* `O`,
valid or not. The classical SE is exact asymptotically even under a wrong adjustment set. Outside
linear-Gaussian this breaks and a sandwich SE is required. The licensed box in `PREREG.md` §8 says so.

## A5. ρ\* is a sensitivity quantity, not a test — no p-value may ever be attached

`d(ρ)` is a max over up to `2^|K| − 1 = 15` members compared against one `z·SE`. Read as a test this
is uncorrected multiplicity. It is **not a test**: it is a worst-case-over-the-ball sensitivity
radius, in the Rosenbaum-Γ tradition. **No p-value, no "significance", ever.** (Vault precedent: the
`rho-breakdown` page already forbids stating this as a continuity bound; this is the same discipline
one level down.)

## A6. ARM 2 coverage at ρ ≥ r is near-definitional — the deliverable is the WIDTH

If the ball around `K_asserted` reaches `K_true`, it contains a valid adjustment set, so the robust
interval covers `tau` almost by construction. Declared in `PREREG.md` §5.4. The two non-trivial
numbers are (i) the **width ratio** — a vacuous interval also covers — and (ii) coverage at
**ρ < r**, where the ball misses the truth. The registered rule is written on the width for exactly
this reason.

## A7. Non-amenable members must not silently shrink the ball

A perturbation that destroys amenability is **loud**: the analyst sees identification fail and does
not report a number. Those members are excluded from `d(ρ)` (there is no estimate to compare) but
are **counted and reported separately** as the `ident` event. Letting them vanish would make the
ball look quieter than it is.

## A8. Selection on `|tau| ≥ 1e-3` interacts with `se` in a way it did not with `sign`

Small-|τ| SCMs are exactly the ones with small t-statistics, i.e. where sampling error dominates.
The pilot's filter therefore biases the `se` arm toward **large-effect** problems, where breakdown
is *harder* (again conservative). Kept for comparability with E1; **sensitivity at 1e-6 reported**.

## A9. What must be reported even if it embarrasses the design

- `censored(n=50)` and `censored(n=200)` (A1).
- The **Spearman against the t-statistic** (`PREREG.md` §3c) at whatever value it takes. If it
  exceeds 0.90 the leg collapses to a sentence and the paper says so.
- The `large`-ensemble number against the registered prediction in `PREREG.md` §4, including the two
  named surprise directions.
- Every gate G1–G7, pass or fail, before any verdict.

## A10. 🔴 Vault-level gap this run does not close, and must not paper over

`wiki/activities/active-claims.md` holds **no** entry for ρ-breakdown — `RANKING.md:352` flags it as
a hard 18/09 item, and `grep -rl 'relevance:.*rho-breakdown' wiki/literature/` returns **zero notes**.
So `check_claims.py` is structurally blind to this paper's nearest competitors
(Taeb–Guo–Henckel 2511.10625, gadjid 2402.08616, Fang 2207.05067, CausalGuard 2605.21928 — the last
is *in* the corpus, `status: to-read` since 07/07, mis-routed).

**This run produces evidence for a claim that is not registered.** Under CLAUDE.md Operating Rule #1
the register entry must exist before the claim is written into a paper. Registering it requires Zé's
own claim wording (three ids: measurement / theorem / conformal §5, per the `kar-vs-kcm-dominance`
post-mortem) — **proposed, not committed by an agent.**
