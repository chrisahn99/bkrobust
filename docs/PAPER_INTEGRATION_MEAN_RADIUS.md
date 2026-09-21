# Paper integration — the mean bias radius as a companion to `r_val`

**Who this is for.** A session integrating the mean-bias-radius results into
`paper/sections/epsilon-radius.tex`. You do not need to re-run anything; every number below
traces to a committed file. §1 is the narrative and the one thing that must not be got wrong,
§4 is the table that goes in the **main body**, §6 is what must not be claimed.

---

## 1. The narrative, which is the whole point

**`r_mu` is an additional reporting column next to `r_val`. It is not a replacement, not a
correction, and not a competitor.** Every sentence written about it must preserve that.

The motivation is a reporting gap the paper already has:

- `r_val` takes **3 distinct values** over the 39 real instances with non-zero bias, and **35 of
  39 are simply `r_val = 1`**.
- The full `r_eps` profile over five tolerances yields **5 distinct signatures**.

So the certified apparatus sorts 39 real analyses into five boxes and puts 90% of them in one.
That is not a flaw in `r_val` — a validity radius answers *whether and when* the adjustment set
breaks, and on this corpus the answer really is "at the first revision, almost everywhere". It is
a limit on what that question can tell a practitioner.

`r_mu` answers a different question on the same geometry: *how much does it cost, on average, as
it breaks?* It yields **23 distinct signatures** over the same 39 instances, **29** with the
ceiling `B(Ĉ)`.

The sentence to aim for: **`r_val` says when the analysis breaks; `r_mu` says how badly, and the
two together separate analyses that `r_val` alone cannot tell apart.**

## 2. Files to read, in priority order

| Tier | File | Why |
|---|---|---|
| 1 | `docs/MEAN_TABLE_FULLSPACE.md` | The result. Table, coverage, caveats. |
| 1 | `results/mean_table_fullspace/instances.jsonl` | Per-instance rows; every number here is in it. |
| 2 | `results/mean_table_fullspace/summary.json` | Denominators and exclusions. |
| 2 | `docs/RESULT_EPSILON_ABOVE_ONE.md` | Why `r_eps` is inert: `beta_up` saturates at `r_val`. |
| 3 | `docs/RESULT_MEAN_SOUNDNESS_GATE.md` | Why the retraction up-set is not used for the mean. |
| 3 | `docs/MEAN_REPORTING.md` | Superseded for reporting; retains methodology and negative results. |

## 3. The findings, with their evidence

### 3.1 The statistic, and the space it lives on

```
mu(d) = mean { B(G) : G in G_Chat, d(G0, G) = d }        r_mu(eps) = min { d : mu(d) > eps }
```

`G_Chat` is the paper's own space — every valid MPDAG refining `Chat`, ordered by model
inclusion, distance the BFS hop count on the covering graph. **This is the space `r_val` is
defined on.** The retraction up-set is *not* used: it is licensed for the worst case by the
proved chain (Anti-Exchange ⟹ Lemma R ⟹ Property S ⟹ Lemma L ⟹ Conjecture 2), verified here
directly (948/948 first crossings agreed), but **no such equivalence holds for a mean**.
→ `results/mean_table_fullspace/`

### 3.2 `r_mu` is an onset, not a threshold

`mu` is **not monotone** on this space: far shells contain heavily oriented states that pin the
effect down and pull the average back. `mediator I→Y` rises to 4.57 at `d = 4` and falls to 1.00
by `d = 7`. So `r_mu(eps)` is the *smallest perturbation size at which expected bias exceeds
`eps`*, not a bound below which you are safe. **14 of 213 finite cells on 7 of 39 instances** have
a crossing that does not persist. This must be stated wherever `r_mu` is defined.

### 3.3 Where the graded behaviour is prominent, and why

Gradedness (distinct finite radii per instance) is driven by the **diameter `D`** of the
perturbation space, not by `|K_G0|`: Pearson **+0.729** against `D`, **+0.556** against `|K_G0|`.

| diameter band | instances | mean distinct radii |
|---|--:|--:|
| `D` ≤ 3 | 9 | **1.00** |
| `D` 4–6 | 14 | 2.29 |
| `D` 7–9 | 9 | 2.78 |
| `D` ≥ 10 | 7 | **4.57** |

Every flat instance has `D ≤ 4` except `mediator I→Y`, whose bias saturates at the first shell.
The interpretation is clean and should go in the prose: **the metric is informative exactly when
the data left enough structure undetermined for the analyst's state to have somewhere to move.**
A CPDAG that is nearly oriented has a shallow space, and nothing graded can happen in it.

**The prominent datasets are `hepar2`, `barley`, `Sebastiani_2005`, `Polzer_2012`, `Schipf_2010`,
`water` and `alarm`.** The flat ones are `Acid_1996`, `Thoemmes_2013` and `Didelez_2010`
(`D` = 2, 2 and 4).

### 3.4 Coverage: 48 of 66, and the exclusions are not random

The space is `3^k` in the `k` undirected edges of `Chat` with a quadratic covering step.
**48 of the 66 informative instances were enumerated exhaustively**; 9 of those have
`r_val = UNREACHED` and zero bias everywhere, leaving **39 in the table**. **18 instances are
excluded** because their space could not be built: `Kampen_2014`, `andes`, `child`, `ecoli70`,
`magic-irri`, `magic-niab`, `paths`, `win95pts`. Nothing was sampled or partially enumerated.

The exclusions are the networks with the **most** undirected edges — which by §3.3 are the ones
where gradedness would be most pronounced. **Say this.** The table under-represents its own best
case, and a reader who notices the missing networks should find the admission already there.

## 4. The table for the MAIN BODY

Restricted to the instances where the graded behaviour is prominent: ≥ 4 distinct finite radii,
`mu` monotone, no `UNREACHED` cell. The full 39-row table goes to the **appendix**.

| Network | Query | `r_val` | `D` | 2% | 10% | 25% | 40% | 60% | 90% | `B(Ĉ)` |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| `hepar2` | `RHepatitis → THepatitis` | 1 | 11 | 1 | 2 | 3 | 4 | 6 | 9 | 1.00 |
| `barley` | `dg25 → s2225` | 3 | 13 | 4 | 5 | 6 | 7 | 8 | 9 | 1.52 |
| `Sebastiani_2005` | `ANXA2.5 → ANXA2.11` | 1 | 11 | 1 | 1 | 2 | 4 | 5 | 8 | 1.00 |
| `Schipf_2010` | `A → U` | 1 | 8 | 1 | 1 | 3 | 4 | 5 | 7 | 1.00 |
| `barley` | `frspdag → dg25` | 2 | 13 | 2 | 3 | 5 | 5 | 6 | 8 | 1.00 |
| `Polzer_2012` | `Caries → Diabetes` | 1 | 12 | 1 | 1 | 2 | 3 | 4 | 5 | 1.00 |
| `alarm` | `ANAPHYLAXIS → BP` | 1 | 5 | 1 | 1 | 2 | 3 | 3 | 4 | 1.00 |
| `alarm` | `ANAPHYLAXIS → HRBP` | 1 | 5 | 1 | 1 | 2 | 3 | 3 | 4 | 1.00 |

LaTeX, ready to paste (the two `alarm` rows are genuinely distinct queries with identical
profiles; cut one for space if needed and say so):

```latex
\begin{tabular}{llrrrrrrrr}
\toprule
Network & Query & $r_{\mathrm{val}}$ & $D$ & \multicolumn{6}{c}{$r_\mu$ at tolerance $\varepsilon$} \\
\cmidrule(lr){5-10}
 & & & & $2\%$ & $10\%$ & $25\%$ & $40\%$ & $60\%$ & $90\%$ \\
\midrule
\texttt{hepar2} & \texttt{RHepatitis}$\to$\texttt{THepatitis} & 1 & 11 & 1 & 2 & 3 & 4 & 6 & 9 \\
\texttt{barley} & \texttt{dg25}$\to$\texttt{s2225} & 3 & 13 & 4 & 5 & 6 & 7 & 8 & 9 \\
\texttt{Sebastiani\_2005} & \texttt{ANXA2.5}$\to$\texttt{ANXA2.11} & 1 & 11 & 1 & 1 & 2 & 4 & 5 & 8 \\
\texttt{Schipf\_2010} & \texttt{A}$\to$\texttt{U} & 1 & 8 & 1 & 1 & 3 & 4 & 5 & 7 \\
\texttt{barley} & \texttt{frspdag}$\to$\texttt{dg25} & 2 & 13 & 2 & 3 & 5 & 5 & 6 & 8 \\
\texttt{Polzer\_2012} & \texttt{Caries}$\to$\texttt{Diabetes} & 1 & 12 & 1 & 1 & 2 & 3 & 4 & 5 \\
\texttt{alarm} & \texttt{ANAPHYLAXIS}$\to$\texttt{BP} & 1 & 5 & 1 & 1 & 2 & 3 & 3 & 4 \\
\texttt{alarm} & \texttt{ANAPHYLAXIS}$\to$\texttt{HRBP} & 1 & 5 & 1 & 1 & 2 & 3 & 3 & 4 \\
\bottomrule
\end{tabular}
```

Caption must carry: `D` is the diameter of the perturbation space; `r_mu` is SEM-conditional
while `r_val` is not; `r_mu` is an onset (§3.2); and the rows are the prominent subset, with all
39 in the appendix.

### The two contrast examples for the prose

These are the strongest argument for the column and belong in the text, not the table:

- `Acid_1996 x1→x10` and `x1→x13` are **identical under both `r_val` and `r_eps`** (both
  `r_val = 1`). Their ceilings are **1.93** and **0.126**: the first can be more than wiped out,
  the second can never move by more than 13% of its own size, however wrong the analyst is.
- `Didelez_2010 Age→HRT` and `Shrier_2008 Coach→IntraGameProprioception` are also both
  `r_val = 1`. The first has expected bias **3.45×** the effect at distance 1; the second never
  exceeds **0.30** at any distance.

Same certificate, opposite practical advice, in both pairs.

## 5. Concrete edits

| Where | What |
|---|---|
| `paper/sections/epsilon-radius.tex`, `sec:epsilon-eval` | Add a `\paragraph{On real graphs.}` carrying §1 and §3.2–3.4, then the §4 table. This section currently contains **no** tabular environment (verified) — introducing one is a deliberate change; keep it to this single table. |
| Appendix | The full 39-row table from `docs/MEAN_TABLE_FULLSPACE.md` §4, plus the coverage statement from §3.4. |
| `docs/PAPER_NARRATIVE.md` | Record that `r_mu` is a reporting companion, never a replacement for `r_val`, and that it is SEM-conditional. |

## 6. What must NOT be claimed

- **Never** present `r_mu` as replacing, improving on, or correcting `r_val`. It answers a
  different question on the same space.
- **Never** state a guarantee against `r_mu`. It is an average and **bounds nothing**;
  `beta_up` remains the quantity every guarantee is stated against, and `mu(d) ≤ beta_up(d)`.
- **Never** call `r_mu` a threshold below which the analysis is safe. It is an onset (§3.2).
- **Never** average `UNREACHED`. It is the sentinel `-1` meaning no distance crosses that
  tolerance — a status, not a large radius. On this statistic it means the tolerance exceeds the
  instance's own ceiling `B(Ĉ)`.
- **Never** paraphrase `r = 1` as "tolerates one error". `r` is the distance to the nearest
  failure; the number of safe moves is `r - 1`, so `r = 1` is **zero** safe revisions.
- **Never** report a per-network median of `r_mu`. Per query only — the instance set with a
  finite radius shrinks as `eps` rises, and a median over a moving denominator has already put
  one false claim in this paper.
- **Never** omit the SEM caveat on an `r_mu` or `B(Ĉ)` number. The corpus ships graph structure
  only; these are computed against a seeded linear-Gaussian SEM on the ground-truth DAG.
  `r_val` is purely graph-theoretic and carries no such assumption.
- **Never** quietly drop the 18 excluded instances. State the denominator, name the networks, and
  note that they are the high-diameter ones (§3.4).
- **Never** claim the result transfers beyond `D_LLM` or beyond this corpus. Nine other
  elicitation conditions sit unused in `results/elicit/knowledge.json`.
- **Not pre-registered.** `docs/EVALUATION_PLAN.md` mentions neither `r_eps` nor `mu`.

## 7. Open decisions, deliberately not taken

1. **Whether the appendix table shows the `*` non-persistence marks.** They are in the data. They
   are honest but visually noisy across 39 rows.
2. **Whether to report `mu`'s peak** (`mu_peak`, `mu_peak_at` in the JSONL) alongside the onset.
   It is well defined even where the onset is not persistent, and it is arguably the more stable
   summary — but it is a second number per row.
3. **Whether to raise coverage** by pushing `--max-undirected` past 10 with a long budget. The
   gain would be `andes` and `magic-niab` (k=10, 5 instances); `Kampen_2014` upward is out of
   reach. Worth it only if a reviewer presses on coverage.
