# The mean bias radius, on the paper's own geometry

The reporting table this session was asked for: the mean bias radius `r_mu`, computed on the
perturbation space the paper's geometry is defined on, reported per query alongside `r_val`.

Every number traces to `results/mean_table_fullspace/`, produced by
`experiments/mean_table_fullspace.py`. Condition `D_LLM`, seed 20260917.

---

## 1. The space, stated first because it is the point

The geometry is

```
G_Chat = { G a valid MPDAG : [G] ⊆ [Chat] },   ordered by model inclusion,
distance = BFS hop count on the covering graph
```

This is the space `r_val` is defined on, and the space `docs/R_EPSILON_THEORY.md` Definition 5
defines `r_eps` on. **Everything below is computed on it.**

The **retraction up-set** `↑G0 = { Meek(Chat, K_G0 \ S) }` is a strictly smaller object. It
entered the project as a computational device, and for the **worst case** it is licensed by a
theorem chain that is now closed — Anti-Exchange is PROVED unconditionally, hence Lemma R,
Property S, Lemma L, and Conjecture 2 as a PROVED THEOREM — so the first crossing of the
retraction profile equals the first crossing of the full-space profile. `r_val` and the
max-based `r_eps` computed over retractions are therefore the *same numbers* the full space
gives. That was also checked directly rather than assumed: over 79 instances the shell maxima
matched the ball maxima on every shell, and the first crossings agreed on all 948 cells
(`results/mean_fullspace/`).

**No such equivalence exists for a mean**, so the up-set is not used here at all. An earlier
version of this statistic averaged over the retraction shell; that is a different quantity, the
two disagree in both directions, and it is superseded for reporting by this table
(`docs/MEAN_REPORTING.md` retains the methodology and the negative results).

The distance column `D` in the table makes the difference concrete. `barley` has `|K_G0| = 7`
but a diameter of 13: the full space reaches states the analyst could hold that are *more*
oriented than `G0`, not only less. Those states are absent from the up-set by construction, and
they are exactly the ones that make the mean behave differently.

## 2. The statistic

```
mu(d) = mean { B(G) : G ∈ G_Chat, d(G0, G) = d }
```

the average bias over the **sphere** at perturbation distance `d` from the analyst's state, and

```
r_mu(eps) = min { d : mu(d) > eps }
```

`B` is unchanged: at each state it is already the worst case over every DAG that state leaves
open. Only the aggregation across a shell is new.

**`mu` is not monotone on this space, and that is a property of the geometry, not a defect of
the statistic.** Far shells contain heavily oriented states that pin the effect down and so carry
little bias, which pulls the average back. `mediator I→Y` rises to 4.57 at `d = 4` and falls to
1.00 by `d = 7`. So `r_mu` is an **onset** — the smallest perturbation size at which the
expected bias exceeds the tolerance — and *not* a threshold beyond which it stays exceeded.
Cells where the crossing does not persist to the outermost shell are marked `*`: 14 of 213
finite cells, on 7 of 39 instances.

Two further things that must travel with any reported number. `mu` is an **average-case**
quantity and bounds nothing; `beta_up` remains what every guarantee is stated against. And every
bias here is **SEM-conditional** — the corpus ships graph structure only, so `B` is evaluated
against a seeded linear-Gaussian SEM attached to the ground-truth DAG. `r_val` carries no such
assumption.

## 3. Coverage, and why it is not the whole corpus

Staying in the real geometry has a price, and it is the reason the retraction acceleration was
invented in the first place. The space is `3^k` in the `k` undirected edges of `Chat`, and
building its covering relation is quadratic in the result. Measured:

| k | network | \|space\| | build |
|--:|---|--:|--:|
| 6 | `water` | 100 | 0.3 s |
| 7 | `Sebastiani_2005` | 891 | 3.6 s |
| 9 | `barley` | 1,980 | 48 s |
| 9 | `hepar2` | 2,997 | 174 s |
| 10 | `magic-niab` | — | exceeded 900 s |
| 11 | `Kampen_2014` | — | `3^k` = 177,147 |
| 25 | `ecoli70` | — | `3^k` = 8.5 × 10¹¹ |

**48 of the 66 informative instances were enumerated exhaustively on the full space.** The other
**18 are excluded** because their space could not be built: `Kampen_2014`, `andes`, `child`,
`ecoli70`, `magic-irri`, `magic-niab`, `paths`, `win95pts`. Nothing was sampled, estimated or
partially enumerated — an instance is either exact or absent.

Of the 48, **9 have `r_val = UNREACHED`** and therefore zero bias everywhere (Theorem A), so
their `r_mu` is `UNREACHED` at every tolerance and they are not shown. The table is the
remaining **39**.

For contrast, the retraction up-set covers all 66 at 49,427 subsets in two minutes. That gap is
the honest cost of the geometry, and it is worth stating in the paper rather than hiding: the
worst case gets the full corpus for free because a theorem licenses the shortcut; the mean does
not.

## 4. The table

`r_val` is graph-theoretic; every `r_mu` column is SEM-conditional. `D` is the diameter — the
largest BFS distance from `G0` in the full space. `UNR` is the `UNREACHED` sentinel: no distance
crosses that tolerance. It is a status, never a large radius, and is never averaged. `*` marks a
crossing that does not persist to the outermost shell. Recall `r = 1` means **zero** safe moves.

| Network | Query | \|K_G0\| | `r_val` | D | 2% | 10% | 25% | 40% | 60% | 90% | peak (at d) | `B(Ĉ)` |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| `barley` | `dg25 → aks_m2` | 7 | 3 | 13 | 5 | 9 | UNR | UNR | UNR | UNR | 0.18 (d=13) | 0.184 |
| `barley` | `dg25 → s2225` | 7 | 3 | 13 | 4 | 5 | 6 | 7 | 8 | 9 | 1.52 (d=12) | 1.525 |
| `water` | `CKNI_12_15 → CBODN_12_45` | 6 | 3 | 12 | 3 | 3 | 3 | 5 | 5 | 7 | 2.25 (d=11) | 2.245 |
| `barley` | `frspdag → dg25` | 7 | 2 | 13 | 2 | 3 | 5 | 5 | 6 | 8 | 1.00 (d=10) | 1.000 |
| `Acid_1996` | `x1 → x10` | 1 | 1 | 2 | 1 | 1 | 1 | 1 | 1 | 1 | 1.93 (d=1) | 1.927 |
| `Acid_1996` | `x1 → x11` | 1 | 1 | 2 | 1 | 1 | 1 | 1 | 1 | 1 | 1.63 (d=1) | 1.629 |
| `Acid_1996` | `x1 → x12` | 1 | 1 | 2 | 1 | 1 | 1 | UNR | UNR | UNR | 0.38 (d=1) | 0.384 |
| `Acid_1996` | `x1 → x13` | 1 | 1 | 2 | 1 | 1 | UNR | UNR | UNR | UNR | 0.13 (d=1) | 0.126 |
| `Didelez_2010` | `Age → HRT` | 3 | 1 | 4 | 1 | 1 | 1 | 1 | 1 | 1 | 3.45 (d=1) | 3.453 |
| `Didelez_2010` | `Age → Occ` | 3 | 1 | 4 | 1 | 1 | 1 | 1 | 1 | 1 | 1.00 (d=1) | 1.000 |
| `Didelez_2010` | `Age → S` | 3 | 1 | 4 | 1 | 1 | 1 | 1 | 1 | 1 | 1.22 (d=1) | 1.223 |
| `Didelez_2010` | `Age → Smo` | 3 | 1 | 4 | 1 | 1 | 1 | 1 | 1 | 1 | 1.00 (d=1) | 1.000 |
| `Didelez_2010` | `Age → TCI` | 3 | 1 | 4 | 1 | 1 | 1 | 1 | 1 | 1 | 1.47 (d=1) | 1.473 |
| `Polzer_2012` | `Caries → Diabetes` | 5 | 1 | 12 | 1 | 1 | 2 | 3 | 4 | 5 | 1.00 (d=8) | 1.000 |
| `Schipf_2010` | `A → PA` | 6 | 1 | 8 | 1 | 1\* | 1\* | 2\* | UNR | UNR | 0.56 (d=5) | 1.000 |
| `Schipf_2010` | `A → S` | 6 | 1 | 8 | 1\* | 1\* | 1\* | 3\* | UNR | UNR | 0.53 (d=5) | 1.000 |
| `Schipf_2010` | `A → T2DM` | 6 | 1 | 8 | 1 | 1\* | 2\* | 3\* | 5\* | UNR | 0.64 (d=5) | 0.791 |
| `Schipf_2010` | `A → TT` | 6 | 1 | 8 | 1 | 1 | 1 | 1 | 1 | 5 | 1.29 (d=7) | 1.422 |
| `Schipf_2010` | `A → U` | 6 | 1 | 8 | 1 | 1 | 3 | 4 | 5 | 7 | 1.00 (d=7) | 1.000 |
| `Sebastiani_2005` | `ANXA2.5 → ANXA2.11` | 4 | 1 | 11 | 1 | 1 | 2 | 4 | 5 | 8 | 1.00 (d=10) | 1.000 |
| `Shrier_2008` | `Coach → IntraGameProprioception` | 3 | 1 | 6 | 1 | 1 | 5 | UNR | UNR | UNR | 0.30 (d=5) | 0.296 |
| `Shrier_2008` | `Coach → TeamMotivation` | 3 | 1 | 6 | 1 | 1 | 1 | 1 | 3 | 5 | 1.00 (d=5) | 1.000 |
| `Thoemmes_2013` | `e0 → s1` | 1 | 1 | 2 | 1 | 1 | 1 | 1 | 1 | 1 | 1.00 (d=1) | 1.000 |
| `Thoemmes_2013` | `e0 → s2` | 1 | 1 | 2 | 1 | 1 | 1 | 1 | 1 | 1 | 1.00 (d=1) | 1.000 |
| `Thoemmes_2013` | `e0 → s3` | 1 | 1 | 2 | 1 | 1 | 1 | 1 | 1 | 1 | 1.00 (d=1) | 1.000 |
| `Thoemmes_2013` | `e0 → x` | 1 | 1 | 2 | 1 | 1 | 1 | 1 | 1 | 1 | 1.00 (d=1) | 1.000 |
| `Thoemmes_2013` | `e0 → y` | 1 | 1 | 2 | 1 | 1 | 1 | 1 | 1 | 1 | 1.00 (d=1) | 1.000 |
| `alarm` | `ANAPHYLAXIS → BP` | 1 | 1 | 5 | 1 | 1 | 2 | 3 | 3 | 4 | 1.00 (d=4) | 1.000 |
| `alarm` | `ANAPHYLAXIS → HRBP` | 1 | 1 | 5 | 1 | 1 | 2 | 3 | 3 | 4 | 1.00 (d=4) | 1.000 |
| `asia` | `asia → dysp` | 3 | 1 | 4 | 1 | 1 | 1 | 2 | 3 | 3 | 1.00 (d=3) | 1.000 |
| `asia` | `asia → either` | 3 | 1 | 4 | 1 | 1 | 1 | 2 | 3 | 3 | 1.00 (d=3) | 1.000 |
| `asia` | `asia → tub` | 3 | 1 | 4 | 1 | 1 | 1 | 2 | 3 | 3 | 1.00 (d=3) | 1.000 |
| `asia` | `asia → xray` | 3 | 1 | 4 | 1 | 1 | 1 | 2 | 3 | 3 | 1.00 (d=3) | 1.000 |
| `asia` | `bronc → dysp` | 3 | 1 | 4 | 1 | 1 | 1 | 1 | 2 | UNR | 0.84 (d=4) | 1.256 |
| `hepar2` | `RHepatitis → THepatitis` | 5 | 1 | 11 | 1 | 2 | 3 | 4 | 6 | 9 | 1.00 (d=10) | 1.000 |
| `mediator` | `I → X` | 4 | 1 | 7 | 1 | 1 | 1 | 2 | 2 | 4 | 1.00 (d=6) | 1.000 |
| `mediator` | `I → Y` | 4 | 1 | 7 | 1 | 1 | 1 | 1 | 1 | 1 | 4.57 (d=4) | 4.967 |
| `mediator` | `I → Z` | 4 | 1 | 7 | 1 | 1 | 2\* | 2\* | 3\* | UNR | 0.75 (d=4) | 1.000 |
| `mediator` | `X → Y` | 4 | 1 | 7 | 1 | 1 | 2 | 2 | 2 | 6 | 1.00 (d=6) | 1.000 |

## 5. What it adds to `r_val`

Over these 39 instances:

| Statistic | distinct signatures |
|---|--:|
| `r_val` alone | **3** — and 35 of 39 are simply `r_val = 1` |
| `r_eps` profile (max, five tolerances) | 5 |
| **`r_mu` profile (six tolerances)** | **23** |
| **`r_mu` + `B(Ĉ)`** | **29** |

`r_val` sorts 39 real analyses into three boxes and puts 90% of them in one. The mean radius
sorts the same 39 into 23. That is the reporting gap this was meant to close.

The concrete reading. `Acid_1996 x1→x10` and `x1→x13` are identical under `r_val` (both 1) and
identical under `r_eps`. Their mean radii separate at 25%: `x1→x10` is still crossing at every
tolerance up to 90% because its ceiling is 1.93, while `x1→x13` has a ceiling of 0.126 — no
perturbation anywhere can move it by more than 13% of its own size. One is fragile in a way the
validity radius cannot express; the other is safe in a way it cannot express either.

Similarly `Didelez_2010 Age→HRT` and `Shrier_2008 Coach→IntraGameProprioception` both have
`r_val = 1`. The first has expected bias already 3.45x the effect at distance 1; the second never
exceeds 30% at any distance. Same certificate, opposite practical advice.

## 6. Limitations

- **Coverage is 48 of 66**, and the 18 exclusions are not random: they are the networks with the
  largest CPDAG undirected-edge counts, which correlates with the queries having more room to
  perturb. The table under-represents exactly the instances where the geometry is richest.
- **`mu` is not monotone**, so `r_mu` is an onset, not a threshold. 7 of 39 instances and 14 of
  213 finite cells are affected. Reporting the peak alongside the onset is what keeps this honest.
- **Average case under a uniform prior** over the sphere. It bounds nothing.
- **SEM-conditional.** A different seed gives a different, equally valid, incomparable panel.
- **One condition** (`D_LLM`), and claims are scoped to this corpus.
- **Not pre-registered.**

## 7. Reproducing

```bash
PYTHONPATH=src .venv/bin/python experiments/mean_table_fullspace.py \
  --budget 900 --max-undirected 10 --out results/mean_table_fullspace
```

35 minutes. `--max-undirected` skips networks whose `3^k` makes the build hopeless rather than
burning the full budget on each; they are recorded as excluded with their `k`.
