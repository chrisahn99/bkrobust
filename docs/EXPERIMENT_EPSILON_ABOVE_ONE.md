# Experiment brief — the ε-radius above a tolerance of 1

**Who this is for.** A session that will re-run `experiments/final_table.py` at bias
tolerances ε ∈ {1, 2, 5, 10, 20} and report the result per query. You do not need to read
the paper. §2 says why this cannot be answered from committed data, §3 bounds the work,
§4 is what to run, §6 is what must not be claimed.

---

## 1. The question, and why it is a fair one

`r_eps` is defined in the paper as

```
r_eps = min{ d(G0, G) : B(G) > ε } = min{ d : β↑(d) > ε }
```

where `β↑(d)` is the worst-case bias over the shell at distance `d`, measured relative to
the effect's own magnitude. **ε is a threshold on that staircase and is not bounded above.**
ε = 1.00 means "the estimate is off by 100% of its own magnitude" — a large error, but not a
maximum. ε = 20 is a perfectly well-posed question: how far must the analyst's knowledge be
revised before the reported effect could be wrong by twenty times its own size?

The committed sweep used the grid `0.01, 0.05, 0.10, 0.25, 0.50, 1.00`. That grid is a
convention of one run, not a property of the construction. Nothing in the theory, and
nothing in `epsilon.certify`, caps ε at 1.

**The scientific point.** On the committed grid, 44 of the 66 informative instances reach
`UNREACHED` by ε = 1.00 — the threshold is never crossed anywhere in the space. The paper
currently reads that as a *ceiling*: no perturbation moves the estimate by its own
magnitude. The open question this experiment settles is whether the remaining 22 instances
**ramp** — whether `r_eps` grows with ε above 1 — or whether they too go straight to
`UNREACHED`. Both answers are publishable and they say different things.

---

## 2. This is NOT derivable from `results/final_table/`

Checked directly, do not re-check:

- `results/final_table/instances.jsonl` has keys `condition, error, len_k, n_knowledge,
  network, r_eps, r_val, seconds, status, theta_z, x, y`. `r_eps` is a dict keyed by the
  **six committed ε only**.
- `results/final_table/summary.json` per-network entries hold `r_eps` keyed by the same six ε,
  each as `{value, n_resolved, n_unreached}`.

**Neither file stores the staircase `β↑(d)`.** Only its thresholded values at six points are
kept. `r_eps` at ε = 2 cannot be recovered from `r_eps` at ε = 1 — the radii are order
statistics of a curve that was not retained. A re-run is required.

---

## 3. The work is bounded: only 22 instances can move

The paper proves `F_ε' ⊆ F_ε` whenever `ε ≤ ε'` (monotonicity of the violating set in the
tolerance), hence **`r_eps` is non-decreasing in ε**. Therefore:

> Any instance already `UNREACHED` at ε = 1.00 is `UNREACHED` at every ε > 1.00.

That is 44 of the 66 informative instances, settled a priori. It does not need computing,
and a re-run that reports it as a fresh finding is wasting the reader's time.

**Exactly 22 instances have a finite `r_eps` at ε = 1.00 and can still move.** They are the
whole target of this experiment:

| Network | Query | `r_val` | `r_eps` @ 1.00 |
|---|---|---|---|
| `paths` | `11 → 14` | 11 | 11 |
| `ecoli70` | `cspG → hupB` | 4 | 4 |
| `barley` | `dg25 → s2225` | 3 | 3 |
| `water` | `CKNI_12_15 → CBODN_12_45` | 3 | 3 |
| `magic-irri` | `G3212 → FT` | 2 | 2 |
| `Acid_1996` | `x1 → x10` | 1 | 1 |
| `Acid_1996` | `x1 → x11` | 1 | 1 |
| `Didelez_2010` | `Age → HRT` | 1 | 1 |
| `Didelez_2010` | `Age → S` | 1 | 1 |
| `Didelez_2010` | `Age → TCI` | 1 | 1 |
| `Kampen_2014` | `AFF → FTW` | 1 | 1 |
| `Kampen_2014` | `AFF → SAN` | 1 | 1 |
| `Schipf_2010` | `A → TT` | 1 | 1 |
| `Sebastiani_2005` | `ANXA2.5 → ANXA2.11` | 1 | 1 |
| `Thoemmes_2013` | `e0 → s2` | 1 | 1 |
| `Thoemmes_2013` | `e0 → s3` | 1 | 1 |
| `andes` | `HORIZ53 → SNode_118` | 1 | 1 |
| `asia` | `asia → either` | 1 | 1 |
| `asia` | `bronc → dysp` | 1 | 1 |
| `child` | `CO2 → DuctFlow` | 1 | 1 |
| `magic-niab` | `G1217 → YLD` | 1 | 1 |
| `mediator` | `I → Y` | 1 | 1 |

Source: `results/final_table/instances.jsonl`, filtered to `status == "ok"`,
`r_val != 0`, `r_eps["1.0"] != -1`.

**Four of these five carry the paper's current Table 1** (`paths 11→14`, `ecoli70
cspG→hupB`, `barley dg25→s2225`, `water CKNI_12_15→CBODN_12_45`). The fifth Table 1 row,
`barley dg25→aks_m2`, is `UNREACHED` from ε = 0.25 and will stay so. If the budget forces a
choice, **run those four first** — they are the rows the paper table will be rebuilt from.

---

## 4. What to run

### 4.1 The sweep

```bash
PYTHONPATH=src .venv/bin/python experiments/final_table.py \
  --condition D_LLM \
  --epsilons 1,2,5,10,20 \
  --queries-per-network 5 \
  --time-limit 120 \
  --seed 20260917 \
  --out results/final_table_eps_gt1
```

**The `--seed 20260917` is not optional.** The corpus ships graph structure only, so `r_eps`
is evaluated against a seeded linear-Gaussian SEM attached to the ground-truth DAG, with the
per-instance seed derived as `sha256(f'{seed}:{network}:{x}:{y}')`. Reusing `20260917`
reproduces the identical SEM per instance, which is the only thing that makes the new radii
comparable to the committed ones. A different seed produces a different — equally valid, but
**incomparable** — panel, and any table mixing the two is wrong.

Write to a **new** `--out` directory. Do not overwrite `results/final_table/`; the paper
currently cites it and the ε ≤ 1 panel must remain reproducible.

### 4.2 Flag spellings — verify before trusting

The flag names above are inferred from the `args` block of the committed
`results/final_table/summary.json` (`condition, epsilons, queries_per_network, networks,
skip, seed, time_limit, max_nodes, out, report_only`), **not** read from the script. The
brief's author did not have read access to `experiments/final_table.py`. Run
`--help` first and correct any spelling that differs. Two values in that block are of
unknown provenance and matter:

- `"skip": "pathfinder"` — is this an argparse default, or was it passed? If a default, the
  32-network corpus already excludes `pathfinder` and the command above reproduces it.
- `"seed": 20260917` — same question. If the default differs, pass it explicitly as above.

Settle both from `--help` / the argparse block and record the answer in the run notes.

### 4.3 Strongly recommended: dump the staircase

This is the second time the ε grid has had to be chosen in advance, and a third re-run is
likely. If `certify` already computes `β↑(d)` per shell — the paper states one traversal
answers every ε — then **persisting that curve per instance makes every future ε free**, and
this whole experiment becomes a post-processing step rather than a 74-minute sweep. Add a
`staircase: {d: beta}` field to each `instances.jsonl` record if the cost is a few lines.
Raise it as a finding even if you do not implement it.

---

## 5. Deliverables

1. `results/final_table_eps_gt1/instances.jsonl` and `summary.json`, same schema as the
   committed panel.
2. A **per-query** table over the 22 movable instances: `network, x, y, |K|, r_val`, then
   `r_eps` at ε = 1, 2, 5, 10, 20. `UNREACHED` rendered as a sentinel, never a number.
3. A rebuilt Table 1 for the paper over the five `r_val > 2` queries, at the new grid.
   The current version lives at the end of `paper/sections/06_epsilon.tex` (label
   `tab:realeps`); match its column layout so the swap is a replacement, not a redesign.
4. One paragraph answering the actual question: **does the radius ramp above ε = 1, or does
   it go to `UNREACHED`?** Give the count of instances in each class.
5. Run notes recording the resolved flags from §4.2, wall-clock, and any new timeouts.

---

## 6. House rules — these are not negotiable

- **Never invent a number.** Every figure must trace to `results/final_table_eps_gt1/`. If a
  number is not in the output, it is a `TODO`, not a guess.
- **`r` vs `r − 1`.** `r` is the distance to the nearest failure; the number of *safe* moves
  is `r − 1`. `r = 1` is **zero** safe retractions. Never paraphrase `r = 1` as "tolerates
  one error".
- **The radius is `r_val`, repo-wide, with no alias.** Do not introduce `r_hop` or any other
  name. When the *unit* rather than the radius is meant, write "hop units" / "claim units" in
  words. See `PAPER_NARRATIVE.md` §11 and commit `06eb9e5`.
- **Report per query, not per-network medians.** See §7 — this is the trap that has already
  produced one false claim in the paper.
- **`UNREACHED` is never averaged.** It is `−1` in the JSON, a sentinel meaning no element of
  the space fails, not a large radius.
- **Degenerate ≠ fragile.** `r_val = 0` means the optimal adjustment set at `G₀` was already
  empty and invalid — the query is ill-posed under that `K`, not close to breaking. Exclude
  from every median; never fold into a "the radius is low" claim.
- **Honest denominators.** 21 of 32 networks are informative; 11 are blank (timeout or
  degenerate) and are **unknown, not robust**. Never silently drop them.
- **`r_eps` is SEM-conditional; `r_val` is not.** Every sentence reporting an `r_eps` number
  on real data must carry the assumption. `r_val` is purely graph-theoretic. Canonical
  wording is in the `assumptions` block of `summary.json`.
- **One condition only.** Everything here is `D_LLM`. Nine other elicitation conditions sit
  unused in `results/elicit/knowledge.json`; claims do not transfer to them.
- **Not pre-registered.** `docs/EVALUATION_PLAN.md` does not mention `r_eps` anywhere
  (verified). Present any result as new, never as the delivery of a planned table.
- **Scope the ε findings to this corpus and this knowledge source.** Say "does not separate
  here", not "does not separate".

---

## 7. The trap that has already cost the paper a false claim

**Do not aggregate `r_eps` into per-network medians across ε.** The set of instances with a
finite `r_eps` *shrinks* as ε rises, so a median taken per ε drifts onto a changing
denominator and produces artefacts in both directions:

- `barley`: median `r_val` = 3, but median `r_eps` = 2.5 at ε = 0.25 and 0.50 — appearing to
  show `r_eps < r_val`, which would contradict the theorem. One instance left the
  denominator.
- `paths`: median `r_val` = 1, median `r_eps` = 11 at ε = 1.00 — appearing to show a gain of
  1 → 11. **There is no gain.** `paths` has three non-degenerate instances at radius 1, 1 and
  11; at that tolerance the two radius-1 instances stop being reached and leave the median,
  which then rests on the instance that was 11 all along. No `paths` instance ever has
  `r_eps > r_val`.

The second artefact was asserted as a real gain in **§3.4 of
`docs/PAPER_INTEGRATION_FINAL_TABLE.md`**, was copied from there into
`paper/sections/06_epsilon.tex`, and has since been retracted in the paper (commit
`a4ef293`). **§3.4 of that brief is still wrong and should be corrected** — it is the
likeliest way for this error to re-enter.

This effect will be *stronger* at ε ∈ {1, 2, 5, 10, 20}, because more instances go
`UNREACHED` at every step. Report per query. If an aggregate is genuinely wanted, fix the
instance set first and state it.

---

## 8. What the answer will look like

Two outcomes, both worth reporting, and the prose is different for each:

- **Ramp.** Some of the 22 carry a strictly larger `r_eps` at ε = 2, 5, 10 or 20. This is the
  first real-network evidence that the ε-refinement buys a *longer runway* and not only a
  ceiling, and it would qualify the paper's current reading. Report which instances, and at
  which ε.
- **Ceiling.** All 22 go from their ε = 1.00 value straight to `UNREACHED`. Then the ceiling
  reading is confirmed above 100% as well, and the honest sentence is that on this corpus, at
  this knowledge source, the bias crosses *any* tolerance essentially at the validity
  boundary or never.

Either way the claim is about **this corpus and `D_LLM`**, and the synthetic sweep — where
the refinement demonstrably does pay — is not contradicted by it.

---

## 9. The paper edit that follows

`paper/sections/06_epsilon.tex`, `\myparagraph{On real graphs.}` at the end of
`sec:epsilon-eval`, and the table `tab:realeps` directly after it. The section reports
results as prose under `\myparagraph` headers; `tab:realeps` is the only tabular environment
in it. Keep it that way — one table, replaced in place.

The paragraph currently states that `r_eps` "does not separate from `r_val` here" on the
ε ≤ 1 grid, and that what the refinement buys is a ceiling. If §8 returns *ramp*, that
paragraph needs amending, not appending: do not leave two paragraphs making opposite claims.
