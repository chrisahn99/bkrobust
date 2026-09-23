# Result — the ε-radius above a tolerance of 1

Answers `docs/EXPERIMENT_EPSILON_ABOVE_ONE.md`. Every number below traces to
`results/final_table_eps_gt1/`; nothing is estimated. Condition is `D_LLM` only.

**`r_eps` is SEM-conditional.** The corpus ships graph structure only, so every `r_eps` here is
evaluated against a seeded linear-Gaussian SEM attached to the ground-truth DAG, seeded per
instance from `--seed 20260917`. `r_val` is purely graph-theoretic and carries no such
assumption. Canonical wording is in the `assumptions` block of `summary.json`.

**Not pre-registered.** `docs/EVALUATION_PLAN.md` does not mention `r_eps`. This is a new
result, not the delivery of a planned table.

---

## 1. The answer: ceiling, not ramp

**No instance ramps. 0 of 22 gain a single certified shell above ε = 1.00.**

Of the 22 instances with a finite `r_eps` at ε = 1.00, every one holds its ε = 1.00 value for
as long as it is reached at all and then goes to `UNREACHED`. Not one carries a strictly larger
`r_eps` at ε = 2, 5, 10 or 20. Splitting by how long they survive:

| Class | Count | Meaning |
|---|--:|---|
| **Ramp** — strictly larger `r_eps` above ε = 1 | **0** | — |
| **Ceiling** — `UNREACHED` at ε = 2 | 18 | the ceiling sits in (1, 2] |
| **Ceiling** — holds its value past ε = 2, then `UNREACHED` | 4 | a higher ceiling, still no gain |

The §8 *ceiling* reading of the brief is confirmed, and confirmed well above 100%: on this
corpus, at this knowledge source, the bias crosses any tolerance essentially at the validity
boundary or never. The synthetic sweep, where the refinement does buy extra shells, is not
contradicted — it is simply not reached here.

The mechanism is visible in the newly persisted staircase: **β↑(d) is flat above `r_val` on 21
of the 22**. The entire bias risk arrives at the validity boundary in one step, so raising ε
cannot buy a shell — it can only stop reaching the step at all. This is the paper's
"step, not ramp" shape (84.1% of synthetic instances) holding at 95% here. The single
exception, `Kampen_2014 AFF → FTW`, rises from 1.6098 at d = 1 to 1.6187 at d = 2 and then
flattens — too small a rise to cross any grid point, so it produces no gain either.

## 2. Per-query table over the 22 movable instances

`r_val` and `r_eps` are radii in hop units. `UNREACHED` (`-1` in the JSON) is a sentinel
meaning no element of the space fails — never a large radius, never averaged. Recall that
`r = 1` is **zero** safe retractions, not "tolerates one error".

| Network | Query | \|K_G0\| | `r_val` | ε=1 | ε=2 | ε=5 | ε=10 | ε=20 | `beta_top` |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|
| `paths` | `11 → 14` | 14 | 11 | 11 | UNREACHED | UNREACHED | UNREACHED | UNREACHED | 1.000 † |
| `ecoli70` | `cspG → hupB` | 13 | 4 | 4 | UNREACHED | UNREACHED | UNREACHED | UNREACHED | 1.000 † |
| `barley` | `dg25 → s2225` | 7 | 3 | 3 | UNREACHED | UNREACHED | UNREACHED | UNREACHED | 1.525 |
| `water` | `CKNI_12_15 → CBODN_12_45` | 6 | 3 | 3 | 3 | UNREACHED | UNREACHED | UNREACHED | 2.245 |
| `magic-irri` | `G3212 → FT` | 8 | 2 | 2 | UNREACHED | UNREACHED | UNREACHED | UNREACHED | 1.225 |
| `Acid_1996` | `x1 → x10` | 1 | 1 | 1 | UNREACHED | UNREACHED | UNREACHED | UNREACHED | 1.927 |
| `Acid_1996` | `x1 → x11` | 1 | 1 | 1 | UNREACHED | UNREACHED | UNREACHED | UNREACHED | 1.629 |
| `Didelez_2010` | `Age → HRT` | 3 | 1 | 1 | 1 | UNREACHED | UNREACHED | UNREACHED | 3.453 |
| `Didelez_2010` | `Age → S` | 3 | 1 | 1 | UNREACHED | UNREACHED | UNREACHED | UNREACHED | 1.223 |
| `Didelez_2010` | `Age → TCI` | 3 | 1 | 1 | UNREACHED | UNREACHED | UNREACHED | UNREACHED | 1.473 |
| `Kampen_2014` | `AFF → FTW` | 10 | 1 | 1 | UNREACHED | UNREACHED | UNREACHED | UNREACHED | 1.619 |
| `Kampen_2014` | `AFF → SAN` | 10 | 1 | 1 | 1 | 1 | 1 | 1 | 58.153 |
| `Schipf_2010` | `A → TT` | 6 | 1 | 1 | UNREACHED | UNREACHED | UNREACHED | UNREACHED | 1.422 |
| `Sebastiani_2005` | `ANXA2.5 → ANXA2.11` | 4 | 1 | 1 | UNREACHED | UNREACHED | UNREACHED | UNREACHED | 1.000 † |
| `Thoemmes_2013` | `e0 → s2` | 1 | 1 | 1 | UNREACHED | UNREACHED | UNREACHED | UNREACHED | 1.000 † |
| `Thoemmes_2013` | `e0 → s3` | 1 | 1 | 1 | UNREACHED | UNREACHED | UNREACHED | UNREACHED | 1.000 † |
| `andes` | `HORIZ53 → SNode_118` | 10 | 1 | 1 | UNREACHED | UNREACHED | UNREACHED | UNREACHED | 1.000 † |
| `asia` | `asia → either` | 3 | 1 | 1 | UNREACHED | UNREACHED | UNREACHED | UNREACHED | 1.000 † |
| `asia` | `bronc → dysp` | 3 | 1 | 1 | UNREACHED | UNREACHED | UNREACHED | UNREACHED | 1.256 |
| `child` | `CO2 → DuctFlow` | 10 | 1 | 1 | UNREACHED | UNREACHED | UNREACHED | UNREACHED | 1.000 † |
| `magic-niab` | `G1217 → YLD` | 5 | 1 | 1 | UNREACHED | UNREACHED | UNREACHED | UNREACHED | 1.127 |
| `mediator` | `I → Y` | 4 | 1 | 1 | 1 | UNREACHED | UNREACHED | UNREACHED | 4.967 |

† `beta_top` equals 1.0 to within 4e-15 — see §3. The ε = 1 entry on these rows is a
floating-point tie, not a crossing.

No aggregate is taken over this table. The instance set with a finite `r_eps` shrinks as ε
rises, so a per-ε median would drift onto a changing denominator (brief §7).

## 3. A finding the brief did not anticipate: ε = 1.00 lands on a mass point

`beta_top` = `B(Ĉ)`, the largest bias any perturbation anywhere can produce, is now persisted
per instance. Its distribution over the 66 informative instances is not smooth:

| `beta_top` | instances |
|---|--:|
| < 0.5 | 18 |
| [0.5, 1) | 1 |
| **= 1.0 (to within 4e-15)** | **33** |
| (1, 2] | 10 |
| (2, 5] | 3 |
| > 20 | 1 |

**Half the informative instances sit at β↑ = 1.0 exactly.** This is structural, not
coincidental: relative bias is `|θ_Z − τ| / |θ_Z|`, so β↑ = 1 is precisely the case where the
worst admissible adjustment set explains the effect away entirely, τ → 0. Verified directly at
`asia asia → either`, whose witness at d = 1 reports `tau_max = 5.03e-16` against
`theta_z = −1.752`.

The consequence is that **ε = 1.00 is the worst possible place to put a grid point.** The
radius test is the strict `β↑(d) > ε`, so the 33 tied instances are split by whether
floating-point subtraction lands them a few ulp above or below 1.0:

- 8 land at 1.0000000000000002–1.0000000000000042 and are reported as **crossing** ε = 1;
- 25 land at 1.0 exactly or a few ulp below and are reported **`UNREACHED`**.

Those 8 are `paths 11→14`, `ecoli70 cspG→hupB`, `Sebastiani_2005 ANXA2.5→ANXA2.11`,
`Thoemmes_2013 e0→s2`, `Thoemmes_2013 e0→s3`, `andes HORIZ53→SNode_118`, `asia asia→either`,
`child CO2→DuctFlow`. They are **8 of the 22** "movable" instances, and they include **2 of
the 5 Table 1 rows**. The committed panel's ε = 1.00 column is, on those rows, an arbitrary
tie-break rather than a measurement — `paths 11→14`'s headline `r_eps = 11` at ε = 1.00 among
them.

This does not change the answer in §1 — those instances are `UNREACHED` from ε = 2 upward
either way, and the ceiling reading is unaffected. It does mean **no claim should rest on the
ε = 1.00 column specifically**, and a future grid should straddle the mass point (0.9 / 1.1)
rather than sit on it.

## 4. Table 1, rebuilt at the new grid

The five `r_val > 2` queries. There is no `tab:realeps` in the tree to replace — see §6.

| Network | Query | \|K_G0\| | `r_val` | ε=1 | ε=2 | ε=5 | ε=10 | ε=20 |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| `paths` | `11 → 14` | 14 | 11 | 11 † | UNREACHED | UNREACHED | UNREACHED | UNREACHED |
| `ecoli70` | `cspG → hupB` | 13 | 4 | 4 † | UNREACHED | UNREACHED | UNREACHED | UNREACHED |
| `barley` | `dg25 → aks_m2` | 7 | 3 | UNREACHED | UNREACHED | UNREACHED | UNREACHED | UNREACHED |
| `barley` | `dg25 → s2225` | 7 | 3 | 3 | UNREACHED | UNREACHED | UNREACHED | UNREACHED |
| `water` | `CKNI_12_15 → CBODN_12_45` | 6 | 3 | 3 | 3 | UNREACHED | UNREACHED | UNREACHED |

† floating-point tie at β↑ = 1.0 (§3).

`r_eps` never exceeds `r_val` on any of these rows, at any ε on the grid.

## 5. Run notes

**Command** (flag spellings verified against `--help` and the argparse block before running):

```
PYTHONPATH=src .venv/bin/python experiments/final_table.py \
  --condition D_LLM --epsilons 1,2,5,10,20 --queries-per-network 5 \
  --time-limit 120 --seed 20260917 --out results/final_table_eps_gt1
```

**§4.2 resolved.** Both values in the committed `summary.json` `args` block are **argparse
defaults**, not values that were passed (`experiments/final_table.py:506–514`;
`DEFAULT_SKIP = "pathfinder"` at line 53):

- `"skip": "pathfinder"` — a default. The 32-network corpus already excludes `pathfinder`, and
  the command above reproduces it without passing `--skip`.
- `"seed": 20260917` — also the default. Passing it explicitly is correct and changes nothing.

Every flag name inferred in the brief was spelled correctly; none needed correction.

**Seed reproduction verified.** `r_val` agrees with the committed panel on all 97 shared `ok`
instances (0 disagreements), and `r_eps` at ε = 1.00 agrees on all 22 movable instances
(0 mismatches). The SEMs are identical, so the new radii are comparable to the committed ones.

**Wall-clock.** 15:43:07–17:03:47 = 80.7 min; 4805 instance-seconds over 160 instances.

**Status counts.** New: 97 ok / 27 error / 36 timeout. Committed: 99 / 27 / 34. The informative
denominator is **unchanged at 66**, and 21 of 32 networks are informative — verified in this
run, not carried over. The other 11 are **unknown, not robust**, and split two ways: 7 resolved
no query at all (`arth150`, `diabetes`, `link`, `munin`, `munin2`, `munin3`, `munin4` — all
timeout or error) and 4 resolved only degenerate ones (`hailfinder`, `insurance`, `munin1`,
`sachs` — `r_val = 0`, ill-posed under that `K`, not close to breaking).

**Two new timeouts**, both on `munin4` (`L_LNLE_ULND5_RD_E → L_LNLE_ADM_MALOSS`,
`L_LNLE_ULN_DIFSLOW → L_LNLE_ADM_DE_REGEN`). Both were `r_val = 0` **degenerate** in the
committed run — ill-posed under that `K`, excluded from every median — so no reported number
changes. They completed in ~21.5 s before and exceeded 120 s here, which is machine-load
variance, not an ε effect: the ε grid does not enter the `r_val` search.

**§4.3 implemented.** Each record now carries `beta_top` and `staircase` `{d: β↑(d)}`, both in
the certificate's reporting units. `Certificate` already retained the staircase in memory, so
this was four lines in `run_instance`. **Every future ε is now a lookup, not a re-run**:
`r_eps(ε)` is the first `d` with `staircase[d] > ε`, and `ε ≥ beta_top` is `UNREACHED`
outright. The staircase covers only the shells the traversal touched (it stops early once it
crosses); `beta_top` is over the whole space and is always exact, which is what makes the
`UNREACHED` half of any future question free.

**On `UNREACHED`.** It is exact here, not a budget artefact. `epsilon_radius` runs an infinity
test first (`profile.py:511`): `top_bias` evaluates `B(Ĉ)` over the entire space, and
`B(Ĉ) ≤ ε` returns `UNREACHED` with `exact=True`. The 200,000 `DEFAULT_SHELL_CAP` is a shell
*width* cap and was not what produced any `UNREACHED` above.

## 6. Two corrections to the briefing documents

**(a) `docs/PAPER_INTEGRATION_FINAL_TABLE.md` §3.4 is wrong on two counts**, as brief §7
predicted for the first:

- "the single visible gain is `paths`, 1 → 11 at ε = 100%" — there is **no gain**. `paths` has
  three non-degenerate instances at radius 1, 1 and 11; at ε = 1.00 the two radius-1 instances
  stop being reached and leave the median, which then rests on the instance that was 11 all
  along. No `paths` instance ever has `r_eps > r_val`, at any ε on either grid.
- "2 of **53** informative instances" — the 2 is right (`Kampen_2014 AFF→CDR` at ε = 0.25 and
  0.50, `Schipf_2010 A→T2DM` at ε = 0.50, verified in the committed panel), but the
  denominator is **66**, not 53.

**(b) Brief §9's paper target does not exist in this tree.** There is no
`paper/sections/06_epsilon.tex` (the file is `paper/sections/epsilon-radius.tex`), no
`\myparagraph{On real graphs.}`, and no `tab:realeps` anywhere in `paper/`. `sec:epsilon-eval`
ends at `\paragraph{Finite samples.}` and contains no tabular environment at all — consistent
with §6.1 of `PAPER_INTEGRATION_FINAL_TABLE.md`, which records introducing that table as an
open decision **deliberately not taken**. The §4 table in §9 of the brief is therefore a *new*
table, not a replacement, and the paragraph it describes amending has not been written yet.
No paper edit was made.
