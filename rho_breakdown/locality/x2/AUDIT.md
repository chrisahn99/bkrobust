# X2 DESIGN AUDIT — written against `PREREG.md`, before implementation

Adversarial read of a design I did not write. E1 was refuted 3/3 on *interpretation* while every
number reproduced; the cheapest way to repeat that is to skip this file. X2's PREREG is a better
document than E1's was — the six-event taxonomy, the four-denominator rule and the §0 disclosure
table are all real improvements. It is also, as written, **not runnable**: its load-bearing gate
fails before the run, its primary arm measures a quantity that is identically zero, its one novel
arm is already computed and sitting on disk, and two of its four registered predictions are cleared
a hundred times over by a nine-second pilot.

Everything below is a number, not an argument. All of it was computed **today, before any X2 code
exists**, from `~/latent-causal/e1prime-se/results/arm1_*.json` and from a bounded pilot in
`/tmp/x2audit/` on betelgeuse that reuses `e1prime-se/code/{graphs,adjust,scm,se}.py` unmodified.
Nothing under `e1prime-se/` or `rho-breakdown-knowledge/` was touched.

---

## A1. 🔴 G4 — the gate the PREREG says decides whether X2 is worth reading — FAILS before the run

§6 G4: *"`L_set(0, ρ=1)` ∈ **[0.05, 0.35]** on `B-licensed`. Outside that band the detector is not
demonstrated to work and **no zero at `dmin ≥ 1` is readable**."*

ARM B is a bit-identical re-run of E1′ at seed 20260819, so `L_set(0, ρ=1)` is already on disk.
Recomputed per-statement at `n = ∞` with the PREREG's own `D_am` denominator:

| ensemble | `D_all` | `D_con` | `D_am` | `N_set` | **`L_set(0, ρ=1)`** | in [0.05, 0.35]? |
|---|---|---|---|---|---|---|
| `original` | 6,640 | 4,224 | 2,721 | 1,007 | **0.3701** | ✗ |
| **`licensed`** | 4,339 | 2,923 | 1,944 | 906 | **0.4660** | ✗ **the registered cell** |
| `large` | 2,111 | 1,715 | 1,026 | 208 | **0.2027** | ✓ |
| `k8` / K4 | 256 | 202 | 91 | 69 | **0.7582** | ✗ |
| `k8` / K8 | 2,789 | 1,458 | 1,079 | 599 | **0.5551** | ✗ |

`B-licensed` = **0.466**, outside the band by 33 %. By the PREREG's own rule, X2 aborts at gate G4
and no verdict on S5, S8 or S9 is readable. Four of the five cells fail.

The band was set by widening E1′'s [0.12, 0.18] "because X2's denominator is `D_am` where E1′'s was
the consistent-member count". The widening went the wrong way: moving from `D_con` to `D_am`
*removes* the non-amenable members from the denominator and therefore **raises** the rate
(licensed: 906/2,923 = 0.310 → 906/1,944 = 0.466). The band was recentred on the number the change
was about to invalidate.

> **Mitigation M1 (binding).** Replace the single band with per-ensemble bands taken from the table
> above, ±0.05 absolute: `original` [0.32, 0.42] · `licensed` [0.42, 0.52] · `large` [0.15, 0.25] ·
> `k8`/K4 [0.71, 0.81] · `k8`/K8 [0.51, 0.61]. G4 then becomes what it was meant to be — a check
> that X2's detector reproduces a known positive — instead of an abort switch aimed at the design's
> own head. State in `RESULTS.md` that these targets were fixed from archived E1′ data in this
> AUDIT, before implementation, and that G4 is consequently a **reproduction** check, not evidence.

---

## A2. 🔴 ARM A — "the primary existence test and the arm no prior experiment ran" — is null by construction

§2.1 makes ARM A primary: `K = ∅`, baseline `G0 = C`, every undirected edge, both orientations.
§0 P2 says the unrun parts of X2 are "ARM A, ARM C, ARM D, `E_gain`". I ran ARM A's exact operator.

`/tmp/x2audit/power2.py`, `/tmp/x2audit/power3.py`, seed 777, `e1prime-se` machinery unmodified,
restricted to the SCMs where `E_set` is even defined (`amen(C) = True`, since `G0 = C`):

| ensemble | `amen(C)` rate | `D_am` at `dmin=0` | `E_set(0)` | `D_am` at `dmin≥1` | `E_set(≥1)` | `E_ident` | `E_meek` |
|---|---|---|---|---|---|---|---|
| `A-licensed` | 960 / 7,386 = **0.130** | 1,250 | **0** | 6,508 | **0** | **0** | **0** |
| `A-large` | 600 / 1,136 = **0.528** | 644 | **0** | 5,740 | **0** | **0** | **0** |
| `A-k8` | 22 / 1,162 = **0.019** | 20 | **0** | 402 | **0** | **0** | **0** |
| `C-deep` | 400 / 2,449 = **0.163** | 380 | **0** | 3,586 | **0** | **0** | **0** |
| **pooled** | | **2,294** | **0** | **16,236** | **0** | **0** | **0** |

**18,530 exhaustive single-statement trials. Zero events of every kind, at every distance,
including `dmin = 0`.** ARM A cannot pass G4 — not marginally, but by two to three orders of
magnitude against `B-licensed`'s 0.466. A search whose detector never fires anywhere reports the
same zero whether locality holds or not; ARM A's `N_set(r ≥ 1) = 0` is therefore **uninformative
about S5**, and pooling it with ARM B's trials (§4.1) dilutes the only denominator that carries
signal.

Three of those zeros are structural, and each of them is its own finding:

- **`E_meek ≡ 0` is a theorem, not a measurement.** Every undirected edge of a CPDAG is oriented
  both ways somewhere in its Markov equivalence class, so a single admissible statement can never
  produce a cycle or a new v-structure. 0 of 18,530 confirms it. So in ARM A `D_all = D_con = D_am`
  identically, and §1.4's four-denominator discipline prints three copies of one number. ARM A also
  cannot speak to S3 ("two thirds pass the consistency check") at all.
- **`E_ident ≡ 0`**: a single statement never destroyed amenability of a bare CPDAG in 18,530 trials,
  while in ARM B it fires 109 / 129 / 24 / 40 / 109 times at `dmin = 1` alone. The difference is the
  base `K`, not the distance.
- **The mechanism for `E_set ≡ 0` is measurable and is the actual result.** Under `amen(C)`, `cn`
  has mean size 1.46–2.82 and the mean number of **undirected edges incident to `cn`** is **0.268**
  (`A-licensed`; only 97 of 960 amenable SCMs have even one). Amenability plus Meek closure has
  already compelled essentially every edge touching the causal-node set, so `pa(cn)` is frozen and
  one further orientation cannot move it. Across all four ensembles only 14 / 22 / 20 trials were
  incident to `cn` at all, and none of them moved `O*`.

That last paragraph is worth more than the arm it kills. It is the shape of the theorem R1 asked
for, at `|K| = 0`: *if `C` is amenable w.r.t. `(X,Y)`, then `O*(M({s})) = O*(C)` for every
admissible statement `s`, at any distance.* 18,530 trials without a counterexample, plus a
mechanism, plus exhaustive verification at `p ≤ 6`, is a lemma. "0 of 210,000" is a table row.

> **Mitigation M2 (binding).** ARM A is **demoted from primary to mechanism arm** and its §4.1 role
> is withdrawn: its trials do **not** enter the pooled `D_am` behind the S5 verdict. It is re-scoped
> to two jobs it can actually do: (i) `E_gain`, which is the only event it produces (A4); (ii) the
> lemma — report `amen(C)` rate, `|cn|`, undirected-edges-incident-to-`cn`, and `E_set = 0 / D_am`
> as evidence for *"a single statement cannot move `O*` from a bare amenable CPDAG"*, with an
> exhaustive `p ∈ {4,5,6}` enumeration over all CPDAGs × all queries × all single statements as the
> honest substitute for a proof. If the exhaustive search finds a counterexample at `p ≤ 6`, that
> instance is the paper.

> **Mitigation M3 (binding).** The primary existence test for S5 becomes **ARM B's operator with an
> exhaustive statement sweep**: base `K` of size 1–4 drawn as in `run_arm1.py:139`, then flip each
> of the `|U(C)|` statements in turn rather than sampling 4. That keeps the dynamic range that makes
> `L_set(0) ≈ 0.2–0.76` and removes the `|K| = 4` sampling that §2.1 rightly objects to. Register
> `L_set(0)` for this operator as G4's positive control per M1.

---

## A3. 🔴 ARM D is not an experiment. It is already computed, and its registered prediction is retrodiction — and half of it is wrong

§0's disclosure table names ARM D as one of the four things X2 does not yet know (P2: *"the unrun
parts are ARM A, ARM C, ARM D"*). But §2.5 defines ARM D as a **re-partition of ARM B's ball**, and
ARM B is a same-seed re-run whose members — `est` at `n = ∞`, `consistent`, `amenable`,
`stmt_hopdist` — are all archived in `arm1_*.json`. I computed `Q_r` and `Cost_r` from the archive
in about ninety seconds, exactly as §2.5 specifies (`B^r_ρ` = members whose flipped statements all
lie within `near_r`; endpoints compared at 1e-9 absolute; `K` itself included in both intervals):

| ensemble | `Q_{−1}` | `Q_0` | `Q_1` | `Q_2` | mean `Cost_0` | mean `Cost_1` | mean `Cost_2` |
|---|---|---|---|---|---|---|---|
| `original` | 0.4907 | 0.1244 | **0.0011** (3/2,637) | 0.0000 | 0.4525 | 0.8268 | 0.8982 |
| `licensed` | 0.4854 | 0.1452 | **0.0014** (3/2,087) | 0.0000 | 0.3557 | 0.8083 | 0.8988 |
| `large` | 0.0970 | **0.0143** (35/2,443) | **0.0000** | 0.0000 | 0.1402 | 0.2808 | 0.4899 |
| `k8`/K4 | 0.5118 | 0.1969 (25/127) | 0.0000 | 0.0000 | 0.3853 | 0.9213 | 0.9801 |
| `k8`/K8 | 0.9010 | 0.5239 (307/586) | 0.0000 | 0.0000 | 0.3076 | 0.9663 | 0.9913 |

Against §4.2's registered predictions:

- *"`Q_0 ≥ 0.01` on `large` and on `k8` (REFUTED at r = 0)"* — **correct** (0.0143, 0.1969, 0.5239),
  and also true on `original` and `licensed`, which the prediction did not name.
- *"predicted `Q_1 ≤ 0.001` on `licensed`"* — **wrong**: 0.0014. `licensed` misses the SUPPORTED bar.
- *"predicted `Q_1 > 0.001` on `large`"* — **wrong, and backwards**: `Q_1 = 0.0000` on `large`, the
  cleanest cell in the table. The §2.3 mechanism (longer causal paths in bigger graphs ⇒ locality
  breaks) is refuted by data that was on disk when the prediction was written.

And the full S8 verdict is therefore **already determined**, before implementation:
**REFUTED at r = 0** (`Q_0 ≥ 0.01` on 4 of 5 cells) · **INCONCLUSIVE at r = 1** (max `Q_1` = 0.0014,
between the 0.001 and 0.01 bars, and mean `Cost_1` = 0.81–0.97 on `original`/`licensed`/`k8` — the
radius that is safe prunes nothing) · **NOT SUPPORTED at any r**, because §4.2 quantifies over
*every* ensemble and only `large` clears both halves anywhere.

> **Mitigation M4 (binding).** Move ARM D out of §4's registered predictions and into §0's
> disclosure table as **P8**, with the full `Q_r` / `Cost_r` matrix above, attributed to this AUDIT.
> §4.2's prediction paragraph is deleted, not amended: a prediction stated after the number exists
> is not a prediction, and the PREREG's own §0 sets that standard. The S8 verdict is then reported
> as a **measurement recovered from archived data**, which is honest and costs nothing — ARM D was
> always free.

> **Mitigation M5 (binding).** `Cost_r ≤ 0.5` has **no aggregator** in §4.2. Mean, median and
> per-SCM fraction give different verdicts: on `original` at r=1 the mean is 0.827 (fails) while
> 28.9 % of SCMs are under 0.5. Fix it now as **the mean over SCMs**, and report the per-SCM
> fraction beside it. Also report the `(Q_r, Cost_r)` **Pareto frontier for r ∈ {−1,0,1,2,3}** — §2.5
> already says the trade-off, not either number, is the reportable object; §4.2's pass/fail
> contradicts it.

---

## A4. 🔴 Two registered bars are FREE — cleared ~100× over by a nine-second pilot

The caller asked whether E1′'s `crit='ident'` failure recurs. It does, twice.

**§4.4, the `E_gain` prediction.** *"PREDICTED: `N_gain(r ≥ 1, ρ = 1) ≥ 1` on `A-large` or on
`C-deep`."* Measured in a 400-SCM pilot per ensemble (`/tmp/x2audit/power.py`, 2–9 s each):

| ensemble | `N_gain` at `dmin ≥ 1` | `N_gain` at `dmin = 0` | per SCM |
|---|---|---|---|
| `A-licensed` | **116** | 302 | 0.29 |
| `A-large` | **112** | 372 | 0.28 |
| `A-k8` | **50** | 28 | 0.16 |
| `C-deep` | **98** | 507 | 0.25 |
| `C-deep8` | **104** | 288 | 0.59 |

The bar is `≥ 1`. It is cleared 116 times in the first 400 SCMs of the first ensemble. Registering
it as a prediction registers nothing.

**§4.3, the S9 refutation bar.** *"S9-as-written REFUTED iff `N_ident + N_gain ≥ 10` at `dmin ≥ 1`,
pooled, with `D_am ≥ 20,000`."* `N_gain` alone is ≈ 0.25 per SCM, so a 4,000-SCM ensemble delivers
≈ 1,000 events — the bar is cleared by a factor of 100 in the first 40 SCMs. **S9-as-written is
refuted before X2 starts**, and P3 already foreclosed the `E_ident` half.

This is not an argument for dropping the bar; S9 *is* false. It is an argument that the bar carries
no information, so the pre-registration must be about the **radius and the rate**, which are open.

> **Mitigation M6 (binding).** Replace both bars with quantities that can come out either way:
> (i) **`r*_gain` and `r*_ident`** — the *maximum* `dmin` at which each fires, per ensemble, with the
> count at each radius 1, 2, 3, ≥4 and `disconnected` printed separately. `r*_elicit` is already
> §1.5's estimand; make it the registered statistic instead of a `≥ 10` threshold.
> (ii) **`N_gain / D_gain`** with its own denominator (A6), split by `true_in_D`.
> (iii) A registered **surprise** with teeth: `r*_gain ≥ 3` on any ensemble, or `N_gain(disconnected) ≥ 1`
> (which would be a bug, see A7).
> And register the direction explicitly: §4.4's stated surprise (`N_gain = 0` over 100,000 trials)
> is now known to be unreachable and must be struck.

Two lesser free bars, for completeness: **`Cost_{−1} = 0.0625` exactly** (1/16, arithmetic), so the
dead switch passes the cost half of the S8 bar with the best possible score; and `Q_2 = 0.0000` on
all five cells, so "S8 SUPPORTED at r = 2" is decided by `Cost_2` alone.

---

## A5. 🔴 The S5 SUPPORTED bar is unreachable with the registered ensembles — a bar set to fail

§4.1: SUPPORTED requires `N_set = 0` **and** pooled `D_am(r ≥ 1, ρ=1) ≥ 150,000`. §3's power table
projects ≈ 210,000 far **trials**. But `D_am` requires `amen(G0) ∧ amen(G')`, and in ARM A/C
`G0 = C`, so `D_am ⊆ {amen(C)}`. Measured `amen(C)` rates: `licensed` 0.130, `large` 0.528,
`k8` 0.019, `C-deep` 0.163, `C-deep8` 0.062. Applying them to §3's own projection:

Rather than chain §3's projections, I measured the realised per-SCM yields directly in the pilot
(`/tmp/x2audit/power.py`, 400 SCMs per ensemble — these per-SCM figures already fold in the `amen(C)`
rate, the Meek survival and the distance split), then multiplied by §3's own registered SCM targets:

| ensemble | §3 target SCMs | `D_am(far)` /SCM | ⇒ `D_am(far)` | `D_s5`/SCM (far ∧ misstatement ∧ connected) | ⇒ **`D_s5`** |
|---|---|---|---|---|---|
| `A-licensed` | 6,000 | 0.995 | 5,970 | 0.453 | **2,718** |
| `A-large` | 4,000 | 4.950 | 19,800 | 1.978 | **7,912** |
| `A-k8` | 2,000 | 0.257 | 514 | 0.129 | **258** |
| `C-deep` | 4,000 | 1.895 | 7,580 | 0.760 | **3,040** |
| `C-deep8` | 2,000 | 0.909 | 1,818 | 0.352 | **704** |
| **pooled** | 18,000 | | **35,682** | | **14,632** |

§3 projects ≈ 210,000 far trials; the realised `D_am` is **35,682**, a factor of **5.9** short of the
150,000 SUPPORTED bar. On the only denominator that can carry an S5 verdict (A6: far ∧ misstatement
∧ connected) it is **14,632 — below §4.1's own `< 20,000` INCONCLUSIVE floor.** So the registered
verdict is fixed before the run, at INCONCLUSIVE or SUPPORTED-WEAK depending on which denominator is
used. That is the mirror image of a free bar and corrupts the pre-registration the same way, because
the only available response after the fact is to re-scale the ensembles — choosing `N` after seeing `n`.

Worse, the composition is upside down: `A-large` (p 15–25) supplies **54 %** of `D_s5`, while the
abstract's licensed box (`A-licensed` + `A-k8`, p ≤ 10) supplies **20 %**. Pooling would license a
claim about small graphs on evidence from big ones — the exact failure §3's "never pooled *only*"
paragraph is trying to prevent.

> **Mitigation M7 (binding).** Set the thresholds from the realised yields, now, per ensemble, not
> pooled: SUPPORTED at `D_s5 ≥ 20,000` **within an ensemble**; SUPPORTED-WEAK at `≥ 3,000`;
> INCONCLUSIVE below. To reach `D_s5 ≥ 20,000` on the licensed box, `A-licensed` needs
> ≈ 20,000 / 0.453 ≈ **44,000 kept SCMs**, not 6,000 — and that is cheap, not a reason to lower the
> bar: the pilot kept 7,386 SCMs in 72 s on 12 workers, so 44,000 ≈ **7 min**. Budget it **before**
> the run. If M2 is accepted and ARM A is demoted, the S5 denominator comes from M3's
> exhaustive-sweep ARM B instead, whose base-amenability rate is 0.52 on `licensed` against ARM A's
> 0.13 — a 4× power gain from the same compute, on the arm that also has dynamic range.

---

## A6. Denominator discipline is stated in §1.4 and then broken in §4.1 — three times

§1.4 is right, and the retired 791 is the right cautionary tale. §4.1's pooled `D_am` then commits
the same class of error three more times.

1. **True statements cannot be misstatements.** ARM A/C assert both orientations (§2.1), so exactly
   50 % of trials are `true_in_D`. For those, `D` extends `M({s})`, hence `O*(M({s}))` is a **valid**
   adjustment set whenever it is defined — `E_bias ≡ 0` by Henckel–Perković–Maathuis, not by
   measurement. S5 is a claim about *misstatements*. Half the denominator behind its rule-of-three
   bound is definitionally incapable of producing the forbidden event.
2. **Disconnected statements cannot reach the query, ever.** Meek R1–R4 act on adjacent triples and
   quadruples; `pcp`, `cn`, `forb`, `pa` are all confined to the query's connected component. A
   statement in another component provably cannot fire `E_set`, `E_ident` or `E_gain`. §1.2's G7
   keeps the sentinel out of the `≥ k` buckets — but §1.4 defines `Omega(r,ρ) = {dmin(K') ≥ r}` and
   `1<<20 ≥ 1`, so the disconnected trials **are** inside the `r ≥ 1` denominator. §1.4 and G7
   contradict each other. The stakes: `large` has 1,665 disconnected statements against 3,863 at
   `dmin ≥ 2` — 30 % of that ensemble's far denominator is a definitional zero.
3. **ARM A's `D_con ≡ D_all`** (A2), so `D_con` there is not a filter and reporting it as one of four
   independent denominators overstates the discipline.

> **Mitigation M8 (binding).** The S5 rule-of-three denominator is
> `D_am ∧ (not true_in_D) ∧ (δ finite)`. Name it **`D_s5`** and print it as a fifth column on every
> S5 line, beside `N`, `D_am`, `D_con`, `D_all`. Report the bound as `3 / D_s5`. Resolve the §1.4 ↔
> G7 contradiction explicitly: `Omega(r,ρ)` is defined on **finite** `dmin` only, and `disconnected`
> is a separate row with its own `N` and `D`, never pooled.

---

## A7. The rule-of-three treats `2|U|` trials per SCM as independent. They are not

ARM A yields 12–26 trials per SCM, all sharing one `(D, C, X, Y)`. "Does a far statement move `O*`"
is a property of the graph, so the independent replicate is the **SCM**, not the trial. A cluster
bound is `3 / n_SCM`, not `3 / D_am`. On the pilot's realised numbers the two differ by an order of
magnitude or more:

| | naive `3/D_am` | cluster `3/n_SCM` | ratio |
|---|---|---|---|
| `A-licensed` (960 amenable SCMs, 6,508 far trials) | 4.6 × 10⁻⁴ | **3.1 × 10⁻³** | 6.8× |
| `A-large` (600 SCMs, 5,740 far trials) | 5.2 × 10⁻⁴ | **5.0 × 10⁻³** | 9.6× |
| §4.1's advertised `D_am ≥ 150,000` | 2.0 × 10⁻⁵ | ≈ 1.7 × 10⁻⁴ at 18,000 SCMs | 8.5× |

The abstract will print whichever number the RESULTS file offers.

> **Mitigation M9 (binding).** Report **both** bounds on every S5 line and put the **cluster** bound
> in any sentence that reaches the abstract. State the design effect explicitly
> (`D_am / n_SCM`, the mean trials per cluster).

---

## A8. 🔴 `dmin` vs `dmax`: the claim is 0 % true under one convention and 100 % true under the other

§1.2 fixes `dmin` and §4.1 relegates the `dmax` reading to a reported line with "no verdict
attached, because it is known or definitional". It is neither. Recomputed exactly, from the stored
`K` and `(x,y)` — a statement is far under `dmax` unless **both** endpoints are the query:

| ensemble | `N_set(dmin ≥ 1)` / `D_am` | `N_set(dmax ≥ 1)` / `D_am` | total `N_set(ρ=1)` | stmts that ARE the X–Y edge |
|---|---|---|---|---|
| `original` | **0** / 2,443 | **1,007** / 5,164 | 1,007 | 1,401 |
| `licensed` | **0** / 2,594 | **906** / 4,538 | 906 | 766 |
| `large` | **0** / 6,350 | **208** / 7,376 | 208 | 403 |
| `k8`/K4 | **0** / 157 | **69** / 248 | 69 | 41 |
| `k8`/K8 | **0** / 837 | **599** / 1,916 | 599 | 368 |

Not "close to the whole population" as §1.2 predicts — **exactly the whole population**. Every
`E_set` event in the corpus involves a statement with one endpoint at the query and the other at
distance ≥ 1. So `N_set(≥1)` is 0 under `min` and is 100 % of all events under `max`, on the same
data, in every ensemble. The sentence in `C-triage.tex:13` says *"no misstatement at graph-distance
at least 1 from the query"*, and "the graph-distance of a statement" is undefined for an edge —
a reviewer choosing the other reading gets the opposite verdict with no arithmetic.

The honest content of the zero is not locality. It is: **`O*` moved only when the misstated edge was
incident to `X` or `Y`** — a sharper, more useful and completely defensible claim.

> **Mitigation M10 (binding).** `dmax` is promoted from a footnote to a **headline row** of the S5
> table, with the numbers above. The paper's sentence is rewritten to name the set rather than a
> distance: *"no misstatement on an edge not incident to X or Y changed `O*` (0 of ⟨D_s5⟩); every
> observed change involved an edge incident to the query (⟨N⟩ of ⟨N⟩)."* §5.4 already notes that
> `dmin ≥ 1 ⟺ not incident to {X,Y}` — say **that**, and drop the word "graph-distance" from the
> claim entirely.

---

## A9. §1.2's justification for skeleton distance is backwards

> *"skeleton distance is a lower bound on any orientation-respecting distance, so `{s : d_skel ≥ 1}`
> is a **superset** of `{s : d_directed ≥ 1}`, and a locality claim proved on the larger set is the
> stronger claim."*

`d_skel ≤ d_dir` because skeleton walks are a superset of orientation-respecting walks. Hence
`d_skel(w) ≥ 1 ⟹ d_dir(w) ≥ 1`, i.e. `{d_skel ≥ 1} ⊆ {d_dir ≥ 1}`. The skeleton-far set is the
**smaller** one, so the skeleton claim is the **weaker** one. Harmless at the registered primary
`r = 1`, where the two sets coincide (`δ = 0 ⟺ w ∈ {X,Y}` under either walk); wrong wherever it
matters, which is `r*_elicit` (§1.5) and ARM D at `r = 2`.

> **Mitigation M11 (binding).** Correct the sentence. State that at `r = 1` the conventions coincide
> and the choice is immaterial; at `r ≥ 2` report skeleton distance and say plainly that an
> orientation-respecting distance would give a **larger** far-set, so `r*_elicit` measured on the
> skeleton is a **lower bound** on the elicitation radius.

---

## A10. 🔴 S3, S4 and S5 are one contingency table, not three findings

The abstract cites *"two thirds of false knowledge sets pass the consistency check"* (S3),
*"14.1 % of the consistent errors go on to bias the estimated effect"* (S4) and *"0 of 791"* (S5) as
separate results. Recomputed at ρ=1, pooled over distance, on `original` — the generic, `p ≤ 8`
arm those figures come from:

```
D_all = 10,548   D_con = 6,776   D_con/D_all = 0.642      <- S3's "two thirds"
N_bias = 1,007   N_bias/D_con   = 0.1486                  <- S4's "14.1%"
of those 1,007 biasing trials, 1,007 are at dmin = 0 and 0 are at dmin >= 1   <- S5's zero
```

(`licensed`: 0.676 and 0.1605. `large`: 0.828 and 0.0257.)

S4's numerator and S5's zero are the two cells of **one** partition of the same 6,776 trials. S4 is
not a property of expert error; it is `P(dmin = 0 | consistent) × L_bias(dmin = 0)` — on `original`,
0.62 × 0.24. Both factors move with `p` and degree, which is why the same quantity is 0.0257 on
`large`. Citing 14.1 % and 0-of-791 as independent support for a locality story double-counts one
measurement and hides that the headline rate is a **composition effect of the ensemble**.

Alongside it, §5.6 asks X2 to recompute the `E_bias / E_set` ratio. It is already known: **`E_set`
and `E_bias` coincide exactly in all five cells** (`ostar_changed ∧ ostar_valid` = 0 of 1,007 / 906 /
208 / 69 / 599). This belongs in §0's disclosure table, not in §5 as an open question — and it is
close to definitional, since a false statement expels `D` from the MPDAG, so a changed `O*` is
generically no longer valid.

> **Mitigation M12 (binding).** `RESULTS.md` prints S3/S4/S5 as a **single 2×2×2 table**
> (near/far × consistent/rejected × biased/not), with the decomposition
> `14.1 % = P(dmin=0 | consistent) × L_bias(dmin=0)` written out and its `p`-dependence shown across
> `original` / `licensed` / `large`. Add `E_bias ≡ E_set` to §0 as **P9**. The abstract must not cite
> S4 and S5 as two findings.

---

## A11. Conditioning masquerading as an effect — and it is large, and it is measurable

The caller flagged this as a live risk for X2. It is, in two distinct places.

**(a) Eligibility for the far stratum is anti-correlated with having a near stratum.** `amen(C)`
requires every proper possibly-causal path to start with a **directed** `X →`, which means the
query's own edges onto causal paths are already compelled and therefore not in `U(C)`. So amenable
CPDAGs systematically lack the near statements that produce all the events. Statement-distance
composition, all SCMs vs `amen(C)` SCMs only:

| ensemble | share at `dmin = 0`, all SCMs | share at `dmin = 0`, `amen(C)` only | ratio |
|---|---|---|---|
| `original` | 6,640 / 10,548 = **0.630** | 138 / 764 = **0.181** | 3.5× |
| `licensed` | 4,339 / 8,348 = **0.520** | 185 / 1,384 = **0.134** | 3.9× |
| `large` | 2,111 / 9,772 = **0.216** | 567 / 5,768 = **0.098** | 2.2× |
| `k8`/K4 | 256 / 508 = **0.504** | 5 / 96 = **0.052** | 9.7× |

ARM A conditions on exactly this population, which is why its positive-control stratum is 20–95 ×
thinner than ARM B's and why G4 is uncomputable there (`A-k8`: 20 near trials in 1,162 draws).

**(b) The near/far contrast in ARM B is computed across different SCM populations.** Splitting SCMs
by whether their `K` contains near statements, far statements, or both:

| ensemble | far-only SCMs | near-only | **both** |
|---|---|---|---|
| `original` | 98 | 461 | **2,078** |
| `licensed` | 235 | 220 | **1,632** |
| `large` | 971 | 10 | **1,462** |
| `k8`/K8 | 22 | 1 | **563** |

`large` is the worst: 971 SCMs contribute only far statements and 10 only near ones, so the
"0.203 near vs 0 far" contrast is partly a comparison between two disjoint sets of graphs. The
matched comparison is available and is not registered.

> **Mitigation M13 (binding).** Report `L_set(0)` and `L_set(≥1)` **restricted to the SCMs that
> contribute both** (`n` = 2,078 / 1,632 / 1,462 / 83 / 563), as a within-SCM matched contrast, on
> the same line as the marginal rates. Report the composition table above so the marginal contrast
> can be read as the mixture it is.

---

## A12. `E_set` and `E_ident` are competing risks, and the denominator differs sharply by stratum

`D_am` excludes trials that destroyed amenability — but destroying amenability is itself an outcome
of the perturbation, and it happens at wildly different rates near and far:

| ensemble | `D_am/D_con` at `dmin=0` | at `dmin=1` | at `dmin≥2` | at `disconnected` |
|---|---|---|---|---|
| `original` | 2,721/4,224 = **0.644** | 1,600/1,709 = 0.936 | 346/346 = **1.000** | 497/497 = **1.000** |
| `licensed` | 1,944/2,923 = **0.665** | 1,805/1,934 = 0.933 | 435/435 = **1.000** | 354/354 = **1.000** |
| `large` | 1,026/1,715 = **0.598** | 1,636/1,660 = 0.986 | 3,186/3,186 = **1.000** | 1,528/1,528 = **1.000** |
| `k8`/K8 | 1,079/1,458 = **0.740** | 794/903 = 0.879 | 30/30 = **1.000** | 13/13 = **1.000** |

At `dmin = 0`, 26–40 % of consistent perturbations are removed from the `E_set` denominator by
`E_ident`; at `dmin ≥ 2`, **not one** ever is. So `L_set` conditions on survival of a competing risk
whose hazard is a function of the very variable under study. The composite "did the analysis change
at all", `P(E_set ∪ E_ident | consistent)`, is the artefact-free quantity: `licensed` gives
1,885/2,923 = **0.645** at `dmin=0`, 129/1,934 = **0.067** at `dmin=1`, **0/435** at `dmin ≥ 2`.

> **Mitigation M14 (binding).** Register `P(E_set ∪ E_ident ∪ E_gain | consistent)` — *"the analysis
> changed in any way visible to the analyst or the oracle"* — as a **co-primary** statistic beside
> `L_set`, per stratum, with `D_con` as its denominator. It is the quantity S9 is actually about,
> and it is immune to the competing-risk conditioning.

---

## A13. `is_amenable` conflates "no causal path" with "not identified", and S9's verdict rests on the conflation

`adjust.py`: `if not paths: return False  # no causal path; excluded upstream`. In a perturbed
graph the closure can destroy every proper possibly-causal path, and the trial is then scored
`E_ident` — *"identification lost"* — when an analyst would read *"no causal path, the effect is
zero, identified by the empty set"*. That is not a loss of identification; it is a different
(wrong) answer. §4.3's S9 verdict is built on `N_ident`, which is 109 / 129 / 24 / 40 / 109 at
`dmin = 1` and is not decomposed. `E_gain` inherits the same ambiguity in reverse.

> **Mitigation M15 (binding).** Split both events at recording time, one extra call per trial:
> `E_ident_nopath` (`pcp(G') = ∅`) vs `E_ident_unamenable` (`pcp(G') ≠ ∅`, some path starts
> undirected); likewise `E_gain_paths_appeared` vs `E_gain_oriented`. Report them separately at every
> radius. If the `dmin = 1` `E_ident` counts turn out to be mostly `nopath`, S9's refutation changes
> character and the paper's sentence changes with it.

---

## A14. Oracle leakage — three channels, one of them new to X2

1. **The query is drawn on the true DAG.** `run_arm1.py:127`:
   `pairs = [... if len(possibly_causal_paths(D, a, b)) > 0]`. Inherited from the pilot, so not X2's
   doing — but §2.3 **repeats it and then maximises over it**: ARM C picks `argmax δ(X,Y)` over an
   oracle-defined candidate set. ARM C is the designated place where S5 "should die", so an
   ensemble defined by the truth is exactly where a reviewer will push.
2. **The `τ` floor is unspecified.** `run_arm1.py:34` sets `TAU_FLOOR = 1e-3` "applied at analysis"
   while `analyse_one` keeps `|τ| ≥ 1e-6`. §3 of the PREREG names neither. `τ` is the oracle total
   effect, so the choice is an oracle-conditioned selection rule on the ensemble, and it changes
   `D_am`. Two implementers will produce different denominators.
3. Declared and fine: `E_bias` and `true_in_D` use `D`; `n = ∞` uses population `Σ`. §7.4's CPDAG
   caveat (E1′ A2) carries over unchanged and is correctly stated.

> **Mitigation M16 (binding).** Define ARM C's candidate set on **`C`** (`possibly_causal_paths(C,a,b)`),
> not `D`, and report the overlap with the `D`-based set as a one-line diagnostic. Fix the floor in
> writing as `|τ| ≥ 1e-6` at draw time with `1e-3` reported as a sensitivity, matching E1′ AUDIT A8,
> and state it in §3's ensemble table rather than leaving it to the code.

---

## A15. Five of the nine integrity gates cannot fail

- **G1** (Meek idempotence) — `meek_closure` runs `while changed` to fixpoint; idempotence is
  guaranteed by the loop.
- **G2** (skeleton invariance) — `apply_background_knowledge` and `meek_closure` only ever zero one
  direction of an existing edge. Structurally guaranteed; §1.2 says as much and then calls it a gate.
- **G3** (replication) — deterministic code, same seed, same filters. It is a copy-integrity check
  and carries no information about the science. (It is also only a **prefix** match: E1′'s `large`
  kept 3,000 SCMs where §3 asks for 4,000, so only the first 3,000 overlap. Say so.)
- **G5** (dead-switch placebo) — `I_prune(−1) = {β(O₀)}`, and `est0 ∈ I_full` always, so
  `Q_{−1} ≠ P(I_full non-degenerate)` is arithmetically impossible. I verified the identity holds
  exactly on all five cells. It is a unit test, not a placebo. The real placebo already exists and
  X2 does not use it: E1′'s `members_frozen` ball (`run_arm1.py:96-98`).
- **G7 / G8** are enforced *by the writer* — they are assertions in the emitter, not tests of the
  world. Correct to have; not gates.

The genuinely informative gates are **G4** (which fails, A1), **G6** (the seed-508 anchor) and
**G9** (abort accounting — my pilot saw **0 aborts** at `p ∈ 12…20`, so `C-deep` is safe; §3.1's
guard is still right to exist for `p → 25`).

> **Mitigation M17 (binding).** Relabel G1/G2/G3/G5 as **unit tests**, run them, print them, and
> state that they cannot fail. Add one gate that can: **G10 — placebo ball.** Re-run the entire ARM D
> analysis on `members_frozen` and require `Q_r = 0` for every `r`; a non-zero there is a
> bookkeeping bug in the prune-set construction, which is the only failure mode ARM D actually has.

---

## A16. Ambiguities that would let two implementers print different numbers

Beyond the `Cost_r` aggregator (M5), the `Omega` / G7 contradiction (M8) and the `τ` floor (M16):

1. **`I_full` when the base is non-amenable.** §2.5 takes min/max "over `K'` in `B_ρ` that are
   consistent+amenable". If `M(K)` itself is non-amenable there is no `est0` and the interval may be
   empty; comparing two empty intervals is undefined. In practice `build_arm` returns early
   (`run_arm1.py:71-72`) so those SCMs are absent — but X2's ARM A has no such early return.
2. **ARM D's `ρ`.** §2.5 says `ρ` generally; §4.2 fixes `ρ = 4`; §3 gives `k8` `ρ = 1–4` while E1′'s
   `k8` ran `max_rho = 3`. My table above uses E1′'s values. Fix `ρ` per ensemble in §3.
3. **`dmin(K')` for ρ ≥ 2 is min-over-flipped**, so the far stratum at ρ ≥ 2 means *all* flips far.
   The complementary "at least one flip far" stratum is never defined, and it is the one S8's pruning
   argument is about. Define both or say which.
4. **"per-SCM roll-up printed beside it"** (§1.6) — no aggregation rule is given (any / min / max
   over statements). E1′ used `min` (`analyse_arm1.py:87`); say so.
5. **`cascade`** is recorded (§2.1) and never used in any decision rule, gate or table.

> **Mitigation M18 (binding).** Resolve 1–4 in the PREREG text before implementation. On 5: `cascade`
> is the mechanism statistic that turns a null into a finding — E1′'s archive shows Meek closure
> **firing** on far statements at 0.27 / 0.31 / 0.45 / 0.79 (`dmin ≥ 2`, fraction with `cascade > 0`,
> `original`/`licensed`/`large`/`k8`). So the honest sentence is *"0 of ⟨D_s5⟩ despite the closure
> orienting additional edges in 45 % of those trials"*, which is far stronger than a bare zero.
> Register `P(cascade > 0 | consistent)` per stratum as a reported statistic.

---

## A17. Under-powered cells whose numbers are uninterpretable

- **`k8` far strata:** `D_am` = 14 (`k8`/K4, `dmin ≥ 2`) and 30 (`k8`/K8). `A-k8` in my pilot
  produced 20 near trials and 402 far trials from 1,162 draws, off a 1.9 % `amen(C)` rate. Any `k8`
  locality statement is a statement about tens of trials. §3 already flags `k8`/K4 as "the weakest
  cell in the entire corpus"; the fix is to drop it from the verdict, not to print it.
- **`Q_1` on `original`/`licensed`: 3 events.** The SUPPORTED bar is 0.001 and the REFUTED bar 0.01;
  3/2,087 = 0.0014 has a 95 % interval of roughly [0.0003, 0.0042], which contains the SUPPORTED bar
  and comes within a factor of 2 of the REFUTED bar. At `n_SCM = 2,000` the S8 verdict at `r = 1`
  turns on whether 2 or 3 SCMs fire. Separating 0.001 from 0.01 needs `n_SCM ≈ 20,000`.
- **`C-deep8`:** 6.2 % `amen(C)`; 176 amenable SCMs from the first 2,800 draws.

> **Mitigation M19 (binding).** State each S8 verdict as a **Clopper–Pearson interval on `Q_r`**, not
> a point estimate against a threshold, and declare INCONCLUSIVE whenever the interval spans a bar.
> Either raise `n_SCM` for `original`/`licensed` ARM D to 20,000 (free — it is a re-partition, and
> the ball is already enumerated) or accept INCONCLUSIVE at `r = 1` in writing. Report `k8` cells but
> exclude any cell with `D_am < 300` from every verdict, and say so in the table.

---

## A18. What must be reported even if it embarrasses the design

- G4's realised `L_set(0)` per ensemble, against A1's table, whatever it is.
- ARM A's event counts **including the zeros at `dmin = 0`** (A2), and the `amen(C)` rate, `|cn|`
  and undirected-edges-incident-to-`cn` that explain them.
- The full `Q_r` / `Cost_r` matrix (A3), with §4.2's two failed predictions named as failed.
- Realised `D_am` and `D_s5` per ensemble against §3's projection (A5), and both bounds (A7).
- The `dmax` row (A8), at headline weight.
- The S3/S4/S5 contingency table and the `E_bias ≡ E_set` identity (A10).
- The matched within-SCM contrast (A11) and the composite change rate (A12).
- Every gate, pass or fail, before any verdict — including the four that cannot fail (A15).

---

## A19. Two things the audit endorses without reservation

§7's non-coverage list is the strongest part of the design and must survive editing intact —
especially **§7.3**, the X1 dependency: a spurious edge changes the skeleton and therefore `δ(·)`
itself, so every distance in X2 is conditional on the skeleton being right. That is a real
dependency and it is stated as one.

And **§8b**: `wiki/activities/active-claims.md:570` records b-LOAD's graph-distance falsifier as
"currently NOT firing" on the strength of `0/791`, a number retired 2026-08-14; the b-LOAD CIKM
camera-ready is **2026-08-23**, four days out. This AUDIT adds that the retirement is not merely
bookkeeping: the *replacement* number is 0 under `dmin` and 906 under `dmax` on the licensed box
(A8), so the register entry cannot be repaired by substituting a new denominator — the claim needs
new wording. Zé's decision, Zé's wording; no agent edits a live claim. E1′ AUDIT A10 (no
ρ-breakdown entry in the register at all, zero literature notes routed to it) remains open and X2
does not close it.

---

## Binding mitigations — the implementer will be told to honour these

| # | Mitigation | Blocks |
|---|---|---|
| **M1** | Per-ensemble G4 bands from A1's measured values ±0.05; G4 relabelled a reproduction check | run start |
| **M2** | ARM A demoted to mechanism/`E_gain` arm; removed from the pooled S5 denominator; add exhaustive `p ∈ {4,5,6}` lemma search | run start |
| **M3** | New primary S5 arm: ARM B operator with **exhaustive** statement sweep (base `K` size 1–4, flip each of `|U(C)|`) | run start |
| **M4** | ARM D moved to §0 as disclosure **P8** with the full `Q_r`/`Cost_r` matrix; §4.2's prediction paragraph deleted | before writing |
| **M5** | `Cost_r` aggregator fixed as the mean over SCMs + per-SCM fraction; report the `(Q_r, Cost_r)` Pareto frontier for `r ∈ {−1,0,1,2,3}` | analysis |
| **M6** | §4.4 / §4.3 bars replaced by `r*_gain`, `r*_ident`, `r*_elicit` and rates with own denominators; §4.4's "surprise" struck | before writing |
| **M7** | S5 thresholds re-set per ensemble from realised `D_s5` (20,000 / 3,000); `A-licensed` budget raised to ≈ 44,000 kept SCMs (≈ 7 min) or replaced by M3 | run start |
| **M8** | S5 denominator = `D_am ∧ ¬true_in_D ∧ finite δ`, named `D_s5`, printed as a 5th column; `Omega(r,ρ)` restricted to finite `dmin`; `disconnected` never pooled | writer |
| **M9** | Report naive **and** cluster (`3/n_SCM`) rule-of-three bounds; cluster bound is the one that may reach the abstract | analysis |
| **M10** | `dmax` row promoted to headline; claim reworded as "not incident to X or Y", dropping "graph-distance" | before writing |
| **M11** | Correct §1.2's superset/subset inversion; `r*_elicit` on the skeleton declared a **lower** bound | PREREG text |
| **M12** | S3/S4/S5 printed as one contingency table with the `14.1 %` decomposition and its `p`-dependence; `E_bias ≡ E_set` added as **P9** | analysis |
| **M13** | Within-SCM matched `L_set(0)` vs `L_set(≥1)` on the both-strata SCMs, beside the marginals | analysis |
| **M14** | `P(E_set ∪ E_ident ∪ E_gain \| consistent)` registered as co-primary, per stratum, on `D_con` | analysis |
| **M15** | `E_ident` and `E_gain` split into `nopath` / `unamenable` variants at record time | implementation |
| **M16** | ARM C candidate set defined on `C` not `D`; `τ` floor fixed in writing (1e-6 draw, 1e-3 sensitivity) | run start |
| **M17** | G1/G2/G3/G5 relabelled unit tests; add **G10**, the `members_frozen` placebo ball with `Q_r = 0` required | implementation |
| **M18** | Resolve the four ambiguities (empty `I_full`, ARM D `ρ` per ensemble, "≥1 flip far" stratum, per-SCM roll-up rule); register `P(cascade > 0)` per stratum | PREREG text |
| **M19** | Clopper–Pearson intervals on `Q_r` with INCONCLUSIVE when a bar is spanned; exclude cells with `D_am < 300` from verdicts | analysis |

---

## Provenance of every number in this file

- E1′ archive, recomputed per statement at `n = ∞`:
  `~/latent-causal/e1prime-se/results/arm1_{original,licensed,large,k8}.json` on betelgeuse, via
  `/tmp/audit_probe2.py`, `/tmp/probe3.py`, `/tmp/dmax.py`. Member index ↔ flip bijection verified:
  **0 length mismatches** over all five ensemble/arm cells.
- ARM A / ARM C pilot: `/tmp/x2audit/power.py` (400 SCMs × 5 ensembles), `power2.py` (2,000-target,
  `A-licensed`, 960 amenable SCMs, 72 s), `power3.py` (`A-large` 600, `A-k8` 22, `C-deep` 400).
  `e1prime-se/code/{graphs,adjust,scm,se}.py` copied unmodified; seeds 20260820 and 777.
- Line references verified against source: `run_arm1.py:34,46,49,71-72,127,139,157`;
  `graphs.py:190,192` (the PREREG cites `graphs.py:179-180` for the non-adjacent `MeekFail`; it is at
  **192**, and 179 is the function header — a citation slip worth fixing);
  `adjust.py:39`; `se.py:23,47`; `analyse_arm1.py:87,217`.
- No file under `~/latent-causal/e1prime-se/` or `~/latent-causal/rho-breakdown-knowledge/` was
  modified. No git commit was made.
