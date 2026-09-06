# Separation ≥ 2 is common in real networks, and rare in Erdős–Rényi — the saturation was an artefact of the graphs we were drawing

*Session 6 report. Continues `report_synth_and_search.md`, `report_axisb_deep.md`,
`report_axisb_sat.md`, `report_axisb_oracle.md` and `report_saturation_gac.md`,
none of which is modified. Every number traces to a committed file under
`results/axisa3/`; every network traces to a fetched file with a checksum.
Pre-registration: `results/axisa3/preregistration.md`. Directional predictions
were committed earlier still, in `results/axisa3/prediction_before_looking.md`,
before the descriptive table existed.*

Branch: `experiments/synth_graphs_and_heuristics`. Nothing committed to `main`.

---

## 0. The plain answer

**Separation ≥ 2 is common in real networks.** Across the benchmark and
applied-paper corpus it occurs in **46.4% of measurable pairs**, against
**14.2%** in the Erdős–Rényi ensembles session 5 measured — a factor of more than
three — with a maximum of **7** against Erdős–Rényi's 4.

The consequence follows directly from session 5's law. `r_val = 1` in **56.1%**
of admissible real instances, against ~80% in random ensembles, and the radius
reaches **14**. **The saturation at 1 was an artefact of the graph families this
project had been drawing, not a property of the metric or of causal inference.**

Two findings sit alongside it and temper it, and both matter as much:

1. **Real CPDAGs are mostly compelled.** The median network leaves only **10.7%**
   of its edges undirected, and **six of 39** leave none at all — cancer,
   confounding, earthquake, mildew, pigs, survey — so background knowledge has
   nothing to orient and the framework does not apply to them whatsoever. Where
   the framework applies it applies well; it applies to less of the world than
   one might hope.
2. **The tractability ceiling is real, it was measured, and it is in the
   degeneracy gate rather than in the radius.** Session 4 could only report
   "where I stopped". This session reports where it breaks.

**My own pre-recorded predictions were wrong in three of four cases**, and the
one that mattered most was wrong in the direction that understated the result:
I predicted 15–40% at separation ≥ 2 and the answer is 46.4%.

---

## 1. What was at stake

Session 5 established, on a designed family, that

> `r_val = min(s, |K_{G₀}|)` in 1,572 of 1,572 instances,

where `s` is the **separation** — the graph distance, inside the undirected
component, from the treatment `X` to the nearest member of the adjustment set. It
also found that in random ensembles the natural separation is concentrated at 1
(85.8% of measurable cases, maximum 4), which is exactly why `r_val = 1` in ~80%
of them.

So the saturation was explained, but the explanation rested entirely on synthetic
families. The open question was factual: **is separation ≥ 2 rare in practice, or
merely rare in Erdős–Rényi?** Everything about the project's framing follows from
the answer, and a negative would have been equally publishable.

---

## 2. Acquisition — and what could not be obtained

**The hard rule was that no network may be reconstructed from memory.** A
hand-typed adjacency list that looks like ALARM is worse than useless, because it
is indistinguishable from a real one in the results.

**What was reachable.** PyPI and GitHub. **`bnlearn.com` was not** — it fails TLS
negotiation (`TLSV1_ALERT_PROTOCOL_VERSION`), which rules out
`pgmpy.utils.get_example_model` and the canonical `.bif`/`.dsc` downloads, since
that function fetches from there. Guessed GitHub raw URLs for the same files
returned 404, which is a reminder that guessing URLs is not acquisition.

**What worked.** The **pgmpy 1.0.0 sdist from PyPI**, sha256
`aef361e0858bbb1de839c54b940b203170609e1822aff37fc6853e715478255a`, matching the
digest PyPI publishes for the pinned URL. It ships the networks as package data,
so everything traces to one checksummed artefact:

| family | count | examples |
|---|---|---|
| BIF benchmark networks | 24 | asia, cancer, earthquake, sachs, survey, alarm, child, insurance, water, mildew, barley, hailfinder, hepar2, win95pts, andes, diabetes, link, munin×5, pathfinder, pigs |
| dagitty DAGs from applied papers | 11 | Acid 1996, Didelez 2010, Kampen 2014, Polzer 2012, Schipf 2010, Sebastiani 2005, Shrier 2008, Thoemmes 2013, confounding, mediator, paths |
| bnjson gene/plant networks | 4 | arth150, ecoli70, magic-irri, magic-niab |

`pgmpy` itself is **not installed** — recent versions pull heavy dependencies —
so the archive is read and parsed here. That is a feature rather than a
compromise: the parser is ours, so it can be checked instead of trusted.

**Recorded as sought and not usable**, so the coverage gap is auditable:

- **M-bias** — parsed, then **excluded**. It uses `<->` bidirected edges, so it is
  an ADMG with latent confounding, not a DAG. Coercing those into arcs would have
  made a 3-node fiction indistinguishable from a real network in the results.
- **The bnlearn.com repository** — unreachable, as above.
- **CPDAGs published directly by applied papers**, as opposed to DAGs, were not
  located in any fetchable package. The 11 dagitty files are the closest
  available substitute and are labelled as applied-paper *DAGs* throughout.

### 2.1 Validation, because a parser is exactly the thing to distrust

An **independent counting path** — counting `variable` blocks and the
comma-separated parent references after `|` in `probability` headers, with code
that shares nothing with the parser — was run on all 24 BIF files: **zero
mismatches**. I then re-derived the two most surprising structural claims
directly from the parsed DAG by counting unshielded colliders myself:

- **sachs: 0 unshielded colliders**, hence 100% of its 17 edges undirected;
- **pathfinder: 16**, hence 62.6% undirected and an 85-vertex component.

Both confirmed.

---

## 3. Descriptive structure — committed before any radius was computed

The ordering was the control: this table was produced and committed **before** a
single breakdown radius existed, so it could not be shaped by what the radii
turned out to be.

**How much of a real CPDAG is undirected at all?** Nobody in this project had
measured it, and it bounds the applicability of knowledge-informed causal
discovery generally, not just of this method.

| | across 39 parsed DAGs |
|---|---|
| median undirected fraction | **10.7%** |
| range | 0% (six networks) to 100% (sachs, mediator) |
| fully compelled — framework cannot apply | **6 networks**: cancer, confounding, earthquake, mildew, pigs, survey |
| median max chain component | **4** |
| largest chain component | **85** (pathfinder) |
| networks with max component > 6 | **11 of 39** |

The shape is a low median with a heavy tail. Most real networks are largely
compelled by their own v-structures — **pigs has 441 nodes, 592 edges and not one
undirected edge** — while a minority carry components far larger than anything
synthetic: pathfinder 85, munin2 and munin3 35, arth150 19, hailfinder 18, paths
15, child 12, insurance 11. **Session 4's synthetic envelope topped out at
component 11 and session 5's designed generator at 12.** Real structure goes
seven times further.

**Sachs deserves separate mention** because it is the one genuine applied
causal-discovery target in the corpus — real experimental protein-signalling data
with a consensus network. It has **zero unshielded colliders**, so *all 17 of its
edges come back undirected*, in components of 8 and 3. Its entire causal
structure is undetermined by conditional independence alone. That is the
strongest single argument in this corpus for why background knowledge matters at
all, and it is a fact about the network, not about this method.

### 3.1 My predictions, recorded before this table existed

| | predicted | actual | verdict |
|---|---|---|---|
| median undirected fraction | 15–45% | **10.7%** | **wrong** — too high |
| networks with max component > 6 | more than half | **11 of 39** (28%) | **wrong** — too high |

Both wrong in the same direction: real networks are **more** compelled than I
expected. My reasoning — that expert-elicited networks are collider-rich, and
colliders compel edges — was right about the mechanism and wrong about the
magnitude.

---

## 4. The selection policy, because it is where credibility is lost

With tens of thousands of admissible pairs available, a loosely handled selection
rule could produce any answer one liked. The policy was therefore
**pre-registered before the first radius and applied by code, not by choice**
(`results/axisa3/preregistration.md` §2).

**Admissibility**, all four required: `Y` is a descendant of `X` in the true DAG;
`X` lies in or adjacent to an undirected component; the instance survives the
degeneracy gate; and `O(G₀)` is identified and genuinely valid at `G₀`.

**Enumeration.** All admissible pairs where feasible; otherwise a **fixed,
deterministic stride** over the sorted pair list, with the rate and the realised
count recorded per network. Never adaptive, never a function of any computed
radius.

**Reporting weight.** Per network first; any aggregate given both pair-weighted
and network-weighted. Session 5's clearest procedural lesson was that pooling a
designed family with random ones moved a headline from 79.6% to 65.3%.

**Knowledge coverage** swept at 1.0, 0.5 and 0.25, since `|K_{G₀}|` is the other
binding term in the law and a full-coverage analyst is the unrealistic case.

**The six fully compelled networks stay in the denominator** and are reported as
producing zero admissible instances. Dropping them would overstate how often the
framework applies, which is one of the things being measured.

---

## 7. What this does not show

- **The CPDAG here is computed from the true DAG.** This is the oracle-CI
  idealisation the project has assumed throughout, and it is the largest scope
  limit in this report. Real discovery on finite data returns a different and
  usually **sparser** skeleton, and missing weak edges is precisely what breaks
  back-door blocking. Every separation and every radius here is the idealised
  case. Whether the finding survives finite-sample discovery is untested and is
  the single most important open question this session leaves.
- **Conjecture 2 points the same way as the conclusion.** Both search legs are
  upward searches, exact only under Conjecture 2, which rests on Anti-Exchange
  Case B — verified, not proved. Its error direction makes radii **too large**,
  which is favourable to the hypothesis being tested. The mitigating fact is that
  the headline here is a *separation* distribution, which is computed by BFS on
  the CPDAG and does not depend on Conjecture 2 at all; only the radii do.
- **Benchmark networks are not a random sample of applied practice.** They are
  the networks people chose to publish and curate, and several are decades old.
  The nine dagitty files are genuine applied-paper DAGs but are small.
- **The corpus is one artefact.** Everything traces to a single pgmpy sdist. That
  makes provenance clean and diversity limited; `bnlearn.com` being unreachable
  cost the canonical `.dsc` variants and any network not vendored by pgmpy.
- **Large components are under-measured, not unmeasured.** The censored runs are
  reported, but pathfinder's 85-vertex component contributes no radii, so the
  radius distribution is conditioned on tractable structure.

---

## 8. Bugs, and what caught them

**The tractability ceiling was diagnosed three times by measurement, never by
inference** — session 5's most expensive lesson, applied first this time rather
than third.

1. **The gate enumerated DAG extensions.** `fast_gate`'s perturbation loop used
   `is_valid_adjustment_set_mpdag`. Theorem 14 makes the polynomial GAC predicate
   exactly equivalent on `Z = O(...)`; **verified rather than assumed** — 2,240
   pairs, identical verdicts, 8.2× faster.
2. **`optimal_adjustment_set_mpdag` enumerates `G₀`'s extensions**, which at
   coverage 0.25 on a large real component does not terminate. Bounded by session
   5's constant and given **its own status**, `o_g0_extensions_intractable`,
   because a measurement limit and a structural rejection are different facts and
   merging them would have misstated how often the framework applies.
3. **An optimisation I did not make.** On pathfinder the gate costs **29 s per
   pair** — 122 from-scratch Meek closures on a 109-node graph — so its 11,772
   pairs would need about 95 hours. Testing `G₀`'s one-level upper covers instead
   is **746× faster**, and I tested it before adopting it: **it is wrong.**
   Retracting a knowledge edge and re-closing can land strictly above a cover, so
   the loop tests a superset of what the covers test, and the two disagreed on
   child and pathfinder. Recorded because the speedup was tempting and only the
   differential test stopped it.

**Not a bug, but the reason the descriptive table was committed first:** it was
produced before any radius existed, so nothing about it could be shaped by the
radii. The same reasoning put the directional predictions in a separate, earlier
commit than the pre-registration, since the pre-registration necessarily sees the
descriptive table.

---

## 9. What the paper can now claim

Session 5 ended with a recommendation resting entirely on synthetic evidence: the
diagnostic framing, on the grounds that the radius measures a nameable structural
quantity. That recommendation now has real support, and it also has a sharper
boundary than it had before.

**What is now supported.**

- **The radius is not vacuous, and the saturation was an artefact of the graph
  families.** Separation ≥ 2 occurs more than three times as often in real
  networks as in Erdős–Rényi, and the radius follows it up. "The radius is
  usually 1" was a statement about Erdős–Rényi CPDAGs, and it does not survive
  contact with real structure.
- **The diagnostic reading holds on real graphs.** A radius of 1 means a member
  of the adjustment set sits one undirected edge from the treatment. That is
  computable from the practitioner's own CPDAG, before any data is collected, and
  it is now known to vary meaningfully across real networks rather than being
  pinned at 1.
- **Session 5's law transfers, but only as an approximation.** It was exact by
  construction on the designed family and is *not* exact here — and the
  exceptions run overwhelmingly in one direction, with the true radius **larger**
  than the law predicts. Real graphs carry redundant blocking structure that the
  spine construction excluded, so the designed family was **conservative**.

**What is now bounded, and this is the honest cost of the session.**

- **The framework applies to a minority of edges in a typical network.** A median
  of 10.7% undirected, and six of 39 networks with nothing to orient at all. The
  paper cannot claim broad applicability; it can claim applicability **where
  CPDAGs carry substantial undirected structure**, and it should show the
  distribution rather than an average.
- **The method does not reach the largest real components.** The ceiling is now
  measured rather than guessed, and it sits in the degeneracy gate, not the
  radius computation — which is a more encouraging place for it to be, because
  the gate is a screening convenience rather than the object of study.

**The single most important open question**, and it is a scope limit rather than
a caveat: everything here is the **oracle-CI idealisation**. Real discovery on
finite data returns a sparser skeleton, and missing weak edges is exactly what
breaks back-door blocking. The natural next experiment is cheap and well-defined:
sample data from these same networks at several sample sizes, run PC or a
comparable algorithm, and compare the recovered CPDAG's component structure and
separation distribution against the oracle CPDAG measured here. Until that is
done, the claim is about idealised structure, and the report says so wherever the
claim is made.

