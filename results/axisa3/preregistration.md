# Pre-registration — session 6 (real graphs)

**Written before the first breakdown radius is computed, and not edited
afterwards.** Directional predictions were committed earlier still, in
`prediction_before_looking.md`, before the descriptive table existed; this file
can see that table, because it must fix the network list.

---

## 1. The network list

39 graphs parsed from the pgmpy 1.0.0 sdist (sha256 in
`results/axisa3/networks/acquisition_manifest.json`), in four families:

- **24 BIF benchmark networks**: alarm, andes, asia, barley, cancer, child,
  diabetes, earthquake, hailfinder, hepar2, insurance, link, mildew, munin,
  munin1, munin2, munin3, munin4, pathfinder, pigs, sachs, survey, water,
  win95pts.
- **11 dagitty DAGs from applied papers**: Acid 1996, Didelez 2010, Kampen 2014,
  Polzer 2012, Schipf 2010, Sebastiani 2005, Shrier 2008, Thoemmes 2013, plus
  confounding, mediator, paths.
- **4 bnjson gene/plant networks**: arth150, ecoli70, magic-irri, magic-niab.

**Sought and not usable, recorded so the gap is auditable:**

- **M-bias** — parsed but **excluded**: it uses `<->` bidirected edges, so it is
  an ADMG with latent confounding, not a DAG. Dropping those edges would make it
  indistinguishable from a real 3-node DAG in the results.
- **The bnlearn.com repository** — unreachable from this machine
  (`TLSV1_ALERT_PROTOCOL_VERSION`), so `pgmpy.utils.get_example_model` and the
  canonical `.bif`/`.dsc` downloads were unavailable. Guessed GitHub raw URLs for
  the same files returned 404. No network was obtained by any route other than
  the checksummed sdist.
- **CPDAGs published directly by applied papers** (as opposed to DAGs) were not
  located in any fetchable package. The 11 dagitty files are the closest
  available substitute and are labelled as applied-paper DAGs, not as published
  CPDAGs.

**Six networks are fully compelled** (0 undirected edges in the CPDAG): cancer,
confounding, earthquake, mildew, pigs, survey. They admit no background-knowledge
perturbation at all. They are **kept in the denominator** of every
network-level count and reported as producing zero admissible instances, because
dropping them would overstate the framework's applicability — which is itself one
of the questions.

## 2. The (X, Y) selection policy — fixed, mechanical, applied without exception

With thousands of admissible pairs it would be trivial to produce any answer one
likes, so the policy is stated here and applied by code, not by choice.

**Admissibility.** An ordered pair `(X, Y)`, `X ≠ Y`, is admissible iff:

1. `Y` is a descendant of `X` in the true DAG (a causal path exists);
2. `X` lies in, or is adjacent to, an undirected component of the CPDAG;
3. the instance survives the degeneracy gate — the five structural checks of
   `synth.runner.gate`, via the `fast_gate` equivalent introduced in session 5
   (differentially tested at 46,800 cases, 16 one-directional disagreements,
   0.034%); and
4. `O(G₀)` is identified and is genuinely valid at `G₀`.

**Enumeration.** For networks where the admissible set is small enough, **all
admissible pairs** are taken. Where it is not, pairs are drawn by a **fixed,
seeded, deterministic sample** over the sorted pair list, with the sampling rate
and the seed recorded per network in the results, and the **realised** count
reported alongside the total admissible count. Sampling is never adaptive and
never depends on any computed radius.

**Domain-documented pairs.** The dagitty files carry `exposure` and `outcome`
annotations. Where present, that pair is reported as a **separate labelled
stratum** — never merged into the main distribution, never the headline.

**Reporting weight.** Every distribution is reported **per network first**, and
any aggregate is given twice: weighted by pair, and weighted by network (each
network contributing equally). Session 5's clearest procedural lesson was that
pooling moved a headline from 79.6% to 65.3%.

## 3. Knowledge coverage

Swept at **1.0, 0.5 and 0.25**, matching session 5's design, since `|K_{G₀}|` is
the other binding term in `r_val = min(s, |K_{G₀}|)` and a full-coverage analyst
is the unrealistic case. Coverage selects which component edges the analyst
asserts, by a deterministic seeded rule.

## 4. Hypotheses

### H11 — the factual question. Is separation ≥ 2 rare in practice?

**Baseline to beat:** Erdős–Rényi gave **85.8%** of measurable separations at
`s = 1`, maximum 4.

**Prediction** (already committed in `prediction_before_looking.md`, before the
descriptive table): **15–40% of measurable pairs at `s ≥ 2`**, and a maximum
separation above 4 in at least one network.

**Falsified if** the fraction at `s ≥ 2` is within a few points of the
Erdős–Rényi rate — i.e. if real structure is no better at separating the
adjustment set from the treatment than a random graph.

**Reported as:** the full separation distribution per network, plus the
undefined-separation rate (no `O(G₀)` member in `X`'s component), which is a
**status and never a number**.

### H12 — generality of the law `r_val = min(s, |K_{G₀}|)`

Session 5 was explicit that this is a law of a designed family, close to true by
construction. Real networks are the first genuine test.

**Prediction:** agreement **above 70% but below 100%**.

**Commitment:** every exception is **characterised, not just counted**. An
exception is more informative than the rule here, because it identifies routes
to failure the spine construction excluded. If agreement is exactly 100% I will
treat that as suspicious and look for a way the measurement could be inheriting
the construction rather than testing it.

### H13 — undirectedness

Already measured in Phase 2 and committed before this file: median undirected
fraction **10.7%**, six fully compelled networks, maximum component **85**. What
remains for Phase 3 is what that implies for admissible-instance counts per
network.

### H14 — tractability on real structure

Session 4's envelope reached component 11 on synthetic graphs and never failed,
so its ceiling was "where I stopped". Real components reach **85**.

**Prediction:** the hybrid fails somewhere between component 15 and 40.

**Commitment:** report which networks and component sizes completed, which hit
the cap, and where the transition happened. Censored instances carry
`wall_until_timeout_s` and are never counted as measurements.

## 5. Analysis commitments

- **Per network first**, aggregates second and labelled with their weighting.
- **Sentinels are statuses.** Undefined separation and `UNREACHED` radii are
  never averaged and never plotted as numbers.
- **Assumption status on every radius.** Both search legs are upward searches,
  exact under Conjecture 2, which rests on Anti-Exchange Case B — verified, not
  proved. Its error direction inflates radii, which is **favourable to this
  project's hypothesis**, and that caveat stays visible.
- **The scope limit, stated in the report's own voice and not buried.** The
  CPDAG here is computed from the true DAG: the oracle-CI idealisation this
  project has assumed throughout. Real discovery on finite data returns a
  different, usually sparser skeleton, and missing weak edges is precisely what
  breaks back-door blocking. This session measures the idealised case only.
- **If separation ≥ 2 turns out rare**, the follow-up is fixed in advance: test
  whether it is *predictable* from properties a practitioner can see — component
  size, degree of `X`, position of `Y`. A radius that is usually 1 but whose
  exceptions are identifiable ahead of time is still a useful diagnostic, and
  that is the salvage path.
