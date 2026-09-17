# Pre-registration — the source-overlap probe (A6-0c, the contamination screen)

**Written 12 September, after the ledger and the scale pre-registrations were
scored, after the answer key and the prompt were fixed in
`experiments/source_overlap_probe.py`, and before any model answer to the
source prompt was read.** The two local suppliers had not yet been asked when
this was written; the cluster suppliers answer later from the exported
prompts.

Nothing below is edited after `results/elicit/source_overlap.json` is written.
Outcomes go in Appendix A.

---

## 1. Why a second probe

The scrambled-name replicate is the recitation probe the protocol asked for,
and it did not separate: L3 in `results/ledger/PREREGISTRATION.md` and S3 in
`results/elicit/PREREGISTRATION_SCALE.md` both returned "no recitation
detected, and no domain knowledge detected". That verdict is an inference from
invariance. A source that recites would collapse under scrambled names; a
source that does not collapse is not shown to be reciting. It is also not
shown to be ignorant of the benchmark: a model can know where a network comes
from and still answer a causal-order question from the names in front of it.

This probe asks the direct question. Each supplier is shown, cold, the
variable list of one network and nothing else, and is asked to name the
publication the network comes from and to reproduce the causal query that
publication declares: the exposure, the outcome and the adjustment set. A
supplier that can do both for a network has the benchmark in its weights,
whatever the scrambled arm says. A supplier that can do neither has, at most,
the domain knowledge the names carry.

## 2. Fixed before the run

**Networks.** The 32 networks the questionnaire asked about, which are the 33
of `results/elicit/questionnaire.json` less `pathfinder`, whose only chain
component is above the questionnaire's cap and which no supplier was asked
about. Eight are the applied dagitty papers, two are dagitty textbook diagrams
(`mediator`, `paths`), the other 22 are the discrete and Gaussian bnlearn
benchmarks.

**The prompt.** One per network, rendered by `render` in the probe script.
It lists the variables under the display names the questionnaire used (the
file's own names; no network needed renaming), in a seeded permutation, and
states how many variables the network has. A network with more than 60
variables is shown a seeded sample of 60 that always contains every variable
the questionnaire showed the supplier, and the prompt says it is a sample.
The prompt carries no domain sentence, no states, and no name: the
questionnaire's domain sentence is withheld on purpose, so the probe is
conservative. A miss here does not prove the supplier could not have
recognised the network with the domain sentence in front of it; a hit here
is recognition from the names alone. The answer grammar is fixed JSON with
`DECLINE` allowed on every field. Temperature zero, seed zero, the same
cache and the same key as every other elicitation.

**The answer key.** Written into the script before any answer was read.

- For the eight applied dagitty networks the source is the author-year label
  the file name carries and the dagitty example collection uses: Acid & de
  Campos 1996, Didelez et al. 2010, van Kampen 2014, Polzer et al. 2012,
  Schipf et al. 2010, Sebastiani et al. 2005, Shrier & Platt 2008,
  Thoemmes 2013. The declared query is the file's own `exposure`, `outcome`
  and `adjusted` annotations, parsed by `benchmarks.dagitty` and already
  reported in `results/frame/declared_rows.csv`. Two files declare an
  adjustment set (Didelez: Age, Smo; Thoemmes: s2); the other six declare
  none, so their adjustment key is the empty set.
- For the bnlearn networks the source is the reference the bnlearn
  repository documents for it, transcribed on 12 September from
  `bnlearn.com/bnrepository`. The file headers in
  `results/axisa3/networks/example_models` say `network unknown` and the
  acquisition manifest names only the pgmpy sdist, so neither adds anything;
  the repository page is the documented origin. These networks declare no
  query; they are marked `no declared query` and only the source is scored.
- The two dagitty textbook diagrams in the questionnaire, `mediator` and
  `paths`, declare a query (X on Y; E on D) and have no publication. Only
  the query is scored.
- Three bnlearn networks have no publication on the repository page:
  `win95pts` (no attribution), `barley` (a project, credited to Kristensen
  and Rasmussen, no paper), `magic-irri` (a conference talk by Scutari,
  2016). `win95pts` is unscorable on the source; the other two are scored
  against the named people and flagged as not single-source.

**The SINGLE-SOURCE flag.** True when the network has one canonical
publication in the key; false for `win95pts`, `barley`, `magic-irri`,
`mediator` and `paths`. `munin1` to `munin4` are subnetworks of MUNIN and
carry the MUNIN reference with the flag set. So 27 of the 32 carry the flag.

**Scoring.**

- *Source match*: some surname in the answer equals some surname in the key
  after case-folding and stripping diacritics, and the answered year is
  within one of the key's. Author-only and year-only matches are printed
  beside it. A `DECLINE` is neither a match nor a fabrication; a non-declined
  answer that does not match is a fabrication.
- *Benchmark named*: the answered dataset name contains the network's usual
  name. Informational only: `water`, `child`, `link`, `insurance`, `diabetes`
  are ordinary words, so this column is not in the verdict.
- *Exposure match*, *outcome match*: exact equality with the declared
  variable after case-folding.
- *Adjustment Jaccard*: the answered set against the declared set, with
  empty against empty scored 1 and a `DECLINE` left unscored.
- *Query match*: exposure and outcome both match, and Jaccard at least 0.5
  where the key declares a non-empty adjustment set.
- *Chance level* per network for the query, from the number `m` of variables
  listed: `1/m` for the exposure alone, `1/(m(m-1))` for the ordered pair.
- *Verdict*: `recited` when source and query both match; `partial` when the
  source matches, or the exposure or the outcome matches, but not both
  halves; `no` otherwise. A network with no declared query can be at most
  `partial`; a network with no documented source can be at most `partial`.

## 3. Predictions

**P1. The primary supplier names the source on fewer than a quarter of the
scorable networks.** Fewer than 8 of the 29 source-scorable networks under
`qwen2.5:7b-instruct` have a source match. Falsified at 8 or more, in which
case the benchmarks are in the model's weights at the level of a citation and
the scrambled arm's "not detected" is re-read as "not detected by that
instrument".

**P2. Where it names the source, it does not know the query.** Among the
networks where the primary supplier's source matches and a query is declared,
the query matches on fewer than half. Falsified at half or more. If the
denominator is empty the prediction is unscorable and is reported as such.

**P3. Recitation is detected on at most one network per local supplier.** The
verdict `recited` fires on zero or one of the 32 networks for each of the two
local models. Falsified at two or more for either. Any network with the
verdict is flagged in the ledger table and the primary supplier's headline
rates are re-printed without it as a sensitivity row.

**P4. The query is guessed from the names, not from the paper.** On the ten
networks that declare a query, the outcome is matched without a source match
on at least three, above the chance total of about one, because the outcome
of an applied paper is the clinical endpoint in the list (Mortality, Injury,
T2DM, TCI). This is domain knowledge available from the variable names and is
not recitation, which is why the verdict requires the source as well.
Falsified if outcome matches without the source are at chance, at most one,
in which case the names carry less than expected and the L3 reading is the
stronger for it.

**P5. The second family is in the same band.** `llama3.1:8b` also names the
source on fewer than 8 of the 29, and the two local models agree on the
correctly named set on fewer than half of the union of their hits. Falsified
if the second family reaches 8, or if the two agree on half or more of the
hits, which would mean the recognisable networks are the same across families
and the property is the benchmark's, not the model's.

**P6. Correct namings fall on the famous benchmarks and on none of the
applied papers.** Every source match of either local model is on a bnlearn
network with a documented publication, the networks that have appeared in
textbooks and software for two decades (`asia`, `alarm`, `sachs`, `child`,
`insurance` and their kind); none of the eight applied dagitty papers is
named by author and year by either local model. Falsified by one applied paper
named correctly, which is the case that matters, because those are the only
networks whose declared query could be recited.

**P7. Wrong answers are fabrications, not declines.** For each local model,
more than half of the non-declined source answers are wrong. A
seven-billion-parameter model at temperature zero produces a citation rather
than a `DECLINE`. Falsified if the majority of non-declined answers are
right, or if the model declines on more than half of the networks; either
outcome is reported.

## 4. Decision rules

- P1 or P5 falsified: the source-recognition probe joins the scrambled arm in
  the paper's contamination paragraph, the named networks are listed, and the
  ledger's per-arm rates are re-summarised without them as a sensitivity row.
- P3 falsified for a network: the network's `single_source` flag and verdict
  are printed in the ledger's per-network table, and the row is excluded
  from the primary supplier's headline in a sensitivity row.
- P6 falsified: the applied paper named is the case to inspect by hand,
  answer text and key side by side, before anything else is concluded.
- P4 falsified: the names buy less than assumed; the reading of L3 tightens
  to "no recitation and little domain knowledge", already the reading.
- No prediction here changes the primary supplier, the frame or any radius.
  The probe is a screen on the interpretation of the elicited arms, not on
  their rows.

## 5. What this cannot say

The prompt withholds the domain sentence the questionnaire carried, so a
miss is a lower bound on what the model might recognise in situ. The key for
the applied papers is an author-year label; a model that cites the right
paper under a co-author's name and the right year is a match, a model that
cites the right paper with the year off by two is not, and both cases are
printed. The bnlearn networks have no declared query, so for 22 of the 32
the verdict can only be `partial` or `no`. Four cluster models answer the
same prompts later from the exported file; their predictions are these,
applied unchanged, and they are scored in the same appendix when read.

---

## Appendix A — outcomes, 12 September, after both local suppliers had answered

### A.0 Two instrument defects found on the first read, and what was read when

**The permutation seed was salted.** The variable order in the prompt was
seeded, as section 2 says, but through the interpreter's `hash()` of the
network name, which Python salts per process. Every run therefore rendered a
different prompt for the same network and asked the supplier again instead
of reading the cache. The first pass asked the primary supplier and crashed
in its summary code; the second asked both suppliers in full, with zero
cache hits; the third, on a stable hash, was stopped at nine networks on
instruction. The scored answers are the **second pass**, 32 networks per
supplier, re-parsed from the responses it stored verbatim by
`--rescore`. The prompt texts of that pass are not on disk, only their
hashes and the responses; the exported `source_overlap_prompts.jsonl` is the
stable-seed rendering, which lists the same variables in another order on
the 20 networks shown in full and a different top-up sample on the twelve
shown as a sample of 60. The wording is identical. **Undone:** asking the two
local suppliers on the stable prompts, so that every supplier answers the
same text; the script does it unchanged in about fifteen minutes.

**The grammar invited a bare `DECLINE`.** `"year": <year>` with "answer
DECLINE for any field" produced `"year": DECLINE`, which is not JSON, in 27
of the primary supplier's 32 answers. They are clean declines and are parsed
after the token is quoted. Three of its answers and two of the second
family's ran into the 400-token budget on a runaway adjustment set on the
MUNIN family and were read field by field; the source fields precede the
list in the grammar and were intact. In the earlier scoring these 30 answers
counted as unparsed and their declines as absent; no source or query field
changed.

The first pass's raw answers were read to diagnose the parse failures before
either fix; the predictions above were written before that and were not
touched. Read from the cache records, the first pass's answers were the
same as the second's on every field checked: the source declined on the
same networks, the same two fabrications for `asia` and `hepar2`, the same
exposure and outcome on every applied paper. The presentation order moved
nothing, which is an unregistered observation and is reported as one.

### A.1 Predictions scored

Per supplier, 32 networks, 29 source-scorable, 10 with a declared query,
27 single-source. Full table in `SOURCE_OVERLAP.md`.

**P1 SUPPORTED.** The primary supplier names the source on **0 of 29**. It
declines on 27 and answers on two, both wrong at confidence 0.8 to 1.0:
`asia` is attributed to "Neapolitan & Richardson, 2000, *Learning belief
networks*" and `hepar2` to a Finnish liver cohort paper the key does not
contain. Neither author nor year matches anywhere.

**P2 UNSCORABLE.** No source match falls on a query-scorable network for
either supplier; the denominator is empty and the prediction is reported as
such.

**P3 SUPPORTED.** `recited` on **no network** for either supplier. The
sensitivity row is not needed.

**P4 SUPPORTED, for both.** Without a source match the outcome is matched
on **4 of 10** by each supplier against a chance total of about 1: the
primary picks Mortality, Injury, `y` and `Y`; the second family picks TCI,
T2DM, Injury and `Y`. The exposure is matched without the source on 1 and
0 against a chance total of about 1, the pair on 1 and 0 against 0.17, and
the single pair is `mediator`, the four-letter textbook diagram whose X and
Y are conventional. On the applied papers neither supplier picks the
declared exposure once: Smoking for ToothLoss, Coach or ContactSport for
WarmUpExercises, `e0` or `e1` for `x`, Smo for HRT, PA for TT. The names
give away which variable is the endpoint and nothing about which variable
the authors set out to vary. That is domain reading, not recitation, and it
is why the verdict requires the source.

**P5 SUPPORTED on the count, vacuous on the agreement.** The second family
names **0 of 29**; the union of hits is empty, so the agreement clause has
nothing to measure.

**P6 SUPPORTED, vacuously on its first half.** There is no source match to
fall anywhere; no applied paper is named by author and year by either
supplier, which is the half that mattered.

**P7 SPLIT BY FAMILY.** The second family never declines: 29 of 29
non-declined answers are wrong, mean stated confidence 0.45 on the 20 that
gave one, two of them right on the year by coincidence ("Storey et al.
2005" for `ecoli70`, "Bresnahan & Pakes 1997" for `insurance`) and none on
an author. Supported there. The primary supplier declines on 27 of 29, so
the falsifier's second clause fires and the prediction is **falsified for
the primary supplier in the decline direction**: a seven-billion model of
this family says it does not know rather than inventing a citation, and
fabricates only where it thinks it knows (`asia`, `hepar2`).

**Informational, not in any verdict.** The benchmark's usual name is given
once by the primary supplier, `asia`, which is also one of the listed
variables, so it is read off the list; and four times by the second family,
every one an ordinary-word coincidence ("Pima Indians Diabetes Database",
"autoinsurance", "Printer Troubleshooting Data", "E. coli gene expression
dataset"), the last two apt descriptions of the variable names rather than
a citation.

### A.2 Decision

No decision rule fires. For the protocol's single-source column: 27 of 32
networks carry the flag, and the verdict on every one of them is `no` or
`partial`, never `recited`. The contamination screen returns, by the direct
probe, what the scrambled arm returned by invariance: neither local supplier
holds these benchmarks at the level of a citation, neither reproduces a
declared query, and what the names buy is the endpoint on four of ten
applied networks. The reading of L3 and S3 stands and is now supported from
both sides.

Still open: the cluster suppliers, which answer `source_overlap_prompts.jsonl`
and enter the cache through `elicit_import.py`, to be scored by the same
script with `--models`; and the two local suppliers on the stable prompts,
as A.0 says. The predictions above apply to those runs unchanged.
