# Pre-registration — the cluster suppliers (scale and a third family)

**Written 12 September, after four batch jobs were submitted on the cluster and
before any of their answers was read.** The jobs answer the same frozen,
hashed questionnaire (540 items) with four checkpoints from the cluster's
shared model store through vLLM, one user message per item through each
model's own chat template, temperature zero, seed zero, 700-token budget:
`Qwen2.5-72B-Instruct` on real names at both presentation orders and on
scrambled names, `Qwen2.5-32B-Instruct`, `Llama-3.3-70B-Instruct` and
`gemma-3-27b-it` on real names at order zero. Their answers enter the same
cache the local suppliers used and go through the same parser, repair and
sweep. The primary supplier of the ledger stays the seven-billion model; the
chance and truthful audit arms stay matched to its edges. These rows are a
scale axis beside it, not a replacement, unless S1 and S6 both hold, in which
case the paper's Table 2 leads with the 72B row and prints the 7B row under it
as the cost-free supplier.

Nothing below is edited after `results/elicit/remote/*.jsonl` is read.

## Predictions

**S1. Scale reduces commission within a family.** The 72B model's false-claim
fraction on the open pairs is below the 7B model's, paired by network, with the
interval of the difference excluding zero. The 32B model lies between them.
Falsified if the paired interval contains zero: size buys nothing on this
substrate and the paper says the supplier's error is not a matter of scale.

**S2. Scale reduces omission.** The 72B model declines fewer open pairs than
the 7B model. Falsified if it declines at least as many.

**S3. Names carry information at scale.** L3 failed at 7B: scrambled and real
names gave the same commission. At 72B the scrambled arm's commission exceeds
the real-name arm's, paired by network, interval excluding zero. Falsified if
they are again within each other's intervals, in which case the recitation
probe is reported as returning "no recitation detected, and no domain
knowledge detected" at both sizes.

**S4. The claim certificate stays dangerous on no row.** Structural; recorded
so its being checked is on the record. Falsified by one row, which is a bug.

**S5. The unit gap does not close with scale.** The hop instrument's dangerous
rate on the 72B arm is strictly positive. The gap is a property of the
operator, not of the supplier. Falsified if the 72B arm has zero hop-dangerous
rows, which would mean a better supplier's errors happen not to sit on
propagating claims.

**S6. Scale reduces the invalid-set rate.** The committed set is invalid at
the truth on a smaller fraction of certifiable rows under the 72B model than
under the 7B model, paired by network, interval excluding zero.

**S7. Certified moves still do not order by truth.** L7 failed at 7B. The 72B
model's certified moves per row are within noise of the truthful arm's. This
is the prediction that the certificate measures the shape of the asserted set
and not its truth, and it is expected to hold at every size.

**S8. The third family is not the outlier.** The gemma model's commission lies
within the interval spanned by the two Llama and Qwen models of comparable or
larger size. Falsified if it lies outside on either side.

**S9. Source variance at scale is smaller than at 8 billion.** The disagreement
on shared edges between the 72B Qwen and the 70B Llama is below the 0.220
measured between their small relatives. Falsified if it is not.

**S10. The compelled-edge control improves with scale but the bridge stays
biased.** The 72B model's error on compelled edges is below the 7B's 0.191,
and it still sits below the 72B model's open-edge commission, so the bridge's
correction is still needed.

## Decision rules

- S4 falsified: halt and inspect before anything else is read.
- S1 and S6 both supported: the 72B row leads Table 2 and the 7B row is
  printed under it; the abstract's supplier sentence names the larger model.
- S1 supported and S6 not: scale buys accuracy and not validity; the paper
  says the certificate is what turns the first into the second.
- S3 supported at 72B and not at 7B: the recitation probe is reported per
  size, and the scrambled arm is the headline where it separates.

## What this cannot say

The CPDAG is still the oracle's on every row. A 72B model is not a domain
expert and the paper does not call it one. Four checkpoints, one run each,
temperature zero; a second seed would measure decoding variance and is not
run.

---

## Appendix A — outcomes, 12 September, after the four answer files were read

Four checkpoints, 540 items each, every one of the 2,160 answers cached under
the same key as the local suppliers'. Six conditions parsed, repaired and
swept beside the four local ones: 8,308 ledger rows on 27 networks, 3,014
certifiable deployment rows. The 72B model on real names asserts 281 of the
369 open pairs and declines 27; the 7B asserted 171 and declined 137. Six of
its 540 answers hit the token budget.

**S1 NOT SUPPORTED as registered, direction right at every size.** Commission
paired by network, 72B minus 7B, **−0.100 [−0.224, +0.017]** on 31 networks;
32B minus 7B −0.056 [−0.183, +0.071]; 72B minus 32B −0.025 [−0.112, +0.040].
Network means 0.325 at 7B, 0.233 at 72B. The order is the registered one and
no interval excludes zero. The second family repeats it: Llama 70B minus
Llama 8B −0.050 [−0.140, +0.050].

**S2 SUPPORTED.** Omission 72B minus 7B **−0.280 [−0.390, −0.164]**; the
larger model answers what the smaller one declined.

**S3 NOT SUPPORTED.** Scrambled minus real at 72B **+0.088 [−0.064, +0.238]**,
against +0.029 at 7B. The gap grows with size and is not resolved on 27
networks. As at 7B, scrambling makes the model decline, 96 against 27, more
than it makes it err. The recitation probe returns the same verdict at both
sizes: no recitation detected, and the names buy a direction, not a resolved
difference.

**S4 SUPPORTED.** The claim certificate is dangerous on no row of the six new
conditions, 3,014 certifiable deployment rows in all.

**S5 SUPPORTED, thinly.** The hop instrument is dangerous on **one** row of
the 72B arm (barley, `frspdag -> dg25`, one false claim, claim radius 1, hop
radius 2), three on its order replicate, three on scrambled names, two on the
32B, one on Llama 70B, and **none** on gemma. The gap needs a false claim that
propagates; a better supplier makes fewer, so the count falls from four to
one, and it does not close except on the one arm with the fewest certifiable
rows per claim.

**S6 NOT SUPPORTED.** The committed set is invalid at the truth on **0.298
[0.171, 0.431]** of the 72B's certifiable rows, 85 of 285, against 0.288
[0.155, 0.435] at 7B, 66 of 229. The 32B sits at 0.201, Llama 70B at 0.247,
gemma at 0.308. Scale bought assertions, and with them certifiable rows, 285
against 229, and a lower false fraction per assertion; it did not buy a lower
rate of invalid sets, because a set is invalid when one wrong claim sits on
the query and more claims give that more chances. This is the third branch of
the decision rules in substance, with the caveat that S1 did not separate:
scale buys accuracy and not validity, and the certificate is what stands
between the two.

**S7 NOT RESOLVED, as predicted.** The 72B certifies 0.846 [0.618, 1.103]
moves per row, the truthful arm 1.086 [0.767, 1.469], overlapping; the
truthful arm is matched to the 7B's edges, so this is a comparison across
edge sets and is reported only as such.

**S8 NOT SUPPORTED.** gemma minus Qwen 72B **+0.112 [+0.029, +0.213]**, gemma
minus Llama 70B +0.073 [−0.052, +0.200]. The third family is the outlier on
the high side against the 72B. It also never reaches 93 of the 369 pairs and
invents four, which the two larger models do not.

**S9 SUPPORTED in direction.** Source disagreement at scale 0.182 [0.079,
0.301] on 236 shared edges, against 0.220 [0.118, 0.333] at eight billion;
overlapping. Instrument disagreement at 72B is **0.047 [0.019, 0.079]**
against 0.160 at 7B: the larger model is far more stable under presentation
order, which is the cleaner statement.

**S10 NOT SUPPORTED, both halves, and the second half fails the useful way.**
The 72B's compelled-edge error is 0.232 [0.121, 0.358], above the 7B's 0.191:
size did not improve it. But control minus open-edge commission at 72B is
**+0.004 [−0.153, +0.159]** against −0.153 [−0.324, +0.009] at 7B: the
truth-free bridge is unbiased on this supplier. L9's failure was the 7B's,
not the bridge's; the bridge is reported as usable at 72B and biased at 7B,
per supplier.

**Decision.** S1 did not separate and S6 failed, so the 7B row stays the
primary supplier of Table 2 and the 72B row is printed as the scale
replicate beside it. The paired instrument ladder holds on every new arm:
against the free rule the claim certificate is as safe and sharper by 0.19
to 0.27, four wins of four; against the hop radius it is safe where the hop is
dangerous on one to three rows and less sharp by 0.03 to 0.12.
