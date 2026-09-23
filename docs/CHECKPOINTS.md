# Checkpoints

Sources, revisions and licences for the audited causal foundation models.

**Nothing here downloads automatically.** `cfm/registry.py` resolves paths and
raises `FileNotFoundError` with instructions when a checkpoint is absent.
`scripts/download_checkpoints.sh` is a stub — fill it in per checkpoint after
reading that checkpoint's licence.

---

## Status

| Key | Model | Source | Revision | Licence | Estimands | Conditioning | Obtained |
|---|---|---|---|---|---|---|---|
| `causalpfn` | CausalPFN | TODO | TODO | TODO | TODO | TODO | ☐ |
| `causalfm` | CausalFM | TODO | TODO | TODO | TODO | TODO | ☐ |
| `dopfn` | Do-PFN | TODO | TODO | TODO | TODO | TODO | ☐ |

Fill in **before** any audit number is reported. Every column is load-bearing:

- **Revision** — these checkpoints get updated. A number produced against "the
  CausalPFN checkpoint" without a revision is not reproducible, and
  `CheckpointSpec` requires one.
- **Licence** — not uniformly permissive. The paper must state terms accurately,
  and some may restrict the comparative evaluation being run.
- **Estimands** — ATE-only checkpoints cannot produce PEHE, and the audit table
  must show a gap rather than a zero.
- **Conditioning** — which injection modes the architecture actually supports.
  Not every checkpoint accepts every mode, and forcing one that is unsupported
  produces numbers about the forcing.

---

## Layout

```
checkpoints/
├── causalpfn/
├── causalfm/
└── dopfn/
```

Gitignored in full, along with `*.pt`, `*.ckpt`, `*.safetensors`.

---

## Conditioning modes

How background knowledge reaches a model that has no explicit graph. See
`cfm/conditioning.py`.

### `none`

No knowledge. The baseline arm.

### `attention_bias`

An additive bias on attention logits between variable tokens — positive for
asserted relations, negative for forbidden ones — following the approach of
*Use What You Know*. Works with any transformer-based checkpoint without
retraining, which is why it is the default.

**The zero-scale invariant.** At `bias_scale = 0.0` this must reproduce
`mode="none"` bit for bit. If it does not, the injection is changing the forward
pass through some route other than the intended bias — a shape bug, a masking
bug, a normalisation that sees the added term — and every downstream number is
measuring that instead of the knowledge. `verify_zero_scale_identity` checks it
and the audit aborts on failure.

### `ancestral_matrix`

An `{-1, 0, 1}` matrix over variable pairs, passed to checkpoints with a
structured input slot. Closer to the model's intended interface than an
attention bias, and correspondingly less widely supported.

### `prompt`

Knowledge serialised into a textual prefix. Most fragile. The exact template
goes into the run manifest — phrasing changes results, and an unrecorded prompt
makes the number unreproducible.

---

## What the audit can and cannot conclude

**Can.** How these models' effect estimates move when conditioned on knowledge
that is consistent with the data and false — and whether that movement tracks
the classical adjustment-set pipeline's, which would indicate a shared
structural failure.

**Cannot.** That attention-bias injection *is* knowledge conditioning in the
sense the classical pipeline means. Imposing knowledge on a CPDAG has a
semantics: the forced orientations are entailed by the rules. An attention bias
has none. The model may follow it, ignore it, or overreact to it, and which of
those it does is an empirical question the audit is designed to answer rather
than assume.

State this in the paper. An audit that reads as "these models are fragile to
background knowledge" when it actually measured "these models respond to
attention biases in ways we did not characterise" would be overclaiming.

---

## Sanity gates

Both run before any sweep, and abort rather than warn.

1. **Zero-scale identity** — above.
2. **True knowledge helps** — the true-knowledge arm must beat the unconditioned
   arm. Without this gate, a flat bias-versus-δ curve is uninterpretable: it
   looks like robustness and is indistinguishable from the model not reading the
   conditioning input at all.

Report both results in the paper regardless of outcome. A checkpoint that
ignores its conditioning input is a finding about that checkpoint, not an
inconvenience to drop from the table.
