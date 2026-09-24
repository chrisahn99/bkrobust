
## Follow-up (coordinator request): separation scoped to the exact 2,422-unit panel

**(1) Pooled, exact panel.** Restricting to the exact 2,422-unit/24-network
panel (status=="ok", radius and AUC_frac both defined -- same cleaning as
`tau_for_stratum`), then further to the subset where `separation` is defined:
**n=1,389, 20 networks** (separation is undefined on 1,033 of the 2,422 units).
tau_b(separation, AUC_frac) = **+0.284 [-0.029, 0.492]**. Paired r_val -
separation on this exact subset: delta = **+0.039 [-0.182, 0.284]**,
P(delta>0)=0.670 -- not significant. (Numerically identical to the earlier
Task 3 result: the internal missing-value cleaning in `tau_for_stratum` /
`paired_cluster_bootstrap` already converges on this same 1,389-unit set
regardless of whether the 2,422-unit filter is applied first; this run
confirms it explicitly against the paper's exact panel.)

**(2) Within-state, matched to the 1,213 separation-defined units.** Of the
92 within-state units (1,973 units, 18 networks) used for the paper's
tau_b=0.42 headline, restricting to the 1,213 units/15 networks/79 states
where separation is also defined:
- r_val within-state tau_b on these same 1,213 units: **+0.678 [0.459, 0.774]**
  (higher than on the full 1,973-unit set, +0.419).
- separation within-state tau_b (same units): **+0.491 [-0.135, 0.598]**
  (as before).
- Paired within-state bootstrap, r_val - separation, network-cluster,
  10,000 resamples: delta = **+0.187 [0.066, 0.895]**, P(delta>0)=0.998 --
  r_val significantly outranks separation within-state on this matched set.
- Agreement: of 1,213 units where both are defined, r_val == separation on
  **679 (56.0%)**.

Full detail: `task5_followup.json`.
