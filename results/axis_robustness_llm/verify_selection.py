"""
verify_selection.py

Audits the LLM-elicited-K survival/ranking sweep in results/axis_robustness_llm/
for selection contamination: entry into the ranking sample requires
status == "ok" (the elicited K Meek-closes to a usable G0 and determines an
optimal adjustment set Z*), and conditions that assert more claims (higher n_k)
plausibly clear that bar more often. Since n_k is itself one of the ranked
predictors, the pooled and even the "balanced panel" strata may be selected on
the predictor being tested.

This script performs five checks against results/axis_robustness_llm/analysis_units.csv
(5400 rows, one per (condition, network, X, Y) unit; n_k is populated for every
row regardless of status, which is what makes checks 1-2 possible):

  1. Quantify the selection: per-condition ok-share vs mean/median n_k, both
     across conditions and within network.
  2. Decompose status by n_k bucket (0, 1-2, 3-5, 6-10, 11+).
  3. Re-derive the fully balanced panel (ok in all 8 real-naming conditions)
     from analysis_units.csv directly, and recompute Kendall tau-b with a
     network-cluster bootstrap (resample networks with replacement, take all
     rows of each drawn network, recompute tau-b; 10000 resamples, seed 0).
  4. Build an intermediate "partially balanced" panel (ok in >= 6 of 8 real
     conditions) and rerun the same taus, to see whether the n_k effect's
     apparent significance tracks the degree of balancing (a diagnostic for
     the selection story) or is stable (evidence it is not just selection).
  5. Paired bootstrap difference in tau (pooled stratum) between radius and
     each baseline predictor, using the same network-cluster resampling, to
     ask whether radius is shown to beat the baselines or merely also excludes
     the n_k_zero cases.

Convention: a radius value of -1 would be the UNREACHED sentinel and must be
excluded from any correlation/averaging. This dataset (results/axis_robustness_llm)
in fact contains no -1 sentinel in `radius` (checked explicitly below) -- that
sentinel appears in the real-network sweep referenced in project memory, not
here -- but the exclusion is still applied defensively.

Run:
  PYTHONPATH=src .venv/bin/python results/axis_robustness_llm/verify_selection.py \
      > results/axis_robustness_llm/VERIFY_SELECTION.txt
"""

import numpy as np
import pandas as pd
from scipy import stats

DATA_PATH = "results/axis_robustness_llm/analysis_units.csv"
N_BOOT = 10000
SEED = 0
PREDICTORS = ["radius", "n_k", "k_g0", "shd_truth"]
ENDPOINT = "AUC_frac"


def hr(title=""):
    print("\n" + "=" * 100)
    if title:
        print(title)
        print("=" * 100)


def load():
    df = pd.read_csv(DATA_PATH)
    n_sentinel = (df["radius"] == -1).sum()
    print(f"[sanity] rows with radius == -1 (UNREACHED sentinel): {n_sentinel}")
    if n_sentinel:
        print("  -> excluding these from all radius correlations/averages, per instructions.")
        df = df[df["radius"] != -1].copy()
    return df


def clean_radius(df):
    """Drop -1 sentinel rows defensively before any radius-touching computation."""
    if "radius" in df.columns:
        return df[df["radius"] != -1].copy()
    return df


# ---------------------------------------------------------------------------
# Task 1
# ---------------------------------------------------------------------------

def task1(df):
    hr("TASK 1: Is ok-share correlated with mean n_k across conditions?")

    def cond_table(sub):
        g = sub.groupby("condition").agg(
            n_units=("status", "size"),
            n_ok=("status", lambda s: (s == "ok").sum()),
            mean_n_k=("n_k", "mean"),
            median_n_k=("n_k", "median"),
        )
        g["ok_share"] = g["n_ok"] / g["n_units"]
        return g

    all10 = cond_table(df).sort_values("ok_share")
    print("Per-condition ok-share and n_k, all 10 elicitation conditions (real+scrambled):")
    print(all10.to_string(float_format=lambda v: f"{v:.4f}"))

    real8 = cond_table(df[df["naming"] == "real"]).sort_values("ok_share")
    print("\nSame, restricted to the 8 real-naming conditions:")
    print(real8.to_string(float_format=lambda v: f"{v:.4f}"))

    for label, table in [("all 10 conditions", all10), ("8 real-naming conditions", real8)]:
        sp = stats.spearmanr(table["ok_share"], table["mean_n_k"])
        kt = stats.kendalltau(table["ok_share"], table["mean_n_k"])
        print(f"\nAcross {label} (n={len(table)}): "
              f"Spearman rho={sp.correlation:.4f} (p={sp.pvalue:.4f}); "
              f"Kendall tau-b={kt.correlation:.4f} (p={kt.pvalue:.4f})")

    hr("TASK 1b: within-network -- does higher n_k predict ok within a fixed network?")
    print("Unit-level: within each network, Kendall tau-b between n_k and is_ok (0/1),")
    print("pooling all conditions' units placed in that network (holds the network fixed,")
    print("so removes between-network confounding of both ok-rate and typical n_k).\n")

    rows = []
    for net, g in df.groupby("network"):
        is_ok = (g["status"] == "ok").astype(int)
        if is_ok.nunique() < 2 or g["n_k"].nunique() < 2:
            rows.append((net, len(g), np.nan, np.nan, g["n_k"].mean(), is_ok.mean()))
            continue
        kt = stats.kendalltau(g["n_k"], is_ok)
        rows.append((net, len(g), kt.correlation, kt.pvalue, g["n_k"].mean(), is_ok.mean()))
    resdf = pd.DataFrame(rows, columns=["network", "n_units", "tau_b", "p_value", "mean_n_k", "ok_share"])
    resdf = resdf.sort_values("tau_b")
    print(resdf.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    valid = resdf.dropna(subset=["tau_b"])
    n_pos = (valid["tau_b"] > 0).sum()
    print(f"\nNetworks with usable tau_b: {len(valid)}/{len(resdf)}")
    print(f"Of those, networks where n_k correlates POSITIVELY with being ok: {n_pos}/{len(valid)}")
    print(f"Median within-network tau_b: {valid['tau_b'].median():.4f}")
    print(f"Mean within-network tau_b:   {valid['tau_b'].mean():.4f}")

    # also condition-level version within network, as an alternate framing
    hr("TASK 1c: within-network, condition-level -- correlate a condition's mean n_k with its ok-share, per network")
    rows2 = []
    for net, g in df.groupby("network"):
        cg = g.groupby("condition").agg(mean_n_k=("n_k", "mean"), ok_share=("status", lambda s: (s == "ok").mean()))
        if cg["mean_n_k"].nunique() < 2 or len(cg) < 3:
            rows2.append((net, len(cg), np.nan, np.nan))
            continue
        kt = stats.kendalltau(cg["mean_n_k"], cg["ok_share"])
        rows2.append((net, len(cg), kt.correlation, kt.pvalue))
    resdf2 = pd.DataFrame(rows2, columns=["network", "n_conditions", "tau_b", "p_value"]).sort_values("tau_b")
    print(resdf2.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    valid2 = resdf2.dropna(subset=["tau_b"])
    print(f"\nNetworks with usable tau_b: {len(valid2)}/{len(resdf2)}; "
          f"median tau_b={valid2['tau_b'].median():.4f}, mean tau_b={valid2['tau_b'].mean():.4f}")


# ---------------------------------------------------------------------------
# Task 2
# ---------------------------------------------------------------------------

def task2(df):
    hr("TASK 2: status by n_k bucket (all 5400 rows, all conditions)")
    d = df.copy()
    bins = [-0.5, 0.5, 2.5, 5.5, 10.5, np.inf]
    labels = ["0", "1-2", "3-5", "6-10", "11+"]
    d["n_k_bucket"] = pd.cut(d["n_k"], bins=bins, labels=labels)

    ct = pd.crosstab(d["n_k_bucket"], d["status"])
    ct = ct[["ok", "optimal_set_undefined", "o_g0_extensions_intractable", "n_k_zero"]]
    print("Counts:")
    print(ct.to_string())

    row_pct = ct.div(ct.sum(axis=1), axis=0) * 100
    print("\nRow percentages:")
    print(row_pct.round(1).to_string())

    print("\nok-share by n_k bucket:")
    print((ct["ok"] / ct.sum(axis=1)).round(4).to_string())

    print("\noptimal_set_undefined share by n_k bucket:")
    print((ct["optimal_set_undefined"] / ct.sum(axis=1)).round(4).to_string())

    chi2, p, dof, _ = stats.chi2_contingency(ct)
    print(f"\nChi2 test of independence (status x n_k bucket): chi2={chi2:.2f}, dof={dof}, p={p:.3e}")

    osu_share = ct["optimal_set_undefined"] / ct.sum(axis=1)
    peak_bucket = osu_share.idxmax()
    print(f"\noptimal_set_undefined peaks at n_k bucket = '{peak_bucket}' "
          f"(share={osu_share.max():.4f}); it is {'HIGHER' if osu_share.index.get_loc(peak_bucket) <= 1 else 'LOWER'}"
          f" at small n_k than at large n_k.")
    ok_share = ct["ok"] / ct.sum(axis=1)
    print(f"ok-share by bucket: {dict(ok_share.round(4))}")
    direction = "RISES" if ok_share.iloc[-1] > ok_share.iloc[0] else "FALLS"
    print(f"ok-share {direction} from the smallest to the largest n_k bucket "
          f"({ok_share.iloc[0]:.4f} -> {ok_share.iloc[-1]:.4f}).")


# ---------------------------------------------------------------------------
# Task 3 & 4 helpers
# ---------------------------------------------------------------------------

def build_panel(df, min_real_conditions_ok):
    """Return (panel_df, triples, n_real_conditions) for units 'ok' in at least
    min_real_conditions_ok of the real-naming conditions, restricted to real-naming
    ok rows belonging to such a triple."""
    real = df[df["naming"] == "real"]
    n_real_conditions = real["condition"].nunique()
    ok = real[real["status"] == "ok"]
    grp = ok.groupby(["network", "x", "y"])["condition"].nunique()
    triples = grp[grp >= min_real_conditions_ok].index
    idx = pd.MultiIndex.from_frame(ok[["network", "x", "y"]])
    mask = idx.isin(triples)
    panel = ok[mask].copy()
    return panel, triples, n_real_conditions


def cluster_bootstrap_tau(panel, predictor, endpoint, n_boot=N_BOOT, seed=SEED):
    """Cluster bootstrap over networks: each resample draws len(distinct networks)
    networks with replacement, takes ALL rows of each drawn network, recomputes
    tau-b. Returns (point_tau, ci_lo, ci_hi, n_valid, n_degenerate, boot_taus_array)."""
    sub = panel.dropna(subset=[predictor, endpoint])
    point = stats.kendalltau(sub[predictor], sub[endpoint]).correlation

    networks = sub["network"].unique()
    by_net = {net: g for net, g in sub.groupby("network")}
    rng = np.random.default_rng(seed)
    n_net = len(networks)

    boot_taus = []
    n_degenerate = 0
    for _ in range(n_boot):
        draw = rng.choice(networks, size=n_net, replace=True)
        parts = [by_net[net] for net in draw]
        resampled = pd.concat(parts, ignore_index=True)
        if resampled[predictor].nunique() < 2 or resampled[endpoint].nunique() < 2 or len(resampled) < 4:
            n_degenerate += 1
            continue
        tau = stats.kendalltau(resampled[predictor], resampled[endpoint]).correlation
        if np.isnan(tau):
            n_degenerate += 1
            continue
        boot_taus.append(tau)

    boot_taus = np.array(boot_taus)
    n_valid = len(boot_taus)
    if n_valid == 0:
        return point, np.nan, np.nan, n_valid, n_degenerate, boot_taus, len(sub), n_net
    ci_lo, ci_hi = np.percentile(boot_taus, [2.5, 97.5])
    return point, ci_lo, ci_hi, n_valid, n_degenerate, boot_taus, len(sub), n_net


def report_panel_taus(panel, label, n_real_conditions_required, n_triples):
    print(f"\n{label}: n_units={len(panel)}, n_triples={n_triples}, "
          f"n_networks={panel['network'].nunique()} "
          f"(ok in >= {n_real_conditions_required} of 8 real-naming conditions)")
    print(f"Networks: {sorted(panel['network'].unique())}")
    results = {}
    for pred in PREDICTORS:
        p = clean_radius(panel) if pred == "radius" else panel
        point, lo, hi, n_valid, n_deg, boot_taus, n_sub, n_net = cluster_bootstrap_tau(p, pred, ENDPOINT)
        results[pred] = (point, lo, hi, n_valid, n_deg, n_sub, n_net)
        print(f"  {pred:12s} tau_b={point:+.4f}  boot CI=[{lo:+.4f}, {hi:+.4f}]  "
              f"(n_units_used={n_sub}, n_networks={n_net}, "
              f"valid_resamples={n_valid}/{N_BOOT}, degenerate={n_deg}/{N_BOOT})")
    return results


def task3(df):
    hr("TASK 3: Re-derive the fully balanced panel and its taus")
    panel, triples, n_real = build_panel(df, min_real_conditions_ok=8)
    print(f"n_real_conditions (real-naming conditions found) = {n_real}")
    print(f"Fully balanced triples (ok in ALL {n_real} real-naming conditions): {len(triples)}")
    print(f"Expected from task brief: 31 triples / 248 units / 7 networks.")
    print(f"Observed: {len(triples)} triples / {len(panel)} units / {panel['network'].nunique()} networks.")
    match = (len(triples) == 31 and len(panel) == 248 and panel['network'].nunique() == 7)
    print(f"MATCH: {match}")

    results = report_panel_taus(panel, "Fully balanced panel (8/8)", 8, len(triples))

    print("\nComparison to the analysis_tau.csv panel_balanced_real_naming point estimates given in the brief:")
    brief = {"radius": 0.499, "n_k": 0.333, "k_g0": 0.380}
    for pred, expected in brief.items():
        point = results[pred][0]
        print(f"  {pred:8s}: recomputed tau_b={point:+.4f} vs brief {expected:+.3f} "
              f"(diff={point-expected:+.4f})")

    n_net = panel["network"].nunique()
    print(f"\nWith only {n_net} network clusters, the bootstrap draws only {n_net} networks with "
          f"replacement each time -- there are at most {n_net}**{n_net} distinct resamples' worth of "
          f"cluster-composition, and a resample that happens to omit a network with unusual (predictor, "
          f"endpoint) values, or that draws a degenerate single-network composition, can swing tau a lot.")
    for pred in PREDICTORS:
        n_deg = results[pred][4]
        print(f"  {pred:8s}: {n_deg}/{N_BOOT} resamples degenerate/NaN "
              f"({100*n_deg/N_BOOT:.2f}%) -- {'a substantial share; CI should be read as approximate.' if n_deg > 0 else 'none.'}")
    return panel, results


def task4(df):
    hr("TASK 4: Intermediate control -- panel balanced in >= 6 of 8 real-naming conditions")
    panel6, triples6, n_real = build_panel(df, min_real_conditions_ok=6)
    results6 = report_panel_taus(panel6, "Partially balanced panel (>=6/8)", 6, len(triples6))

    hr("TASK 4b: trend across degrees of balancing (8/8 -> >=7/8 -> >=6/8 -> pooled/no balancing)")
    panel8, triples8, _ = build_panel(df, min_real_conditions_ok=8)
    panel7, triples7, _ = build_panel(df, min_real_conditions_ok=7)
    pooled_real = df[(df["naming"] == "real") & (df["status"] == "ok")]

    rows = []
    for name, panel in [
        (">=8/8 (fully balanced)", panel8),
        (">=7/8", panel7),
        (">=6/8", panel6),
        ("pooled (no balancing, all real-naming ok units)", pooled_real),
    ]:
        for pred in PREDICTORS:
            p = clean_radius(panel) if pred == "radius" else panel
            sub = p.dropna(subset=[pred, ENDPOINT])
            if len(sub) < 4 or sub[pred].nunique() < 2:
                rows.append((name, pred, len(panel), panel["network"].nunique(), np.nan, np.nan))
                continue
            tau = stats.kendalltau(sub[pred], sub[ENDPOINT]).correlation
            rows.append((name, pred, len(panel), panel["network"].nunique(), tau, len(sub)))
    trend = pd.DataFrame(rows, columns=["panel", "predictor", "n_units", "n_networks", "tau_b", "n_used"])
    print(trend.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    print("\nFocus on n_k's point estimate and whether its bootstrap CI (from task 3/4's cluster "
          "bootstrap) crosses zero, across balancing levels:")
    for name, panel in [(">=8/8", panel8), (">=7/8", panel7), (">=6/8", panel6)]:
        p = panel
        point, lo, hi, n_valid, n_deg, *_ = cluster_bootstrap_tau(p, "n_k", ENDPOINT)
        crosses = "CI crosses 0" if lo < 0 < hi else "CI excludes 0"
        print(f"  n_k @ {name:8s} (n_units={len(panel):4d}, n_networks={panel['network'].nunique()}): "
              f"tau_b={point:+.4f}, CI=[{lo:+.4f},{hi:+.4f}] -> {crosses}")

    print("\nVERDICT (report honestly whichever way it goes): as the panel is loosened from fully "
          "balanced (8/8) toward no balancing (pooled), see whether n_k's tau_b and CI stability move "
          "monotonically -- if tighter balancing (which most directly matches the selection-on-n_k "
          "story, since only conditions/triples that all clear the bar survive) produces a WEAKER or "
          "more uncertain n_k effect, that is evidence the pooled n_k correlation is inflated by "
          "selection; if it stays similar or strengthens, that is evidence against the selection story "
          "dominating.")


# ---------------------------------------------------------------------------
# Task 5
# ---------------------------------------------------------------------------

def task5(df):
    hr("TASK 5: Does radius beat the baselines in the pooled stratum, or merely also exclude n_k=0?")
    pooled = df[(df["naming"] == "real") & (df["status"] == "ok")].copy()
    print(f"Pooled stratum (real-naming, status==ok): n_units={len(pooled)}, "
          f"n_networks={pooled['network'].nunique()}")

    networks = pooled["network"].unique()
    n_net = len(networks)
    by_net = {net: g for net, g in pooled.groupby("network")}
    rng = np.random.default_rng(SEED)

    baselines = ["n_k", "k_g0", "shd_truth"]

    # point estimates (radius sentinel-cleaned defensively; none present here)
    pooled_r = clean_radius(pooled)
    point = {}
    for pred in ["radius"] + baselines:
        p = pooled_r if pred == "radius" else pooled
        sub = p.dropna(subset=[pred, ENDPOINT])
        point[pred] = stats.kendalltau(sub[pred], sub[ENDPOINT]).correlation
    print("Point tau_b (pooled, real-naming, ok):")
    for k, v in point.items():
        print(f"  {k:10s}: {v:+.4f}")

    for baseline in baselines:
        diffs = []
        n_degenerate = 0
        for _ in range(N_BOOT):
            draw = rng.choice(networks, size=n_net, replace=True)
            parts = [by_net[net] for net in draw]
            resampled = pd.concat(parts, ignore_index=True)
            resampled_r = clean_radius(resampled)

            sub_r = resampled_r.dropna(subset=["radius", ENDPOINT])
            sub_b = resampled.dropna(subset=[baseline, ENDPOINT])
            if (sub_r["radius"].nunique() < 2 or sub_r[ENDPOINT].nunique() < 2 or len(sub_r) < 4 or
                    sub_b[baseline].nunique() < 2 or sub_b[ENDPOINT].nunique() < 2 or len(sub_b) < 4):
                n_degenerate += 1
                continue
            tau_r = stats.kendalltau(sub_r["radius"], sub_r[ENDPOINT]).correlation
            tau_b = stats.kendalltau(sub_b[baseline], sub_b[ENDPOINT]).correlation
            if np.isnan(tau_r) or np.isnan(tau_b):
                n_degenerate += 1
                continue
            diffs.append(tau_r - tau_b)

        diffs = np.array(diffs)
        n_valid = len(diffs)
        point_diff = point["radius"] - point[baseline]
        if n_valid == 0:
            print(f"\nradius - {baseline}: point diff={point_diff:+.4f}; bootstrap fully degenerate, no CI available.")
            continue
        lo, hi = np.percentile(diffs, [2.5, 97.5])
        crosses = lo < 0 < hi
        verdict = ("CI INCLUDES ZERO -> radius is NOT shown to beat this baseline"
                   if crosses else
                   "CI EXCLUDES ZERO -> radius tau_b is significantly higher than this baseline's")
        print(f"\nradius - {baseline}: point diff={point_diff:+.4f}, "
              f"cluster-bootstrap CI=[{lo:+.4f}, {hi:+.4f}] "
              f"(valid={n_valid}/{N_BOOT}, degenerate={n_degenerate}/{N_BOOT})")
        print(f"  -> {verdict}")


def main():
    df = load()
    task1(df)
    task2(df)
    _, _ = task3(df)
    task4(df)
    task5(df)
    hr("DONE")


if __name__ == "__main__":
    main()
