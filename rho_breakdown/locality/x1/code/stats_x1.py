"""
X1 inferential tools.  Three, and no more.

  wilson              -- every rate in RESULTS carries one.  No fraction without
                         its denominator.
  cluster_bootstrap_* -- AUDIT M1: the arms are PAIRED (same SCMs), so
                         independent-sample CI disjointness is the wrong
                         instrument for a DIFFERENCE.  Resample SCMs, not
                         statements: the 4 statements of a rho=1 ball share one
                         SCM (M13).
  cliffs_delta        -- AUDIT: the ONLY two-sample test in the run, and it is on
                         a bias MAGNITUDE distribution.  No p-value on rho*, ever.
"""
import numpy as np

Z95 = 1.959963984540054


def wilson(k, n, z=Z95):
    """Wilson score interval.  (0.0, 0.0) for n = 0 -- an empty stratum is
    printed as empty, never suppressed."""
    if n <= 0:
        return (0.0, 0.0)
    ph = k / n
    d = 1.0 + z * z / n
    c = ph + z * z / (2 * n)
    h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n))
    return (float(max(0.0, (c - h) / d)), float(min(1.0, (c + h) / d)))


def rate(k, n):
    return dict(k=int(k), n=int(n), rate=(float(k) / n if n else None),
                ci=list(wilson(k, n)))


def _boot_idx(m, B, seed, chunk=500):
    rng = np.random.default_rng(seed)
    done = 0
    while done < B:
        b = min(chunk, B - done)
        yield rng.integers(0, m, size=(b, m))
        done += b


def cluster_bootstrap_rr(num_s, den_s, num_r, den_r, B=10000, seed=20260819,
                         alpha=0.05):
    """Paired SCM-level cluster bootstrap for the RISK RATIO
        RR = (sum num_s / sum den_s) / (sum num_r / sum den_r).
    The SAME resampled SCM index set is used for both arms -- that is what makes
    it paired.  Returns (point, lo, hi)."""
    num_s = np.asarray(num_s, float); den_s = np.asarray(den_s, float)
    num_r = np.asarray(num_r, float); den_r = np.asarray(den_r, float)
    m = len(num_s)
    ps = num_s.sum() / den_s.sum() if den_s.sum() else np.nan
    pr = num_r.sum() / den_r.sum() if den_r.sum() else np.nan
    point = ps / pr if pr else np.nan
    out = []
    for idx in _boot_idx(m, B, seed):
        a = num_s[idx].sum(1); b = den_s[idx].sum(1)
        c = num_r[idx].sum(1); d = den_r[idx].sum(1)
        with np.errstate(divide="ignore", invalid="ignore"):
            out.append((a / b) / (c / d))
    out = np.concatenate(out)
    out = out[np.isfinite(out)]
    if out.size == 0:
        return float(point), float("nan"), float("nan")
    lo, hi = np.quantile(out, [alpha / 2, 1 - alpha / 2])
    return float(point), float(lo), float(hi)


def cluster_bootstrap_diff(num_s, den_s, num_r, den_r, B=10000, seed=20260819,
                           alpha=0.05):
    """Same machinery for the DIFFERENCE of two rates."""
    num_s = np.asarray(num_s, float); den_s = np.asarray(den_s, float)
    num_r = np.asarray(num_r, float); den_r = np.asarray(den_r, float)
    m = len(num_s)
    point = num_s.sum() / den_s.sum() - num_r.sum() / den_r.sum()
    out = []
    for idx in _boot_idx(m, B, seed):
        a = num_s[idx].sum(1); b = den_s[idx].sum(1)
        c = num_r[idx].sum(1); d = den_r[idx].sum(1)
        with np.errstate(divide="ignore", invalid="ignore"):
            out.append(a / b - c / d)
    out = np.concatenate(out)
    out = out[np.isfinite(out)]
    lo, hi = np.quantile(out, [alpha / 2, 1 - alpha / 2])
    return float(point), float(lo), float(hi)


def design_effect(counts, m0):
    """AUDIT M13: variance of the per-SCM event count against the binomial it
    would have if the m0 statements of a ball were independent.  > 1.25 forces
    every CI to the cluster bootstrap."""
    c = np.asarray(counts, float)
    if c.size < 2:
        return float("nan")
    p = c.sum() / (m0 * c.size)
    var_bin = m0 * p * (1 - p)
    if var_bin <= 0:
        return float("nan")
    return float(c.var(ddof=1) / var_bin)


def cliffs_delta(a, b, B=10000, seed=20260819, alpha=0.05):
    """Cliff's delta = P(A > B) - P(A < B), with a bootstrap CI.
    delta > 0 means A stochastically dominates B."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    if a.size == 0 or b.size == 0:
        return dict(delta=None, ci=[None, None], n_a=int(a.size), n_b=int(b.size))

    def d(u, v):
        u = np.sort(u)
        gt = np.searchsorted(u, v, side="left").sum()      # #(u < v)
        ge = np.searchsorted(u, v, side="right").sum()     # #(u <= v)
        n = u.size * v.size
        lt = n - ge                                        # #(u > v)
        return (lt - gt) / n

    point = d(a, b)
    rng = np.random.default_rng(seed)
    boots = np.empty(B)
    for i in range(B):
        boots[i] = d(a[rng.integers(0, a.size, a.size)],
                     b[rng.integers(0, b.size, b.size)])
    lo, hi = np.quantile(boots, [alpha / 2, 1 - alpha / 2])
    return dict(delta=float(point), ci=[float(lo), float(hi)],
                n_a=int(a.size), n_b=int(b.size),
                median_a=float(np.median(a)), median_b=float(np.median(b)),
                iqr_a=[float(np.quantile(a, .25)), float(np.quantile(a, .75))],
                iqr_b=[float(np.quantile(b, .25)), float(np.quantile(b, .75))])


def standardise(strata_s, strata_ref):
    """AUDIT M9.3: direct standardisation of arm S's rate onto arm R's hop
    distribution.  strata_* : {stratum: (k, n)}.  Weights are R's shares over the
    strata where S has data; the weights are returned so they can be printed."""
    keys = [k for k in strata_ref if strata_ref[k][1] > 0 and
            strata_s.get(k, (0, 0))[1] > 0]
    tot = sum(strata_ref[k][1] for k in keys)
    if tot == 0:
        return None
    w = {k: strata_ref[k][1] / tot for k in keys}
    val = sum(w[k] * strata_s[k][0] / strata_s[k][1] for k in keys)
    cov = sum(strata_s.get(k, (0, 0))[1] for k in keys)
    return dict(value=float(val), weights={str(k): float(w[k]) for k in keys},
                n_covered=int(cov))
