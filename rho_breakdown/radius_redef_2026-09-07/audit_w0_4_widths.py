import sys; sys.path.insert(0,'.')
from audit_w0_rk import *

Z95 = 1.959963985

def ols(data, idx, x, zs, y):
    cols = [idx[x]] + [idx[v] for v in sorted(zs)]
    A = np.column_stack([np.ones(len(data))] + [data[:, c] for c in cols])
    b, *_ = np.linalg.lstsq(A, data[:, idx[y]], rcond=None)
    resid = data[:, idx[y]] - A @ b
    s2 = resid @ resid / (len(data) - A.shape[1])
    XtXi = np.linalg.inv(A.T @ A)
    return float(b[1]), float(np.sqrt(s2 * XtXi[1, 1]))

def methods_for(pr, rng, n=2000):
    dag, cp, x, y, Z, K = pr['dag'], pr['cpdag'], pr['x'], pr['y'], pr['z'], pr['K']
    sem = random_sem(dag, rng)
    truth = sem.true_total_effect(x, y)
    sigma = sem.covariance()
    nodes = list(dag.nodes); idx = {v: i for i, v in enumerate(nodes)}
    data = rng.multivariate_normal(np.zeros(len(nodes)), sigma, size=n)
    cache = {}
    def est(zs):
        k = frozenset(zs)
        if k not in cache: cache[k] = ols(data, idx, x, k, y)
        return cache[k]
    out = {}
    e, s = est(Z); out['point_Z'] = (e, e, s)                       # elicited set
    fibre = pr['fibre']
    ms = {d: m_of(d, K) for d in fibre}
    for t in (0, 1, 2, 99):
        S = [d for d in fibre if ms[d] <= t]
        if not S: S = [d for d in fibre if ms[d] == min(ms.values())]
        vals = [est(optimal_adjustment_set_dag(d, x, y)) for d in S]
        lo = min(v[0] for v in vals); hi = max(v[0] for v in vals)
        se = max(v[1] for v in vals)
        out['blanket' if t == 99 else f'S_t={t}'] = (lo, hi, se)
    return out, truth

def calibrate(rows, key):
    """smallest scale c giving >=0.95 coverage; then mean width"""
    lo = np.array([r[0][key][0] for r in rows]); hi = np.array([r[0][key][1] for r in rows])
    se = np.array([r[0][key][2] for r in rows]); tr = np.array([r[1] for r in rows])
    def cov(c): return np.mean((tr >= lo - c*Z95*se) & (tr <= hi + c*Z95*se))
    if cov(0) >= 0.95: a = 0.0
    else:
        a, b = 0.0, 1.0
        while cov(b) < 0.95 and b < 1e6: b *= 2
        if cov(b) < 0.95: return cov(1.0), float('inf'), float('inf')
        for _ in range(60):
            mid = (a+b)/2
            if cov(mid) >= 0.95: b = mid
            else: a = mid
        a = b
    return cov(1.0), a, float(np.mean(hi - lo + 2*a*Z95*se))

for nf in (0, 1, 2):
    rng = np.random.default_rng(777 + nf)
    probs, _ = gen(rng, 400, n=7, p=0.40, k_claims=3, n_false=nf)
    rng2 = np.random.default_rng(999 + nf)
    rows = [methods_for(pr, rng2) for pr in probs]
    print(f"\n=== truly-false claims in K: f={nf}   (n={len(rows)} problems, n_obs=2000) ===")
    print(f"{'method':10s} {'raw cov(c=1)':>13s} {'calib factor':>13s} {'MEAN WIDTH @95%':>17s}")
    for k in ('point_Z','S_t=0','S_t=1','S_t=2','blanket'):
        c1, c, w = calibrate(rows, k)
        print(f"{k:10s} {c1:13.3f} {c:13.2f} {w:17.4f}")
