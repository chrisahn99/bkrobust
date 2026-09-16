"""
iSCM parameterisation (Ormaniec, Sussex, Lorch, Rothfuss, Krause,
"Standardizing Structural Causal Models", arXiv:2406.11601).

Each variable is standardized *internally*, during ancestral sampling, using the
statistics of a large reference sample; the standardization constants (mu_j,
sd_j) then become FROZEN parts of the SCM.  Consequence for the linear case: the
frozen SCM is again a linear SEM with coefficients w_ij/sd_j and noise variance
sigma_j^2/sd_j^2, so Sigma and the true total effect are available in closed
form, with no Monte-Carlo error.

Control arm: the naive linear SEM sampler ("naive"), which is the thing
Reisach/Tami (arXiv:2303.18211) show leaks the causal order through marginal
variance.
"""
import numpy as np

from graphs import topological_order

N_REF = 100_000  # reference sample used to fix the iSCM standardization constants


def _draw_weights(D, rng, lo=0.5, hi=2.0):
    W = np.zeros(D.shape)
    idx = np.argwhere(D == 1)
    signs = rng.choice([-1.0, 1.0], size=len(idx))
    mags = rng.uniform(lo, hi, size=len(idx))
    for (i, j), s, m in zip(idx, signs, mags):
        W[i, j] = s * m
    return W


def make_linear_iscm(D, rng, noise_scale=(0.5, 1.5)):
    """
    Returns dict with the FROZEN linear SEM: x = A^T x + c + e.
    A[i,j] = w_ij / sd_j ; e_j ~ N(0, (sigma_j/sd_j)^2).
    Intercepts are irrelevant for effects/covariances and are dropped.
    """
    p = D.shape[0]
    W = _draw_weights(D, rng)
    sigma = rng.uniform(*noise_scale, size=p)
    order = topological_order(D)

    # reference ancestral sample -> standardization constants (as in iSCM)
    Xref = np.zeros((N_REF, p))
    mu = np.zeros(p)
    sd = np.ones(p)
    for j in order:
        pa = np.flatnonzero(D[:, j])
        raw = Xref[:, pa] @ W[pa, j] + rng.normal(0, sigma[j], size=N_REF)
        mu[j] = raw.mean()
        sd[j] = raw.std()
        Xref[:, j] = (raw - mu[j]) / sd[j]

    A = W / sd[None, :]          # effective coefficients
    omega = (sigma / sd) ** 2    # effective noise variances
    return dict(kind="linear_iscm", D=D, W=W, sigma=sigma, mu=mu, sd=sd,
                A=A, omega=omega, order=order)


def make_linear_naive(D, rng, noise_scale=(0.5, 1.5)):
    """Control: no standardization (the var-sortability-leaking sampler)."""
    p = D.shape[0]
    W = _draw_weights(D, rng)
    sigma = rng.uniform(*noise_scale, size=p)
    return dict(kind="linear_naive", D=D, W=W, sigma=sigma,
                mu=np.zeros(p), sd=np.ones(p), A=W.copy(), omega=sigma ** 2,
                order=topological_order(D))


def sample_linear(scm, n, rng):
    p = scm["D"].shape[0]
    X = np.zeros((n, p))
    for j in scm["order"]:
        pa = np.flatnonzero(scm["D"][:, j])
        raw = X[:, pa] @ scm["W"][pa, j] + rng.normal(0, scm["sigma"][j], size=n)
        X[:, j] = (raw - scm["mu"][j]) / scm["sd"][j]
    return X


# ------------------------------------------------------------------ nonlinear iSCM
_FUNCS = [
    lambda z: z,
    np.tanh,
    np.sin,
    lambda z: z ** 2 - 1.0,
    lambda z: np.log1p(np.exp(z)) - 0.7,
]


def make_nonlinear_iscm(D, rng, noise_scale=(0.5, 1.5)):
    p = D.shape[0]
    W = _draw_weights(D, rng)
    fidx = rng.integers(0, len(_FUNCS), size=(p, p))
    sigma = rng.uniform(*noise_scale, size=p)
    order = topological_order(D)

    Xref = np.zeros((N_REF, p))
    mu = np.zeros(p)
    sd = np.ones(p)
    for j in order:
        pa = np.flatnonzero(D[:, j])
        raw = np.zeros(N_REF)
        for i in pa:
            raw += W[i, j] * _FUNCS[fidx[i, j]](Xref[:, i])
        raw += rng.normal(0, sigma[j], size=N_REF)
        mu[j] = raw.mean()
        sd[j] = raw.std()
        Xref[:, j] = (raw - mu[j]) / sd[j]

    return dict(kind="nonlinear_iscm", D=D, W=W, fidx=fidx, sigma=sigma,
                mu=mu, sd=sd, order=order)


def sample_nonlinear(scm, n, rng, intervene=None):
    """intervene: dict {node: value} for do()."""
    D, W, fidx = scm["D"], scm["W"], scm["fidx"]
    p = D.shape[0]
    X = np.zeros((n, p))
    for j in scm["order"]:
        if intervene is not None and j in intervene:
            X[:, j] = intervene[j]
            continue
        pa = np.flatnonzero(D[:, j])
        raw = np.zeros(n)
        for i in pa:
            raw += W[i, j] * _FUNCS[fidx[i, j]](X[:, i])
        raw += rng.normal(0, scm["sigma"][j], size=n)
        X[:, j] = (raw - scm["mu"][j]) / scm["sd"][j]
    return X


def sample_linear_do(scm, n, rng, intervene):
    D, W = scm["D"], scm["W"]
    p = D.shape[0]
    X = np.zeros((n, p))
    for j in scm["order"]:
        if j in intervene:
            X[:, j] = intervene[j]
            continue
        pa = np.flatnonzero(D[:, j])
        raw = X[:, pa] @ W[pa, j] + rng.normal(0, scm["sigma"][j], size=n)
        X[:, j] = (raw - scm["mu"][j]) / scm["sd"][j]
    return X


# ------------------------------------------------------------------ sortability
def var_sortability(X, D):
    """Reisach et al. 2303.18211: fraction of directed paths whose variance
    increases along the path. 0.5 = uninformative."""
    p = D.shape[0]
    v = X.var(axis=0)
    inc = tot = 0
    reach = (D > 0).astype(int)
    # all directed paths of any length via transitive closure over path counts
    M = np.eye(p, dtype=float)
    A = (D > 0).astype(float)
    for _ in range(p - 1):
        M = M @ A
        for i in range(p):
            for j in range(p):
                if M[i, j] > 0:
                    tot += M[i, j]
                    if v[j] > v[i]:
                        inc += M[i, j]
                    elif v[j] == v[i]:
                        inc += 0.5 * M[i, j]
    del reach
    return 0.5 if tot == 0 else inc / tot


def r2_sortability(X, D):
    """Reisach et al. 2303.18211 (R^2-sortability): same, with R^2 of each node
    regressed on all others in place of the variance."""
    p = D.shape[0]
    r2 = np.zeros(p)
    S = np.cov(X, rowvar=False)
    for j in range(p):
        others = [k for k in range(p) if k != j]
        if not others:
            r2[j] = 0.0
            continue
        Sxx = S[np.ix_(others, others)]
        Sxy = S[np.ix_(others, [j])]
        beta = np.linalg.solve(Sxx + 1e-10 * np.eye(len(others)), Sxy)
        resid = S[j, j] - float((Sxy.T @ beta).item())
        r2[j] = 1.0 - resid / S[j, j]
    inc = tot = 0
    M = np.eye(p, dtype=float)
    A = (D > 0).astype(float)
    for _ in range(p - 1):
        M = M @ A
        for i in range(p):
            for j in range(p):
                if M[i, j] > 0:
                    tot += M[i, j]
                    if r2[j] > r2[i]:
                        inc += M[i, j]
                    elif r2[j] == r2[i]:
                        inc += 0.5 * M[i, j]
    return 0.5 if tot == 0 else inc / tot
