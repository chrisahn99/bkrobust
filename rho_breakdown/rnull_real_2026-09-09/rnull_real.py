#!/usr/bin/env python3
"""r_0 sobre estrutura REAL, calculado APENAS de (C-hat, K, X, Y, Sigma-hat).

A pergunta do Ze: podemos calcular o raio de nulidade nos dados reais?
Este script responde executando. Nenhuma funcao no caminho de r_0 recebe o DAG
verdadeiro; ele so aparece (a) para gerar a amostra, como a natureza faria, e
(b) no fim, para AUDITAR, nunca para computar.
"""
import sys, itertools, json, time
from pathlib import Path
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
import numpy as np
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag, is_valid_adjustment_set_mpdag
from bkrobust.demo.example import knowledge_to_recover
from bkrobust.demo.graph import undirected_components


def select_knowledge(dag, cpdag, coverage):
    """Copia verbatim de bkrobust.benchmarks.measure (evita importar ortools)."""
    k_true = sorted(knowledge_to_recover(dag, cpdag))
    if coverage >= 1.0:
        return k_true
    keep = round(len(k_true) * coverage)
    if keep <= 0:
        return []
    step = len(k_true) / keep
    return [k_true[min(len(k_true) - 1, int(i * step))] for i in range(keep)]


def component_of(graph, node):
    for comp in undirected_components(graph):
        if node in comp:
            return comp
    return None
from load_real import load, sem_of


# ---------- 1. o estimando, calculado de uma matriz de covariancia -----------
def beta_se_from_sigma(sigma, idx, x, y, z, n):
    """(beta, erro-padrao) de x na regressao de y sobre {x}+z, a partir de Sigma.

    se = sqrt( sigma^2_{y|x,z} / (n * sigma^2_{x|z}) ), a forma Frisch-Waugh-Lovell
    que o `asymptotic_variance` do Chris usa; aqui de Sigma-hat em vez do SEM.
    """
    reg = [x] + sorted(set(z) - {x, y})
    ri = [idx[a] for a in reg]; yi = idx[y]
    A = sigma[np.ix_(ri, ri)]
    b = sigma[np.ix_(ri, [yi])].ravel()
    try:
        coef = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        coef = np.linalg.lstsq(A, b, rcond=None)[0]
    beta = float(coef[0])
    if n is None:
        return beta, 0.0
    s_y = float(sigma[yi, yi] - b @ coef)                       # var residual de y | x,z
    zi = ri[1:]
    if zi:
        Az = sigma[np.ix_(zi, zi)]; bz = sigma[np.ix_(zi, [idx[x]])].ravel()
        try:
            cz = np.linalg.solve(Az, bz)
        except np.linalg.LinAlgError:
            cz = np.linalg.lstsq(Az, bz, rcond=None)[0]
        s_x = float(sigma[idx[x], idx[x]] - bz @ cz)            # var residual de x | z
    else:
        s_x = float(sigma[idx[x], idx[x]])
    if s_x <= 0 or s_y <= 0:
        return beta, float("inf")
    return beta, float(np.sqrt(s_y / (n * s_x)))


def beta_from_sigma(sigma, idx, x, y, z):
    """Coeficiente de x na regressao de y sobre {x} + z, direto de Sigma.

    Vale tanto para Sigma populacional quanto para Sigma-hat amostral: e a
    mesma formula. E' esse o ponto -- o estimando nao precisa da verdade.
    """
    reg = [x] + sorted(set(z) - {x, y})
    ri = [idx[n] for n in reg]
    A = sigma[np.ix_(ri, ri)]
    b = sigma[np.ix_(ri, [idx[y]])].ravel()
    try:
        return float(np.linalg.solve(A, b)[0])
    except np.linalg.LinAlgError:
        return float(np.linalg.lstsq(A, b, rcond=None)[0][0])


# ---------- 2. o casco de efeitos sobre as extensoes -------------------------
def sub_component(g, x):
    """Sub-MPDAG induzido na componente nao-orientada de x (a 'chain component')."""
    comp = component_of(g, x)
    if comp is None or len(comp) < 2:
        return None
    ns = sorted(comp)
    dire = [(a, b) for (a, b) in g.directed_edges if a in comp and b in comp]
    und = [(a, b) for (a, b) in g.undirected_edges if a in comp and b in comp]
    return MPDAG(nodes=ns, directed=dire, undirected=und), comp


def theta(g, x, y, sigma, idx, cache, n=None, ext_cap=4096):
    """Theta(g) = {efeito ajustando por pa_D(x)} sobre as extensoes D de g.

    So a componente de x importa: pa_D(x) nao depende da orientacao das outras.
    """
    k = g.key()
    if k in cache:
        pasets = cache[k]
    else:
        fixed = {a for (a, b) in g.directed_edges if b == x}
        sc = sub_component(g, x)
        if sc is None:
            pasets = {frozenset(fixed)}
        else:
            sub, comp = sc
            exts = enumerate_dag_extensions(sub, limit=ext_cap)
            if not exts:
                return None
            pasets = {frozenset(fixed | {a for (a, b) in d.directed_edges if b == x})
                      for d in exts}
        cache[k] = pasets
    return [beta_se_from_sigma(sigma, idx, x, y, z, n) for z in pasets]


TOL = 1e-9
Z975 = 1.959963984540054


def hull(vals, band=False):
    """(lo, hi) do casco. Com band=True, alarga cada ponto por +-1.96 se."""
    if not vals:
        return None
    if band:
        return (min(b - Z975 * s for b, s in vals), max(b + Z975 * s for b, s in vals))
    return (min(b for b, _ in vals), max(b for b, _ in vals))


# ---------- 3. os dois raios -------------------------------------------------
def radii(cpdag, K, x, y, Z, sigma, idx, jmax=4, n=None):
    """Para cada orcamento j: o casco sobre a bola de retracao, e a validade de Z.

    Entradas: cpdag, K, x, y, Z, Sigma. O DAG verdadeiro NAO entra aqui.
    """
    cache, out = {}, {}
    for j in range(jmax + 1):
        vals, invalid = [], False
        for r in range(j + 1):
            for S in itertools.combinations(range(len(K)), r):
                Kp = [K[i] for i in range(len(K)) if i not in S]
                g = apply_orientations(cpdag, Kp) if Kp else cpdag
                if g is None:
                    continue
                t = theta(g, x, y, sigma, idx, cache, n)
                if t is None:
                    continue
                vals += t
                # r_val NAO e' calculado aqui: is_valid_adjustment_set_mpdag nao
                # retornou em 120 s numa rede de 46 nos (ver TETO.md). Chris ja o
                # mediu em estrutura real; o que falta e' r_0, e r_0 nao precisa dele.
        h, hb = hull(vals), hull(vals, band=True)
        if h is None:
            continue
        structural = any(abs(b) <= TOL for b, _ in vals)
        out[j] = dict(lo=h[0], hi=h[1],
                      has_zero=bool(h[0] <= TOL and h[1] >= -TOL),
                      band_lo=hb[0], band_hi=hb[1],
                      band_zero=bool(hb[0] <= 0.0 <= hb[1]),
                      structural_zero=structural,
                      n_eff=len(set(round(b, 9) for b, _ in vals)))
    return out


def first(d, key):
    for j in sorted(d):
        if d[j][key]:
            return j
    return None


# ---------- 4. amostragem: a natureza, com os coeficientes reais -------------
def sample(dag, W, var, n, rng):
    """n observacoes do SEM linear-gaussiano ajustado aos dados reais."""
    order, seen = [], set()
    def visit(v):
        if v in seen: return
        seen.add(v)
        for p in dag.parents(v): visit(p)
        order.append(v)
    for v in dag.nodes: visit(v)
    idx = {v: i for i, v in enumerate(dag.nodes)}
    X = np.zeros((n, len(dag.nodes)))
    for v in order:
        col = rng.normal(0.0, np.sqrt(var[v]), size=n)
        for p in dag.parents(v):
            col = col + W[(p, v)] * X[:, idx[p]]
        X[:, idx[v]] = col
    return X, idx


def pairs_for(dag, cpdag, max_pairs, min_comp=3, max_comp=9, rng=None):
    """(X, Y) com X numa componente nao-orientada de tamanho tratavel e Y descendente."""
    out = []
    for x in dag.nodes:
        comp = component_of(cpdag, x)
        if comp is None or not (min_comp <= len(comp) <= max_comp):
            continue
        for y in sorted(dag.descendants(x)):
            if y in comp:
                continue
            out.append((x, y))
    if rng is not None:
        rng.shuffle(out)
    return out[:max_pairs]


def run(name="ecoli70", n_obs=(500, 5000), coverage=1.0, max_pairs=40,
        seed=20260909, jmax=4):
    rng = np.random.default_rng(seed)
    dag, W, var = load(name)
    sem = sem_of(dag, W, var)
    cpdag = dag_to_cpdag(dag)
    sigma_pop = sem.covariance()
    idx_pop = {v: i for i, v in enumerate(dag.nodes)}
    prs = pairs_for(dag, cpdag, max_pairs, rng=rng)
    rows = []
    samples = {n: sample(dag, W, var, n, rng) for n in n_obs}
    for (x, y) in prs:
        K = sorted(select_knowledge(dag, cpdag, coverage))
        g0 = apply_orientations(cpdag, K) if K else cpdag
        if g0 is None:
            continue
        Z = optimal_adjustment_set_mpdag(g0, x, y)
        if Z is None:
            continue
        Kloc = [e for e in K if component_of(cpdag, e[0]) == component_of(cpdag, x)]
        if not Kloc:
            continue
        rec = dict(network=name, x=x, y=y, k_local=len(Kloc), Z=sorted(Z),
                   comp=len(component_of(cpdag, x)))
        t0 = time.time()
        # -- oraculo: Sigma populacional (nao disponivel na pratica) --
        d = radii(cpdag, Kloc, x, y, Z, sigma_pop, idx_pop, jmax)
        rec["oracle"] = dict(r0=first(d, "has_zero"), rval=None,
                             hull0=[d[0]["lo"], d[0]["hi"]] if 0 in d else None,
                             n_eff0=d[0]["n_eff"] if 0 in d else None)
        # -- o caso real: Sigma-hat de n observacoes --
        for n, (Xs, idxs) in samples.items():
            sh = np.cov(Xs, rowvar=False)
            dn = radii(cpdag, Kloc, x, y, Z, sh, idxs, jmax)
            rec[f"n{n}"] = dict(r0=first(dn, "has_zero"), rval=None,
                                hull0=[dn[0]["lo"], dn[0]["hi"]] if 0 in dn else None)
        rec["secs"] = round(time.time() - t0, 2)
        rows.append(rec)
        print(f"  [{len(rows):2d}/{len(prs)}] {x}->{y} |C|={rec['comp']} "
              f"|K|={rec['k_local']} r0={rec['oracle']['r0']} "
              f"rval={rec['oracle']['rval']} {rec['secs']}s", flush=True)
    return rows


if __name__ == "__main__":
    name = sys.argv[1] if len(sys.argv) > 1 else "ecoli70"
    rows = run(name)
    print(f"\n{len(rows)} instancias em {name}\n")
    hdr = f"{'X':>7} {'Y':>7} {'|C|':>3} {'|K|':>3} {'r0 orac':>7} {'rval':>4} " \
          f"{'r0 n=500':>8} {'r0 n=5000':>9} {'casco (oraculo)':>26}"
    print(hdr); print("-" * len(hdr))
    for r in rows:
        h = r["oracle"]["hull0"]
        hs = f"[{h[0]:+.3f}, {h[1]:+.3f}]" if h else "--"
        f = lambda v: "nunca" if v is None else str(v)
        print(f"{r['x']:>7} {r['y']:>7} {r['comp']:>3} {r['k_local']:>3} "
              f"{f(r['oracle']['r0']):>7} {f(r['oracle']['rval']):>4} "
              f"{f(r['n500']['r0']):>8} {f(r['n5000']['r0']):>9} {hs:>26}")
    json.dump(rows, open(f"resultados_{name}.json", "w"), indent=1)
