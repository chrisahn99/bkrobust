"""Sortability probes on the Gaussian tier: the substrate check for Table 2.

Var-sortability (Reisach, Seiler and Weichwald, NeurIPS 2021, arXiv:2102.13647,
section 3.1) and R2-sortability (Reisach, Tami, Seiler, Chambaz and Weichwald,
NeurIPS 2023, arXiv:2303.18211, equation 2) of ancestral samples drawn from the
four networks that ship with fitted linear-Gaussian parameters: ``arth150``,
``ecoli70``, ``magic-irri`` and ``magic-niab``. Both criteria are implemented
here from the published definitions; nothing is imported from CausalDisco.

The shared functional. For a DAG with adjacency ``E`` and a per-node score
``tau``, sum over path lengths ``k = 1 .. d-1`` and over the ordered pairs
``(s, t)`` joined by at least one directed path of length exactly ``k`` the
value ``incr(tau_s, tau_t)``, which is 1 when the score rises along the path,
1/2 on a tie and 0 when it falls, and divide by the number of such ``(k, s, t)``
terms. Var-sortability is that functional with the marginal variance as the
score; R2-sortability with the coefficient of determination of a least-squares
regression of each variable on all the others. Both papers define the count per
path *length*, not per path: a pair joined by two paths of the same length is
one term.

Every other network in the corpus has no continuous parameters (conditional
probability tables for the BIF files, structure only for the dagitty diagrams)
and is written out as ``UNDEFINED`` with its reason rather than omitted, so the
substrate table can print the gap instead of hiding it.

Beside the sampled values, the population value is printed: the exact
covariance ``(I - B)^-1 Omega (I - B)^-T`` gives the marginal variances and, via
the precision of the correlation matrix, the population R2 of the best linear
predictor of each variable from the others. The seed dispersion is then read
against a fixed target rather than against itself, and a sampler bug would show
up as a gap between the two.

    python experiments/sortability_probe.py --n 5000 --seeds 10
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for entry in (ROOT, ROOT / "src"):
    if str(entry) not in sys.path:
        sys.path.insert(0, str(entry))

from bkrobust.benchmarks.describe import parse_file  # noqa: E402
from experiments.ledger_sweep import TIER  # noqa: E402

MODELS = ROOT / "results" / "axisa3" / "networks" / "example_models"
OUT = ROOT / "results" / "substrate" / "sortability.json"

DEFINITIONS = {
    "var_sortability": (
        "Reisach, Seiler and Weichwald, 'Beware of the Simulated DAG!', NeurIPS 2021, "
        "arXiv:2102.13647, section 3.1: fraction of (path length, source, target) terms whose "
        "marginal variance increases along the path; ties count one half."
    ),
    "r2_sortability": (
        "Reisach, Tami, Seiler, Chambaz and Weichwald, 'A Scale-Invariant Sorting Criterion "
        "to Find a Causal Order in Additive Noise Models', NeurIPS 2023, arXiv:2303.18211, "
        "equation 2 with tau = R2 of a least-squares regression of each variable on all others."
    ),
}


class GaussianNetwork:
    """A linear-Gaussian network read from a bnlearn-style JSON file.

    Args:
        name: Short network name.
        nodes: Node names in file order.
        arcs: Directed arcs as ``(parent, child)``.
        intercept: Per-node intercept.
        coef: Per-node map from parent to fitted coefficient.
        noise_var: Per-node residual variance.
    """

    def __init__(
        self,
        name: str,
        nodes: list[str],
        arcs: list[tuple[str, str]],
        intercept: dict[str, float],
        coef: dict[str, dict[str, float]],
        noise_var: dict[str, float],
    ) -> None:
        self.name = name
        self.nodes = nodes
        self.arcs = arcs
        self.intercept = intercept
        self.coef = coef
        self.noise_var = noise_var
        self.index = {v: i for i, v in enumerate(nodes)}

    @property
    def adjacency(self) -> np.ndarray:
        """Boolean adjacency, ``adj[s, t]`` true for the arc ``s -> t``."""
        adj = np.zeros((len(self.nodes), len(self.nodes)), dtype=bool)
        for s, t in self.arcs:
            adj[self.index[s], self.index[t]] = True
        return adj

    def topological_order(self) -> list[str]:
        """Kahn's algorithm over the arcs; raises if the file is not acyclic."""
        indeg = dict.fromkeys(self.nodes, 0)
        children: dict[str, list[str]] = {v: [] for v in self.nodes}
        for s, t in self.arcs:
            indeg[t] += 1
            children[s].append(t)
        ready = [v for v in self.nodes if indeg[v] == 0]
        order: list[str] = []
        while ready:
            v = ready.pop(0)
            order.append(v)
            for c in children[v]:
                indeg[c] -= 1
                if indeg[c] == 0:
                    ready.append(c)
        if len(order) != len(self.nodes):
            raise ValueError(f"{self.name}: the arcs contain a cycle")
        return order

    def population_covariance(self) -> np.ndarray:
        """Exact covariance ``(I - B)^-1 Omega (I - B)^-T`` with ``B[child, parent]``."""
        d = len(self.nodes)
        b = np.zeros((d, d))
        for t, parents in self.coef.items():
            for s, w in parents.items():
                b[self.index[t], self.index[s]] = w
        omega = np.diag([self.noise_var[v] for v in self.nodes])
        inv = np.linalg.inv(np.eye(d) - b)
        return inv @ omega @ inv.T

    def sample(self, n: int, rng: np.random.Generator) -> np.ndarray:
        """Ancestral sampling: each node is its fitted equation plus Gaussian noise."""
        x = np.zeros((n, len(self.nodes)))
        for v in self.topological_order():
            j = self.index[v]
            mean = np.full(n, self.intercept[v])
            for s, w in self.coef[v].items():
                mean += w * x[:, self.index[s]]
            x[:, j] = mean + np.sqrt(self.noise_var[v]) * rng.standard_normal(n)
        return x


def load_gaussian(path: Path) -> GaussianNetwork:
    """Read a ``{nodes, arcs, cpds}`` file and check the parameters cover the arcs.

    Args:
        path: The JSON file.

    Returns:
        The network with its fitted coefficients and residual variances.

    Raises:
        ValueError: If any node lacks a parameter block, if the coefficient keys
            disagree with the arcs, or if a residual variance is not positive.
    """
    payload = json.loads(path.read_text(encoding="utf-8"))
    nodes = list(payload["nodes"])
    arcs = [(a, b) for a, b in payload["arcs"]]
    cpds = payload["cpds"]
    parents: dict[str, set[str]] = {v: set() for v in nodes}
    for s, t in arcs:
        parents[t].add(s)
    intercept: dict[str, float] = {}
    coef: dict[str, dict[str, float]] = {}
    noise_var: dict[str, float] = {}
    for v in nodes:
        block = cpds.get(v)
        if block is None:
            raise ValueError(f"{path.name}: node {v!r} has no parameter block")
        raw = block["coefficients"]
        intercept[v] = float(raw.get("(Intercept)", [0.0])[0])
        coef[v] = {p: float(c[0]) for p, c in raw.items() if p != "(Intercept)"}
        if set(coef[v]) != parents[v]:
            raise ValueError(f"{path.name}: coefficients of {v!r} do not match its arcs")
        variance = float(block["variance"][0])
        if not variance > 0:
            raise ValueError(f"{path.name}: residual variance of {v!r} is {variance}")
        noise_var[v] = variance
    return GaussianNetwork(path.stem, nodes, arcs, intercept, coef, noise_var)


def order_alignment(adj: np.ndarray, score: np.ndarray) -> tuple[float, int]:
    """The tau-sortability functional (Reisach et al. 2023, equation 2).

    Args:
        adj: Boolean adjacency of the DAG, ``adj[s, t]`` for ``s -> t``.
        score: One value per node, in adjacency order.

    Returns:
        The sortability in ``[0, 1]`` and the number of ``(k, s, t)`` terms it
        averages over. A graph with no arcs has no terms and returns ``nan``.
    """
    d = adj.shape[0]
    reach = adj.copy()
    numerator = 0.0
    terms = 0
    for _ in range(d - 1):
        if not reach.any():
            break
        s, t = np.nonzero(reach)
        rises = score[t] > score[s]
        ties = score[t] == score[s]
        numerator += rises.sum() + 0.5 * ties.sum()
        terms += len(s)
        reach = (reach.astype(np.int64) @ adj.astype(np.int64)) > 0
    return (numerator / terms if terms else float("nan")), terms


def r2_on_all_others(x: np.ndarray) -> np.ndarray:
    """Coefficient of determination of each column regressed on all the others.

    Least squares with an intercept, exactly the fit the 2023 paper prescribes:
    ``R2_t = 1 - Var(X_t - fit_t) / Var(X_t)``.

    Args:
        x: Samples, one column per variable.

    Returns:
        One R2 per column.
    """
    n, d = x.shape
    out = np.empty(d)
    ones = np.ones((n, 1))
    for t in range(d):
        others = np.delete(x, t, axis=1)
        design = np.hstack([ones, others])
        beta, *_ = np.linalg.lstsq(design, x[:, t], rcond=None)
        resid = x[:, t] - design @ beta
        out[t] = 1.0 - resid.var() / x[:, t].var()
    return out


def population_r2(sigma: np.ndarray) -> np.ndarray:
    """Population R2 of the best linear predictor of each variable from the others.

    From the precision of the correlation matrix: ``R2_t = 1 - 1 / (C^-1)_tt``.

    Args:
        sigma: Population covariance.

    Returns:
        One R2 per variable.
    """
    sd = np.sqrt(np.diag(sigma))
    corr = sigma / np.outer(sd, sd)
    return 1.0 - 1.0 / np.diag(np.linalg.inv(corr))


def summarise(values: list[float]) -> dict[str, float | list[float]]:
    """Mean, seed dispersion and range of a per-seed list."""
    arr = np.asarray(values)
    return {
        "mean": round(float(arr.mean()), 4),
        "sd": round(float(arr.std(ddof=1)) if len(arr) > 1 else 0.0, 4),
        "min": round(float(arr.min()), 4),
        "max": round(float(arr.max()), 4),
        "per_seed": [round(float(v), 4) for v in arr],
    }


def probe(net: GaussianNetwork, n: int, seeds: list[int]) -> dict[str, object]:
    """Both sortabilities on ``len(seeds)`` ancestral samples of size ``n``."""
    adj = net.adjacency
    sigma = net.population_covariance()
    pop_var, terms = order_alignment(adj, np.diag(sigma))
    pop_r2, _ = order_alignment(adj, population_r2(sigma))
    var_values: list[float] = []
    r2_values: list[float] = []
    worst_var_error = 0.0
    for seed in seeds:
        x = net.sample(n, np.random.default_rng(seed))
        sample_var = x.var(axis=0, ddof=1)
        worst_var_error = max(
            worst_var_error, float(np.max(np.abs(sample_var / np.diag(sigma) - 1.0)))
        )
        var_values.append(order_alignment(adj, sample_var)[0])
        r2_values.append(order_alignment(adj, r2_on_all_others(x))[0])
    return {
        "status": "defined",
        "n_nodes": len(net.nodes),
        "n_edges": len(net.arcs),
        "n_path_terms": terms,
        "var_sortability": {**summarise(var_values), "population": round(float(pop_var), 4)},
        "r2_sortability": {**summarise(r2_values), "population": round(float(pop_r2), 4)},
        "sampler_check_max_rel_var_error": round(worst_var_error, 4),
    }


def undefined_reason(path: Path) -> str:
    """Why a network file carries no sortability."""
    if path.name.endswith((".bif", ".bif.gz")):
        return "discrete conditional probability tables (BIF); no continuous parameters"
    parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
    if parsed.bidirected:
        return "bidirected edges (latent confounding); not a DAG and no parameters"
    return "structure only (dagitty diagram); no parameters"


def main(argv: list[str] | None = None) -> None:
    """Run the probes and write ``results/substrate/sortability.json``."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--n", type=int, default=5000, help="sample size per seed")
    ap.add_argument("--seeds", type=int, default=10, help="seeds 0 .. seeds-1")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    seeds = list(range(args.seeds))

    t0 = time.perf_counter()
    networks: dict[str, dict[str, object]] = {}
    for path in sorted(MODELS.iterdir()):
        if not path.is_file() or path.name.startswith("_"):
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        base = {"source_file": path.name, "sha256": digest, "tier": TIER.get(path.stem, "T1")}
        if path.suffix == ".json":
            net = load_gaussian(path)
            networks[net.name] = {**base, **probe(net, args.n, seeds)}
            row = networks[net.name]
            var = row["var_sortability"]
            r2 = row["r2_sortability"]
            assert isinstance(var, dict) and isinstance(r2, dict)
            print(
                f"{net.name:12s} d={len(net.nodes):4d} e={len(net.arcs):4d} "
                f"terms={row['n_path_terms']:6} "
                f"var={var['mean']:.3f}±{var['sd']:.3f} (pop {var['population']:.3f}) "
                f"r2={r2['mean']:.3f}±{r2['sd']:.3f} (pop {r2['population']:.3f}) "
                f"varerr={row['sampler_check_max_rel_var_error']}"
            )
        else:
            networks[path.stem if not path.name.endswith(".bif.gz") else path.name[:-7]] = {
                **base,
                "status": "UNDEFINED",
                "reason": undefined_reason(path),
            }

    payload = {
        "definitions": DEFINITIONS,
        "n": args.n,
        "seeds": seeds,
        "seconds": round(time.perf_counter() - t0, 1),
        "n_defined": sum(1 for r in networks.values() if r["status"] == "defined"),
        "n_undefined": sum(1 for r in networks.values() if r["status"] == "UNDEFINED"),
        "networks": networks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    print(
        f"defined {payload['n_defined']}, undefined {payload['n_undefined']}, "
        f"{payload['seconds']}s -> {args.out.relative_to(ROOT)}"
    )


if __name__ == "__main__":
    main()
