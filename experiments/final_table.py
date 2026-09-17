"""The paper's final table: r_val and r_eps on every real network, at LLM knowledge.

For each real benchmark network in ``results/axisa3/networks/example_models/``,
at its LLM-elicited background-knowledge set ``K`` (from
``results/elicit/knowledge.json``), this script certifies a handful of causal
queries (from ``results/frame/frame.jsonl``) and reports the validity radius
``r_val`` and the epsilon-bias radius ``r_eps`` at a grid of epsilons.

The benchmark corpus ships graph structure only -- no continuous observational
data -- but ``r_eps`` needs a covariance. We supply one by attaching a seeded
linear-Gaussian SEM to the ground-truth DAG (:func:`bkrobust.demo.evaluate.
random_sem`); see ``ASSUMPTIONS`` below, which is also written into the JSON
output and printed as a footnote under the table.

Usage::

    PYTHONPATH=src .venv/bin/python experiments/final_table.py [options]

Writes ``<out>/instances.jsonl`` (one record per network/query) and
``<out>/summary.json`` (per-network aggregates), and prints a summary table.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import signal
import statistics
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from bkrobust.benchmarks.describe import parse_file, to_mpdag  # noqa: E402
from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.demo.evaluate import random_sem  # noqa: E402
from bkrobust.demo.example import dag_to_cpdag  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.epsilon.certify import certify  # noqa: E402

NETWORKS_DIR = ROOT / "results" / "axisa3" / "networks" / "example_models"
KNOWLEDGE_PATH = ROOT / "results" / "elicit" / "knowledge.json"
FRAME_PATH = ROOT / "results" / "frame" / "frame.jsonl"

DEFAULT_EPSILONS = "0.01,0.05,0.10,0.25,0.50,1.00"
DEFAULT_SKIP = "pathfinder"

#: The covariance assumption behind r_eps on this corpus. Recorded verbatim in
#: summary.json and printed as a table footnote -- see the module docstring.
ASSUMPTIONS: dict[str, Any] = {
    "covariance": (
        "The example_models corpus ships graph structure only (no continuous "
        "observational data), but r_eps needs a covariance matrix to evaluate "
        "bias. We attach a seeded linear-Gaussian SEM to the ground-truth DAG "
        "via bkrobust.demo.evaluate.random_sem(dag, rng) and pass it to "
        "certify(..., sem=sem); certify derives the covariance, and the "
        "reported estimate theta_z, from that SEM. The seed for each "
        "(network, x, y) instance is derived deterministically from --seed "
        "(sha256 of f'{seed}:{network}:{x}:{y}'), so the sweep is exactly "
        "reproducible from the CLI arguments alone. The SEM is a synthetic "
        "stand-in for a data-generating process these benchmark networks do "
        "not ship -- it is never treated as the network's real parameters."
    ),
    "random_sem_defaults": {"weight_range": [0.3, 1.5], "noise_range": [0.5, 1.5]},
}


class InstanceTimeout(Exception):  # noqa: N818 - a control-flow signal, not an error condition
    """Raised when one (network, x, y) instance exceeds --time-limit."""


def _alarm_handler(signum: int, frame: Any) -> None:
    """Turn SIGALRM into an :class:`InstanceTimeout`."""
    raise InstanceTimeout("instance exceeded --time-limit")


def run_with_timeout(time_limit: float, func: Any, *args: Any, **kwargs: Any) -> Any:
    """Call ``func(*args, **kwargs)`` under a SIGALRM wall-clock cap.

    Args:
        time_limit: Seconds before raising :class:`InstanceTimeout`. ``<= 0``
            disables the cap.
        func: Callable to invoke.
        *args: Positional arguments to ``func``.
        **kwargs: Keyword arguments to ``func``.

    Returns:
        Whatever ``func`` returns.

    Raises:
        InstanceTimeout: If ``func`` has not returned within ``time_limit``
            seconds.
    """
    if time_limit and time_limit > 0 and hasattr(signal, "SIGALRM"):
        old_handler = signal.signal(signal.SIGALRM, _alarm_handler)
        signal.alarm(max(1, round(time_limit)))
        try:
            return func(*args, **kwargs)
        finally:
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old_handler)
    return func(*args, **kwargs)


def derive_seed(base_seed: int, network: str, x: str, y: str) -> int:
    """Deterministically derive a per-instance RNG seed from the base seed.

    Uses sha256 rather than Python's hash() so the derived seed does not
    depend on PYTHONHASHSEED and is stable across processes and runs.

    Args:
        base_seed: The CLI ``--seed``.
        network: Network name.
        x: Treatment node.
        y: Outcome node.

    Returns:
        A seed in ``[0, 2**32)``.
    """
    payload = f"{base_seed}:{network}:{x}:{y}".encode()
    digest = hashlib.sha256(payload).digest()
    return int.from_bytes(digest[:8], "big") % (2**32)


def load_networks(
    names: set[str] | None, skip: set[str], max_nodes: int
) -> dict[str, dict[str, MPDAG]]:
    """Parse real networks and build their CPDAGs, but only where needed.

    Two passes: a cheap ``parse_file`` learns each file's network name,
    bidirectedness, and node count; the expensive step -- ``to_mpdag`` plus
    the Meek closure in ``dag_to_cpdag`` -- is then paid only for names that
    survive ``names``/``skip``/``max_nodes`` filtering. This matters on this
    corpus: several real networks have 700-1000+ nodes, and building their
    CPDAGs is the dominant cost of a small ``--networks`` run if it is not
    skipped for everything else.

    Args:
        names: If given, only these network names are built past the parse
            step (still subject to ``skip``/``max_nodes``). ``None`` means no
            name filter -- every non-skipped, non-oversize network is built.
        skip: Network names to exclude outright.
        max_nodes: If positive, drop networks with more nodes than this.

    Returns:
        ``{network_name: {"dag": MPDAG, "cpdag": MPDAG}}``.
    """
    out: dict[str, dict[str, MPDAG]] = {}
    for path in sorted(NETWORKS_DIR.iterdir()):
        if not path.is_file() or path.name.startswith("."):
            continue
        try:
            parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
        except Exception as exc:
            print(f"[final_table] parse failed for {path.name}: {exc}", file=sys.stderr)
            continue
        if parsed.bidirected:
            continue  # an ADMG (e.g. M-bias.txt), not a DAG
        if parsed.name in skip:
            continue
        if names is not None and parsed.name not in names:
            continue
        if max_nodes and len(parsed.nodes) > max_nodes:
            continue
        try:
            dag = to_mpdag(parsed)
            cpdag = dag_to_cpdag(dag)
        except Exception as exc:
            print(f"[final_table] cpdag build failed for {parsed.name}: {exc}", file=sys.stderr)
            continue
        out[parsed.name] = {"dag": dag, "cpdag": cpdag}
    return out


def load_knowledge(condition: str) -> dict[str, list[tuple[str, str]]]:
    """Load the LLM-elicited K for every network under one condition.

    Args:
        condition: A key of ``results/elicit/knowledge.json`` (e.g. ``D_LLM``).

    Returns:
        ``{network_name: [(tail, head), ...]}``.

    Raises:
        SystemExit: If ``condition`` is not present in the file.
    """
    data = json.loads(KNOWLEDGE_PATH.read_text())
    if condition not in data:
        raise SystemExit(f"unknown --condition {condition!r}; choices: {sorted(data)}")
    nets = data[condition]["networks"]
    return {net: [tuple(e) for e in info["k"]] for net, info in nets.items()}


def load_queries(queries_per_network: int) -> dict[str, list[tuple[str, str]]]:
    """Load the first ``queries_per_network`` (x, y) pairs per network, in file order.

    Args:
        queries_per_network: How many rows to keep per network.

    Returns:
        ``{network_name: [(x, y), ...]}``.
    """
    out: dict[str, list[tuple[str, str]]] = {}
    with FRAME_PATH.open() as f:
        for line in f:
            row = json.loads(line)
            net = row["network"]
            bucket = out.setdefault(net, [])
            if len(bucket) < queries_per_network:
                bucket.append((row["x"], row["y"]))
    return out


def run_instance(
    condition: str,
    net: str,
    x: str,
    y: str,
    dag: MPDAG,
    cpdag: MPDAG,
    k: list[tuple[str, str]],
    epsilons: tuple[float, ...],
    seed: int,
    time_limit: float,
    search_budget: int,
) -> dict[str, Any]:
    """Certify one (network, x, y) instance, catching and recording any failure.

    Args:
        condition: The knowledge condition name, recorded on the record.
        net: Network name.
        x: Treatment node.
        y: Outcome node.
        dag: The ground-truth DAG (for the SEM).
        cpdag: The estimated CPDAG passed to certify.
        k: The asserted background knowledge.
        epsilons: The epsilon grid, in relative units.
        seed: Per-instance RNG seed (see :func:`derive_seed`).
        time_limit: Wall-clock cap in seconds; ``<= 0`` disables it.
        search_budget: Depth budget forwarded to ``certify``.

    Returns:
        A JSON-serialisable record: network, x, y, condition, len_k,
        n_knowledge, r_val, r_eps, theta_z, status, seconds, error.
    """
    record: dict[str, Any] = {
        "network": net,
        "n_nodes": len(cpdag.nodes),
        "x": x,
        "y": y,
        "condition": condition,
        "len_k": len(k),
        "n_knowledge": None,
        "r_val": None,
        "r_eps": None,
        "theta_z": None,
        "status": None,
        "seconds": None,
        "error": None,
    }
    t0 = time.perf_counter()
    try:
        rng = np.random.default_rng(seed)
        sem = random_sem(dag, rng)
        cert = run_with_timeout(
            time_limit,
            certify,
            cpdag,
            k,
            x,
            y,
            sem=sem,
            epsilons=epsilons,
            units="relative",
            method="semilocal",
            search_budget=search_budget,
        )
        record["n_knowledge"] = cert.n_knowledge
        record["r_val"] = cert.r_val
        record["r_eps"] = {str(eps): radius for eps, radius in cert.r_eps.items()}
        record["theta_z"] = cert.theta_z
        record["status"] = "ok"
    except InstanceTimeout as exc:
        record["status"] = "timeout"
        record["error"] = str(exc)
    except Exception as exc:
        record["status"] = "error"
        record["error"] = f"{type(exc).__name__}: {exc}"
    record["seconds"] = round(time.perf_counter() - t0, 4)
    return record


def _median_finite(values: list[Any]) -> float | None:
    """Median of the *finite* radii in ``values``.

    :data:`UNREACHED` is a sentinel, not a radius, so it is dropped before the
    median is taken -- averaging it in would produce a negative "radius", which
    is what an earlier version of this function did. ``None`` is returned when
    nothing finite survives, which the caller must distinguish from "no data".
    """
    finite = [v for v in values if v is not None and v != UNREACHED]
    return statistics.median(finite) if finite else None


def _cell(values: list[Any]) -> dict[str, Any]:
    """Summarise one column's radii for one network into a renderable cell.

    Three outcomes are kept apart, because collapsing them loses the meaning:
    a finite median, "every resolved instance was UNREACHED", and "no instance
    resolved at all" (every query timed out or errored).
    """
    resolved = [v for v in values if v is not None]
    n_unr = sum(1 for v in resolved if v == UNREACHED)
    return {
        "value": _median_finite(resolved),
        "n_resolved": len(resolved),
        "n_unreached": n_unr,
    }


def summarise_network(
    n_nodes: int,
    net_instances: list[dict[str, Any]],
    epsilons: tuple[float, ...],
) -> dict[str, Any]:
    """Aggregate one network's instance records into the table row's numbers.

    Args:
        n_nodes: ``|V|`` for the network.
        net_instances: This network's instance records.
        epsilons: The epsilon grid.

    Returns:
        Aggregate dict used both for summary.json and the printed table.
    """
    ok = [r for r in net_instances if r["status"] == "ok"]
    len_k = net_instances[0]["len_k"] if net_instances else 0
    eps_cells = {
        str(eps): _cell(
            [r["r_eps"][str(eps)] for r in ok if r["r_eps"] is not None]
        )
        for eps in epsilons
    }
    # |K_G0| is a count of closure-oriented edges, never a radius: it has no
    # UNREACHED case. Take it from any solved instance; it is graph-determined
    # and so identical across this network's queries.
    n_knowledge = next((r["n_knowledge"] for r in ok if r["n_knowledge"] is not None), None)
    return {
        "n_nodes": n_nodes,
        "len_k": len_k,
        "n_knowledge": n_knowledge,
        "n_queries": len(net_instances),
        "n_ok": len(ok),
        "n_error": sum(1 for r in net_instances if r["status"] == "error"),
        "n_timeout": sum(1 for r in net_instances if r["status"] == "timeout"),
        "r_val": _cell([r["r_val"] for r in ok]),
        "r_eps": eps_cells,
    }


def _fmt_count(value: Any) -> str:
    """Render a plain count, or an em dash when it is unknown."""
    return "-" if value is None else str(int(value))


def _fmt_cell(cell: dict[str, Any]) -> str:
    """Render a :func:`_cell` summary for a fixed-width table column.

    ``-`` means no query resolved (all timed out or errored) -- unknown, not
    unreachable. ``UNR`` means every resolved query was genuinely UNREACHED. A
    trailing ``*`` marks a finite median that co-exists with some UNREACHED
    queries, so the reader knows the cell rests on a subset.
    """
    if cell["n_resolved"] == 0:
        return "-"
    if cell["value"] is None:
        return "UNR"
    value = cell["value"]
    text = str(int(value)) if float(value).is_integer() else f"{value:.1f}"
    return text + ("*" if cell["n_unreached"] else "")


def render_table(
    universe: list[str],
    summaries: dict[str, dict[str, Any]],
    all_ok: list[dict[str, Any]],
    epsilons: tuple[float, ...],
    condition: str,
) -> str:
    """Render the fixed-width stdout summary table.

    Args:
        universe: Network names, in the order to print them.
        summaries: Per-network aggregates from :func:`summarise_network`.
        all_ok: Every ``status == "ok"`` instance across all networks, for the
            pooled TOTAL/median row.
        epsilons: The epsilon grid (used for column headers).
        condition: The knowledge condition, printed above the table.

    Returns:
        The table as one multi-line string.
    """
    net_w = max([len("network"), len("TOTAL/median")] + [len(n) for n in universe]) + 1
    eps_headers = [f"eps={e:g}" for e in epsilons]
    col_w = max([9] + [len(h) + 1 for h in eps_headers])
    header = (
        f"{'network':<{net_w}}"
        f"{'|V|':>6}{'K':>6}{'|K_G0|':>8}{'ok/to/er':>10}{'r_val':>{col_w}}"
        + "".join(f"{h:>{col_w}}" for h in eps_headers)
    )
    rule = "-" * len(header)
    lines = [f"final_table  condition={condition}", rule, header, rule]

    for net in universe:
        s = summaries[net]
        support = f"{s['n_ok']}/{s['n_timeout']}/{s['n_error']}"
        row = (
            f"{net:<{net_w}}"
            f"{s['n_nodes']:>6}{s['len_k']:>6}{_fmt_count(s['n_knowledge']):>8}"
            f"{support:>10}{_fmt_cell(s['r_val']):>{col_w}}"
        )
        for eps in epsilons:
            row += f"{_fmt_cell(s['r_eps'][str(eps)]):>{col_w}}"
        lines.append(row)

    lines.append(rule)
    pooled_support = (
        f"{sum(summaries[n]['n_ok'] for n in universe)}"
        f"/{sum(summaries[n]['n_timeout'] for n in universe)}"
        f"/{sum(summaries[n]['n_error'] for n in universe)}"
    )
    total_row = (
        f"{'POOLED median':<{net_w}}"
        f"{'':>6}{'':>6}{'':>8}{pooled_support:>10}"
        f"{_fmt_cell(_cell([r['r_val'] for r in all_ok])):>{col_w}}"
    )
    for eps in epsilons:
        key = str(eps)
        pooled = _cell([r["r_eps"][key] for r in all_ok if r["r_eps"] is not None])
        total_row += f"{_fmt_cell(pooled):>{col_w}}"
    lines.append(total_row)
    lines.append(rule)
    lines.append(
        "Assumption: no real network here ships observational data, so r_eps is "
        "computed against a seeded linear-Gaussian SEM attached to the "
        "ground-truth DAG (bkrobust.demo.evaluate.random_sem), seeded "
        "deterministically per (network, x, y) from --seed; see 'assumptions' "
        "in summary.json for the full statement."
    )
    lines.append(
        "Legend: a number is the median over queries that resolved to a finite "
        "radius. UNR = every resolved query was UNREACHED (no perturbation on "
        "the enumerated space crosses that threshold) -- a result, not a "
        "radius. '-' = no query resolved (all timed out or errored) -- unknown, "
        "NOT unreachable. '*' = a finite median that co-exists with some "
        "UNREACHED queries. ok/to/er = queries that solved / timed out / "
        "errored. epsilon is relative: a fraction of the reported estimate "
        "theta_z."
    )
    return "\n".join(lines)


def build_arg_parser() -> argparse.ArgumentParser:
    """Build the CLI argument parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--condition", default="D_LLM", help="knowledge.json condition key")
    parser.add_argument("--epsilons", default=DEFAULT_EPSILONS, help="comma list, relative units")
    parser.add_argument("--queries-per-network", type=int, default=5)
    parser.add_argument("--networks", default="", help="comma list; default all available")
    parser.add_argument("--skip", default=DEFAULT_SKIP, help="comma list of networks to exclude")
    parser.add_argument("--seed", type=int, default=20260917)
    parser.add_argument("--time-limit", type=float, default=120.0, help="seconds per instance")
    parser.add_argument("--max-nodes", type=int, default=0, help="0 = no cap")
    parser.add_argument("--out", default="results/final_table")
    parser.add_argument(
        "--report-only",
        action="store_true",
        help="re-render the table and summary.json from an existing "
        "<out>/instances.jsonl, without recomputing any radius",
    )
    return parser


def _report_only(out_dir: Path, epsilons: tuple[float, ...], args: argparse.Namespace) -> None:
    """Re-aggregate and re-render from an existing ``instances.jsonl``.

    The sweep is expensive and its per-instance records are the ground truth,
    so a presentation-layer fix must not require recomputing any radius.
    """
    instances_path = out_dir / "instances.jsonl"
    instances = [json.loads(line) for line in instances_path.read_text().splitlines() if line]
    # |V| is a property of the network, not of an instance; older sweeps did not
    # record it per instance, so fall back to the previous summary.json.
    n_nodes = {r["network"]: r.get("n_nodes") for r in instances if r.get("n_nodes")}
    summary_path_in = out_dir / "summary.json"
    if summary_path_in.exists():
        prior = json.loads(summary_path_in.read_text()).get("networks", {})
        for net, row in prior.items():
            n_nodes.setdefault(net, row.get("n_nodes"))
    per_network: dict[str, list[dict[str, Any]]] = {}
    for record in instances:
        per_network.setdefault(record["network"], []).append(record)
    universe = sorted(per_network)
    summaries = {
        net: summarise_network(n_nodes.get(net) or 0, per_network[net], epsilons)
        for net in universe
    }
    all_ok = [r for r in instances if r["status"] == "ok"]
    summary = {
        "args": vars(args),
        "epsilons": list(epsilons),
        "assumptions": ASSUMPTIONS,
        "networks": summaries,
    }
    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    print(render_table(universe, summaries, all_ok, epsilons, args.condition))
    print(f"\n[final_table] re-rendered from {instances_path}", file=sys.stderr)


def main() -> None:
    """Run the full sweep: parse args, load inputs, certify, write, print."""
    args = build_arg_parser().parse_args()

    epsilons = tuple(float(e) for e in args.epsilons.split(",") if e.strip())
    skip = {s.strip() for s in args.skip.split(",") if s.strip()}
    requested = (
        {n.strip() for n in args.networks.split(",") if n.strip()} if args.networks else None
    )

    out_dir = ROOT / args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.report_only:
        _report_only(out_dir, epsilons, args)
        return

    knowledge = load_knowledge(args.condition)
    queries = load_queries(args.queries_per_network)
    parsed = load_networks(names=requested, skip=skip, max_nodes=args.max_nodes)

    universe = sorted(set(knowledge) & set(queries) & set(parsed) - skip)
    if requested is not None:
        missing = requested - set(universe)
        if missing:
            print(f"[final_table] requested but unavailable: {sorted(missing)}", file=sys.stderr)
        universe = [n for n in universe if n in requested]

    print(
        f"[final_table] condition={args.condition} networks={len(universe)} "
        f"epsilons={epsilons} seed={args.seed} time_limit={args.time_limit}",
        file=sys.stderr,
    )

    instances: list[dict[str, Any]] = []
    per_network: dict[str, list[dict[str, Any]]] = {}
    t_start = time.perf_counter()
    for net in universe:
        dag = parsed[net]["dag"]
        cpdag = parsed[net]["cpdag"]
        k = knowledge[net]
        net_queries = queries[net]
        print(
            f"[final_table] {net}: |V|={len(dag.nodes)} |K|={len(k)} "
            f"n_queries={len(net_queries)} elapsed={time.perf_counter() - t_start:.1f}s",
            file=sys.stderr,
        )
        bucket: list[dict[str, Any]] = []
        for x, y in net_queries:
            seed = derive_seed(args.seed, net, x, y)
            record = run_instance(
                args.condition, net, x, y, dag, cpdag, k, epsilons, seed, args.time_limit, 3
            )
            bucket.append(record)
            instances.append(record)
            print(
                f"[final_table]   {net} {x}->{y} status={record['status']} "
                f"r_val={record['r_val']} seconds={record['seconds']:.2f}"
                + (f" error={record['error']}" if record["error"] else ""),
                file=sys.stderr,
            )
        per_network[net] = bucket

    summaries = {
        net: summarise_network(len(parsed[net]["dag"].nodes), per_network[net], epsilons)
        for net in universe
    }
    all_ok = [r for r in instances if r["status"] == "ok"]

    instances_path = out_dir / "instances.jsonl"
    with instances_path.open("w") as f:
        for record in instances:
            f.write(json.dumps(record) + "\n")

    summary = {
        "args": vars(args),
        "epsilons": list(epsilons),
        "assumptions": ASSUMPTIONS,
        "networks": summaries,
    }
    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))

    table = render_table(universe, summaries, all_ok, epsilons, args.condition)
    print(table)
    print(f"\n[final_table] wrote {instances_path} and {summary_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
