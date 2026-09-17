"""Stage 3: the questionnaire, derived from the frozen frame and hashed before anyone sees it.

One question per chain component that contains a candidate treatment of the
frozen frame. The unit of elicitation is the component, asked for a causal
ORDER over its variables rather than for individual edge directions, because
pairwise questions to an imperfect source produce cycles at scale and an order
cannot. The response grammar per adjacent pair is ``a->b``, ``b->a`` or
``DECLINE``, with a source-reported confidence that is retained and quarantined.

Each component is rendered under two naming conditions and two presentation
orders. Real names carry the data dictionary a practising analyst would have: a
domain sentence that never names the benchmark, and the variable's states where
the source file has them. Scrambled names replace every variable by a token and
the domain by a neutral sentence, so a source that reasons from the presented
structure gives the same answer and a source that recites a published graph
does not. The presentation order of nodes and pairs is a seeded permutation, so
two order seeds give two distinct prompts for the same component, which is what
lets instrument variance be separated from source variance.

A compelled-edge control block is attached to every network: a seeded sample of
edges the data already orient, rendered exactly like the open questions. A wrong
answer there is pure commission, detectable without any truth, and costs the
deployment path nothing.

Nothing here reads the true DAG. The CPDAG is the one the frame was drawn from
and carries that frame's provenance stamp.

Writes ``results/elicit/questionnaire.json`` and prints its hash.

    python experiments/questionnaire_build.py
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from bkrobust.benchmarks.describe import parse_file, to_mpdag  # noqa: E402
from bkrobust.demo.example import dag_to_cpdag  # noqa: E402
from bkrobust.demo.graph import MPDAG, undirected_components  # noqa: E402

MODELS = ROOT / "results" / "axisa3" / "networks" / "example_models"
FRAME = ROOT / "results" / "frame" / "frame.jsonl"
OUT = ROOT / "results" / "elicit"

MAX_COMPONENT = 20
CONTROL_EDGES_PER_NETWORK = 8
ORDER_SEEDS = (0, 1)

#: The domain sentence an analyst would have. It never names the benchmark, a
#: paper or an author, because the point of the scrambled arm is to detect a
#: source reciting a published graph, and naming it would hand that over.
DOMAIN: dict[str, str] = {
    "alarm": "anaesthesia monitoring during surgery: ventilation, circulation, blood gases and "
    "the alarms derived from them",
    "andes": "a student's problem-solving steps and knowledge in an introductory physics tutor",
    "asia": "patients at a chest clinic: travel history, smoking, tuberculosis, lung cancer, "
    "bronchitis, chest X-ray and shortness of breath",
    "barley": "a field trial of malting barley: soil, fertiliser, sowing, weather, pests, grain "
    "yield and grain quality",
    "child": "newborn infants referred with suspected congenital heart disease: history, "
    "examination findings, tests and diagnosis",
    "diabetes": "the physiology of a diabetic patient over one day: meals, insulin, glucose and "
    "related quantities measured hourly",
    "hailfinder": "severe summer weather forecasting on the high plains: atmospheric conditions, "
    "scenarios and forecasts",
    "hepar2": "liver disorders: risk factors, symptoms, laboratory tests and diagnoses",
    "insurance": "car insurance risk: driver characteristics, vehicle, accidents and costs",
    "link": "genetic linkage in a pedigree: genotypes and phenotypes",
    "munin": "electromyography for neuromuscular disorders: muscle and nerve findings",
    "munin1": "electromyography for neuromuscular disorders: muscle and nerve findings",
    "munin2": "electromyography for neuromuscular disorders: muscle and nerve findings",
    "munin3": "electromyography for neuromuscular disorders: muscle and nerve findings",
    "munin4": "electromyography for neuromuscular disorders: muscle and nerve findings",
    "sachs": "signalling in human immune cells: phosphorylation states of signalling proteins "
    "under perturbations",
    "water": "a biological wastewater treatment plant: inflow, process states and outflow",
    "win95pts": "troubleshooting printing from a personal computer: application, spooler, "
    "network, driver and printer",
    "ecoli70": "gene expression in the bacterium E. coli: transcript levels of named genes",
    "arth150": "gene expression in the plant Arabidopsis thaliana: transcript levels of named "
    "genes",
    "magic-irri": "rice: genetic markers and agronomic phenotypes in a multi-parent population",
    "magic-niab": "wheat: genetic markers and agronomic phenotypes in a multi-parent population",
    "Acid_1996": "an applied epidemiological study; the variables are given by their published "
    "codes",
    "Didelez_2010": "hormone replacement therapy and venous thromboembolism in women: age, "
    "smoking, occupation and selection into the study",
    "Kampen_2014": "an applied clinical study; the variables are given by their published codes",
    "Polzer_2012": "tooth loss and mortality: alcohol, smoking, diabetes and socio-economic "
    "factors",
    "Schipf_2010": "type 2 diabetes in men: testosterone, physical activity, waist "
    "circumference, smoking and age",
    "Sebastiani_2005": "stroke in sickle-cell anaemia: gene polymorphisms and clinical outcome",
    "Shrier_2008": "sports injury: warm-up, fitness, connective tissue, genetics, coaching and "
    "team motivation",
    "Thoemmes_2013": "a structural model with coded variables",
    "mediator": "a textbook causal diagram with abstract variables",
    "paths": "a textbook causal diagram with abstract variables",
    "confounding": "a textbook causal diagram with abstract variables",
    "M-structure": "a textbook causal diagram with abstract variables",
}
NEUTRAL_DOMAIN = "an observational dataset; the variables are given by neutral codes"

_BIF_VAR = re.compile(
    r"variable\s+(\S+)\s*\{\s*type\s+discrete\s*\[\s*\d+\s*\]\s*\{([^}]*)\}", re.S
)


def states_from_bif(path: Path) -> dict[str, list[str]]:
    """Variable states from a BIF file, which is the nearest thing it has to a glossary."""
    text = gzip.open(path, "rt").read() if path.suffix == ".gz" else path.read_text()
    out: dict[str, list[str]] = {}
    for m in _BIF_VAR.finditer(text):
        out[m.group(1)] = [s.strip() for s in m.group(2).split(",") if s.strip()]
    return out


def render(
    domain: str,
    names: dict[str, str],
    nodes: list[str],
    pairs: list[tuple[str, str]],
    states: dict[str, list[str]],
) -> str:
    """The prompt. Order of ``nodes`` and ``pairs`` is the caller's, already permuted."""
    lines = [
        "You are asked about causal direction among variables from " + domain + ".",
        "Some of these variables are directly related, but the data alone cannot tell which "
        "of each related pair causes the other.",
        "",
        "Variables in this group:",
    ]
    for n in nodes:
        label = names[n]
        st = states.get(n)
        lines.append(f"  {label}" + (f" (values: {', '.join(st)})" if st else ""))
    lines += ["", "Directly related pairs whose direction is unknown:"]
    for a, b in pairs:
        lines.append(f"  {names[a]} -- {names[b]}")
    lines += [
        "",
        "Give a causal ORDER of the variables in this group, earliest cause first, so that "
        "every related pair is oriented from the earlier variable to the later one. Omit a "
        "variable if you cannot place it. Then, for each listed pair, state the direction "
        "your order implies, or DECLINE if you cannot judge it, with a confidence from 0 to 1.",
        "",
        'Answer ONLY with JSON: {"order": ["...", "..."], "pairs": [{"a": "<first name as '
        'listed>", "b": "<second name as listed>", "direction": "a->b" | "b->a" | "DECLINE", '
        '"confidence": 0.0}, ...]}',
    ]
    return "\n".join(lines)


def stable_seed(name: str) -> int:
    """A per-network seed that does not depend on the interpreter's hash salt."""
    return int(hashlib.sha256(name.encode()).hexdigest()[:8], 16)


def load_cpdag_source(path: Path) -> dict[str, MPDAG]:
    """Estimated CPDAGs, one per network, from the JSON the estimated panel writes.

    The file is ``{"meta": {...}, "networks": {name: {"nodes", "directed",
    "undirected", ...}}}``; only the three graph fields are read here.
    """
    data = json.loads(path.read_text())
    return {
        name: MPDAG(
            e["nodes"],
            directed=[tuple(x) for x in e["directed"]],
            undirected=[tuple(x) for x in e["undirected"]],
        )
        for name, e in data["networks"].items()
    }


def build_network_items(
    name: str,
    cpdag: MPDAG,
    raw: dict[str, str],
    states: dict[str, list[str]],
    xs: set[str],
    seed: int,
    net_seed: int,
) -> tuple[list[dict], dict]:
    """Every rendered item for one network, and its per-network summary row.

    ``xs`` are the frame's treatments on this network, ``raw`` maps normalised
    node names back to the file's, ``states`` the variable states where the file
    has them. ``net_seed`` salts every per-network permutation and is the
    caller's to fix: the frozen questionnaire used the interpreter's string hash,
    the estimated panel uses :func:`stable_seed`.
    """
    real = {n: raw.get(n, n) for n in cpdag.nodes}
    rng_scr = np.random.default_rng([seed, net_seed])
    perm = rng_scr.permutation(len(cpdag.nodes))
    scrambled = {n: f"V{int(perm[i]) + 1}" for i, n in enumerate(sorted(cpdag.nodes))}

    comps = [
        c
        for c in undirected_components(cpdag)
        if len(c) >= 2 and (c & xs or any(cpdag.has_edge(x, v) for x in xs for v in c))
    ]
    kept = [c for c in comps if len(c) <= MAX_COMPONENT]
    dropped = [len(c) for c in comps if len(c) > MAX_COMPONENT]

    # compelled-edge control block: directed edges of the CPDAG touching a kept component
    touching = set().union(*kept) if kept else set()
    compelled = sorted((a, b) for (a, b) in cpdag.directed_edges if a in touching or b in touching)
    rng_ctl = np.random.default_rng([seed, 7, net_seed])
    if len(compelled) > CONTROL_EDGES_PER_NETWORK:
        idx = sorted(rng_ctl.choice(len(compelled), CONTROL_EDGES_PER_NETWORK, replace=False))
        compelled = [compelled[int(i)] for i in idx]

    items: list[dict] = []
    n_items = 0
    for ci, comp in enumerate(sorted(kept, key=lambda c: sorted(c))):
        nodes = sorted(comp)
        pairs = sorted(tuple(sorted((a, b))) for (a, b) in cpdag.undirected_edges if a in comp)
        for naming, names, domain in (
            ("real", real, DOMAIN.get(name, "an applied study")),
            ("scrambled", scrambled, NEUTRAL_DOMAIN),
        ):
            first_nodes: list[str] | None = None
            for order_seed in ORDER_SEEDS:
                rng = np.random.default_rng([seed, order_seed, ci, net_seed])
                node_order = [nodes[int(i)] for i in rng.permutation(len(nodes))]
                pair_order = [pairs[int(i)] for i in rng.permutation(len(pairs))]
                # A small component has few orderings, so a second seed can land
                # on the first seed's permutation. The two prompts must differ,
                # or instrument variance is measured against a copy of itself.
                if first_nodes is None:
                    first_nodes = node_order
                elif node_order == first_nodes:
                    node_order = list(reversed(node_order))
                    pair_order = list(reversed(pair_order))
                prompt = render(
                    domain, names, node_order, pair_order, states if naming == "real" else {}
                )
                items.append(
                    {
                        "network": name,
                        "component": ci,
                        "kind": "open",
                        "naming": naming,
                        "order_seed": order_seed,
                        "nodes": nodes,
                        "pairs": pairs,
                        "label_of": {n: names[n] for n in nodes},
                        "prompt": prompt,
                        "prompt_sha256": hashlib.sha256(prompt.encode()).hexdigest(),
                    }
                )
                n_items += 1
    # control block, real names, order seed 0, rendered as if the edges were open
    if compelled:
        cnodes = sorted({v for e in compelled for v in e})
        cpairs = [tuple(sorted(e)) for e in compelled]
        for naming, names, domain in (
            ("real", real, DOMAIN.get(name, "an applied study")),
            ("scrambled", scrambled, NEUTRAL_DOMAIN),
        ):
            rng = np.random.default_rng([seed, 99, net_seed])
            node_order = [cnodes[int(i)] for i in rng.permutation(len(cnodes))]
            pair_order = [cpairs[int(i)] for i in rng.permutation(len(cpairs))]
            prompt = render(
                domain, names, node_order, pair_order, states if naming == "real" else {}
            )
            items.append(
                {
                    "network": name,
                    "component": -1,
                    "kind": "compelled_control",
                    "naming": naming,
                    "order_seed": 0,
                    "nodes": cnodes,
                    "pairs": cpairs,
                    "compelled_direction": {f"{a}|{b}": [a, b] for (a, b) in compelled},
                    "label_of": {n: names[n] for n in cnodes},
                    "prompt": prompt,
                    "prompt_sha256": hashlib.sha256(prompt.encode()).hexdigest(),
                }
            )
    per_network = {
        "network": name,
        "components_with_candidate": len(comps),
        "components_kept": len(kept),
        "components_dropped_for_size": dropped,
        "open_questions": sum(
            len(c)
            for c in kept and [[e for e in cpdag.undirected_edges if e[0] in c] for c in kept]
        ),
        "control_edges": len(compelled),
        "items": n_items,
    }
    return items, per_network


def main(argv: list[str] | None = None) -> None:
    """Build and freeze the questionnaire."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=20260913)
    ap.add_argument("--frame", default=str(FRAME), help="frozen frame the treatments come from")
    ap.add_argument(
        "--cpdag-source",
        default="",
        help="JSON of estimated CPDAGs per network; by default the true DAG's CPDAG",
    )
    ap.add_argument("--out", default=str(OUT), help="output directory")
    args = ap.parse_args(argv)
    out = Path(args.out)

    frame = [json.loads(line) for line in Path(args.frame).open()]
    cands: dict[str, set[str]] = {}
    for r in frame:
        cands.setdefault(r["network"], set()).add(r["x"])
    estimated = load_cpdag_source(Path(args.cpdag_source)) if args.cpdag_source else None

    items: list[dict] = []
    per_network: list[dict] = []
    for path in sorted(MODELS.iterdir()):
        parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
        if parsed.bidirected or parsed.name not in cands:
            continue
        if estimated is None:
            cpdag: MPDAG = dag_to_cpdag(to_mpdag(parsed))
            # The frozen questionnaire of 2026-09-12 was built with Python's salted hash() here,
            # so it cannot be regenerated byte for byte; the frozen file stays authoritative.
            net_seed = stable_seed(parsed.name)
        elif parsed.name in estimated:
            cpdag = estimated[parsed.name]
            net_seed = stable_seed(parsed.name)
        else:
            continue
        raw = dict(parsed.name_map)
        states = states_from_bif(path) if parsed.source_format == "bif" else {}
        states = {n: states[raw.get(n, n)] for n in cpdag.nodes if raw.get(n, n) in states}
        net_items, per = build_network_items(
            parsed.name, cpdag, raw, states, cands[parsed.name], args.seed, net_seed
        )
        items += net_items
        per_network.append(per)

    # the startup assertion the protocol asks for: two order seeds, two distinct keys
    by_key: dict[tuple, set[str]] = {}
    for it in items:
        by_key.setdefault((it["network"], it["component"], it["naming"]), set()).add(
            it["prompt_sha256"]
        )
    collisions = [k for k, v in by_key.items() if len(v) < len(ORDER_SEEDS) and k[1] >= 0]
    assert not collisions, f"order seeds did not produce distinct prompts: {collisions[:3]}"

    out.mkdir(parents=True, exist_ok=True)
    payload = {
        "seed": args.seed,
        "max_component": MAX_COMPONENT,
        "order_seeds": list(ORDER_SEEDS),
        "control_edges_per_network": CONTROL_EDGES_PER_NETWORK,
        "cpdag_source": (
            f"estimated CPDAG from {Path(args.cpdag_source).name}; true_dag_on_path = false"
            if args.cpdag_source
            else "dag_to_cpdag(true DAG); true_dag_on_path = via_Chat_only"
        ),
        "per_network": per_network,
        "items": items,
    }
    text = json.dumps(payload, indent=1, sort_keys=True)
    digest = hashlib.sha256(text.encode()).hexdigest()
    frozen = out / "questionnaire.sha256"
    if frozen.exists() and frozen.read_text().strip() != digest:
        # Never overwrite a frozen questionnaire that suppliers have answered. The one frozen on
        # 2026-09-12 used a process-salted seed for the scrambled permutation and the control
        # sample, so a rebuild differs from it by construction; write the rebuild beside it.
        (out / "questionnaire.rebuilt.json").write_text(text)
        (out / "questionnaire.rebuilt.sha256").write_text(digest + "\n")
        print(f"frozen questionnaire kept; rebuild written beside it ({digest[:12]})")
        return
    (out / "questionnaire.json").write_text(text)
    frozen.write_text(digest + "\n")
    opens = [i for i in items if i["kind"] == "open"]
    print(
        f"{len(items)} rendered items over {len(per_network)} networks: "
        f"{len({(i['network'], i['component']) for i in opens})} components x "
        f"{len(ORDER_SEEDS)} orders x 2 namings = {len(opens)} open items, "
        f"{len(items) - len(opens)} control items; "
        f"{sum(len(i['pairs']) for i in opens) // (2 * len(ORDER_SEEDS))} open pairs"
    )
    print("questionnaire sha256", digest)


if __name__ == "__main__":
    main()
