"""A6-0c, the source-overlap probe: can a supplier name a network's source, cold?

The scrambled-name replicate tests recitation by invariance and did not
separate from real names. This probe asks the direct question. Each supplier
is shown the variable list of one network, under the display names the
questionnaire used and nothing else, and is asked for the publication the
network comes from and for the causal query that publication declares: the
exposure, the outcome and the adjustment set. The answer key is fixed here
before any answer is read: the author-year label the dagitty files carry and
their own ``exposure`` / ``outcome`` / ``adjusted`` annotations for the eight
applied papers; the reference the bnlearn repository documents for each
benchmark, which declares no query and is scored on the source alone.

Every call goes through ``elicit_run.ask``: temperature zero, seed zero,
cached under the hash of model and prompt. The prompts are also exported in
the format of ``elicit_export.py`` so the cluster suppliers can answer them
later and ``elicit_import.py`` can put the answers into the same cache. If the
GPU host is unreachable the export is written and the script stops.

Writes ``results/elicit/source_overlap_prompts.jsonl``,
``results/elicit/source_overlap.json`` and ``results/elicit/SOURCE_OVERLAP.md``.

    python experiments/source_overlap_probe.py
    python experiments/source_overlap_probe.py --models qwen2.5:7b-instruct,Qwen2.5-72B-Instruct
    python experiments/source_overlap_probe.py --export-only
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import time
import unicodedata
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

from bkrobust.benchmarks.describe import parse_file  # noqa: E402
from elicit_run import CACHE, CLUSTER_MODELS, ELICIT, HOST, MODELS, ask  # noqa: E402

SEED = 20260913
MAX_LISTED = 60
LOCAL_MODELS = ("qwen2.5:7b-instruct", "llama3.1:8b")

#: The answer key for the source. Applied dagitty papers: the author-year label
#: the file name carries and the dagitty example collection uses. bnlearn
#: networks: the reference the repository documents (bnlearn.com/bnrepository,
#: read 12 September 2026); the BIF headers say ``network unknown`` and the
#: acquisition manifest names only the pgmpy sdist. ``single_source`` is true
#: when the network has one canonical publication. ``names`` are the usual
#: names of the benchmark, used for an informational column only.
SOURCE_KEY: dict[str, dict] = {
    "Acid_1996": {
        "authors": ["Acid", "de Campos"],
        "year": 1996,
        "single_source": True,
        "key_from": "file name; dagitty example label 'Acid & de Campos, 1996'",
        "names": [],
    },
    "Didelez_2010": {
        "authors": ["Didelez", "Kreiner", "Keiding"],
        "year": 2010,
        "single_source": True,
        "key_from": "file name; dagitty example label 'Didelez et al, 2010'",
        "names": [],
    },
    "Kampen_2014": {
        "authors": ["van Kampen", "Kampen"],
        "year": 2014,
        "single_source": True,
        "key_from": "file name; dagitty example label 'van Kampen, 2014'",
        "names": [],
    },
    "Polzer_2012": {
        "authors": ["Polzer"],
        "year": 2012,
        "single_source": True,
        "key_from": "file name; dagitty example label 'Polzer et al., 2012'",
        "names": [],
    },
    "Schipf_2010": {
        "authors": ["Schipf"],
        "year": 2010,
        "single_source": True,
        "key_from": "file name; dagitty example label 'Schipf et al., 2010'",
        "names": [],
    },
    "Sebastiani_2005": {
        "authors": ["Sebastiani", "Ramoni", "Nolan", "Baldwin", "Steinberg"],
        "year": 2005,
        "single_source": True,
        "key_from": "file name; dagitty example label 'Sebastiani et al., 2005'",
        "names": [],
    },
    "Shrier_2008": {
        "authors": ["Shrier", "Platt"],
        "year": 2008,
        "single_source": True,
        "key_from": "file name; dagitty example label 'Shrier & Platt, 2008'",
        "names": [],
    },
    "Thoemmes_2013": {
        "authors": ["Thoemmes"],
        "year": 2013,
        "single_source": True,
        "key_from": "file name; dagitty example label 'Thoemmes, 2013'",
        "names": [],
    },
    "mediator": {
        "authors": [],
        "year": None,
        "single_source": False,
        "key_from": "dagitty textbook diagram 'Small model with mediator'; no publication",
        "names": [],
    },
    "paths": {
        "authors": [],
        "year": None,
        "single_source": False,
        "key_from": "dagitty textbook diagram 'Many variables but few paths'; no publication",
        "names": [],
    },
    "asia": {
        "authors": ["Lauritzen", "Spiegelhalter"],
        "year": 1988,
        "single_source": True,
        "key_from": "bnlearn: Lauritzen & Spiegelhalter, JRSS-B 50(2), 1988",
        "names": ["asia", "chest clinic"],
    },
    "alarm": {
        "authors": ["Beinlich", "Suermondt", "Chavez", "Cooper"],
        "year": 1989,
        "single_source": True,
        "key_from": "bnlearn: Beinlich, Suermondt, Chavez & Cooper, AIME 1989",
        "names": ["alarm"],
    },
    "barley": {
        "authors": ["Kristensen", "Rasmussen"],
        "year": 2002,
        "single_source": False,
        "key_from": "bnlearn credits the project (Kristensen, Rasmussen) and no paper; "
        "the year is that of Kristensen & Rasmussen's later publication",
        "names": ["barley"],
    },
    "child": {
        "authors": ["Spiegelhalter", "Cowell"],
        "year": 1992,
        "single_source": True,
        "key_from": "bnlearn: Spiegelhalter & Cowell, Bayesian Statistics 4, 1992",
        "names": ["child"],
    },
    "insurance": {
        "authors": ["Binder", "Koller", "Russell", "Kanazawa"],
        "year": 1997,
        "single_source": True,
        "key_from": "bnlearn: Binder, Koller, Russell & Kanazawa, Machine Learning 29, 1997",
        "names": ["insurance"],
    },
    "water": {
        "authors": ["Jensen", "Kjaerulff", "Olesen", "Pedersen"],
        "year": 1989,
        "single_source": True,
        "key_from": "bnlearn: Jensen, Kjaerulff, Olesen & Pedersen, technical report, 1989",
        "names": ["water"],
    },
    "hailfinder": {
        "authors": ["Abramson", "Brown", "Edwards", "Murphy", "Winkler"],
        "year": 1996,
        "single_source": True,
        "key_from": "bnlearn: Abramson et al., Int. J. Forecasting 12(1), 1996",
        "names": ["hailfinder"],
    },
    "hepar2": {
        "authors": ["Onisko"],
        "year": 2003,
        "single_source": True,
        "key_from": "bnlearn: Onisko, PhD dissertation, 2003",
        "names": ["hepar"],
    },
    "win95pts": {
        "authors": [],
        "year": None,
        "single_source": False,
        "key_from": "bnlearn gives no reference or attribution",
        "names": ["win95pts", "windows 95", "printer"],
    },
    "andes": {
        "authors": ["Conati", "Gertner", "VanLehn", "Druzdzel"],
        "year": 1997,
        "single_source": True,
        "key_from": "bnlearn: Conati, Gertner, VanLehn & Druzdzel, User Modeling 1997",
        "names": ["andes"],
    },
    "diabetes": {
        "authors": ["Andreassen", "Hovorka", "Benn", "Olesen", "Carson"],
        "year": 1991,
        "single_source": True,
        "key_from": "bnlearn: Andreassen, Hovorka, Benn, Olesen & Carson, AIME 1991",
        "names": ["diabetes"],
    },
    "link": {
        "authors": ["Jensen", "Kong"],
        "year": 1999,
        "single_source": True,
        "key_from": "bnlearn: Jensen & Kong, Am. J. Human Genetics, 1999",
        "names": ["link"],
    },
    "munin": {
        "authors": ["Andreassen", "Jensen", "Andersen", "Falck", "Kjaerulff", "Woldbye"],
        "year": 1989,
        "single_source": True,
        "key_from": "bnlearn: Andreassen et al., MUNIN, 1989",
        "names": ["munin"],
    },
    "munin1": {
        "authors": ["Andreassen", "Jensen", "Andersen", "Falck", "Kjaerulff", "Woldbye"],
        "year": 1989,
        "single_source": True,
        "key_from": "subnetwork of MUNIN; bnlearn: Andreassen et al., 1989",
        "names": ["munin"],
    },
    "munin2": {
        "authors": ["Andreassen", "Jensen", "Andersen", "Falck", "Kjaerulff", "Woldbye"],
        "year": 1989,
        "single_source": True,
        "key_from": "subnetwork of MUNIN; bnlearn: Andreassen et al., 1989",
        "names": ["munin"],
    },
    "munin3": {
        "authors": ["Andreassen", "Jensen", "Andersen", "Falck", "Kjaerulff", "Woldbye"],
        "year": 1989,
        "single_source": True,
        "key_from": "subnetwork of MUNIN; bnlearn: Andreassen et al., 1989",
        "names": ["munin"],
    },
    "munin4": {
        "authors": ["Andreassen", "Jensen", "Andersen", "Falck", "Kjaerulff", "Woldbye"],
        "year": 1989,
        "single_source": True,
        "key_from": "subnetwork of MUNIN; bnlearn: Andreassen et al., 1989",
        "names": ["munin"],
    },
    "pathfinder": {
        "authors": ["Heckerman", "Horvitz", "Nathwani"],
        "year": 1992,
        "single_source": True,
        "key_from": "bnlearn: Heckerman, Horvitz & Nathwani, Methods Inf. Med., 1992",
        "names": ["pathfinder"],
    },
    "sachs": {
        "authors": ["Sachs", "Perez", "Pe'er", "Lauffenburger", "Nolan"],
        "year": 2005,
        "single_source": True,
        "key_from": "bnlearn: Sachs, Perez, Pe'er, Lauffenburger & Nolan, Science 308, 2005",
        "names": ["sachs"],
    },
    "ecoli70": {
        "authors": ["Schaefer", "Schafer", "Strimmer"],
        "year": 2005,
        "single_source": True,
        "key_from": "bnlearn: Schaefer & Strimmer, Stat. Appl. Genet. Mol. Biol. 4, 2005",
        "names": ["ecoli", "e. coli", "e.coli"],
    },
    "arth150": {
        "authors": ["Opgen-Rhein", "Opgen", "Strimmer"],
        "year": 2007,
        "single_source": True,
        "key_from": "bnlearn: Opgen-Rhein & Strimmer, BMC Systems Biology 1(37), 2007",
        "names": ["arth", "arabidopsis"],
    },
    "magic-niab": {
        "authors": ["Scutari", "Howell", "Balding", "Mackay"],
        "year": 2014,
        "single_source": True,
        "key_from": "bnlearn: Scutari, Howell, Balding & Mackay, Genetics 198(1), 2014",
        "names": ["magic", "niab", "wheat"],
    },
    "magic-irri": {
        "authors": ["Scutari"],
        "year": 2016,
        "single_source": False,
        "key_from": "bnlearn credits a 2016 conference talk by Scutari and no paper",
        "names": ["magic", "irri", "rice"],
    },
}

_FENCE = re.compile(r"^```(?:json)?\s*|\s*```$", re.M)
_YEAR = re.compile(r"\b(1[89]\d\d|20\d\d)\b")
#: The grammar writes ``"year": <year>``; a model that declines writes a bare
#: ``DECLINE`` there, which is not JSON. It is quoted before parsing.
_BARE_DECLINE = re.compile(r'(?<!["\w])DECLINE(?!["\w])')
#: Field-wise recovery when the completion budget truncated the JSON: the
#: source fields come first in the grammar and survive; the list that ran
#: away is read as far as it got.
_FIELD_STR = re.compile(r'"(dataset_name|title|exposure|outcome)"\s*:\s*"([^"]*)"')
_FIELD_YEAR = re.compile(r'"year"\s*:\s*"?([^,"}\s]+)')
_FIELD_LIST = re.compile(r'"(authors|adjustment_set)"\s*:\s*\[([^\]]*)')
_STOP = {"de", "van", "der", "von", "den", "and", "et", "al", "the", "of"}


def fold(text: str) -> str:
    """Case-fold and strip diacritics and punctuation, for surname comparison."""
    text = unicodedata.normalize("NFKD", text)
    text = "".join(ch for ch in text if not unicodedata.combining(ch))
    return re.sub(r"[^a-z0-9 ]", " ", text.casefold())


def surname_tokens(names: list[str]) -> set[str]:
    """The comparable tokens of a list of author strings."""
    out: set[str] = set()
    for name in names:
        for tok in fold(name).split():
            if len(tok) >= 3 and tok not in _STOP:
                out.add(tok)
    return out


def render(n_total: int, listed: list[str]) -> str:
    """The cold prompt: a variable list, a count, and the question."""
    sampled = len(listed) < n_total
    lines = [
        "You are shown the variables of a dataset whose causal structure has been published "
        "as a directed acyclic graph.",
        f"The dataset has {n_total} variables. "
        + (
            f"A sample of {len(listed)} of them is listed below, in no particular order:"
            if sampled
            else "All of them are listed below, in no particular order:"
        ),
    ]
    lines += [f"  {v}" for v in listed]
    lines += [
        "",
        "From this list alone, identify the dataset and the publication it comes from: the "
        "authors' surnames, the year and the title. Then state the causal query that "
        "publication declares for it: the exposure (treatment) variable, the outcome "
        "variable, and the adjustment set it uses, each given by its name exactly as listed "
        "above. Answer DECLINE for any field you cannot fill; DECLINE is a valid answer.",
        "",
        'Answer ONLY with JSON: {"dataset_name": "<name>" | "DECLINE", "source": {"authors": '
        '["<surname>", ...], "year": <year>, "title": "<title>"} | "DECLINE", "exposure": '
        '"<name as listed>" | "DECLINE", "outcome": "<name as listed>" | "DECLINE", '
        '"adjustment_set": ["<name as listed>", ...] | "DECLINE", "confidence": 0.0}',
    ]
    return "\n".join(lines)


def build_prompts() -> list[dict]:
    """One prompt per network the questionnaire asked about, with its key."""
    q = json.loads((ELICIT / "questionnaire.json").read_text())
    asked = {p["network"] for p in q["per_network"] if p["items"] > 0}
    shown: dict[str, set[str]] = {}
    for it in q["items"]:
        if it["naming"] == "real":
            shown.setdefault(it["network"], set()).update(it["label_of"].values())

    out = []
    for path in sorted(MODELS.iterdir()):
        parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
        if parsed.name not in asked:
            continue
        raw = dict(parsed.name_map)
        display = {n: raw.get(n, n) for n in parsed.nodes}
        names = sorted(display.values())
        # a stable hash: the interpreter's is salted per process, and a prompt
        # rendered differently on each run is asked again on each run
        stable = int(hashlib.sha256(parsed.name.encode()).hexdigest()[:8], 16)
        rng = np.random.default_rng([SEED, 3, stable])
        if len(names) > MAX_LISTED:
            must = sorted(shown.get(parsed.name, set()))
            rest = [v for v in names if v not in set(must)]
            extra = rng.choice(len(rest), MAX_LISTED - len(must), replace=False)
            listed = must + [rest[int(i)] for i in sorted(extra)]
        else:
            listed = list(names)
        listed = [listed[int(i)] for i in rng.permutation(len(listed))]
        prompt = render(len(names), listed)

        roles = parsed.roles
        exposures = [display[n] for n, r in roles.items() if "exposure" in r]
        outcomes = [display[n] for n, r in roles.items() if "outcome" in r]
        adjusted = sorted(display[n] for n, r in roles.items() if "adjusted" in r)
        declared = bool(exposures and outcomes)
        key = SOURCE_KEY[parsed.name]
        m = len(listed)
        out.append(
            {
                "network": parsed.name,
                "n_variables": len(names),
                "n_listed": m,
                "listed": listed,
                "prompt": prompt,
                "prompt_sha256": hashlib.sha256(prompt.encode()).hexdigest(),
                "key": {
                    **key,
                    "source_scorable": bool(key["authors"]),
                    "query_scorable": declared,
                    "exposure": exposures[0] if declared else None,
                    "outcome": outcomes[0] if declared else None,
                    "adjusted": adjusted if declared else None,
                    "chance_exposure": 1.0 / m,
                    "chance_pair": 1.0 / (m * (m - 1)),
                },
            }
        )
    return out


def export(items: list[dict]) -> Path:
    """The prompts in ``elicit_export.py``'s format, for the cluster suppliers."""
    out = ELICIT / "source_overlap_prompts.jsonl"
    with out.open("w") as fh:
        for it in items:
            fh.write(
                json.dumps(
                    {
                        "prompt_sha256": it["prompt_sha256"],
                        "network": it["network"],
                        "kind": "source_overlap",
                        "naming": "real",
                        "order_seed": 0,
                        "prompt": it["prompt"],
                    }
                )
                + "\n"
            )
    return out


def _as_list(value: object) -> list[str] | None:
    """A JSON field that should be a list of names, tolerating a string."""
    if value is None:
        return None
    if isinstance(value, str):
        if value.strip().upper() == "DECLINE":
            return None
        items = [s.strip() for s in re.split(r"[,;]|\band\b|&", value) if s.strip()]
    elif isinstance(value, list):
        items = [str(v).strip() for v in value if str(v).strip()]
    else:
        items = [str(value)]
    kept = [s for s in items if s.upper() != "DECLINE"]
    # ``[]`` is an answer (the empty set); ``["DECLINE"]`` is a decline
    return None if items and not kept else kept


def _as_name(value: object) -> str | None:
    """A JSON field that should be one variable name."""
    if value is None or isinstance(value, (list, dict)):
        return None
    s = str(value).strip()
    return None if not s or s.upper() == "DECLINE" else s


def parse(response: str) -> dict:
    """Read the model's JSON into the five answer fields, ``None`` for a decline."""
    out: dict = {
        "parse_ok": False,
        "parse_mode": "failed",
        "truncated": False,
        "dataset_name": None,
        "authors": None,
        "year": None,
        "title": None,
        "source_declined": True,
        "exposure": None,
        "outcome": None,
        "adjustment_set": None,
        "confidence": None,
    }
    text = _FENCE.sub("", response.strip())
    start, end = text.find("{"), text.rfind("}")
    data: object = None
    if start >= 0 and end > start:
        candidate = text[start : end + 1]
        for mode, body in (
            ("json", candidate),
            ("json_quoted_decline", _BARE_DECLINE.sub('"DECLINE"', candidate)),
        ):
            try:
                data = json.loads(body)
            except json.JSONDecodeError:
                continue
            out["parse_mode"] = mode
            break
    if not isinstance(data, dict):
        body = _BARE_DECLINE.sub('"DECLINE"', text)
        fields = {m.group(1): m.group(2) for m in _FIELD_STR.finditer(body)}
        lists = {
            m.group(1): [s.strip().strip('"') for s in m.group(2).split(",") if s.strip('" ')]
            for m in _FIELD_LIST.finditer(body)
        }
        year = _FIELD_YEAR.search(body)
        if not fields and not lists and year is None:
            return out
        data = {
            "dataset_name": fields.get("dataset_name"),
            "source": {
                "authors": lists.get("authors"),
                "year": year.group(1) if year else None,
                "title": fields.get("title"),
            },
            "exposure": fields.get("exposure"),
            "outcome": fields.get("outcome"),
            "adjustment_set": lists.get("adjustment_set"),
        }
        out["parse_mode"] = "regex_fallback"
        out["truncated"] = True
    out["parse_ok"] = True
    out["dataset_name"] = _as_name(data.get("dataset_name"))
    src = data.get("source")
    if isinstance(src, dict):
        out["authors"] = _as_list(src.get("authors"))
        m = _YEAR.search(str(src.get("year", "")))
        out["year"] = int(m.group(1)) if m else None
        out["title"] = _as_name(src.get("title"))
        out["source_declined"] = not (out["authors"] or out["year"] or out["title"])
    elif isinstance(src, str) and src.strip().upper() != "DECLINE" and src.strip():
        m = _YEAR.search(src)
        out["year"] = int(m.group(1)) if m else None
        out["authors"] = [src]
        out["title"] = src
        out["source_declined"] = False
    out["exposure"] = _as_name(data.get("exposure"))
    out["outcome"] = _as_name(data.get("outcome"))
    out["adjustment_set"] = _as_list(data.get("adjustment_set"))
    try:
        out["confidence"] = float(data.get("confidence")) if "confidence" in data else None
    except (TypeError, ValueError):
        out["confidence"] = None
    return out


def score(item: dict, ans: dict) -> dict:
    """Source match, query match and the verdict for one network and one answer."""
    key = item["key"]
    listed = {fold(v).strip(): v for v in item["listed"]}
    rec: dict = {
        "source_scorable": key["source_scorable"],
        "source_declined": ans["source_declined"],
        "author_match": None,
        "year_match": None,
        "source_match": None,
        "fabricated": None,
        "benchmark_named": None,
        "query_scorable": key["query_scorable"],
        "exposure_match": None,
        "outcome_match": None,
        "exposure_listed": None,
        "outcome_listed": None,
        "adjustment_jaccard": None,
        "query_match": None,
        "verdict": "no",
    }
    if key["source_scorable"] and not ans["source_declined"]:
        got = surname_tokens(ans["authors"] or [])
        want = surname_tokens(key["authors"])
        rec["author_match"] = bool(got & want)
        rec["year_match"] = (
            ans["year"] is not None
            and key["year"] is not None
            and abs(ans["year"] - key["year"]) <= 1
        )
        rec["source_match"] = bool(rec["author_match"] and rec["year_match"])
        rec["fabricated"] = not rec["source_match"]
    elif key["source_scorable"]:
        rec["source_match"] = False
        rec["fabricated"] = False
    if key["names"] and ans["dataset_name"]:
        nm = fold(ans["dataset_name"])
        rec["benchmark_named"] = any(fold(n).strip() in nm for n in key["names"])
    elif key["names"]:
        rec["benchmark_named"] = False

    if key["query_scorable"]:
        for field in ("exposure", "outcome"):
            got_name = ans[field]
            if got_name is None:
                rec[f"{field}_match"] = False
                continue
            g = fold(got_name).strip()
            rec[f"{field}_listed"] = g in listed
            rec[f"{field}_match"] = g == fold(key[field]).strip()
        want_set = {fold(v).strip() for v in key["adjusted"]}
        if ans["adjustment_set"] is not None:
            got_set = {fold(v).strip() for v in ans["adjustment_set"]}
            union = got_set | want_set
            rec["adjustment_jaccard"] = 1.0 if not union else len(got_set & want_set) / len(union)
        need_adj = bool(want_set)
        rec["query_match"] = bool(
            rec["exposure_match"]
            and rec["outcome_match"]
            and (not need_adj or (rec["adjustment_jaccard"] or 0.0) >= 0.5)
        )

    if rec["source_match"] and rec["query_match"]:
        rec["verdict"] = "recited"
    elif rec["source_match"] or rec["exposure_match"] or rec["outcome_match"]:
        rec["verdict"] = "partial"
    return rec


def summarise(rows: list[dict]) -> dict:
    """The per-supplier summary the pre-registration scores."""
    src = [r for r in rows if r["source_scorable"]]
    named = [r for r in src if r["source_match"]]
    answered = [r for r in src if not r["source_declined"]]
    qry = [r for r in rows if r["query_scorable"]]
    qry_named = [r for r in qry if r["source_match"]]
    qry_unnamed = [r for r in qry if not r["source_match"]]
    return {
        "networks_probed": len(rows),
        "answered": sum(1 for r in rows if r["answered"]),
        "parse_ok": sum(1 for r in rows if r["answer"]["parse_ok"]),
        "parsed_after_quoting_decline": sum(
            1 for r in rows if r["answer"]["parse_mode"] == "json_quoted_decline"
        ),
        "truncated_by_budget": sum(1 for r in rows if r["answer"]["truncated"]),
        "source_scorable": len(src),
        "source_declined": sum(1 for r in src if r["source_declined"]),
        "source_answered": len(answered),
        "source_match": len(named),
        "source_match_networks": [r["network"] for r in named],
        "author_match_only": sum(1 for r in answered if r["author_match"] and not r["year_match"]),
        "year_match_only": sum(1 for r in answered if r["year_match"] and not r["author_match"]),
        "fabricated": sum(1 for r in answered if r["fabricated"]),
        "fabricated_fraction_of_answered": (
            round(sum(1 for r in answered if r["fabricated"]) / len(answered), 3)
            if answered
            else None
        ),
        "benchmark_named": sum(1 for r in rows if r["benchmark_named"]),
        "benchmark_named_networks": [r["network"] for r in rows if r["benchmark_named"]],
        "query_scorable": len(qry),
        "query_match": sum(1 for r in qry if r["query_match"]),
        "query_match_networks": [r["network"] for r in qry if r["query_match"]],
        "exposure_match": sum(1 for r in qry if r["exposure_match"]),
        "outcome_match": sum(1 for r in qry if r["outcome_match"]),
        "source_match_on_query_scorable": len(qry_named),
        "query_match_given_source_match": sum(1 for r in qry_named if r["query_match"]),
        "outcome_match_without_source": sum(1 for r in qry_unnamed if r["outcome_match"]),
        "exposure_match_without_source": sum(1 for r in qry_unnamed if r["exposure_match"]),
        "pair_match_without_source": sum(
            1 for r in qry_unnamed if r["exposure_match"] and r["outcome_match"]
        ),
        "chance_exposure_total": round(sum(r["chance_exposure"] for r in qry_unnamed), 3),
        "chance_pair_total": round(sum(r["chance_pair"] for r in qry_unnamed), 3),
        "verdicts": {
            v: sum(1 for r in rows if r["verdict"] == v) for v in ("recited", "partial", "no")
        },
        "recited_networks": [r["network"] for r in rows if r["verdict"] == "recited"],
        "single_source": sum(1 for r in rows if r["single_source"]),
        "cache_hits": sum(1 for r in rows if r.get("cache_hit")),
    }


def host_reachable() -> bool:
    """Whether the GPU host answers an ssh in batch mode."""
    try:
        proc = subprocess.run(
            ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=8", HOST, "true"],
            capture_output=True,
            timeout=20,
        )
    except (subprocess.TimeoutExpired, OSError):
        return False
    return proc.returncode == 0


def cached(model: str, prompt: str) -> dict | None:
    """The cache record for a model and prompt, or ``None``."""
    key = hashlib.sha256(f"{model}|{prompt}".encode()).hexdigest()
    path = CACHE / f"{key}.json"
    if not path.exists():
        return None
    rec = json.loads(path.read_text())
    rec["cache_hit"] = True
    return rec


def _yn(value: object) -> str:
    return "" if value is None else ("yes" if value else "no")


def _short_source(row: dict) -> str:
    if row["source_declined"]:
        return "DECLINE"
    authors = row["answer"]["authors"] or []
    first = authors[0] if authors else "?"
    year = row["answer"]["year"]
    return f"{first[:22]} {year if year else '?'}"


def write_markdown(payload: dict) -> Path:
    """The per-network tables and per-supplier summaries."""
    lines = [
        "# Source-overlap probe (A6-0c)",
        "",
        "One cold prompt per network: the variable list under the questionnaire's display "
        "names, no domain sentence, no states. The supplier names the publication and the "
        "declared exposure, outcome and adjustment set. Key and scoring rules: "
        "`PREREGISTRATION_SOURCE_OVERLAP.md`; generated by "
        "`experiments/source_overlap_probe.py`.",
        "",
        f"Networks probed: {payload['n_networks']}; source-scorable: "
        f"{payload['n_source_scorable']}; with a declared query: "
        f"{payload['n_query_scorable']}; single-source: {payload['n_single_source']}.",
        "",
    ]
    for model, block in payload["suppliers"].items():
        s = block["summary"]
        lines += [f"## {model}", ""]
        if not block["rows"] or s["answered"] == 0:
            lines += ["No answers in the cache for this model.", ""]
            continue
        lines += [
            f"- parsed {s['parse_ok']} of {s['answered']} answers "
            f"({s['parsed_after_quoting_decline']} after quoting a bare DECLINE, "
            f"{s['truncated_by_budget']} read field-wise after the token budget cut the JSON)",
            f"- source named (author and year within one): **{s['source_match']} of "
            f"{s['source_scorable']}** scorable networks"
            + (f" ({', '.join(s['source_match_networks'])})" if s["source_match"] else ""),
            f"- source declined on {s['source_declined']}; answered on {s['source_answered']}, "
            f"of which {s['fabricated']} do not match the key "
            f"(fabricated fraction {s['fabricated_fraction_of_answered']})",
            f"- author-only matches {s['author_match_only']}, year-only {s['year_match_only']}",
            f"- benchmark's usual name given (informational): {s['benchmark_named']}"
            + (f" ({', '.join(s['benchmark_named_networks'])})" if s["benchmark_named"] else ""),
            f"- declared query matched on {s['query_match']} of {s['query_scorable']}"
            + (f" ({', '.join(s['query_match_networks'])})" if s["query_match"] else "")
            + f"; exposure {s['exposure_match']}, outcome {s['outcome_match']}",
            f"- source named on {s['source_match_on_query_scorable']} of the query-scorable "
            f"networks; query matched on {s['query_match_given_source_match']} of those",
            f"- without the source: outcome matched on {s['outcome_match_without_source']}, "
            f"exposure on {s['exposure_match_without_source']}, both on "
            f"{s['pair_match_without_source']}; chance totals "
            f"{s['chance_exposure_total']} and {s['chance_pair_total']}",
            f"- verdicts: recited {s['verdicts']['recited']}, partial {s['verdicts']['partial']}, "
            f"no {s['verdicts']['no']}"
            + (f"; recited: {', '.join(s['recited_networks'])}" if s["recited_networks"] else ""),
            "",
            "| network | n | listed | single-source | key source | answered source | src | "
            "declared query | answered E / O | E | O | J | verdict |",
            "|---|---:|---:|---|---|---|---|---|---|---|---|---:|---|",
        ]
        for r in block["rows"]:
            k = r["key"]
            ksrc = f"{k['authors'][0]} {k['year']}" if k["authors"] else "none"
            kq = f"{k['exposure']} on {k['outcome']}" if k["query_scorable"] else "none"
            if k["query_scorable"] and k["adjusted"]:
                kq += " | {" + ", ".join(k["adjusted"]) + "}"
            aq = (
                f"{r['answer']['exposure'] or 'DECLINE'} / {r['answer']['outcome'] or 'DECLINE'}"
                if r["answered"]
                else ""
            )
            j = r["adjustment_jaccard"]
            lines.append(
                f"| {r['network']} | {r['n_variables']} | {r['n_listed']} | "
                f"{_yn(r['single_source'])} | {ksrc} | "
                f"{_short_source(r) if r['answered'] else 'not answered'} | "
                f"{_yn(r['source_match'])} | {kq} | {aq} | {_yn(r['exposure_match'])} | "
                f"{_yn(r['outcome_match'])} | {'' if j is None else f'{j:.2f}'} | "
                f"{r['verdict']} |"
            )
        lines.append("")
    out = ELICIT / "SOURCE_OVERLAP.md"
    out.write_text("\n".join(lines))
    return out


def main(argv: list[str] | None = None) -> None:
    """Render, export, ask the local suppliers, score, and write the tables."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", default=",".join(LOCAL_MODELS))
    ap.add_argument("--export-only", action="store_true")
    ap.add_argument(
        "--rescore",
        default="",
        help="re-parse and re-score the responses stored in an earlier source_overlap.json "
        "instead of asking anyone; for a parser change",
    )
    args = ap.parse_args(argv)
    models = [m for m in args.models.split(",") if m]

    items = build_prompts()
    out = export(items)
    print(f"{len(items)} prompts -> {out}", flush=True)
    if args.export_only:
        return
    stored: dict[str, dict[str, dict]] = {}
    if args.rescore:
        earlier = json.loads(Path(args.rescore).read_text())
        stored = {
            model: {r["network"]: r for r in block["rows"] if r.get("answered")}
            for model, block in earlier["suppliers"].items()
        }
        models = list(stored)
    local = [m for m in models if m not in CLUSTER_MODELS]
    reachable = host_reachable() if local and not stored else True
    if local and not reachable:
        print(f"{HOST} unreachable: prompts exported, no supplier asked", flush=True)
        return

    suppliers: dict[str, dict] = {}
    for model in models:
        rows = []
        t0 = time.perf_counter()
        for it in items:
            if stored:
                old = stored[model].get(it["network"])
                rec = (
                    {
                        "key": old["cache_key"],
                        "cache_hit": True,
                        "response": old["response"],
                        "prompt_sha256": old["prompt_sha256"],
                    }
                    if old
                    else None
                )
            else:
                rec = cached(model, it["prompt"]) if model in CLUSTER_MODELS else None
            if rec is None and (model in CLUSTER_MODELS or stored):
                ans = parse("")
                rows.append(
                    {
                        "network": it["network"],
                        "n_variables": it["n_variables"],
                        "n_listed": it["n_listed"],
                        "single_source": it["key"]["single_source"],
                        "chance_exposure": it["key"]["chance_exposure"],
                        "chance_pair": it["key"]["chance_pair"],
                        "key": it["key"],
                        "prompt_sha256": it["prompt_sha256"],
                        "answered": False,
                        "answer": ans,
                        **score(it, ans),
                    }
                )
                continue
            if rec is None:
                rec = ask(model, it["prompt"], num_predict=400)
            ans = parse(rec["response"])
            rows.append(
                {
                    "network": it["network"],
                    "n_variables": it["n_variables"],
                    "n_listed": it["n_listed"],
                    "single_source": it["key"]["single_source"],
                    "chance_exposure": it["key"]["chance_exposure"],
                    "chance_pair": it["key"]["chance_pair"],
                    "key": it["key"],
                    "prompt_sha256": rec.get("prompt_sha256") or it["prompt_sha256"],
                    "cache_key": rec["key"],
                    "cache_hit": rec["cache_hit"],
                    "answered": True,
                    "answer": ans,
                    "response": rec["response"],
                    **score(it, ans),
                }
            )
            print(
                f"  {model} {it['network']:16s} src={_yn(rows[-1]['source_match']) or '-':3s} "
                f"verdict={rows[-1]['verdict']:8s} {time.perf_counter() - t0:.0f} s",
                flush=True,
            )
        suppliers[model] = {"summary": summarise(rows), "rows": rows}
        print(f"== {model}: {json.dumps(suppliers[model]['summary'])}", flush=True)

    payload = {
        "seed": SEED,
        "max_listed": MAX_LISTED,
        "n_networks": len(items),
        "n_source_scorable": sum(1 for it in items if it["key"]["source_scorable"]),
        "n_query_scorable": sum(1 for it in items if it["key"]["query_scorable"]),
        "n_single_source": sum(1 for it in items if it["key"]["single_source"]),
        "prompts_sha256": hashlib.sha256(
            "".join(it["prompt_sha256"] for it in items).encode()
        ).hexdigest(),
        "rescored_from": Path(args.rescore).name if args.rescore else None,
        "suppliers": suppliers,
    }
    (ELICIT / "source_overlap.json").write_text(json.dumps(payload, indent=1))
    md = write_markdown(payload)
    print(f"-> {ELICIT / 'source_overlap.json'}\n-> {md}")


if __name__ == "__main__":
    main()
