"""Export the frozen questionnaire's prompts for a supplier that runs elsewhere.

One line per rendered item: the prompt, its SHA-256, and the item's frame
fields, so that a batch runner on a cluster can answer every item without the
questionnaire file and the answers can be imported back into the cache that
``elicit_run.py`` reads. The prompt text is exactly what the local suppliers
saw; the hash is the same one the cache records carry.

    python experiments/elicit_export.py            # -> results/elicit/prompts.jsonl
    python experiments/elicit_export.py --questionnaire results/pest/questionnaire.json \
        --out results/pest/prompts.jsonl
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ELICIT = ROOT / "results" / "elicit"


def export(questionnaire: Path, out: Path) -> int:
    """Write every rendered item's prompt with its hash; return the item count."""
    q = json.loads(questionnaire.read_text())
    with out.open("w") as fh:
        for it in q["items"]:
            sha = hashlib.sha256(it["prompt"].encode()).hexdigest()
            assert sha == it["prompt_sha256"], it["network"]
            fh.write(
                json.dumps(
                    {
                        "prompt_sha256": sha,
                        "network": it["network"],
                        "kind": it["kind"],
                        "naming": it["naming"],
                        "order_seed": it["order_seed"],
                        "prompt": it["prompt"],
                    }
                )
                + "\n"
            )
    return len(q["items"])


def main(argv: list[str] | None = None) -> None:
    """Export the questionnaire named on the command line, the ledger's by default."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--questionnaire", default=str(ELICIT / "questionnaire.json"))
    ap.add_argument("--out", default=str(ELICIT / "prompts.jsonl"))
    args = ap.parse_args(argv)
    n = export(Path(args.questionnaire), Path(args.out))
    print(f"{n} prompts -> {args.out}")


if __name__ == "__main__":
    main()
