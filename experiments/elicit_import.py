"""Import a cluster supplier's answers into the cache ``elicit_run.py`` reads.

Each answer line carries the prompt hash and the model id; the cache key is
the SHA-256 of ``model | prompt``, which needs the prompt text, so the
questionnaire is read to recover it. A record already in the cache is never
overwritten. Nothing here calls a model.

    python experiments/elicit_import.py results/elicit/remote/Qwen2.5-72B-Instruct.jsonl
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ELICIT = ROOT / "results" / "elicit"
CACHE = ELICIT / "cache"


def main(paths: list[str]) -> None:
    """Write one cache record per answer, keyed like the local supplier's."""
    q = json.loads((ELICIT / "questionnaire.json").read_text())
    prompt_of = {it["prompt_sha256"]: it["prompt"] for it in q["items"]}
    # the source-overlap probe exports its prompts in the same format; a cluster
    # answer to one of them enters the same cache under the same key
    extra = ELICIT / "source_overlap_prompts.jsonl"
    if extra.exists():
        for line in extra.open():
            it = json.loads(line)
            prompt_of.setdefault(it["prompt_sha256"], it["prompt"])
    CACHE.mkdir(parents=True, exist_ok=True)
    for p in paths:
        written = kept = truncated = 0
        for line in Path(p).open():
            a = json.loads(line)
            prompt = prompt_of[a["prompt_sha256"]]
            key = hashlib.sha256(f"{a['model']}|{prompt}".encode()).hexdigest()
            path = CACHE / f"{key}.json"
            truncated += a.get("finish_reason") == "length"
            if path.exists():
                kept += 1
                continue
            path.write_text(
                json.dumps(
                    {
                        "key": key,
                        "model": a["model"],
                        "prompt_sha256": a["prompt_sha256"],
                        "response": a["response"],
                        "eval_count": a.get("eval_count"),
                        "eval_ms": None,
                        "wall_s": None,
                        "cache_hit": False,
                        "source": "cluster_vllm",
                        "finish_reason": a.get("finish_reason"),
                    },
                    indent=1,
                )
            )
            written += 1
        print(f"{p}: {written} written, {kept} already cached, {truncated} hit the token budget")


if __name__ == "__main__":
    main(sys.argv[1:])
