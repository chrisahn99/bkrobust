"""Answer the exported prompts with a local checkpoint through vLLM, on a cluster.

Every prompt is sent as one user message through the model's own chat
template, at temperature zero with a fixed seed and the same completion
budget the local suppliers had, so a record produced here is the same kind
of record ``elicit_run.py`` caches. Thinking is turned off where a template
offers it; the questionnaire asks for an order, not a derivation.

    python elicit_vllm.py --model <checkpoint> --model-id Qwen2.5-72B-Instruct
        --prompts prompts.jsonl --out answers.jsonl --tp 4
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path


def run_one(
    model: str,
    model_id: str,
    items: list[dict],
    out: Path,
    tp: int,
    max_tokens: int,
    no_think: bool,
    eager: bool = False,
) -> None:
    """Answer every item with one checkpoint and write the file atomically."""
    from vllm import LLM, SamplingParams

    extra: dict = {}
    if "gemma" in model_id.lower():
        extra["limit_mm_per_prompt"] = {"image": 0}  # text only; the questionnaire has no images
    llm = LLM(
        model=model,
        tensor_parallel_size=tp,
        dtype="bfloat16",
        seed=0,
        max_model_len=8192,
        gpu_memory_utilization=0.90,
        enforce_eager=eager,  # no compile cache on a quota-bound home directory
        **extra,
    )
    params = SamplingParams(temperature=0.0, seed=0, max_tokens=max_tokens)
    messages = [[{"role": "user", "content": it["prompt"]}] for it in items]
    kwargs = {"chat_template_kwargs": {"enable_thinking": False}} if no_think else {}
    t = time.perf_counter()
    outs = llm.chat(messages, params, use_tqdm=True, **kwargs)
    wall = time.perf_counter() - t
    tmp = out.with_suffix(".partial")
    with tmp.open("w") as fh:
        for it, o in zip(items, outs, strict=True):
            fh.write(
                json.dumps(
                    {
                        "prompt_sha256": it["prompt_sha256"],
                        "model": model_id,
                        "response": o.outputs[0].text,
                        "eval_count": len(o.outputs[0].token_ids),
                        "finish_reason": o.outputs[0].finish_reason,
                    }
                )
                + "\n"
            )
    tmp.rename(out)
    print(f"{len(items)} prompts, {wall:.0f} s, model {model_id}", flush=True)
    del llm
    import gc

    import torch

    gc.collect()
    torch.cuda.empty_cache()


def main() -> None:
    """Run the requested checkpoint, then any queued one whose output is missing."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--model-id", required=True, help="the name the cache key is built from")
    ap.add_argument("--prompts", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tp", type=int, default=4)
    ap.add_argument("--max-tokens", type=int, default=700)
    ap.add_argument("--no-think", action="store_true", help="pass enable_thinking=False")
    ap.add_argument("--eager", action="store_true", help="enforce_eager, no compile cache")
    ap.add_argument(
        "--queue",
        default="queue.txt",
        help="optional file of '<checkpoint> <model id> [--no-think]' lines to run after",
    )
    args = ap.parse_args()
    items = [json.loads(line) for line in Path(args.prompts).open()]
    out = Path(args.out)
    if out.exists():
        print(f"{out} exists, skipping {args.model_id}", flush=True)
    else:
        run_one(
            args.model,
            args.model_id,
            items,
            out,
            args.tp,
            args.max_tokens,
            args.no_think,
            args.eager,
        )
    q = Path(args.queue)
    if q.exists():
        for line in q.read_text().splitlines():
            parts = line.split()
            if len(parts) < 2:
                continue
            path, model_id = parts[0], parts[1]
            target = out.parent / f"{model_id}.jsonl"
            if target.exists():
                print(f"{target} exists, skipping {model_id}", flush=True)
                continue
            run_one(
                path,
                model_id,
                items,
                target,
                args.tp,
                args.max_tokens,
                "--no-think" in parts,
                args.eager,
            )


if __name__ == "__main__":
    main()
