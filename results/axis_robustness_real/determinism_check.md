# Determinism Check

Snippet run twice, once with `PYTHONHASHSEED=0` and once with `PYTHONHASHSEED=12345`, both times with `PYTHONPATH=src /usr/bin/python3` from `/Users/ahn/Documents/Research/iclr27/bkrobust`.

## Output with `PYTHONHASHSEED=0`

```
81a5abbdc461082b7367713031a4ae207252a3c0d6407b13f8f999edd5d86b71
[["asia", [["asia", "tub"], ["smoke", "bronc"], ["smoke", "lung"]], [["bronc", "smoke"], ["smoke", "lung"], ["tub", "asia"]], [["asia", "tub"], ["bronc", "smoke"], ["lung", "smoke"]], [["bronc", "dysp"], ["bronc", "smoke"], ["either", "dysp"], ["either", "xray"], ["lung", "either"], ["smoke", "lung"], ["tub", "asia"], ["tub", "either"]]], ["sachs", [["PIP3", "PIP2"], ["PKA", "Jnk"], ["PKA", "P38"]
```
(truncated to first 400 chars of the JSON payload by the snippet itself, as specified)

## Output with `PYTHONHASHSEED=12345`

```
81a5abbdc461082b7367713031a4ae207252a3c0d6407b13f8f999edd5d86b71
[["asia", [["asia", "tub"], ["smoke", "bronc"], ["smoke", "lung"]], [["bronc", "smoke"], ["smoke", "lung"], ["tub", "asia"]], [["asia", "tub"], ["bronc", "smoke"], ["lung", "smoke"]], [["bronc", "dysp"], ["bronc", "smoke"], ["either", "dysp"], ["either", "xray"], ["lung", "either"], ["smoke", "lung"], ["tub", "asia"], ["tub", "either"]]], ["sachs", [["PIP3", "PIP2"], ["PKA", "Jnk"], ["PKA", "P38"]
```

## Result

Both runs produced the identical SHA-256 hash:
```
81a5abbdc461082b7367713031a4ae207252a3c0d6407b13f8f999edd5d86b71
```

`diff` between the two full captured stdout files (`det_seed0.txt` vs `det_seed12345.txt`) produced **no output** (exit code 0), confirming the two outputs are byte-identical.

**Determinism holds**: the two hashes match exactly. Output is independent of `PYTHONHASHSEED` for this snippet on this interpreter/library set.
