import json, sys
from collections import Counter
import numpy as np
from scipy.stats import spearmanr

def summarise(tag):
    d = json.load(open(f"sweep_{tag}.json")); R = d["rows"]
    n = len(R); m = np.array([r["m"] for r in R]); c = np.array([r["c"] for r in R])
    kap = np.array([r["kappa"] for r in R]); phi1 = np.array([r["phi1"] for r in R])
    print(f"--- sweep {tag}: n={n} secs={d['secs']} tries={d['tries']}")
    print(f"  |K*| mean={m.mean():.2f} min={m.min()} max={m.max()}")
    print(f"  c: distinct={len(set(c.tolist()))} dist={sorted(Counter(c.tolist()).items())}")
    print(f"     mean={c.mean():.3f} P(c=0)={np.mean(c==0):.4f} P(c=|K*|)={np.mean(c==m):.4f} modal_mass={max(Counter(c.tolist()).values())/n:.4f}")
    print(f"  kappa: distinct={len(set(np.round(kap,6).tolist()))} mean={kap.mean():.3f} P(k=1)={np.mean(kap==1):.4f} P(k=0)={np.mean(kap==0):.4f}")
    print(f"  c == phi1 (singleton-only unsafe): {np.mean(c==phi1):.4f}")
    mm = np.array([r["max_min_size"] for r in R])
    print(f"  P(some minimal unsafe set size>=2) = {np.mean(mm>=2):.4f}")
    ch = np.array([r["n_checks"] for r in R])
    print(f"  validity checks/problem: mean={ch.mean():.2f} max={ch.max()}  vs 2^m mean={np.mean(2.0**m):.2f}  3^m mean={np.mean(3.0**m):.2f}")
    rc = np.array([r.get("retract_changed", False) for r in R])
    print(f"  retraction arm changed c on {rc.sum()}/{n}")
    thm = [r["thm"] for r in R if r["thm"] is not None]
    print(f"  c=0 theorem (Z valid in bare CPDAG): {sum(thm)}/{len(thm)}")
    gp = np.array([r["gap"] for r in R])
    print(f"  conditional-vs-joint gap fires on {np.mean(gp):.4f}")
    uq = np.array([r["unique_core"] for r in R])
    print(f"  core unique on {np.mean(uq):.4f}")
    W = np.array([r["W"] for r in R], float); M = np.array([r["M"] for r in R], float)
    pw = np.array([r["pow_m"] for r in R], float); pc = np.array([r["pow_c"] for r in R], float)
    print(f"  hedge SUPPORT: full 2^m={pw.mean():.2f} admissible |W|={W.mean():.2f} must-hedge |M|={M.mean():.2f} claimed 2^c={pc.mean():.2f}")
    print(f"     honest support ratio |W|/|M| mean={np.mean(W/M):.3f}  claimed 2^m/2^c mean={np.mean(pw/pc):.3f}")
    if np.std(kap)>0 and np.std(phi1/m)>0:
        print(f"  spearman(kappa, phi1/|K*|) = {spearmanr(kap, phi1/m).correlation:.3f}")
    pad = [r for r in R if "pad_c" in r]
    if pad:
        pc0 = np.array([r["c"] for r in pad]); pc1 = np.array([r["pad_c"] for r in pad])
        pk0 = np.array([r["kappa"] for r in pad]); pk1 = np.array([r["pad_kappa"] for r in pad])
        same = np.array([r["pad_same_G0"] for r in pad])
        print(f"  PADDING (n={len(pad)}, G0 identical on {np.mean(same):.4f}):")
        print(f"     c changed on {np.mean(pc1!=pc0):.4f}; kappa rose on {np.mean(pk1>pk0):.4f}, fell on {np.mean(pk1<pk0):.4f}, mean delta={np.mean(pk1-pk0):+.4f} max={np.max(pk1-pk0):+.4f}")
        print(f"     mean |K| {np.mean([r['m'] for r in pad]):.2f} -> {np.mean([r['pad_m'] for r in pad]):.2f}; mean kappa {pk0.mean():.3f} -> {pk1.mean():.3f}")

for t in sys.argv[1:]:
    summarise(t)
