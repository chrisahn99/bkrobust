"""Density/size sweep summary: does O*-locality survive dense graphs?"""
import json, sys
import numpy as np

RNG = np.random.default_rng(13)

def boot(v, B=2000):
    v = np.asarray(v, float)
    if not len(v): return (np.nan,)*3
    bs = v[RNG.integers(0, len(v), size=(B, len(v)))].mean(axis=1)
    return float(v.mean()), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))

def main():
    cells = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "../results/density_sweep.json"))
    out = {}
    rows = []
    for key in sorted(cells, key=lambda k: (int(k.split("_")[0][1:]), float(k.split("deg")[1]))):
        rs = [r for r in cells[key] if "generic" in r["arms"]]
        if len(rs) < 20: continue
        ch = inv = cons = tot = 0
        wn, ex = [], []
        for r in rs:
            arm = r["arms"]["generic"]
            for m in arm["members"]:
                if m["rho"] != 1: continue
                tot += 1
                if not m["consistent"]: continue
                cons += 1
                if not m["amenable"]: ch += 1; continue
                if m["ostar_changed"]:
                    ch += 1
                    if not m["ostar_valid"]: inv += 1
            ests = [arm["est0"]] + [m["est"] for m in arm["members"]
                                    if m["rho"] <= 1 and m["consistent"] and m.get("amenable")]
            wn.append((max(ests) - min(ests)) / abs(r["tau"]))
            ex.append(float(min(ests) > 0 or max(ests) < 0))
        p, d = key.split("_"); p = int(p[1:]); d = float(d[3:])
        w = boot(wn); e = boot(ex)
        rows.append(dict(p=p, deg=d, n_scm=len(rs),
                         n_undirected=float(np.mean([r["n_undirected"] for r in rs])),
                         frac_meek_inconsistent=round(1 - cons / tot, 4),
                         NUMBER1_changed=round(ch / cons, 4),
                         NUMBER1_invalid=round(inv / cons, 4),
                         width_mean=round(w[0], 3), width_ci=[round(w[1],3), round(w[2],3)],
                         width_median=round(float(np.median(wn)), 3),
                         excl0=round(e[0], 4), excl0_ci=[round(e[1],4), round(e[2],4)]))
    print(f"{'p':>3} {'deg':>4} {'n':>4} {'|U|':>5} {'inconsist':>10} {'NUM1 chg':>9} {'NUM1 inv':>9} {'width_med':>9} {'width_mean':>10} {'excl0':>7}")
    for r in rows:
        print(f"{r['p']:>3} {r['deg']:>4} {r['n_scm']:>4} {r['n_undirected']:>5.1f} "
              f"{r['frac_meek_inconsistent']:>10.3f} {r['NUMBER1_changed']:>9.3f} "
              f"{r['NUMBER1_invalid']:>9.3f} {r['width_median']:>9.3f} {r['width_mean']:>10.3f} {r['excl0']:>7.3f}")
    json.dump(rows, open("../results/summary_sweep.json", "w"), indent=1)

if __name__ == "__main__":
    main()
