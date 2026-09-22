import csv, collections, statistics
from scipy.stats import binomtest, wilcoxon
rows=[u for u in csv.DictReader(open('results/axis_robustness_real/analysis_units.csv'))
      if u['arm']=='flip' and u['status']=='ok' and u['base_wrongness'] not in ('','None') and float(u['base_wrongness'])==0.0]
covs=collections.defaultdict(dict)
for u in rows: covs[(u['network'],u['x'],u['y'])][float(u['coverage'])]=u
matched={t:d for t,d in covs.items() if set(d)=={0.25,0.5,1.0}}
print("MATCHED-POPULATION COVERAGE CONTRAST ([RE-6]) -- 105 triples on 8 networks")
print("Conditioning on a fixed instance set removes the composition effect that makes")
print("the unmatched 543 -> 182 -> 106 comparison meaningless. It does NOT remove the")
print("non-nesting of the coverage sweep ([RE-4]), which is reported, not fixed.\n")
same=diff=0; dd=collections.Counter()
for t,d in matched.items():
    r=[int(float(d[c]['radius'])) for c in (0.25,0.5,1.0)]
    (same:=same+1) if len(set(r))==1 else (diff:=diff+1, dd.__setitem__(tuple(r), dd[tuple(r)]+1))
print(f"radius identical at all three coverages: {same} of {len(matched)} triples")
print(f"radius changes: {diff}; patterns (r@0.25, r@0.5, r@1.0): {dd.most_common()}")
for c in (0.25,0.5,1.0):
    a=[float(d[c]['AUC_frac_usable']) for d in matched.values()]
    print(f"coverage {c}: median AUC {statistics.median(a):.4f}  mean {statistics.mean(a):.4f}")
pairs=[(float(d[0.25]['AUC_frac_usable']),float(d[1.0]['AUC_frac_usable'])) for d in matched.values()]
up=sum(1 for a,b in pairs if b>a); dn=sum(1 for a,b in pairs if b<a)
print(f"\nper-triple AUC 1.0 vs 0.25: higher {up}, lower {dn} (n={len(pairs)})")
print(f"  triple-level sign test p={binomtest(up,up+dn).pvalue:.3g}, Wilcoxon p={wilcoxon([b-a for a,b in pairs]).pvalue:.3g}")
print("  NOTE: the triples come from only 8 networks, so a triple-level test is")
print("  anti-conservative. The network-level test below is the honest one.\n")
bynet=collections.defaultdict(list)
for t,d in matched.items(): bynet[t[0]].append(d)
print(f"{'network':18s} {'n':>4s} {'medAUC@0.25':>12s} {'medAUC@1.0':>11s} {'rises':>6s} {'medR@0.25':>10s} {'medR@1.0':>9s}")
nup=nr=0
for n in sorted(bynet):
    ds=bynet[n]
    a=statistics.median(float(d[0.25]['AUC_frac_usable']) for d in ds)
    b=statistics.median(float(d[1.0]['AUC_frac_usable']) for d in ds)
    ra=statistics.median(int(float(d[0.25]['radius'])) for d in ds)
    rb=statistics.median(int(float(d[1.0]['radius'])) for d in ds)
    nup+= b>a; nr += rb>ra
    print(f"{n:18s} {len(ds):4d} {a:12.3f} {b:11.3f} {('yes' if b>a else 'no'):>6s} {ra:10.1f} {rb:9.1f}")
print(f"\nnetworks whose median survival rises with coverage: {nup} of {len(bynet)}  (sign test p={binomtest(nup,len(bynet)).pvalue:.3g})")
print(f"networks whose median radius rises with coverage:   {nr} of {len(bynet)}")
