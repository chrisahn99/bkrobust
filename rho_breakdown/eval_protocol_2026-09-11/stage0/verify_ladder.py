"""Recompute the free-rule (B2) agreement, the direction of its errors, and the law's sign, from instances.jsonl."""
import json, collections, os, pathlib, math
BK = pathlib.Path(os.path.expanduser("~/bkrobust"))
rows = [json.loads(l) for l in open(BK/"results/axisa3/instances.jsonl")]
adm = [r for r in rows if r.get("admissible")]
print("admissible", len(adm))
T0 = {"confounding","mediator","paths","M-bias"}
def b2(r):
    return r["separation"] if r["separation_status"]=="measured" else 1
exact = sum(1 for r in adm if b2(r)==r["radius"])
within1 = sum(1 for r in adm if abs(b2(r)-r["radius"])<=1)
print("B2 exact", exact, "within1", within1)
# direction of B2 errors
over = [r for r in adm if b2(r) > r["radius"]]   # free rule over-certifies (dangerous)
under = [r for r in adm if b2(r) < r["radius"]]  # free rule under-certifies (conservative)
print("B2 > r (over-certify, DANGEROUS):", len(over), collections.Counter(r["network"] for r in over))
print("B2 < r (under-certify, conservative):", len(under), "by network", collections.Counter(r["network"] for r in under).most_common(8))
print("  of which in undefined stratum:", sum(1 for r in under if r["separation_status"]!="measured"))
# undefined stratum: r distribution
und = [r for r in adm if r["separation_status"]!="measured"]
print("undefined stratum n", len(und), "r dist", collections.Counter(r["radius"] for r in und))
print("undefined stratum status reasons", collections.Counter(r["separation_status"] for r in und))
meas = [r for r in adm if r["separation_status"]=="measured"]
print("measurable n", len(meas))
law_eq = sum(1 for r in meas if r["radius"]==min(r["separation"], r["k_g0"]))
r_gt = sum(1 for r in meas if r["radius"]>min(r["separation"], r["k_g0"]))
r_lt = sum(1 for r in meas if r["radius"]<min(r["separation"], r["k_g0"]))
print("law: eq", law_eq, "r>min", r_gt, "r<min", r_lt, " => r>=min(s,k) holds in", law_eq+r_gt, "; r<=min(s,k) holds in", law_eq+r_lt)
print("r > k_g0 anywhere?", sum(1 for r in adm if r["radius"]>r["k_g0"]), "r == k_g0", sum(1 for r in adm if r["radius"]==r["k_g0"]))
print("min(s,k)==s in measurable:", sum(1 for r in meas if min(r["separation"],r["k_g0"])==r["separation"]))
# the 16 r<s cases
lt = [r for r in meas if r["radius"]<r["separation"]]
print("r<s cases:", [(r["network"],r["x"],r["y"],r["coverage"],r["radius"],r["separation"],r["k_g0"],r["method"]) for r in lt])
# |diff|>=2 set
d2 = [r for r in adm if abs(b2(r)-r["radius"])>=2]
print("|B2-r|>=2:", len(d2), collections.Counter(r["network"] for r in d2))
# outside T0
nt0 = [r for r in adm if r["network"] not in T0]
print("outside T0 n", len(nt0), "max r", max(r["radius"] for r in nt0), "r dist", sorted(collections.Counter(r["radius"] for r in nt0).items()))
print("outside T0 B2 exact", sum(1 for r in nt0 if b2(r)==r["radius"]), "B2>r", sum(1 for r in nt0 if b2(r)>r["radius"]))
# method vs radius
print("method x radius", collections.Counter((r["method"], r["radius"]) for r in adm))
# per coverage
for c in (1.0,0.5,0.25):
    s=[r for r in adm if r["coverage"]==c]
    print("coverage",c,"n",len(s),"r=1",sum(1 for r in s if r["radius"]==1),"B2 exact",sum(1 for r in s if b2(r)==r["radius"]),"B2>r",sum(1 for r in s if b2(r)>r["radius"]),"undefined",sum(1 for r in s if r["separation_status"]!="measured"), "median k_g0", sorted(r["k_g0"] for r in s)[len(s)//2])
# k_g0 vs coverage: how many with k_g0 <= 3
print("k_g0 dist", sorted(collections.Counter(r["k_g0"] for r in adm).items())[:15])
# distinct pairs and coverage triples
pairs = collections.defaultdict(dict)
for r in adm: pairs[(r["network"],r["x"],r["y"])][r["coverage"]] = r["radius"]
print("distinct admissible pairs", len(pairs), "admissible at all three coverages", sum(1 for v in pairs.values() if len(v)==3))
same3 = [v for v in pairs.values() if len(v)==3]
print("  among those: r identical across coverages", sum(1 for v in same3 if len(set(v.values()))==1), "r rises as coverage falls", sum(1 for v in same3 if v[0.25]>v[1.0]), "r falls", sum(1 for v in same3 if v[0.25]<v[1.0]))
# ICC quick: one-way ANOVA on r==1
groups = collections.defaultdict(list)
for r in adm: groups[r["network"]].append(1.0 if r["radius"]==1 else 0.0)
k=len(groups); N=len(adm); gm=sum(sum(g) for g in groups.values())/N
msb=sum(len(g)*(sum(g)/len(g)-gm)**2 for g in groups.values())/(k-1)
msw=sum(sum((x-sum(g)/len(g))**2 for x in g) for g in groups.values())/(N-k)
m0=(N-sum(len(g)**2 for g in groups.values())/N)/(k-1)
icc=(msb-msw)/(msb+(m0-1)*msw)
print("ICC r==1", round(icc,4), "k",k,"m0",round(m0,2))
