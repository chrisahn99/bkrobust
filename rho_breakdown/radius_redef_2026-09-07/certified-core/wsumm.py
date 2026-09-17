import json, sys
import numpy as np
def s(tag):
    d=json.load(open(f"width_{tag}.json")); R=d["rows"]; n=len(R)
    m=np.array([r["m"] for r in R]); c=np.array([r["c"] for r in R])
    wM=np.array([r["wM"] for r in R]); wW=np.array([r["wW"] for r in R]); wC=np.array([r["wC"] for r in R])
    print(f"--- width {tag}: n={n} |K*| mean={m.mean():.2f}  c mean={c.mean():.3f} P(c=0)={np.mean(c==0):.4f}")
    print(f"  blanket-over-class width wC mean={wC.mean():.4f}   admissible-reversal width wW mean={wW.mean():.4f}  (equal on {np.mean(np.abs(wC-wW)<1e-9):.4f})")
    print(f"  certified-core hedge width wM mean={wM.mean():.4f}")
    print(f"  mean width ratio wW/wM (where wM>0) = {np.mean((wW/np.maximum(wM,1e-12))[wM>1e-12]):.3f}")
    print(f"  mean width SAVED  wW-wM = {np.mean(wW-wM):.4f}; wM/wW mean (wW>0) = {np.mean((wM/np.maximum(wW,1e-12))[wW>1e-12]):.4f}")
    print(f"  P(wM==0) = {np.mean(wM<1e-12):.4f}   P(wW==0) = {np.mean(wW<1e-12):.4f}")
    z=c==0
    if z.sum(): print(f"  on c=0 rows (n={z.sum()}): P(wM==0)={np.mean(wM[z]<1e-12):.4f}  P(point exact)={np.mean([r['covPoint'] for r,k in zip(R,z) if k]):.4f}  mean wW there={wW[z].mean():.4f}")
    nz=c>0
    if nz.sum(): print(f"  on c>0 rows (n={nz.sum()}): mean wM={wM[nz].mean():.4f} mean wW={wW[nz].mean():.4f} ratio of means={wW[nz].mean()/max(wM[nz].mean(),1e-12):.3f}")
    print(f"  STRUCTURAL COVERAGE of true effect: certified-core hedge {np.mean([r['covM'] for r in R]):.4f}  blanket {np.mean([r['covW'] for r in R]):.4f}  bare point estimate {np.mean([r['covPoint'] for r in R]):.4f}")
    nW=np.array([r["nW"] for r in R],float); nM=np.array([r["nM"] for r in R],float)
    print(f"  support |W|={nW.mean():.2f} |M|={nM.mean():.2f} ratio mean={np.mean(nW/nM):.3f}")
for t in sys.argv[1:]: s(t)
