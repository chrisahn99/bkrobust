import sys; sys.path.insert(0,'.')
from audit_w0_rk import *
from audit_w0_4_widths import methods_for, calibrate
# realistic: the analyst does not know how many of her claims are false.
rows=[]
for nf,share in ((0,160),(1,160),(2,80)):
    rng=np.random.default_rng(6000+nf); rng2=np.random.default_rng(7000+nf)
    ps,_=gen(rng, share, n=7, p=0.40, k_claims=3, n_false=nf)
    rows += [methods_for(pr,rng2) for pr in ps]
print(f"MIXTURE over f in {{0,1,2}}  (n={len(rows)})")
print(f"{'method':10s} {'raw cov(c=1)':>13s} {'calib factor':>13s} {'MEAN WIDTH @95%':>17s}")
best=None
for k in ('point_Z','S_t=0','S_t=1','S_t=2','blanket'):
    c1,c,w = calibrate(rows,k)
    print(f"{k:10s} {c1:13.3f} {c:13.2f} {w:17.4f}")
