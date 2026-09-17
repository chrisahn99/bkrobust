import sys; sys.path.insert(0,"/Users/josecosta/bkrobust/src")
import numpy as np
from rnull_real import *
dag,W,var = load("ecoli70"); sem=sem_of(dag,W,var)
sig=sem.covariance(); idx={v:i for i,v in enumerate(dag.nodes)}
rng=np.random.default_rng(7)
X,ix = sample(dag,W,var,200000,rng)
sh=np.cov(X,rowvar=False)
print("ordem igual:", ix==idx)
print("erro relativo max Sigma-hat vs Sigma:", float(np.max(np.abs(sh-sig))/np.max(np.abs(sig))))
print("diag pop :", np.round(np.diag(sig)[:6],4))
print("diag hat :", np.round(np.diag(sh)[:6],4))
