# Does `r_eps` survive an estimated covariance?

120 instances, 120 bootstrap replicates, alpha = 0.05.

`r_val` is not in this table because it has no sampling distribution: it reads
the CPDAG and the asserted knowledge, and no numbers. Only `r_eps` is estimated.

The column that matters is **too large**. A radius that is too large claims more
certified shells than the truth supports, which overstates robustness; a radius
that is too small merely wastes some. The conservative construction exists to
drive the first rate below the nominal alpha, and is only useful if it does.

| cell | estimator | exact | too small (safe) | too large (overstates) | rate too large |
|---|---|---|---|---|---|
| n=100, eps=0.05 | plugin | 115/120 | 2 | 3 | 2.5% |
| n=100, eps=0.05 | conservative | 117/120 | 3 | 0 | 0.0% |
| n=100, eps=0.25 | plugin | 106/120 | 8 | 6 | 5.0% |
| n=100, eps=0.25 | conservative | 100/120 | 20 | 0 | 0.0% |
| n=500, eps=0.05 | plugin | 118/120 | 1 | 1 | 0.8% |
| n=500, eps=0.05 | conservative | 117/120 | 3 | 0 | 0.0% |
| n=500, eps=0.25 | plugin | 116/120 | 1 | 3 | 2.5% |
| n=500, eps=0.25 | conservative | 106/120 | 14 | 0 | 0.0% |
| n=2000, eps=0.05 | plugin | 118/120 | 0 | 2 | 1.7% |
| n=2000, eps=0.05 | conservative | 118/120 | 2 | 0 | 0.0% |
| n=2000, eps=0.25 | plugin | 119/120 | 0 | 1 | 0.8% |
| n=2000, eps=0.25 | conservative | 114/120 | 6 | 0 | 0.0% |
| n=10000, eps=0.05 | plugin | 119/120 | 0 | 1 | 0.8% |
| n=10000, eps=0.05 | conservative | 118/120 | 2 | 0 | 0.0% |
| n=10000, eps=0.25 | plugin | 120/120 | 0 | 0 | 0.0% |
| n=10000, eps=0.25 | conservative | 119/120 | 1 | 0 | 0.0% |