# Hellinger distance enrichment for categorical composition data.

Compare subject-level categorical compositions (cluster counts,
cell-type proportions) across experimental groups using
Hellinger-distance enrichment ratios with permutation or conjugate
Bayesian posterior-predictive inference.

The enrichment ratio, Jeffreys softening, label-permutation null, and
collapse/subset contrasts follow methods developed by Paul Edlefsen. The
nested Dirichlet Bayes option implements his planned extension that
propagates uncertainty in estimated simplex compositions.

## References

Edlefsen, P. Original Hellinger between/within enrichment procedure for
subject-level categorical compositions (label permutation;
Jeffreys-softened Hellinger distances), with Bayesian composition
uncertainty as a subsequent step.
