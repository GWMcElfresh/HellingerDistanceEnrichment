# Compare categorical compositions across groups.

Tests whether between-group Hellinger distances exceed within-group
distances using a permutation null or conjugate Dirichlet
posterior-predictive Bayes. Computes an omnibus statistic plus all
pairwise group contrasts and optional custom contrasts (subsets or
collapses).

## Usage

``` r
CompareGroupCompositions(
  x,
  method = c("permutation", "bayes"),
  contrasts = NULL,
  nPermutations = 1000,
  nPosterior = 100,
  nCores = 1,
  seed = NULL,
  priorPseudocounts = 1/2,
  pAdjustMethod = "holm",
  ...
)
```

## Arguments

- x:

  A CategoryComposition object or long counts table.

- method:

  Inference method: `"permutation"` or `"bayes"`.

- contrasts:

  Optional named list of custom contrast specifications. Each element
  may be a length-2 character vector of group levels (pairwise subset)
  or a list with `collapse` (named remap of group levels) and/or
  `groups` (subset to these levels before testing).

- nPermutations:

  Number of label permutations (permutation method only).

- nPosterior:

  Number of posterior draws (bayes method).

- nCores:

  Number of parallel workers (default 1).

- seed:

  Random seed for reproducibility.

- priorPseudocounts:

  Jeffreys prior increment per category (default 1/2).

- pAdjustMethod:

  Multiple-testing adjustment for permutation contrasts (`"holm"`,
  `"BH"`, `"none"`). Ignored for Bayes.

- ...:

  Ignored.

## Value

A HellingerEnrichmentResult with omnibus summary, pairwise `contrasts`,
and (Bayes only) long-format posterior `draws` for plotting.

## Details

The between/within Hellinger ratio, Jeffreys softening,
observed-inclusive label-permutation p-value, and collapse/subset
contrasts follow Paul Edlefsen's original procedure. `method = "bayes"`
instead draws subject compositions from a conjugate Dirichlet posterior,
rebuilds the Hellinger ratio \\R\\ on each draw, and reports the
posterior probability that \\R \> 1\\ (`PPGT1`).

## References

Edlefsen, P. Hellinger between/within enrichment for subject-level
categorical compositions (Jeffreys softening, ratio statistic, label
permutation, collapse/subset contrasts).
