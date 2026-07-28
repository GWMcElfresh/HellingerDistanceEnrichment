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

  Number of label permutations (permutation method).

- nPosterior:

  Number of posterior draws (bayes method).

- nCores:

  Number of parallel workers (default 1).

- seed:

  Random seed for reproducibility.

- priorPseudocounts:

  Jeffreys prior increment per category (default 1/2).

- pAdjustMethod:

  Multiple-testing adjustment for contrasts (`"holm"`, `"BH"`,
  `"none"`).

- ...:

  Ignored.

## Value

A HellingerEnrichmentResult with omnibus summary, pairwise `contrasts`,
and (Bayes only) long-format posterior `draws` for plotting.
