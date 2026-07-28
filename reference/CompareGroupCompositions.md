# Compare compositions across groups

Tests between-versus-within Hellinger distance ratios with omnibus and
contrast statistics.

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

  `CategoryComposition` or long counts table.

- method:

  `"permutation"` or `"bayes"`.

- contrasts:

  Optional named custom contrasts.

- nPermutations:

  Permutation count.

- nPosterior:

  Posterior draw count for Bayes method.

- nCores:

  Parallel workers.

- seed:

  Random seed.

- priorPseudocounts:

  Jeffreys prior increment per category.

- pAdjustMethod:

  Contrast p-value adjustment method.

- ...:

  Unused.

## Value

A `HellingerEnrichmentResult` object.
