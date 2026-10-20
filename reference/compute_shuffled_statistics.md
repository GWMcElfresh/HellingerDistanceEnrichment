# Ratio and energy distance after shuffling group labels on one draw's matrix.

Runs the same [`sample()`](https://rdrr.io/r/base/sample.html) calls as
compute_shuffled_ratios, so a Bayes draw consumes the same random
numbers whichever statistics it needs, and builds two group summaries
per shuffle: plain Hellinger distances for \\R\\ and squared distances
for \\E\\. The permutation test keeps the ratio-only helper, which costs
half as much per shuffle.

## Usage

``` r
compute_shuffled_statistics(contrast_group, dist_matrix, n_shuffles)
```

## Arguments

- contrast_group:

  Named factor aligned to `dist_matrix` rownames.

- dist_matrix:

  Symmetric subject Hellinger distance matrix.

- n_shuffles:

  Number of label shuffles.

## Value

Matrix with `n_shuffles` rows and columns `ratio` and `energy`.
