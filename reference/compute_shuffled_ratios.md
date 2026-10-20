# Between/within ratios after shuffling group labels on a fixed distance matrix.

Serves the permutation null on one fixed matrix from softened counts;
the Bayes NOPEs use compute_shuffled_statistics, which makes the same
[`sample()`](https://rdrr.io/r/base/sample.html) calls and adds \\E\\.
Shuffling reassigns labels among subjects, so group sizes and the
distance matrix stay fixed while the shuffle breaks any link between
label and composition.

## Usage

``` r
compute_shuffled_ratios(contrast_group, dist_matrix, n_shuffles)
```

## Arguments

- contrast_group:

  Named factor aligned to `dist_matrix` rownames.

- dist_matrix:

  Symmetric subject Hellinger distance matrix.

- n_shuffles:

  Number of label shuffles.

## Value

Numeric vector of `n_shuffles` shuffled ratios.
