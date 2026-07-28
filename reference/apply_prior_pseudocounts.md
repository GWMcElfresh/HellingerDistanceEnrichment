# Apply Jeffreys pseudocounts to a count matrix.

Apply Jeffreys pseudocounts to a count matrix.

## Usage

``` r
apply_prior_pseudocounts(counts, prior_pseudocounts = 1/2)
```

## Arguments

- counts:

  Subjects-by-categories count matrix.

- prior_pseudocounts:

  Scalar or per-category pseudocount (default 1/2).

## Value

Softened count matrix.
