# Sample one Dirichlet draw per subject from conjugate posteriors.

Sample one Dirichlet draw per subject from conjugate posteriors.

## Usage

``` r
sample_dirichlet_compositions(counts, prior_pseudocounts = 1/2)
```

## Arguments

- counts:

  Subjects-by-categories count matrix.

- prior_pseudocounts:

  Jeffreys or user-specified Dirichlet alpha increment.

## Value

Matrix of frequency draws (subjects by categories).
