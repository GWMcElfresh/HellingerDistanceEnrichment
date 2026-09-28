# Execute Bayes enrichment from Dirichlet posterior draws of R.

For each contrast, returns both a one-row posterior summary and the long
draw table used by half-eye plots. Draw storage scales as
`n_posterior * n_contrasts` effect sizes (omnibus excluded from the
returned `draws`, matching `contrasts`).

## Usage

``` r
run_bayes_enrichment(
  composition,
  contrast_specs,
  n_posterior,
  n_cores,
  prior_pseudocounts
)
```
