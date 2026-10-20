# Execute Bayes enrichment from Dirichlet posterior draws of R and E.

For each contrast, returns both a one-row posterior summary and the long
draw table used by half-eye plots. Draw storage scales as
`n_posterior * n_contrasts` effect sizes, including omnibus, so the
posterior of the global ratio can be inspected alongside pairwise draws.

## Usage

``` r
run_bayes_enrichment(
  composition,
  contrast_specs,
  n_posterior,
  n_cores,
  prior_pseudocounts,
  credible_interval = "hpdi",
  nope = "shuffle",
  energy_nope = "shuffle",
  n_nope_shuffles = 200,
  nope_quantile = 0.95
)
```

## Details

`R_PPGTnope` compares each draw's \\R\\ with that same draw's \\R\\
edge, and `E_PPGTnope` does the same for \\E\\; `R_nopeHigh` and
`E_nopeHigh` in the summary are the medians of the per-draw edges.
`R_hdiAboveNope` is `TRUE` when the 95% HPDI of \\R\\ sits entirely
above `R_nopeHigh`; `E_hdiAboveNope` is the same comparison for \\E\\.
Shuffles consume the random stream between Dirichlet draws, so a fixed
seed gives different draws when `nope` or `energy_nope` is `"shuffle"`
than when neither is.
