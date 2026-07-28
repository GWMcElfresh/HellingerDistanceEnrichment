# Plot Hellinger enrichment contrasts with uncertainty.

Takes a HellingerEnrichmentResult, maps each pairwise contrast onto the
between/within Hellinger ratio axis, and returns a ggplot. Bayes results
are drawn from the retained posterior effect-size draws (`nPosterior`
rows per contrast) via
[`ggdist::stat_halfeye`](https://mjskay.github.io/ggdist/reference/stat_halfeye.html),
because a density slab needs samples rather than summarized
`effectCiLow`/`effectCiHigh` intervals. Permutation results have no
draws, so they render as points colored by adjusted significance. When
`showOmnibus` is TRUE, the omnibus ratio is overlaid as a dashed
reference with its p-value labeled beside the line.

## Usage

``` r
PlotCompositionContrasts(
  result,
  theme = egg::theme_article(),
  showOmnibus = TRUE,
  ...
)
```

## Arguments

- result:

  A HellingerEnrichmentResult from CompareGroupCompositions.

- theme:

  ggplot2 theme (default egg::theme_article()).

- showOmnibus:

  If TRUE, annotate the omnibus effect size as a dashed line.

- ...:

  Passed to ggdist::stat_halfeye for Bayes plots.

## Value

A ggplot object.
