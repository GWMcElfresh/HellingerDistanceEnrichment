# Plot composition contrast effect sizes

Visualizes contrast-level enrichment ratios with ggdist uncertainty
intervals.

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

  `HellingerEnrichmentResult`.

- theme:

  ggplot2 theme (default
  [`egg::theme_article()`](https://rdrr.io/pkg/egg/man/theme_article.html)).

- showOmnibus:

  Annotate omnibus effect size.

- ...:

  Passed to ggdist geoms.

## Value

A ggplot object.
