# Plot Hellinger enrichment contrasts with uncertainty.

Takes a HellingerEnrichmentResult and returns a ggplot with one row per
pairwise contrast. By default the plot stacks two panels that share
those rows: the between/within Hellinger ratio \\R\\ on top and the
energy distance \\E\\ below, each with its own x axis, so both read left
to right. \\R\\ carries the p-values, `PPGT1`, and `R_PPGTnope`; \\E\\
shows the size of the shift in squared Hellinger units. A dotted line
marks each statistic's no-difference value (\\R = 1\\, \\E = 0\\). Set
`showEnergyDistance = FALSE` for the ratio alone.

## Usage

``` r
PlotCompositionContrasts(
  result,
  theme = egg::theme_article(),
  showOmnibus = TRUE,
  showEnergyDistance = TRUE,
  ...
)
```

## Arguments

- result:

  A HellingerEnrichmentResult from CompareGroupCompositions.

- theme:

  ggplot2 theme (default egg::theme_article()). With the default, the
  legend sits above the panels; a theme you supply keeps its own legend
  placement.

- showOmnibus:

  If TRUE, mark the omnibus value in each panel with a dashed line.

- showEnergyDistance:

  If TRUE (default), stack an energy distance panel below the ratio
  panel.

- ...:

  Passed to ggdist::stat_halfeye for Bayes plots.

## Value

A ggplot object.

## Details

Bayes results draw from the retained posterior draws (`nPosterior` rows
per contrast) via
[`ggdist::stat_halfeye`](https://mjskay.github.io/ggdist/reference/stat_halfeye.html),
because a density slab needs samples rather than summarized intervals.
Half-eye nested intervals follow `result$settings$credibleInterval`
(HPDI via
[`ggdist::mean_hdi`](https://mjskay.github.io/ggdist/reference/point_interval.html),
CPI via
[`ggdist::mean_qi`](https://mjskay.github.io/ggdist/reference/point_interval.html))
unless `...` overrides `point_interval`. When a Bayes result carries a
NOPE, permutation-calibrated or fixed, a gray band behind each half-eye
spans that statistic's interval: \\R\\ from 1 to `R_nopeHigh` in the top
panel, \\E\\ from 0 to `E_nopeHigh` in the bottom panel. You can then
see how much of each posterior clears its edge. The legend names each
band: "NOPE (permutation)" for the default permutation-calibrated bound,
whose band shows the median of the per-draw bounds, or "NOPE (fixed)"
for a bound you supplied. Permutation results have no draws, so they
render as points colored by adjusted significance. When `showOmnibus` is
TRUE, a dashed line in each panel marks the omnibus value; the \\R\\
panel labels it with its p-value (permutation) or `PPGT1` and
`R_PPGTnope` (Bayes), and the \\E\\ panel labels it with \\E\\ and,
under Bayes, `E_PPGTnope`. Omnibus draws stay in `result$omnibus$draws`
rather than as a plot row.
