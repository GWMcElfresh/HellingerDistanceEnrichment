# Build the shaded NOPE band for one statistic in each Bayes contrast row.

Each band runs from the statistic's no-difference value (1 for \\R\\, 0
for \\E\\) to the contrast's `R_nopeHigh` or `E_nopeHigh` and covers the
height of that row's half-eye slab, which
[`ggdist::stat_halfeye`](https://mjskay.github.io/ggdist/reference/stat_halfeye.html)
draws upward from the row by `slab_scale` (0.9 by default). For the
NOPE, the edge is the median of the per-draw edges. Contrasts whose NOPE
probability is `NA` (groups too small for the NOPE) get no band, because
the edge there carries no information.

## Usage

``` r
build_nope_band(
  plot_data,
  contrast_levels,
  statistic_prefix,
  nope_setting,
  slab_scale = NULL
)
```

## Arguments

- plot_data:

  Contrast summary rows with `contrastId` as a factor.

- contrast_levels:

  Factor levels in plotting order.

- statistic_prefix:

  `"R"` or `"E"`, matching the column prefixes.

- nope_setting:

  `result$settings$nope` or `result$settings$energyNope`.

- slab_scale:

  Half-eye `scale` argument, or `NULL` for the default.

## Value

data.frame of rectangles with a `statistic` column, or `NULL` when the
result has no NOPE or fixed edge on that statistic.
