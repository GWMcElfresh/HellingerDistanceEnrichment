# Upper NOPE edge for one statistic in one posterior draw.

Upper NOPE edge for one statistic in one posterior draw.

## Usage

``` r
resolve_nope_edge(
  nope_setting,
  shuffled_statistics,
  statistic_column,
  floor_value,
  nope_quantile
)
```

## Arguments

- nope_setting:

  `"shuffle"`, `"none"`, or a fixed numeric edge.

- shuffled_statistics:

  Matrix from compute_shuffled_statistics, or `NULL` when no shuffles
  ran.

- statistic_column:

  `"ratio"` or `"energy"`.

- floor_value:

  The statistic's no-difference value (1 for \\R\\, 0 for \\E\\); the
  NOPE edge never falls below it.

- nope_quantile:

  Quantile of the shuffled statistic used as the edge.

## Value

Scalar edge, or `NA` when `nope_setting` is `"none"`.
