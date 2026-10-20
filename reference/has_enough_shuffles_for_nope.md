# Check whether a contrast has enough distinct label shuffles for a NOPE.

Each NOPE edge is the `nope_quantile` of its statistic (\\R\\ or \\E\\)
across label shuffles, so it only works as a bar when the true grouping
is a small fraction of all shuffles. With few subjects there are few
distinct ways to split them; the true split then sets the quantile
itself and neither statistic can clear it, the same limit that keeps a 3
versus 3 permutation p-value at or above 0.1.

## Usage

``` r
has_enough_shuffles_for_nope(contrast_group, nope_quantile = 0.95)
```

## Arguments

- contrast_group:

  Factor of group labels for one contrast.

- nope_quantile:

  Quantile of each shuffled statistic used as its NOPE edge.

## Value

`TRUE` when \\R\\ and \\E\\ can clear their NOPEs, `FALSE` otherwise.

## Details

\\R\\ and \\E\\ depend only on which subjects share a group, so
relabelings that swap two equal-sized groups count as one split. The
number of distinct splits is `n! / (prod(n_g!) * prod(m_k!))`, where
`n_g` are group sizes and `m_k` counts groups of each size (10 for 3
versus 3, 35 for 4 versus 4). The contrast passes when the true split is
under `1 - nope_quantile` of all shuffles, i.e. more than 20 distinct
splits at the default 0.95. Depends only on group sizes, so it uses no
random numbers and runs once per contrast.
