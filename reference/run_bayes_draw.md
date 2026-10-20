# One Bayes draw: sample compositions, rebuild \\R\\ and \\E\\, and set this draw's NOPE upper edges.

When `nope` or `energy_nope` is `"shuffle"`, the function shuffles
labels on this draw's own distance matrix, so each edge carries the same
count noise as this draw's statistic rather than the noise of one fixed
point-estimate matrix. One set of shuffles serves both statistics. Each
NOPE edge is the `nope_quantile` of its shuffled statistic, and the
function floors it at the no-difference value (1 for \\R\\, 0 for
\\E\\): a quantile below that value would invert the interval and let a
draw below no difference count as clearing it. A numeric setting is a
fixed edge and needs no shuffles.

## Usage

``` r
run_bayes_draw(
  counts,
  group,
  spec,
  prior_pseudocounts,
  nope = "shuffle",
  energy_nope = "shuffle",
  n_nope_shuffles = 200,
  nope_quantile = 0.95
)
```
