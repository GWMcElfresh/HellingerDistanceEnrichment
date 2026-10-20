# 95% posterior interval on finite draws: shortest HPDI or equal-tailed CPI.

HPDI uses the Chen-Shao shortest window covering `ceil(prob * n)` sorted
samples. Quantile intervals spend equal tail mass, which is a poor
summary when the posterior of \\R\\ is asymmetric and sits just above 1;
the shortest covering window is the HPD for a continuous unimodal
sample.

## Usage

``` r
compute_posterior_interval(x, prob = 0.95, type = c("hpdi", "cpi"))
```
