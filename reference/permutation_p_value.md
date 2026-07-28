# Permutation p-value for one contrast (includes observed in numerator).

Matches Edlefsen's workbook form:
`mean(c(permuted_ratios, observed_ratio) >= observed_ratio)`.

## Usage

``` r
permutation_p_value(observed_ratio, permuted_ratios)
```
