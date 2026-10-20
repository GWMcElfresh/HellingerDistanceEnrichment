# Print one omnibus NOPE line for the ratio or the energy distance.

Print one omnibus NOPE line for the ratio or the energy distance.

## Usage

``` r
print_nope_line(
  nope_probability,
  nope_high,
  nope_setting,
  statistic_label,
  no_difference_value,
  omnibus
)
```

## Arguments

- nope_probability:

  `R_PPGTnope` or `E_PPGTnope`, or `NULL` when the result carries no
  NOPE or fixed edge on that statistic.

- nope_high:

  Matching `R_nopeHigh` or `E_nopeHigh`.

- nope_setting:

  `x$settings$nope` or `x$settings$energyNope`.

- statistic_label:

  `"R"` or `"E"`.

- no_difference_value:

  Lower edge of the interval (1 for R, 0 for E).

- omnibus:

  The result's omnibus list, for `enoughShufflesForNope` and
  `R_hdiAboveNope` or `E_hdiAboveNope`.

## Value

`NULL`, invisibly; called for its printed line.
