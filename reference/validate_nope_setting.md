# Validate a `nope` or `energyNope` argument.

Accepts `"shuffle"` (a NOPE interval), `"none"`, or one finite number
above the statistic's no-difference value: above 1 for the ratio \\R\\,
above 0 for the energy distance \\E\\. An edge at or below that value
would let draws consistent with no difference count as clearing the
interval, so the function rejects it.

## Usage

``` r
validate_nope_setting(nope_setting, scale = c("ratio", "energy"))
```

## Arguments

- nope_setting:

  Value passed to CompareGroupCompositions.

- scale:

  `"ratio"` for `nope`, `"energy"` for `energyNope`.

## Value

`"shuffle"`, `"none"`, or the numeric edge.
