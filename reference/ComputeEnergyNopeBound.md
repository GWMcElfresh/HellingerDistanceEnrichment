# Energy distance NOPE bound for a composition shift.

Translates the smallest shift in mean composition that you would call
meaningful into a fixed NOPE bound on the energy distance \\E\\, ready
to pass as `energyNope` to CompareGroupCompositions. You give the
reference composition (e.g. control) and the shifted composition you
have in mind; the function normalizes both to frequencies and returns
\\2(1 - \mathrm{BC})\\, where \\\mathrm{BC} = \sum_k \sqrt{p_k q_k}\\ is
the Bhattacharyya coefficient between them, in squared Hellinger units.

## Usage

``` r
ComputeEnergyNopeBound(referenceComposition, shiftedComposition)
```

## Arguments

- referenceComposition:

  Non-negative numeric vector of category frequencies or counts for the
  reference group. The function normalizes it to sum to 1.

- shiftedComposition:

  Non-negative numeric vector for the shifted group, over the same
  categories. When both vectors are named, the function matches
  categories by name; otherwise it matches them by position.

## Value

Scalar energy distance in squared Hellinger units, above 0.

## Details

The bound equals \\E\\ exactly when every subject sits at its group's
mean composition. With within-group spread, \\E\\ is the squared
distance between group centroids in square-root composition space, so
treat the bound as a close guide rather than an exact translation. No
such helper exists for the ratio \\R\\, because \\R\\ also divides by
within-group spread; a fixed `nope` bound on \\R\\ needs a simulation at
your own depths and group sizes.

## See also

[`CompareGroupCompositions()`](https://gwmcelfresh.github.io/HellingerDistanceEnrichment/reference/CompareGroupCompositions.md),
whose `energyNope` argument takes the bound.

## Examples

``` r
reference <- c(TypeA = 0.50, TypeB = 0.30, TypeC = 0.20)
shifted <- c(TypeA = 0.45, TypeB = 0.35, TypeC = 0.20)
ComputeEnergyNopeBound(reference, shifted)
#> [1] 0.003242632

# Counts work as well; the function normalizes them.
ComputeEnergyNopeBound(c(500, 300, 200), c(450, 350, 200))
#> [1] 0.003242632
```
