# Energy distance between groups from a group summary matrix.

For two groups, \\E\\ is twice the mean between-group Hellinger distance
minus each group's mean within-group distance, on the same distances the
ratio uses. With more groups, the function averages \\E\\ over every
pair of groups, the same pairing the ratio uses. Each within mean then
enters \\G - 1\\ pairs, so the average reduces to
`2 * (mean between - mean within)`.

## Usage

``` r
compute_energy_distance(group_summary_matrix)
```

## Arguments

- group_summary_matrix:

  Square matrix from build_group_distance_summary_matrix on squared
  Hellinger distances.

## Value

Scalar energy distance in squared Hellinger units.

## Details

Callers pass a summary built from squared Hellinger distances. Then
\\E\\ equals the squared distance between group centroids in square-root
composition space, minus the sampling spread each centroid carries, so
it measures a shift in mean composition with no division by within-group
spread. Independent count noise adds equally to between- and
within-group squared distances, so it cancels from \\E\\ in expectation;
Jeffreys softening and the curvature of the square root still shrink
\\E\\ at shallow depth. The sample estimate excludes each subject's
distance to itself, so it can fall slightly below 0 when groups do not
differ.
