# Group-level mean pairwise distance summary matrix.

Group-level mean pairwise distance summary matrix.

## Usage

``` r
build_group_distance_summary_matrix(group, dist_matrix)
```

## Arguments

- group:

  Named factor aligned to dist_matrix rownames.

- dist_matrix:

  Symmetric subject distance matrix.

## Value

Square matrix of mean distances (diagonal = within, off-diagonal =
between).
