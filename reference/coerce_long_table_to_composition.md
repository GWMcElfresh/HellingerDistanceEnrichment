# Coerce a long counts table to CategoryComposition.

Coerce a long counts table to CategoryComposition.

## Usage

``` r
coerce_long_table_to_composition(
  longTable,
  subjectCol = "subjectId",
  categoryCol = "category",
  groupCol = "group",
  countCol = "n"
)
```

## Arguments

- longTable:

  data.frame with subject, category, group, and count columns.

- subjectCol:

  Column name for subject identifiers.

- categoryCol:

  Column name for category labels.

- groupCol:

  Column name for group labels (one per subject).

- countCol:

  Column name for raw counts (default "n").

## Value

A CategoryComposition object.
