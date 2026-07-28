# Build a CategoryComposition object from validated components.

Build a CategoryComposition object from validated components.

## Usage

``` r
CategoryComposition(
  counts,
  group,
  categoryLevels = colnames(counts),
  subjectIds = rownames(counts),
  provenance = list(source = "manual")
)
```

## Arguments

- counts:

  Subjects-by-categories matrix of raw counts (non-negative).

- group:

  Named factor of group labels, one per subject row in counts.

- categoryLevels:

  Character vector of category column order.

- subjectIds:

  Character vector of subject identifiers (rownames).

- provenance:

  List describing source type and column mapping.

## Value

A CategoryComposition object.
