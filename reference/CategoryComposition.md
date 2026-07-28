# Subject-level categorical composition container

Constructs a validated composition object used by comparison and
plotting functions.

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

  Subjects-by-categories matrix of raw counts.

- group:

  Named factor of group labels.

- categoryLevels:

  Category column order.

- subjectIds:

  Subject identifiers.

- provenance:

  Source metadata list.

## Value

A `CategoryComposition` S3 object.
