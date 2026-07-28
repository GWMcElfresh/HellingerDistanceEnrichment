# Extract subject-by-category counts from metadata

Aggregates cell-level metadata to subject-level category counts.

## Usage

``` r
ExtractClusterComposition(
  object,
  subjectCol = "subjectId",
  categoryCol = "category",
  groupCol = "group",
  countCol = "n",
  ...
)
```

## Arguments

- object:

  A data.frame, Seurat object, or anndata object.

- subjectCol:

  Subject identifier column.

- categoryCol:

  Category label column.

- groupCol:

  Group label column (constant within subject).

- countCol:

  Count column for long-table input.

- ...:

  Unused.

## Value

A `CategoryComposition` object.
