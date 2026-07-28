# Extract subject-by-category counts from cell-level metadata.

Aggregates cell-level categorical labels (clusters, cell types) to
subject-level count matrices. Missing subject-category combinations are
filled with zero.

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

  A data.frame, Seurat object, or anndata AnnData object.

- subjectCol:

  Metadata column identifying the subject (donor, sample).

- categoryCol:

  Metadata column with category labels (cluster, cell type).

- groupCol:

  Metadata column with experimental group per cell; must be constant
  within each subject.

- countCol:

  For long-table input, the raw count column (default "n").

- ...:

  Passed to method-specific helpers.

## Value

A CategoryComposition object.
