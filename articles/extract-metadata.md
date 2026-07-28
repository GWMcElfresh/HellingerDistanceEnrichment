# Extracting compositions from Seurat and anndata metadata

Single-cell workflows store the annotations needed for composition
analysis at the cell level, while Hellinger enrichment operates on
subject-level category counts. Bridging those scales depends on metadata
columns that map each cell to a subject, a category (cluster or cell
type), and an experimental group. This vignette builds tiny synthetic
Seurat and anndata objects in memory and shows that both routes produce
the same `CategoryComposition` input used by
[`CompareGroupCompositions()`](https://gwmcelfresh.github.io/HellingerDistanceEnrichment/reference/CompareGroupCompositions.md).

## Seurat metadata path

For a Seurat object,
[`ExtractClusterComposition()`](https://gwmcelfresh.github.io/HellingerDistanceEnrichment/reference/ExtractClusterComposition.md)
reads the requested columns from `meta.data`, aggregates cells to
subjects, and returns the composition object. Expression counts are
irrelevant to the enrichment step; only the metadata mapping matters.
The example below creates a minimal object so the extraction path can be
run without a study-sized dataset.

``` r

library(HellingerDistanceEnrichment)

if (requireNamespace("Seurat", quietly = TRUE)) {
  n_cells <- 120
  meta <- data.frame(
    subjectId = rep(paste0("Donor", 1:6), each = 20),
    category = sample(paste0("Cluster", 0:3), n_cells, replace = TRUE),
    group = rep(c("Control", "Control", "Control", "Case", "Case", "Case"), each = 20),
    stringsAsFactors = FALSE
  )

  counts <- matrix(
    stats::rpois(100 * n_cells, lambda = 2),
    nrow = 100,
    ncol = n_cells
  )
  colnames(counts) <- paste0("cell", seq_len(n_cells))
  rownames(counts) <- paste0("gene", seq_len(100))

  seu <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta)

  composition <- ExtractClusterComposition(
    seu,
    subjectCol = "subjectId",
    categoryCol = "category",
    groupCol = "group"
  )

  print(composition)

  result <- CompareGroupCompositions(
    composition,
    method = "permutation",
    nPermutations = 50,
    seed = 1
  )
  result$omnibus
} else {
  message("Seurat not installed; skipping Seurat example.")
}
#> Warning: Data is of class matrix. Coercing to dgCMatrix.
#> CategoryComposition with 6 subjects, 4 categories, and 2 groups
#> $effectSize
#> [1] 0.9485571
#> 
#> $pValue
#> [1] 0.4705882
```

## anndata obs path

The anndata route is mechanically the same: cell-level fields live in
`obs`, and the extractor aggregates them to subjects. Column names may
differ across pipelines, so the `subjectCol`, `categoryCol`, and
`groupCol` arguments must match the object at hand. Once those columns
are identified, the resulting composition object is interchangeable with
the Seurat-derived one.

``` r

if (requireNamespace("anndata", quietly = TRUE)) {
  obs <- data.frame(
    subjectId = rep(paste0("S", 1:4), each = 15),
    category = sample(paste0("Type", 0:2), 60, replace = TRUE),
    group = rep(c("A", "A", "B", "B"), each = 15),
    stringsAsFactors = FALSE
  )
  X <- matrix(stats::rnorm(60 * 10), nrow = 60)
  adata <- anndata::AnnData(X = X, obs = obs)

  composition <- ExtractClusterComposition(
    adata,
    subjectCol = "subjectId",
    categoryCol = "category",
    groupCol = "group"
  )

  CompareGroupCompositions(
    composition,
    nPermutations = 30,
    seed = 2
  )$omnibus
} else {
  message("anndata not installed; skipping anndata example.")
}
#> $effectSize
#> [1] 0.997142
#> 
#> $pValue
#> [1] 1
```

The practical consequence is that object format need not dictate the
enrichment workflow. Seurat and anndata are alternative containers for
the same subject-category-group mapping. For HPC deployments, bind-mount
study metadata at runtime rather than baking cohort files into container
images; the extraction step then remains portable across environments
while the cohort annotations stay outside the image layer.
