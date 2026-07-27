#' Extract subject-by-category counts from cell-level metadata.
#'
#' Aggregates cell-level categorical labels (clusters, cell types) to subject-level
#' count matrices. Missing subject-category combinations are filled with zero.
#'
#' @param object A data.frame, Seurat object, or anndata AnnData object.
#' @param subjectCol Metadata column identifying the subject (donor, sample).
#' @param categoryCol Metadata column with category labels (cluster, cell type).
#' @param groupCol Metadata column with experimental group per cell; must be
#'   constant within each subject.
#' @param countCol For long-table input, the raw count column (default "n").
#' @param ... Passed to method-specific helpers.
#' @return A CategoryComposition object.
#' @export
ExtractClusterComposition <- function(object,
                                      subjectCol = "subjectId",
                                      categoryCol = "category",
                                      groupCol = "group",
                                      countCol = "n",
                                      ...) {
    UseMethod("ExtractClusterComposition")
}

#' @export
ExtractClusterComposition.data.frame <- function(object,
                                                 subjectCol = "subjectId",
                                                 categoryCol = "category",
                                                 groupCol = "group",
                                                 countCol = "n",
                                                 ...) {
    if (all(c(subjectCol, categoryCol, groupCol, countCol) %in% colnames(object))) {
        return(coerce_long_table_to_composition(
            longTable = object,
            subjectCol = subjectCol,
            categoryCol = categoryCol,
            groupCol = groupCol,
            countCol = countCol
        ))
    }

    # Cell-level metadata: one row per cell, aggregate to counts.
    required_cols <- c(subjectCol, categoryCol, groupCol)
    missing_cols <- setdiff(required_cols, colnames(object))
    if (length(missing_cols) > 0) {
        stop(sprintf(
            "metadata is missing required columns: %s",
            paste(missing_cols, collapse = ", ")
        ))
    }

    meta <- object[, required_cols, drop = FALSE]
    colnames(meta) <- c("subjectId", "category", "group")
    meta$subjectId <- as.character(meta$subjectId)
    meta$category <- as.character(meta$category)
    meta$group <- as.character(meta$group)

    subject_group_map <- unique(meta[, c("subjectId", "group")])
    duplicated_subjects <- subject_group_map$subjectId[
        duplicated(subject_group_map$subjectId)
    ]
    if (length(duplicated_subjects) > 0) {
        stop(sprintf(
            "subject(s) map to multiple groups: %s",
            paste(unique(duplicated_subjects), collapse = ", ")
        ))
    }

    category_levels <- sort(unique(meta$category))
    subject_ids <- sort(unique(meta$subjectId))

    counts <- matrix(
        0,
        nrow = length(subject_ids),
        ncol = length(category_levels),
        dimnames = list(subject_ids, category_levels)
    )

    tab <- table(meta$subjectId, meta$category)
    counts[rownames(tab), colnames(tab)] <- as.numeric(tab)

    group <- stats::setNames(
        factor(subject_group_map$group, levels = sort(unique(subject_group_map$group))),
        subject_group_map$subjectId
    )
    group <- group[subject_ids]

    CategoryComposition(
        counts = counts,
        group = group,
        categoryLevels = category_levels,
        subjectIds = subject_ids,
        provenance = list(
            source = "data.frame",
            subjectCol = subjectCol,
            categoryCol = categoryCol,
            groupCol = groupCol
        )
    )
}

#' @export
ExtractClusterComposition.Seurat <- function(object,
                                             subjectCol = "subjectId",
                                             categoryCol = "category",
                                             groupCol = "group",
                                             countCol = "n",
                                             ...) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Seurat must be installed to extract from Seurat objects")
    }

    meta <- object@meta.data
  if (!all(c(subjectCol, categoryCol, groupCol) %in% colnames(meta))) {
    stop(sprintf(
      "Seurat meta.data is missing required columns among: %s, %s, %s",
      subjectCol, categoryCol, groupCol
    ))
  }

    ExtractClusterComposition.data.frame(
        object = meta,
        subjectCol = subjectCol,
        categoryCol = categoryCol,
        groupCol = groupCol,
        countCol = countCol,
        ...
    )
}

#' @export
ExtractClusterComposition.AnnData <- function(object,
                                              subjectCol = "subjectId",
                                              categoryCol = "category",
                                              groupCol = "group",
                                              countCol = "n",
                                              ...) {
    if (!requireNamespace("anndata", quietly = TRUE)) {
        stop("anndata must be installed to extract from AnnData objects")
    }

    obs <- object$obs
    obs_df <- as.data.frame(obs)

    ExtractClusterComposition.data.frame(
        object = obs_df,
        subjectCol = subjectCol,
        categoryCol = categoryCol,
        groupCol = groupCol,
        countCol = countCol,
        ...
    )
}

#' @export
ExtractClusterComposition.AnnDataR6 <- function(object,
                                                subjectCol = "subjectId",
                                                categoryCol = "category",
                                                groupCol = "group",
                                                countCol = "n",
                                                ...) {
    if (!requireNamespace("anndata", quietly = TRUE)) {
        stop("anndata must be installed to extract from AnnData objects")
    }

    obs <- object$obs
    obs_df <- as.data.frame(obs)

    ExtractClusterComposition.data.frame(
        object = obs_df,
        subjectCol = subjectCol,
        categoryCol = categoryCol,
        groupCol = groupCol,
        countCol = countCol,
        ...
    )
}

#' @export
ExtractClusterComposition.default <- function(object, ...) {
    if (inherits(object, "python.builtin.object")) {
        return(ExtractClusterComposition.AnnData(object, ...))
    }
    stop(sprintf(
        "no ExtractClusterComposition method for class '%s'",
        paste(class(object), collapse = ", ")
    ))
}
