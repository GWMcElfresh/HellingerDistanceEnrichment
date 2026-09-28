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

anndata_py_object <- function(object) {
    if (inherits(object, "AnnDataR6")) {
        return(object$.get_py_object())
    }
    if (inherits(object, "python.builtin.object")) {
        return(object)
    }
    stop(sprintf(
        "no Python AnnData object found for class '%s'",
        paste(class(object), collapse = ", ")
    ))
}

python_string_list_to_character <- function(values) {
    if (is.character(values)) {
        return(values)
    }
    if (is.null(values)) {
        return(character())
    }
    vapply(values, function(value) {
        if (is.null(value) || length(value) != 1L || is.na(value)) {
            NA_character_
        } else {
            as.character(value)
        }
    }, character(1))
}

# Pull obs column names and string values in Python. Do not use object$obs:
# reticulate::py_convert_pandas_df fails on the numpy 2 / pandas dtypes in
# discvr-base ("Not compatible with STRSXP: [type=environment]").
anndata_obs_payload <- function(object, columns) {
    py_ad <- anndata_py_object(object)
    py_env <- reticulate::py_run_string(
        "
def _hde_obs_columns(adata, columns):
    obs = adata.obs
    colnames = [str(name) for name in list(obs.columns)]
    wanted = set(columns)
    values = {}
    for name in colnames:
        if name not in wanted:
            continue
        series = obs[name]
        missing_mask = series.isna().tolist()
        labels = series.astype(str).tolist()
        values[name] = [
            None if is_missing else label
            for is_missing, label in zip(missing_mask, labels)
        ]
    return {'columns': colnames, 'values': values}
",
        local = TRUE,
        convert = FALSE
    )
    raw <- py_env$`_hde_obs_columns`(py_ad, as.list(columns))
    if (inherits(raw, "python.builtin.object")) {
        reticulate::py_to_r(raw)
    } else {
        raw
    }
}

obs_frame_from_payload <- function(payload, columns) {
    obs_names <- as.character(payload$columns)
    missing_cols <- setdiff(columns, obs_names)
    if (length(missing_cols) > 0L) {
        stop(sprintf(
            "metadata is missing required columns: %s",
            paste(missing_cols, collapse = ", ")
        ))
    }

    values <- payload$values
    obs_df <- as.data.frame(
        lapply(columns, function(name) {
            python_string_list_to_character(values[[name]])
        }),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    names(obs_df) <- columns
    obs_df
}

extract_composition_from_anndata <- function(object,
                                             subject_col,
                                             category_col,
                                             group_col,
                                             count_col,
                                             ...) {
    required_cols <- c(subject_col, category_col, group_col)
    # Ask for the count column too, but only keep it when obs actually has it.
    # A cell-level obs table has no count column and must not be reported as
    # missing one.
    probe_cols <- unique(c(required_cols, count_col))
    payload <- anndata_obs_payload(object, probe_cols)
    obs_names <- as.character(payload$columns)
    missing_cols <- setdiff(required_cols, obs_names)
    if (length(missing_cols) > 0L) {
        stop(sprintf(
            "metadata is missing required columns: %s",
            paste(missing_cols, collapse = ", ")
        ))
    }

    fetch_cols <- required_cols
    if (count_col %in% obs_names) {
        fetch_cols <- unique(c(fetch_cols, count_col))
    }
    obs_df <- obs_frame_from_payload(payload, fetch_cols)
    ExtractClusterComposition.data.frame(
        object = obs_df,
        subjectCol = subject_col,
        categoryCol = category_col,
        groupCol = group_col,
        countCol = count_col,
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

    extract_composition_from_anndata(
        object = object,
        subject_col = subjectCol,
        category_col = categoryCol,
        group_col = groupCol,
        count_col = countCol,
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

    extract_composition_from_anndata(
        object = object,
        subject_col = subjectCol,
        category_col = categoryCol,
        group_col = groupCol,
        count_col = countCol,
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
