#' Build a CategoryComposition object from validated components.
#'
#' @param counts Subjects-by-categories matrix of raw counts (non-negative).
#' @param group Named factor of group labels, one per subject row in counts.
#' @param categoryLevels Character vector of category column order.
#' @param subjectIds Character vector of subject identifiers (rownames).
#' @param provenance List describing source type and column mapping.
#' @return A CategoryComposition object.
#' @export
CategoryComposition <- function(counts,
                                  group,
                                  categoryLevels = colnames(counts),
                                  subjectIds = rownames(counts),
                                  provenance = list(source = "manual")) {
    counts <- as.matrix(counts)
    storage.mode(counts) <- "numeric"

    if (is.null(rownames(counts)) && !is.null(subjectIds)) {
        rownames(counts) <- subjectIds
    }
    if (is.null(colnames(counts)) && !is.null(categoryLevels)) {
        colnames(counts) <- categoryLevels
    }

    subject_ids <- rownames(counts)
    if (is.null(subject_ids) || any(subject_ids == "")) {
        stop("counts must have named subject rows (subjectIds)")
    }

    group <- stats::setNames(as.factor(group), names(group))
    if (!all(subject_ids %in% names(group))) {
        stop("group must be named for every subject in counts")
    }
    group <- group[subject_ids]

    if (anyDuplicated(subject_ids) > 0) {
        stop("duplicate subjectIds are not allowed")
    }
    if (any(counts < 0, na.rm = TRUE)) {
        stop("counts must be non-negative")
    }
    if (any(is.na(counts))) {
        stop("counts must not contain NA; missing subject-category cells are 0")
    }

    structure(
        list(
            counts = counts,
            group = group,
            categoryLevels = colnames(counts),
            subjectIds = subject_ids,
            provenance = provenance
        ),
        class = "CategoryComposition"
    )
}

#' @export
print.CategoryComposition <- function(x, ...) {
    cat(
        "CategoryComposition with",
        nrow(x$counts), "subjects,",
        ncol(x$counts), "categories, and",
        length(levels(x$group)), "groups\n"
    )
    invisible(x)
}

#' Coerce a long counts table to CategoryComposition.
#'
#' @param longTable data.frame with subject, category, group, and count columns.
#' @param subjectCol Column name for subject identifiers.
#' @param categoryCol Column name for category labels.
#' @param groupCol Column name for group labels (one per subject).
#' @param countCol Column name for raw counts (default "n").
#' @return A CategoryComposition object.
#' @keywords internal
coerce_long_table_to_composition <- function(longTable,
                                             subjectCol = "subjectId",
                                             categoryCol = "category",
                                             groupCol = "group",
                                             countCol = "n") {
    required_cols <- c(subjectCol, categoryCol, groupCol, countCol)
    missing_cols <- setdiff(required_cols, colnames(longTable))
    if (length(missing_cols) > 0) {
        stop(sprintf(
            "long table is missing required columns: %s",
            paste(missing_cols, collapse = ", ")
        ))
    }

    long_table <- longTable[, required_cols, drop = FALSE]
    colnames(long_table) <- c("subjectId", "category", "group", "n")
    long_table$subjectId <- as.character(long_table$subjectId)
    long_table$category <- as.character(long_table$category)
    long_table$group <- as.character(long_table$group)
    long_table$n <- as.numeric(long_table$n)

    if (any(is.na(long_table$n))) {
        stop("count column must not contain NA")
    }
    if (any(long_table$n < 0)) {
        stop("counts must be non-negative")
    }

    # One group per subject; aggregation would hide conflicting labels.
    subject_group_map <- unique(long_table[, c("subjectId", "group")])
    duplicated_subjects <- subject_group_map$subjectId[
        duplicated(subject_group_map$subjectId)
    ]
    if (length(duplicated_subjects) > 0) {
        stop(sprintf(
            "subject(s) map to multiple groups: %s",
            paste(unique(duplicated_subjects), collapse = ", ")
        ))
    }

    category_levels <- sort(unique(long_table$category))
    subject_ids <- sort(unique(long_table$subjectId))

    counts <- matrix(
        0,
        nrow = length(subject_ids),
        ncol = length(category_levels),
        dimnames = list(subject_ids, category_levels)
    )

    for (row_idx in seq_len(nrow(long_table))) {
        row <- long_table[row_idx, ]
        counts[row$subjectId, row$category] <- counts[row$subjectId, row$category] + row$n
    }

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
            source = "long_table",
            subjectCol = subjectCol,
            categoryCol = categoryCol,
            groupCol = groupCol,
            countCol = countCol
        )
    )
}
