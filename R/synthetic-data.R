#' Build a long counts table with exchangeable multinomial compositions per group.
#'
#' Each subject receives an independent symmetric multinomial draw over categories
#' within its group. Before planting, all groups share the same generative process,
#' so group labels are exchangeable under the null.
#'
#' @param n_subjects_per_group Number of subjects per group level (>= 1).
#' @param n_categories Number of category levels (Category0, Category1, ...).
#' @param groups Character vector of distinct group labels (>= 2 levels).
#' @param seed Random seed for reproducibility.
#' @return data.frame with columns subjectId, category, group, n.
#' @keywords internal
build_synthetic_long_table <- function(n_subjects_per_group = 5,
                                       n_categories = 4,
                                       groups = c("Control", "Treatment"),
                                       seed = 1) {
    if (n_subjects_per_group < 1) {
        stop("n_subjects_per_group must be at least 1")
    }
    if (n_categories < 1) {
        stop("n_categories must be at least 1")
    }
    if (length(groups) < 2) {
        stop("groups must contain at least two levels")
    }
    if (anyDuplicated(groups)) {
        stop("groups must be unique")
    }

    set.seed(seed)
    rows <- list()
    # Total count per subject sets the multinomial precision; symmetric probs keep
    # groups exchangeable until planting targets one arm.
    multinomial_size <- 100

    for (group_label in groups) {
        for (subject_idx in seq_len(n_subjects_per_group)) {
            subject_id <- sprintf("%s_S%02d", group_label, subject_idx)
            category_counts <- stats::rmultinom(
                1,
                size = multinomial_size,
                prob = rep(1, n_categories)
            )
            for (cat_idx in seq_len(n_categories)) {
                rows[[length(rows) + 1]] <- data.frame(
                    subjectId = subject_id,
                    category = paste0("Category", cat_idx - 1),
                    group = group_label,
                    n = category_counts[cat_idx],
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    do.call(rbind, rows)
}

#' Plant a composition shift in one group-category cell.
#'
#' Multiplies counts at (target_group, target_category) by boost. Other group
#' levels are unchanged, so a multi-group design can plant enrichment in one arm
#' (for example Condition1) while leaving another arm exchangeable with Control
#' (for example Condition2).
#'
#' @param long_table Output from build_synthetic_long_table.
#' @param target_group Group label that receives elevated counts.
#' @param target_category Category label to enrich.
#' @param boost Positive multiplier applied only to matching rows.
#' @return Modified long table with the planted shift.
#' @keywords internal
plant_group_structure <- function(long_table,
                                  target_group = "Treatment",
                                  target_category = "Category0",
                                  boost = 5) {
    if (boost <= 0) {
        stop("boost must be positive")
    }
    if (!target_group %in% long_table$group) {
        stop(sprintf("target_group '%s' not found in long_table", target_group))
    }
    if (!target_category %in% long_table$category) {
        stop(sprintf("target_category '%s' not found in long_table", target_category))
    }

    is_target <- long_table$group == target_group & long_table$category == target_category
    long_table$n[is_target] <- long_table$n[is_target] * boost
    long_table
}
