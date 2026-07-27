#' Synthetic long-table builders for tests and vignettes.
#'
#' @param n_subjects_per_group Number of subjects in each group.
#' @param n_categories Number of categorical levels.
#' @param groups Character vector of group names.
#' @param seed Random seed.
#' @return data.frame with subjectId, category, group, n columns.
#' @keywords internal
build_synthetic_long_table <- function(n_subjects_per_group = 5,
                                       n_categories = 4,
                                       groups = c("Control", "Treatment"),
                                       seed = 1) {
    set.seed(seed)
    rows <- list()

    for (group_label in groups) {
        for (subject_idx in seq_len(n_subjects_per_group)) {
            subject_id <- sprintf("%s_S%02d", group_label, subject_idx)
            category_counts <- stats::rmultinom(1, size = 100, prob = rep(1, n_categories))
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

#' Plant group-specific composition structure into a long table.
#'
#' @param long_table Output from build_synthetic_long_table.
#' @param target_group Group that receives elevated counts in target_category.
#' @param target_category Category label to enrich.
#' @param boost Multiplier for target category counts.
#' @return Modified long table.
#' @keywords internal
plant_group_structure <- function(long_table,
                                  target_group = "Treatment",
                                  target_category = "Category0",
                                  boost = 5) {
    is_target <- long_table$group == target_group & long_table$category == target_category
    long_table$n[is_target] <- long_table$n[is_target] * boost
    long_table
}
