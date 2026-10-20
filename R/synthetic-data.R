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

#' Simulate a two-group cohort with controlled depth, heterogeneity, and shift.
#'
#' Each subject draws a latent composition from a Dirichlet centered on its
#' group's mean composition, then draws `cells_per_subject` counts from a
#' multinomial on that composition. Three knobs separate the sources of
#' variation the sampling vignette studies: `cells_per_subject` sets count noise
#' within a subject, `concentration` sets how far subjects scatter around their
#' group mean (larger is more homogeneous), and `shift` moves composition mass
#' between two categories in the second group only.
#'
#' The baseline mean composition decays linearly over categories (e.g. 10
#' categories run from about 0.18 down to 0.02), so the shift moves `shift`
#' mass from the first, most abundant category to the last, rarest one.
#'
#' @param n_subjects_per_group Subjects per group; length 1 or 2.
#' @param n_categories Number of categories (>= 2).
#' @param cells_per_subject Multinomial count total per subject, typically cells
#'   (or events).
#' @param concentration Dirichlet concentration per group; length 1 or 2.
#' @param shift Composition mass moved in the second group, in [0, first
#'   category's baseline share).
#' @param groups Two group labels.
#' @param seed Random seed, or `NULL` to use the current stream.
#' @return List with `long_table` (subjectId, category, group, n),
#'   `latent_counts` (subjects by categories latent compositions scaled to
#'   `cells_per_subject`, for noise-free reference statistics), and `group`
#'   (named factor).
#' @keywords internal
simulate_composition_cohort <- function(n_subjects_per_group = 8,
                                        n_categories = 10,
                                        cells_per_subject = 2000,
                                        concentration = 50,
                                        shift = 0,
                                        groups = c("Control", "Treatment"),
                                        seed = NULL) {
    if (length(groups) != 2 || anyDuplicated(groups)) {
        stop("groups must be two distinct labels")
    }
    if (n_categories < 2) {
        stop("n_categories must be at least 2")
    }
    n_subjects_per_group <- rep_len(n_subjects_per_group, 2)
    concentration <- rep_len(concentration, 2)
    if (any(n_subjects_per_group < 2) || any(concentration <= 0) || cells_per_subject < 1) {
        stop("need at least 2 subjects per group, positive concentration, and at least 1 cell")
    }

    baseline_weights <- rev(seq_len(n_categories))
    baseline_mean <- baseline_weights / sum(baseline_weights)
    if (shift < 0 || shift >= baseline_mean[1]) {
        stop(sprintf("shift must be in [0, %.3f) for %d categories", baseline_mean[1], n_categories))
    }
    shifted_mean <- baseline_mean
    shifted_mean[1] <- shifted_mean[1] - shift
    shifted_mean[n_categories] <- shifted_mean[n_categories] + shift
    group_means <- list(baseline_mean, shifted_mean)

    if (!is.null(seed)) {
        set.seed(seed)
    }

    category_labels <- paste0("Category", seq_len(n_categories) - 1)
    subject_ids <- character()
    subject_groups <- character()
    latent_rows <- list()
    count_rows <- list()

    for (group_idx in 1:2) {
        alpha <- concentration[group_idx] * group_means[[group_idx]]
        for (subject_idx in seq_len(n_subjects_per_group[group_idx])) {
            gamma_draws <- stats::rgamma(n_categories, shape = alpha, rate = 1)
            latent_composition <- gamma_draws / sum(gamma_draws)
            subject_counts <- as.vector(stats::rmultinom(1, cells_per_subject, latent_composition))

            subject_ids <- c(subject_ids, sprintf("%s_S%02d", groups[group_idx], subject_idx))
            subject_groups <- c(subject_groups, groups[group_idx])
            latent_rows[[length(latent_rows) + 1]] <- latent_composition * cells_per_subject
            count_rows[[length(count_rows) + 1]] <- subject_counts
        }
    }

    latent_counts <- do.call(rbind, latent_rows)
    count_matrix <- do.call(rbind, count_rows)
    dimnames(latent_counts) <- list(subject_ids, category_labels)
    dimnames(count_matrix) <- list(subject_ids, category_labels)

    long_table <- data.frame(
        subjectId = rep(subject_ids, times = n_categories),
        category = rep(category_labels, each = length(subject_ids)),
        group = rep(subject_groups, times = n_categories),
        n = as.vector(count_matrix),
        stringsAsFactors = FALSE
    )
    group <- factor(subject_groups, levels = groups)
    names(group) <- subject_ids

    list(long_table = long_table, latent_counts = latent_counts, group = group)
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
