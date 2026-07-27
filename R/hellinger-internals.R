#' Hellinger distance between two frequency vectors on a simplex.
#'
#' @param freqs1 Numeric frequency vector (sums to 1).
#' @param freqs2 Numeric frequency vector (sums to 1).
#' @return Scalar Hellinger distance in [0, 1].
#' @keywords internal
compute_hellinger_distance_freqs <- function(freqs1, freqs2) {
    diff_sqrt <- sqrt(freqs1) - sqrt(freqs2)
    sqrt(sum(diff_sqrt^2)) / sqrt(2)
}

#' Hellinger distance between two count vectors (frequencies derived internally).
#'
#' @param counts1 Numeric count vector.
#' @param counts2 Numeric count vector.
#' @return Scalar Hellinger distance.
#' @keywords internal
compute_hellinger_distance_counts <- function(counts1, counts2) {
    compute_hellinger_distance_freqs(
        counts1 / sum(counts1),
        counts2 / sum(counts2)
    )
}

#' Pairwise Hellinger distance matrix from a count matrix (rows = subjects).
#'
#' @param count_matrix Numeric matrix, subjects by categories.
#' @return Symmetric distance matrix with dimnames from count_matrix.
#' @keywords internal
build_hellinger_distance_matrix <- function(count_matrix) {
    n_subjects <- nrow(count_matrix)
    dist_matrix <- matrix(
        0,
        nrow = n_subjects,
        ncol = n_subjects,
        dimnames = list(rownames(count_matrix), rownames(count_matrix))
    )

    for (i in seq_len(n_subjects)) {
        for (j in i:n_subjects) {
            distance <- compute_hellinger_distance_counts(
                count_matrix[i, ],
                count_matrix[j, ]
            )
            dist_matrix[i, j] <- distance
            dist_matrix[j, i] <- distance
        }
    }

    dist_matrix
}

#' Mean pairwise distance within a set of subjects.
#'
#' @param subject_ids Character vector of subject rownames.
#' @param dist_matrix Symmetric distance matrix.
#' @return Scalar mean of lower-triangle distances.
#' @keywords internal
mean_pairwise_distance <- function(subject_ids, dist_matrix) {
    sub_matrix <- dist_matrix[subject_ids, subject_ids, drop = FALSE]
    mean(stats::as.dist(sub_matrix))
}

#' Mean pairwise distance between two disjoint subject sets.
#'
#' @param subject_ids1 Character vector of subjects in set 1.
#' @param subject_ids2 Character vector of subjects in set 2.
#' @param dist_matrix Symmetric distance matrix.
#' @return Scalar mean cross-set distance.
#' @keywords internal
mean_pairwise_distance_between <- function(subject_ids1, subject_ids2, dist_matrix) {
    if (length(intersect(subject_ids1, subject_ids2)) > 0) {
        stop("between-group distance requires disjoint subject sets")
    }
    sub_matrix <- dist_matrix[subject_ids1, subject_ids2, drop = FALSE]
    mean(sub_matrix)
}

#' Group-level mean pairwise distance summary matrix.
#'
#' @param group Named factor aligned to dist_matrix rownames.
#' @param dist_matrix Symmetric subject distance matrix.
#' @return Square matrix of mean distances (diagonal = within, off-diagonal = between).
#' @keywords internal
build_group_distance_summary_matrix <- function(group, dist_matrix) {
    group_levels <- levels(group)
    n_levels <- length(group_levels)
    summary_matrix <- matrix(
        NA_real_,
        nrow = n_levels,
        ncol = n_levels,
        dimnames = list(group_levels, group_levels)
    )

    subject_ids <- names(group)
    for (level_i in seq_len(n_levels)) {
        subjects_i <- subject_ids[group == group_levels[level_i]]
        summary_matrix[level_i, level_i] <- mean_pairwise_distance(subjects_i, dist_matrix)

        if (level_i < n_levels) {
            for (level_j in seq(level_i + 1, n_levels)) {
                subjects_j <- subject_ids[group == group_levels[level_j]]
                between_distance <- mean_pairwise_distance_between(
                    subjects_i,
                    subjects_j,
                    dist_matrix
                )
                summary_matrix[level_i, level_j] <- between_distance
                summary_matrix[level_j, level_i] <- between_distance
            }
        }
    }

    summary_matrix
}

#' Between-to-within Hellinger distance ratio from a group summary matrix.
#'
#' @param group_summary_matrix Square matrix from build_group_distance_summary_matrix.
#' @return Scalar ratio (mean between / mean within).
#' @keywords internal
compute_between_within_ratio <- function(group_summary_matrix) {
    within_mean <- mean(diag(group_summary_matrix), na.rm = FALSE)
    between_values <- group_summary_matrix[upper.tri(group_summary_matrix)]
    between_mean <- mean(between_values, na.rm = FALSE)

    if (is.na(within_mean) || is.na(between_mean)) {
        stop("cannot compute ratio: insufficient subjects for pairwise distances in a contrast")
    }

    if (within_mean == 0) {
        if (between_mean == 0) {
            return(1)
        }
        return(Inf)
    }

    between_mean / within_mean
}

#' Apply Jeffreys pseudocounts to a count matrix.
#'
#' @param counts Subjects-by-categories count matrix.
#' @param prior_pseudocounts Scalar or per-category pseudocount (default 1/2).
#' @return Softened count matrix.
#' @keywords internal
apply_prior_pseudocounts <- function(counts, prior_pseudocounts = 1/2) {
    counts + prior_pseudocounts
}

#' Sample one Dirichlet draw per subject from conjugate posteriors.
#'
#' @param counts Subjects-by-categories count matrix.
#' @param prior_pseudocounts Jeffreys or user-specified Dirichlet alpha increment.
#' @return Matrix of frequency draws (subjects by categories).
#' @keywords internal
sample_dirichlet_compositions <- function(counts, prior_pseudocounts = 1/2) {
    n_categories <- ncol(counts)
    alpha_increment <- if (length(prior_pseudocounts) == 1) {
        rep(prior_pseudocounts, n_categories)
    } else {
        prior_pseudocounts
    }

    draw_matrix <- matrix(
        NA_real_,
        nrow = nrow(counts),
        ncol = n_categories,
        dimnames = dimnames(counts)
    )

    for (subject_idx in seq_len(nrow(counts))) {
        alpha <- counts[subject_idx, ] + alpha_increment
        gamma_draws <- stats::rgamma(length(alpha), shape = alpha, rate = 1)
        draw_matrix[subject_idx, ] <- gamma_draws / sum(gamma_draws)
    }

    draw_matrix
}

#' Build distance matrix from one posterior draw of compositions.
#'
#' @param draw_matrix Subjects-by-categories frequency matrix (one draw).
#' @return Symmetric Hellinger distance matrix.
#' @keywords internal
build_distance_matrix_from_draw <- function(draw_matrix) {
    build_hellinger_distance_matrix(draw_matrix)
}
