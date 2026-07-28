#' Compare categorical compositions across groups.
#'
#' Tests whether between-group Hellinger distances exceed within-group distances
#' using a permutation null or conjugate Dirichlet posterior-predictive Bayes.
#' Computes an omnibus statistic plus all pairwise group contrasts and optional
#' custom contrasts (subsets or collapses).
#'
#' @param x A CategoryComposition object or long counts table.
#' @param method Inference method: `"permutation"` or `"bayes"`.
#' @param contrasts Optional named list of custom contrast specifications. Each
#'   element may be a length-2 character vector of group levels (pairwise subset)
#'   or a list with `collapse` (named remap of group levels) and/or `groups`
#'   (subset to these levels before testing).
#' @param nPermutations Number of label permutations (permutation method).
#' @param nPosterior Number of posterior draws (bayes method).
#' @param nCores Number of parallel workers (default 1).
#' @param seed Random seed for reproducibility.
#' @param priorPseudocounts Jeffreys prior increment per category (default 1/2).
#' @param pAdjustMethod Multiple-testing adjustment for contrasts (`"holm"`, `"BH"`, `"none"`).
#' @param ... Ignored.
#' @return A HellingerEnrichmentResult object.
#' @export
CompareGroupCompositions <- function(x,
                                     method = c("permutation", "bayes"),
                                     contrasts = NULL,
                                     nPermutations = 1000,
                                     nPosterior = 100,
                                     nCores = 1,
                                     seed = NULL,
                                     priorPseudocounts = 1/2,
                                     pAdjustMethod = "holm",
                                     ...) {
    method <- match.arg(method)

    composition <- if (inherits(x, "CategoryComposition")) {
        x
    } else if (is.data.frame(x)) {
        coerce_long_table_to_composition(x)
    } else {
        stop("x must be a CategoryComposition or long counts data.frame")
    }

    validate_composition_for_inference(composition)

    if (!is.null(seed)) {
        set.seed(seed)
    }

    contrast_specs <- build_contrast_specifications(
        group = composition$group,
        custom_contrasts = contrasts
    )

    if (method == "permutation") {
        result_body <- run_permutation_enrichment(
            composition = composition,
            contrast_specs = contrast_specs,
            n_permutations = nPermutations,
            n_cores = nCores,
            prior_pseudocounts = priorPseudocounts
        )
    } else {
        result_body <- run_bayes_enrichment(
            composition = composition,
            contrast_specs = contrast_specs,
            n_posterior = nPosterior,
            n_permutations = nPermutations,
            n_cores = nCores,
            prior_pseudocounts = priorPseudocounts
        )
    }

    contrasts_df <- result_body$contrasts
    if (nrow(contrasts_df) > 0 && pAdjustMethod != "none") {
        contrasts_df$pAdj <- stats::p.adjust(contrasts_df$pValue, method = pAdjustMethod)
    } else {
        contrasts_df$pAdj <- contrasts_df$pValue
    }

    structure(
        list(
            omnibus = result_body$omnibus,
            contrasts = contrasts_df,
            draws = result_body$draws,
            method = method,
            settings = list(
                nPermutations = nPermutations,
                nPosterior = nPosterior,
                nCores = nCores,
                seed = seed,
                priorPseudocounts = priorPseudocounts,
                pAdjustMethod = pAdjustMethod
            ),
            diagnostics = result_body$diagnostics
        ),
        class = "HellingerEnrichmentResult"
    )
}

#' Validate that groups have enough subjects for pairwise distances.
#' @keywords internal
validate_composition_for_inference <- function(composition) {
    group <- composition$group
    group_sizes <- table(group)
    singleton_groups <- names(group_sizes)[group_sizes < 2]

    if (length(singleton_groups) > 0) {
        stop(sprintf(
            "each group needs at least 2 subjects for pairwise distances; undersized: %s",
            paste(singleton_groups, collapse = ", ")
        ))
    }

    invisible(TRUE)
}

#' @export
print.HellingerEnrichmentResult <- function(x, ...) {
    cat(sprintf(
        "HellingerEnrichmentResult (%s): omnibus effectSize=%.4f, pValue=%.4f\n",
        x$method,
        x$omnibus$effectSize,
        x$omnibus$pValue
    ))
    cat(sprintf("  %d contrast(s) tested\n", nrow(x$contrasts)))
    invisible(x)
}

#' Build contrast specifications from group levels and custom overrides.
#' @keywords internal
build_contrast_specifications <- function(group, custom_contrasts = NULL) {
    group_levels <- levels(group)
    specs <- list()

    specs[["omnibus"]] <- list(
        contrastId = "omnibus",
        groups = group_levels,
        collapse = NULL
    )

    pair_indices <- utils::combn(group_levels, 2, simplify = FALSE)
    for (pair in pair_indices) {
        contrast_id <- paste(pair, collapse = "_vs_")
        specs[[contrast_id]] <- list(
            contrastId = contrast_id,
            groups = pair,
            collapse = NULL
        )
    }

    if (!is.null(custom_contrasts)) {
        if (!is.list(custom_contrasts) || is.null(names(custom_contrasts))) {
            stop("contrasts must be a named list")
        }
        for (contrast_name in names(custom_contrasts)) {
            spec <- normalize_contrast_spec(custom_contrasts[[contrast_name]], contrast_name, group_levels)
            specs[[contrast_name]] <- spec
        }
    }

    specs
}

#' @keywords internal
normalize_contrast_spec <- function(spec, contrast_id, all_group_levels) {
    if (is.character(spec) && length(spec) == 2) {
        return(list(
            contrastId = contrast_id,
            groups = spec,
            collapse = NULL
        ))
    }

    if (!is.list(spec)) {
        stop(sprintf("contrast '%s' must be a length-2 character vector or list", contrast_id))
    }

    groups <- spec$groups
    collapse <- spec$collapse

    if (is.null(groups) && is.null(collapse)) {
        stop(sprintf("contrast '%s' must specify groups and/or collapse", contrast_id))
    }

    list(
        contrastId = contrast_id,
        groups = groups,
        collapse = collapse
    )
}

#' Apply subset and collapse transforms to a group factor for one contrast.
#' @keywords internal
prepare_contrast_group <- function(group, spec) {
    working_group <- group

    if (!is.null(spec$collapse)) {
        collapse_map <- spec$collapse
        collapsed_values <- as.character(working_group)
        for (old_level in names(collapse_map)) {
            collapsed_values[working_group == old_level] <- collapse_map[[old_level]]
        }
        working_group <- factor(collapsed_values)
    }

    if (!is.null(spec$groups)) {
        keep_subjects <- names(working_group)[working_group %in% spec$groups]
        working_group <- working_group[keep_subjects]
        working_group <- factor(working_group, levels = spec$groups)
    }

    working_group
}

#' Subset count matrix and distance inputs for a contrast.
#' @keywords internal
subset_for_contrast <- function(counts, group, spec) {
    contrast_group <- prepare_contrast_group(group, spec)
    subject_ids <- names(contrast_group)
    list(
        counts = counts[subject_ids, , drop = FALSE],
        group = contrast_group
    )
}

#' Observed effect size and ratio for one contrast at fixed counts.
#' @keywords internal
compute_observed_contrast <- function(counts, group, spec, prior_pseudocounts) {
    subset_data <- subset_for_contrast(counts, group, spec)
    softened <- apply_prior_pseudocounts(subset_data$counts, prior_pseudocounts)
    dist_matrix <- build_hellinger_distance_matrix(softened)
    summary_matrix <- build_group_distance_summary_matrix(subset_data$group, dist_matrix)
    effect_size <- compute_between_within_ratio(summary_matrix)

    list(
        effectSize = effect_size,
        dist_matrix = dist_matrix,
        group = subset_data$group,
        summary_matrix = summary_matrix
    )
}

#' Permutation p-value for one contrast (includes observed in numerator).
#' @keywords internal
permutation_p_value <- function(observed_ratio, permuted_ratios) {
    mean(c(permuted_ratios, observed_ratio) >= observed_ratio)
}

#' Run permutation null for one contrast.
#' @keywords internal
run_contrast_permutation <- function(counts, group, spec, n_permutations, prior_pseudocounts) {
    observed <- compute_observed_contrast(counts, group, spec, prior_pseudocounts)
    subset_data <- subset_for_contrast(counts, group, spec)
    softened <- apply_prior_pseudocounts(subset_data$counts, prior_pseudocounts)
    dist_matrix <- build_hellinger_distance_matrix(softened)
    contrast_group <- subset_data$group

    permuted_ratios <- vapply(seq_len(n_permutations), function(perm_idx) {
        shuffled_group <- sample(contrast_group)
        names(shuffled_group) <- names(contrast_group)
        summary_matrix <- build_group_distance_summary_matrix(shuffled_group, dist_matrix)
        compute_between_within_ratio(summary_matrix)
    }, FUN.VALUE = numeric(1))

    p_value <- permutation_p_value(observed$effectSize, permuted_ratios)

    list(
        effectSize = observed$effectSize,
        pValue = p_value,
        permutedRatios = permuted_ratios
    )
}

#' Parallel worker wrapper for permutation contrasts.
#' @keywords internal
run_permutation_worker <- function(spec, counts, group, n_permutations, prior_pseudocounts) {
    run_contrast_permutation(counts, group, spec, n_permutations, prior_pseudocounts)
}

#' Execute permutation enrichment across all contrasts.
#' @keywords internal
run_permutation_enrichment <- function(composition, contrast_specs, n_permutations, n_cores, prior_pseudocounts) {
    counts <- composition$counts
    group <- composition$group

    contrast_names <- names(contrast_specs)
    run_one <- function(contrast_name) {
        spec <- contrast_specs[[contrast_name]]
        result <- run_contrast_permutation(counts, group, spec, n_permutations, prior_pseudocounts)
        data.frame(
            contrastId = spec$contrastId,
            effectSize = result$effectSize,
            pValue = result$pValue,
            stringsAsFactors = FALSE
        )
    }

    contrast_rows <- if (n_cores > 1 && length(contrast_names) > 1) {
        parallel::mclapply(contrast_names, run_one, mc.cores = min(n_cores, length(contrast_names)))
    } else {
        lapply(contrast_names, run_one)
    }

    contrasts_df <- do.call(rbind, contrast_rows)
    rownames(contrasts_df) <- NULL

    omnibus_row <- contrasts_df[contrasts_df$contrastId == "omnibus", , drop = FALSE]
    contrast_only <- contrasts_df[contrasts_df$contrastId != "omnibus", , drop = FALSE]

    list(
        omnibus = list(
            effectSize = omnibus_row$effectSize[1],
            pValue = omnibus_row$pValue[1]
        ),
        contrasts = contrast_only,
        draws = NULL,
        diagnostics = NULL
    )
}

#' One Bayes draw: sample compositions, rebuild distances, nested permutation null.
#' @keywords internal
run_bayes_draw <- function(counts, group, spec, n_permutations, prior_pseudocounts) {
    subset_data <- subset_for_contrast(counts, group, spec)
    draw_matrix <- sample_dirichlet_compositions(subset_data$counts, prior_pseudocounts)
    dist_matrix <- build_distance_matrix_from_draw(draw_matrix)
    contrast_group <- subset_data$group

    summary_matrix <- build_group_distance_summary_matrix(contrast_group, dist_matrix)
    observed_ratio <- compute_between_within_ratio(summary_matrix)

    permuted_ratios <- vapply(seq_len(n_permutations), function(perm_idx) {
        shuffled_group <- sample(contrast_group)
        names(shuffled_group) <- names(contrast_group)
        perm_summary <- build_group_distance_summary_matrix(shuffled_group, dist_matrix)
        compute_between_within_ratio(perm_summary)
    }, FUN.VALUE = numeric(1))

    list(
        effectSize = observed_ratio,
        pValue = permutation_p_value(observed_ratio, permuted_ratios)
    )
}

#' Execute Bayes enrichment with full nesting per posterior draw.
#' @keywords internal
run_bayes_enrichment <- function(composition, contrast_specs, n_posterior, n_permutations, n_cores, prior_pseudocounts) {
    counts <- composition$counts
    group <- composition$group
    contrast_names <- names(contrast_specs)

    run_contrast_bayes <- function(contrast_name) {
        spec <- contrast_specs[[contrast_name]]
        draw_results <- vector("list", n_posterior)

        draw_indices <- seq_len(n_posterior)
        compute_draw <- function(draw_idx) {
            run_bayes_draw(counts, group, spec, n_permutations, prior_pseudocounts)
        }

        if (n_cores > 1) {
            draw_results <- parallel::mclapply(draw_indices, compute_draw, mc.cores = min(n_cores, n_posterior))
        } else {
            draw_results <- lapply(draw_indices, compute_draw)
        }

        effect_sizes <- vapply(draw_results, function(x) x$effectSize, numeric(1))
        p_values <- vapply(draw_results, function(x) x$pValue, numeric(1))

        ci_quantiles <- stats::quantile(
            effect_sizes,
            probs = c(0.025, 0.975),
            na.rm = TRUE,
            names = FALSE
        )

        summary_row <- data.frame(
            contrastId = spec$contrastId,
            effectSize = mean(effect_sizes),
            pValue = mean(p_values),
            effectCiLow = ci_quantiles[1],
            effectCiHigh = ci_quantiles[2],
            stringsAsFactors = FALSE
        )

        draws_df <- data.frame(
            contrastId = spec$contrastId,
            draw = draw_indices,
            effectSize = effect_sizes,
            stringsAsFactors = FALSE
        )

        list(summary = summary_row, draws = draws_df)
    }

    contrast_results <- if (n_cores > 1 && length(contrast_names) > 1) {
        parallel::mclapply(contrast_names, run_contrast_bayes, mc.cores = min(n_cores, length(contrast_names)))
    } else {
        lapply(contrast_names, run_contrast_bayes)
    }

    contrasts_df <- do.call(rbind, lapply(contrast_results, function(x) x$summary))
    rownames(contrasts_df) <- NULL

    draws_df <- do.call(rbind, lapply(contrast_results, function(x) x$draws))
    rownames(draws_df) <- NULL

    omnibus_row <- contrasts_df[contrasts_df$contrastId == "omnibus", , drop = FALSE]
    contrast_only <- contrasts_df[contrasts_df$contrastId != "omnibus", , drop = FALSE]
    draws_only <- draws_df[draws_df$contrastId != "omnibus", , drop = FALSE]

    list(
        omnibus = list(
            effectSize = omnibus_row$effectSize[1],
            pValue = omnibus_row$pValue[1],
            effectCiLow = omnibus_row$effectCiLow[1],
            effectCiHigh = omnibus_row$effectCiHigh[1]
        ),
        contrasts = contrast_only,
        draws = draws_only,
        diagnostics = list(nPosterior = n_posterior, nPermutations = n_permutations)
    )
}
