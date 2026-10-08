#' Compare categorical compositions across groups.
#'
#' Tests whether between-group Hellinger distances exceed within-group distances
#' using a permutation null or conjugate Dirichlet posterior-predictive Bayes.
#' Always reports an omnibus statistic. When `contrasts` is omitted, also tests
#' every pairwise group comparison; when `contrasts` is supplied, those
#' specifications are the contrast family (auto pairwises are omitted unless
#' `includePairwise = TRUE`).
#'
#' The between/within Hellinger ratio, Jeffreys softening, observed-inclusive
#' label-permutation p-value, and collapse/subset contrasts follow Paul
#' Edlefsen's original procedure. `method = "bayes"` instead draws subject
#' compositions from a conjugate Dirichlet posterior, rebuilds the Hellinger
#' ratio \(R\) on each draw, and reports the posterior probability that
#' \(R > 1\) (`PPGT1`) with a 95% highest posterior density interval on \(R\)
#' by default.
#'
#' @param x A CategoryComposition object or long counts table.
#' @param method Inference method: `"permutation"` or `"bayes"`.
#' @param contrasts Planned comparisons, or `NULL` for every pairwise contrast.
#'   A length-2 character vector is one pairwise comparison. A list may mix
#'   length-2 character vectors and lists with `collapse` (named remap of group
#'   levels) and/or `groups` (subset to these levels before testing). Names
#'   become `contrastId` values; unnamed elements get `A_vs_B` identifiers.
#' @param nPermutations Number of label permutations (permutation method only).
#' @param nPosterior Number of posterior draws (bayes method).
#' @param nCores Number of parallel workers (default 1).
#' @param seed Random seed for reproducibility.
#' @param priorPseudocounts Jeffreys prior increment per category (default 1/2).
#' @param pAdjustMethod Multiple-testing adjustment for permutation contrasts (`"holm"`, `"BH"`, `"none"`). Ignored for Bayes.
#' @param credibleInterval Bayes interval type: `"hpdi"` (default, shortest 95%
#'   interval on the draws) or `"cpi"` (equal-tailed 2.5%/97.5% quantiles).
#'   Ignored for permutation.
#' @param includePairwise Include every pairwise comparison in addition to
#'   `contrasts`. Default `NULL` includes pairwises when `contrasts` is omitted
#'   (the usual all-pairs analysis) and skips them when `contrasts` is supplied
#'   so Holm is paid only on the planned family. `FALSE` skips auto pairwises
#'   even when `contrasts` is omitted (omnibus only). `TRUE` keeps all pairwises
#'   even when `contrasts` is supplied (the old additive overlay).
#' @param ... Ignored.
#' @return A HellingerEnrichmentResult with omnibus summary, contrast table, and
#'   (Bayes only) long-format posterior `draws` including the omnibus draws.
#' @references
#' Edlefsen, P. Hellinger between/within enrichment for subject-level
#' categorical compositions (Jeffreys softening, ratio statistic, label
#' permutation, collapse/subset contrasts).
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
                                     credibleInterval = c("hpdi", "cpi"),
                                     includePairwise = NULL,
                                     ...) {
    method <- match.arg(method)
    credibleInterval <- match.arg(credibleInterval)

    if (method == "bayes") {
        unused <- character()
        if (nPermutations != 1000) {
            unused <- c(unused, "nPermutations")
        }
        if (pAdjustMethod != "holm") {
            unused <- c(unused, "pAdjustMethod")
        }
        if (length(unused) > 0) {
            warning(
                sprintf(
                    "%s unused for method = \"bayes\"; ignored",
                    paste(unused, collapse = " and ")
                ),
                call. = FALSE
            )
        }
    } else if (!identical(credibleInterval, "hpdi")) {
        warning("credibleInterval unused for method = \"permutation\"; ignored", call. = FALSE)
    }

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

    include_pairwise <- resolve_include_pairwise(includePairwise, !is.null(contrasts))
    if (!include_pairwise && is.null(contrasts)) {
        warning(
            "includePairwise is FALSE and no contrasts were supplied; only the omnibus test will run",
            call. = FALSE
        )
    }

    contrast_specs <- build_contrast_specifications(
        group = composition$group,
        custom_contrasts = contrasts,
        include_pairwise = include_pairwise
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
            n_cores = nCores,
            prior_pseudocounts = priorPseudocounts,
            credible_interval = credibleInterval
        )
    }

    contrasts_df <- result_body$contrasts
    if (method == "permutation") {
        if (nrow(contrasts_df) > 0 && pAdjustMethod != "none") {
            contrasts_df$pAdj <- stats::p.adjust(contrasts_df$pValue, method = pAdjustMethod)
        } else {
            contrasts_df$pAdj <- contrasts_df$pValue
        }
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
                pAdjustMethod = pAdjustMethod,
                credibleInterval = credibleInterval,
                includePairwise = include_pairwise
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
    if (identical(x$method, "bayes")) {
        cat(sprintf(
            "HellingerEnrichmentResult (%s): omnibus effectSize=%.4f, PPGT1=%.4f\n",
            x$method,
            x$omnibus$effectSize,
            x$omnibus$PPGT1
        ))
    } else {
        cat(sprintf(
            "HellingerEnrichmentResult (%s): omnibus effectSize=%.4f, pValue=%.4f\n",
            x$method,
            x$omnibus$effectSize,
            x$omnibus$pValue
        ))
    }
    cat(sprintf("  %d contrast(s) tested\n", nrow(x$contrasts)))
    invisible(x)
}

#' 95% posterior interval on finite draws: shortest HPDI or equal-tailed CPI.
#'
#' HPDI uses the Chen-Shao shortest window covering `ceil(prob * n)` sorted
#' samples. Quantile intervals spend equal tail mass, which is a poor summary
#' when the posterior of \(R\) is asymmetric and sits just above 1; the
#' shortest covering window is the HPD for a continuous unimodal sample.
#' @keywords internal
compute_posterior_interval <- function(x, prob = 0.95, type = c("hpdi", "cpi")) {
    type <- match.arg(type)
    x <- x[is.finite(x)]
    n <- length(x)
    if (n < 2) {
        return(c(NA_real_, NA_real_))
    }

    if (type == "cpi") {
        tail_mass <- (1 - prob) / 2
        return(stats::quantile(
            x,
            probs = c(tail_mass, 1 - tail_mass),
            na.rm = TRUE,
            names = FALSE
        ))
    }

    x_sorted <- sort(x)
    n_cover <- min(n, ceiling(prob * n))
    n_windows <- n - n_cover + 1
    widths <- x_sorted[seq.int(n_cover, n)] - x_sorted[seq_len(n_windows)]
    best <- which.min(widths)[1]
    c(x_sorted[best], x_sorted[best + n_cover - 1])
}

#' Resolve whether auto pairwises run given planned contrasts.
#'
#' NULL follows the planned-comparison rule: all pairwises when `contrasts` is
#' omitted, none when a family is supplied. FALSE is omnibus-only if no family
#' is supplied; TRUE restores the old additive overlay.
#' @keywords internal
resolve_include_pairwise <- function(include_pairwise, has_custom_contrasts) {
    if (is.null(include_pairwise)) {
        return(!has_custom_contrasts)
    }
    if (!is.logical(include_pairwise) || length(include_pairwise) != 1 || is.na(include_pairwise)) {
        stop("includePairwise must be TRUE, FALSE, or NULL")
    }
    include_pairwise
}

#' Coerce a length-2 vector or (un)named list into a named contrast list.
#' @keywords internal
coerce_contrast_list <- function(contrasts) {
    if (is.null(contrasts)) {
        return(NULL)
    }

    if (is.character(contrasts) && length(contrasts) == 2) {
        contrast_id <- paste(contrasts, collapse = "_vs_")
        out <- list(contrasts)
        names(out) <- contrast_id
        return(out)
    }

    if (!is.list(contrasts)) {
        stop("contrasts must be a length-2 character vector or a list of contrast specifications")
    }
    if (length(contrasts) == 0) {
        stop("contrasts must contain at least one specification")
    }

    contrast_names <- names(contrasts)
    if (is.null(contrast_names)) {
        contrast_names <- rep("", length(contrasts))
    }

    auto_ids <- vapply(seq_along(contrasts), function(idx) {
        if (nzchar(contrast_names[[idx]])) {
            return(contrast_names[[idx]])
        }
        spec <- contrasts[[idx]]
        if (is.character(spec) && length(spec) == 2) {
            return(paste(spec, collapse = "_vs_"))
        }
        if (is.list(spec) && !is.null(spec$groups) && length(spec$groups) > 0) {
            return(paste(spec$groups, collapse = "_vs_"))
        }
        paste0("contrast_", idx)
    }, character(1))

    names(contrasts) <- make.unique(auto_ids, sep = "_")
    contrasts
}

#' Build contrast specifications from group levels and planned comparisons.
#' @keywords internal
build_contrast_specifications <- function(group, custom_contrasts = NULL, include_pairwise = TRUE) {
    group_levels <- levels(group)
    specs <- list()

    specs[["omnibus"]] <- list(
        contrastId = "omnibus",
        groups = group_levels,
        collapse = NULL
    )

    if (isTRUE(include_pairwise)) {
        pair_indices <- utils::combn(group_levels, 2, simplify = FALSE)
        for (pair in pair_indices) {
            contrast_id <- paste(pair, collapse = "_vs_")
            specs[[contrast_id]] <- list(
                contrastId = contrast_id,
                groups = pair,
                collapse = NULL
            )
        }
    }

    custom_contrasts <- coerce_contrast_list(custom_contrasts)
    if (!is.null(custom_contrasts)) {
        for (contrast_name in names(custom_contrasts)) {
            spec <- normalize_contrast_spec(
                custom_contrasts[[contrast_name]],
                contrast_name,
                group_levels
            )
            specs[[contrast_name]] <- spec
        }
    }

    specs
}

#' @keywords internal
normalize_contrast_spec <- function(spec, contrast_id, all_group_levels) {
    if (is.character(spec) && length(spec) == 2) {
        unknown_groups <- setdiff(spec, all_group_levels)
        if (length(unknown_groups) > 0) {
            stop(sprintf(
                "contrast '%s' refers to unknown group(s): %s",
                contrast_id,
                paste(unknown_groups, collapse = ", ")
            ))
        }
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

    if (!is.null(collapse)) {
        unknown_from <- setdiff(names(collapse), all_group_levels)
        if (length(unknown_from) > 0) {
            stop(sprintf(
                "contrast '%s' collapse refers to unknown group(s): %s",
                contrast_id,
                paste(unknown_from, collapse = ", ")
            ))
        }
    }

    if (!is.null(groups) && is.null(collapse)) {
        unknown_groups <- setdiff(groups, all_group_levels)
        if (length(unknown_groups) > 0) {
            stop(sprintf(
                "contrast '%s' refers to unknown group(s): %s",
                contrast_id,
                paste(unknown_groups, collapse = ", ")
            ))
        }
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
#'
#' Matches Edlefsen's workbook form:
#' `mean(c(permuted_ratios, observed_ratio) >= observed_ratio)`.
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

#' One Bayes draw: sample compositions and rebuild the Hellinger ratio.
#' @keywords internal
run_bayes_draw <- function(counts, group, spec, prior_pseudocounts) {
    subset_data <- subset_for_contrast(counts, group, spec)
    draw_matrix <- sample_dirichlet_compositions(subset_data$counts, prior_pseudocounts)
    dist_matrix <- build_distance_matrix_from_draw(draw_matrix)
    contrast_group <- subset_data$group

    summary_matrix <- build_group_distance_summary_matrix(contrast_group, dist_matrix)
    observed_ratio <- compute_between_within_ratio(summary_matrix)

    list(effectSize = observed_ratio)
}

#' Execute Bayes enrichment from Dirichlet posterior draws of R.
#'
#' For each contrast, returns both a one-row posterior summary and the long
#' draw table used by half-eye plots. Draw storage scales as
#' `n_posterior * n_contrasts` effect sizes, including omnibus, so the
#' posterior of the global ratio can be inspected alongside pairwise draws.
#' @keywords internal
run_bayes_enrichment <- function(composition,
                                 contrast_specs,
                                 n_posterior,
                                 n_cores,
                                 prior_pseudocounts,
                                 credible_interval = "hpdi") {
    counts <- composition$counts
    group <- composition$group
    contrast_names <- names(contrast_specs)

    run_contrast_bayes <- function(contrast_name) {
        spec <- contrast_specs[[contrast_name]]
        draw_indices <- seq_len(n_posterior)
        compute_draw <- function(draw_idx) {
            run_bayes_draw(counts, group, spec, prior_pseudocounts)
        }

        draw_results <- if (n_cores > 1) {
            parallel::mclapply(draw_indices, compute_draw, mc.cores = min(n_cores, n_posterior))
        } else {
            lapply(draw_indices, compute_draw)
        }

        effect_sizes <- vapply(draw_results, function(x) x$effectSize, numeric(1))
        ci_bounds <- compute_posterior_interval(
            effect_sizes,
            prob = 0.95,
            type = credible_interval
        )

        summary_row <- data.frame(
            contrastId = spec$contrastId,
            effectSize = mean(effect_sizes),
            PPGT1 = mean(effect_sizes > 1),
            effectCiLow = ci_bounds[1],
            effectCiHigh = ci_bounds[2],
            stringsAsFactors = FALSE
        )

        # Keep per-draw effect sizes for PlotCompositionContrasts half-eyes;
        # summarizing to CI alone would force geom_interval and lose the slab.
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
    omnibus_draws <- draws_df$effectSize[draws_df$contrastId == "omnibus"]

    list(
        omnibus = list(
            effectSize = omnibus_row$effectSize[1],
            PPGT1 = omnibus_row$PPGT1[1],
            effectCiLow = omnibus_row$effectCiLow[1],
            effectCiHigh = omnibus_row$effectCiHigh[1],
            draws = omnibus_draws
        ),
        contrasts = contrast_only,
        draws = draws_df,
        diagnostics = list(nPosterior = n_posterior)
    )
}
