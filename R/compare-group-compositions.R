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
#' ratio \eqn{R} on each draw, and reports the posterior probability that
#' \eqn{R > 1} (`PPGT1`) with a 95% highest posterior density interval on \eqn{R}
#' by default.
#'
#' `PPGT1` alone overstates enrichment when counts are deep. The Dirichlet
#' posterior only describes count noise within each subject, so as depth grows
#' it concentrates on the sample \eqn{R}, which is rarely exactly 1 even when
#' groups do not differ. With `nope = "shuffle"` (the default), the package
#' reports a Null-Oriented Permutation Equivalence (NOPE) interval on \eqn{R}:
#' a permutation-calibrated null reference interval for conditional effects. Within each posterior draw, the package shuffles group
#' labels on that draw's own distance matrix, so the interval is conditional on
#' the observed design (group sizes, depths, within-group spread) and on that
#' draw's compositions. The goal parallels a permutation test, but NOPE
#' reports a region instead of a p-value. The NOPE interval for \eqn{R} runs from 1 to the 95th
#' percentile of the shuffled \eqn{R}; `R_nopeHigh` is the median of that upper
#' bound across draws, and `R_PPGTnope` is the posterior probability that
#' \eqn{R} clears it. `energyNope = "shuffle"` (also the default) builds the
#' matching interval for the energy distance \eqn{E} from the same shuffles,
#' from 0 to `E_nopeHigh`, and reports `E_PPGTnope`.
#'
#' NOPE borrows the equivalence logic of Kruschke's region of practical
#' equivalence (ROPE) and changes two things. In Kruschke's procedure the
#' analyst fixes the region before seeing the data, from domain knowledge of
#' which effects are negligible, and the decision compares the 95% highest
#' density interval with that region. Few analysts have that intuition for a
#' ratio of Hellinger distances or an energy distance, so NOPE takes its upper
#' bound from the data: the largest value of the statistic these subjects
#' routinely produce when labels carry no information. The bound also differs
#' from draw to draw, and each NOPE probability counts the draws that clear
#' their own bound. NOPE therefore acts as a null reference calibrated by
#' permutation, closer to a permutation test than to a statement of practical
#' significance. When you know the smallest effect that matters, pass it as a
#' fixed NOPE bound, which plays the role of Kruschke's preset region:
#' `nope = <number above 1>` on \eqn{R}, or `energyNope = <number above 0>` on
#' \eqn{E} in squared Hellinger units. A fixed bound needs no shuffles, so it
#' also works for groups too small for a NOPE interval.
#'
#' To set a fixed bound on \eqn{E} from a composition shift, pass the reference
#' and shifted compositions to ComputeEnergyNopeBound and hand its result to
#' `energyNope`. It returns \eqn{2(1 - \mathrm{BC})}{2(1 - BC)}, with BC the Bhattacharyya
#' coefficient between the two compositions, which equals \eqn{E} when every
#' subject sits at its group's mean and approximates it otherwise.
#'
#' Because shuffling keeps the observed group sizes, depths, and within-group
#' spread, the NOPE widens or narrows with the data rather than relying on a
#' fixed number. Each draw's edge never falls below the no-difference value (1
#' for \eqn{R}, 0 for \eqn{E}), so `R_PPGTnope <= PPGT1` always holds, and
#' `E_PPGTnope` never exceeds the posterior probability that \eqn{E > 0}.
#' Reading a NOPE probability above 0.95 asks that nearly every draw clear
#' its own 95th shuffle percentile, which runs well below a 5% false-positive
#' rate and costs power; a looser threshold such as 0.5 behaves more like the
#' permutation test. The sampling-depth vignette measures both. Very small
#' groups limit the NOPE the same way they limit the permutation p-value.
#' With 3 versus 3 subjects there are only 10 distinct splits, so the true
#' split alone is 10% of all shuffles and neither statistic can clear its 95th
#' percentile. Whenever the true split is at least 5% of shuffles (e.g. 2
#' versus 2, 3 versus 3, 2 versus 4), the function sets each NOPE
#' probability to `NA` with a warning rather than a misleading 0, and
#' `enoughShufflesForNope` is `FALSE` so the reason travels with the result;
#' 4 versus 4 (35 splits) is the smallest balanced design that supports it.
#'
#' Both methods also report the energy distance between groups
#' (`energyDistance`): twice the mean between-group squared Hellinger distance
#' minus each group's mean within-group squared distance, averaged over group
#' pairs for the omnibus. Squaring the existing distance matrix costs nothing
#' extra. \eqn{R} measures how far apart the groups are relative to how spread
#' out they are; \eqn{E} measures the shift in mean composition itself, as the
#' squared distance between group centroids in square-root composition space,
#' with no division. Independent count noise adds equally to between- and
#' within-group squared distances and cancels from \eqn{E} in expectation, but
#' Jeffreys softening and the curvature of the square root still shrink \eqn{E}
#' at shallow depth, though less than they shrink the excess of \eqn{R} over 1
#' (see the sampling-depth vignette). \eqn{E} responds mainly to shifts in mean
#' composition; a difference in spread alone moves it only slightly, likely
#' through that same curvature. The sample estimate excludes each subject's
#' distance to itself, so it can dip slightly below 0 when groups do not
#' differ. The p-values and `PPGT1` stay on \eqn{R}; under Bayes the package
#' reports a posterior interval on \eqn{E} (`energyCiLow`, `energyCiHigh`) and,
#' with `energyNope`, `E_PPGTnope`.
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
#' @param nope Bayes NOPE interval on \eqn{R}. `"shuffle"` (default) sets a
#'   permutation-calibrated null reference interval whose upper bound comes from label
#'   shuffles within each posterior draw. A single number above 1 sets a fixed
#'   NOPE bound on \eqn{R} that you choose in advance.
#'   Either adds `R_PPGTnope` and `R_nopeHigh` to the output. `"none"` drops
#'   them. Ignored for permutation.
#' @param energyNope Bayes NOPE interval on the energy distance \eqn{E}, in
#'   squared Hellinger units. `"shuffle"` (default) sets a NOPE interval from
#'   the same label shuffles as `nope`. A single number above 0 sets a fixed
#'   NOPE bound on \eqn{E}; [ComputeEnergyNopeBound()] turns a composition shift
#'   into one. Either adds `E_PPGTnope` and
#'   `E_nopeHigh` to the output. `"none"` drops them. Ignored for permutation.
#' @param nNopeShuffles Label shuffles per posterior draw used to set the
#'   NOPEs (default 200). Each shuffle builds one group summary per
#'   shuffled statistic, so cost scales as `nPosterior * nNopeShuffles`
#'   summaries per contrast for one NOPE and twice that for both. Used only
#'   when `method = "bayes"` and `nope` or
#'   `energyNope` is `"shuffle"`.
#' @param ... Ignored.
#' @return A HellingerEnrichmentResult with omnibus summary, contrast table, and
#'   (Bayes only) long-format posterior `draws` including the omnibus draws.
#'   Bayes summaries carry `R_PPGTnope` and `R_nopeHigh` unless
#'   `nope = "none"`, `E_PPGTnope` and `E_nopeHigh` unless
#'   `energyNope = "none"` (each edge column is the median across draws of the
#'   per-draw upper edge), and `enoughShufflesForNope` when `nope` or `energyNope` is
#'   `"shuffle"` (`FALSE` when groups are too small for the NOPE, which is
#'   why the NOPE probabilities are then `NA`); `draws` carries the per-draw
#'   `R_nopeHigh` and `E_nopeHigh`. Both methods carry `energyDistance` beside
#'   the ratio `effectSize`; Bayes adds `energyCiLow`/`energyCiHigh` and
#'   per-draw `energyDistance` in `draws`.
#' @references
#' Edlefsen, P. Hellinger between/within enrichment for subject-level
#' categorical compositions (Jeffreys softening, ratio statistic, label
#' permutation, collapse/subset contrasts).
#'
#' Kruschke, J. K. (2018). Rejecting or accepting parameter values in Bayesian
#' estimation. Advances in Methods and Practices in Psychological Science,
#' 1(2), 270-280. \doi{10.1177/2515245918771304}. Article page:
#' \url{https://www.psychologicalscience.org/journals/ampps/2515245918771304/}; open-access record:
#' \url{https://scholarworks.iu.edu/dspace/items/bada026d-fe1b-4f53-ad5b-db6e8bf99711}
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
                                     nope = "shuffle",
                                     energyNope = "shuffle",
                                     nNopeShuffles = 200,
                                     ...) {
    method <- match.arg(method)
    credibleInterval <- match.arg(credibleInterval)
    nope <- validate_nope_setting(nope, scale = "ratio")
    energyNope <- validate_nope_setting(energyNope, scale = "energy")

    if (!is.numeric(nNopeShuffles) || length(nNopeShuffles) != 1 ||
        is.na(nNopeShuffles) || nNopeShuffles < 1) {
        stop("nNopeShuffles must be a single positive number")
    }

    unused <- character()
    if (method == "bayes") {
        if (nPermutations != 1000) {
            unused <- c(unused, "nPermutations")
        }
        if (pAdjustMethod != "holm") {
            unused <- c(unused, "pAdjustMethod")
        }
        if (!identical(nope, "shuffle") && !identical(energyNope, "shuffle") && nNopeShuffles != 200) {
            warning("nNopeShuffles unused unless nope or energyNope is \"shuffle\"; ignored", call. = FALSE)
        }
    } else {
        if (!identical(credibleInterval, "hpdi")) {
            unused <- c(unused, "credibleInterval")
        }
        if (!identical(nope, "shuffle")) {
            unused <- c(unused, "nope")
        }
        if (!identical(energyNope, "shuffle")) {
            unused <- c(unused, "energyNope")
        }
        if (nNopeShuffles != 200) {
            unused <- c(unused, "nNopeShuffles")
        }
    }
    if (length(unused) > 0) {
        warning(
            sprintf(
                "%s unused for method = \"%s\"; ignored",
                paste(unused, collapse = " and "),
                method
            ),
            call. = FALSE
        )
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
            credible_interval = credibleInterval,
            nope = nope,
            energy_nope = energyNope,
            n_nope_shuffles = nNopeShuffles
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
                includePairwise = include_pairwise,
                nope = nope,
                energyNope = energyNope,
                nNopeShuffles = nNopeShuffles
            ),
            diagnostics = result_body$diagnostics
        ),
        class = "HellingerEnrichmentResult"
    )
}

#' Print one omnibus NOPE line for the ratio or the energy distance.
#'
#' @param nope_probability `R_PPGTnope` or `E_PPGTnope`, or `NULL` when the
#'   result carries no NOPE or fixed edge on that statistic.
#' @param nope_high Matching `R_nopeHigh` or `E_nopeHigh`.
#' @param nope_setting `x$settings$nope` or `x$settings$energyNope`.
#' @param statistic_label `"R"` or `"E"`.
#' @param no_difference_value Lower edge of the interval (1 for R, 0 for E).
#' @param omnibus The result's omnibus list, for `enoughShufflesForNope`.
#' @return `NULL`, invisibly; called for its printed line.
#' @keywords internal
print_nope_line <- function(nope_probability, nope_high, nope_setting, statistic_label, no_difference_value, omnibus) {
    if (is.null(nope_probability)) {
        return(invisible(NULL))
    }
    if (identical(nope_setting, "shuffle") && isFALSE(omnibus$enoughShufflesForNope)) {
        cat(sprintf(
            "  omnibus %s_PPGTnope=NA (groups too small: too few distinct label shuffles for the NOPE)\n",
            statistic_label
        ))
        return(invisible(NULL))
    }
    nope_label <- if (identical(nope_setting, "shuffle")) {
        "NOPE, permutation-calibrated null reference interval"
    } else {
        "fixed NOPE bound"
    }
    cat(sprintf(
        "  omnibus %s_PPGTnope=%.4f (%s on %s: [%g, %.4g])\n",
        statistic_label,
        nope_probability,
        nope_label,
        statistic_label,
        no_difference_value,
        nope_high
    ))
    invisible(NULL)
}

#' Validate a `nope` or `energyNope` argument.
#'
#' Accepts `"shuffle"` (a NOPE interval), `"none"`, or one finite number
#' above the statistic's no-difference value: above 1 for the ratio \eqn{R},
#' above 0 for the energy distance \eqn{E}. An edge at or below that value
#' would let draws consistent with no difference count as clearing the interval,
#' so the function rejects it.
#'
#' @param nope_setting Value passed to CompareGroupCompositions.
#' @param scale `"ratio"` for `nope`, `"energy"` for `energyNope`.
#' @return `"shuffle"`, `"none"`, or the numeric edge.
#' @keywords internal
validate_nope_setting <- function(nope_setting, scale = c("ratio", "energy")) {
    scale <- match.arg(scale)
    no_difference_value <- if (identical(scale, "ratio")) 1 else 0
    if (is.character(nope_setting) && length(nope_setting) == 1 && nope_setting %in% c("shuffle", "none")) {
        return(nope_setting)
    }
    if (is.numeric(nope_setting) && length(nope_setting) == 1 && is.finite(nope_setting) &&
        nope_setting > no_difference_value) {
        return(as.numeric(nope_setting))
    }
    if (identical(scale, "ratio")) {
        stop("nope must be \"shuffle\", \"none\", or a single number above 1 (a fixed upper edge on R)")
    }
    stop("energyNope must be \"shuffle\", \"none\", or a single number above 0 (a fixed upper edge on E)")
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
        print_nope_line(x$omnibus$R_PPGTnope, x$omnibus$R_nopeHigh, x$settings$nope, "R", 1, x$omnibus)
        print_nope_line(x$omnibus$E_PPGTnope, x$omnibus$E_nopeHigh, x$settings$energyNope, "E", 0, x$omnibus)
    } else {
        cat(sprintf(
            "HellingerEnrichmentResult (%s): omnibus effectSize=%.4f, pValue=%.4f\n",
            x$method,
            x$omnibus$effectSize,
            x$omnibus$pValue
        ))
    }
    if (!is.null(x$omnibus$energyDistance)) {
        cat(sprintf("  omnibus energyDistance=%.4f (squared Hellinger units)\n", x$omnibus$energyDistance))
    }
    cat(sprintf("  %d contrast(s) tested\n", nrow(x$contrasts)))
    invisible(x)
}

#' 95% posterior interval on finite draws: shortest HPDI or equal-tailed CPI.
#'
#' HPDI uses the Chen-Shao shortest window covering `ceil(prob * n)` sorted
#' samples. Quantile intervals spend equal tail mass, which is a poor summary
#' when the posterior of \eqn{R} is asymmetric and sits just above 1; the
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
        energyDistance = compute_energy_distance(
            build_group_distance_summary_matrix(subset_data$group, dist_matrix^2)
        ),
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

#' Between/within ratios after shuffling group labels on a fixed distance matrix.
#'
#' Serves the permutation null on one fixed matrix from softened counts; the
#' Bayes NOPEs use compute_shuffled_statistics, which makes the same
#' `sample()` calls and adds \eqn{E}. Shuffling reassigns labels among subjects, so group sizes and the distance matrix stay fixed while the
#' shuffle breaks any link between label and composition.
#'
#' @param contrast_group Named factor aligned to `dist_matrix` rownames.
#' @param dist_matrix Symmetric subject Hellinger distance matrix.
#' @param n_shuffles Number of label shuffles.
#' @return Numeric vector of `n_shuffles` shuffled ratios.
#' @keywords internal
compute_shuffled_ratios <- function(contrast_group, dist_matrix, n_shuffles) {
    vapply(seq_len(n_shuffles), function(shuffle_idx) {
        shuffled_group <- sample(contrast_group)
        names(shuffled_group) <- names(contrast_group)
        summary_matrix <- build_group_distance_summary_matrix(shuffled_group, dist_matrix)
        compute_between_within_ratio(summary_matrix)
    }, FUN.VALUE = numeric(1))
}

#' Ratio and energy distance after shuffling group labels on one draw's matrix.
#'
#' Runs the same `sample()` calls as compute_shuffled_ratios, so a Bayes draw
#' consumes the same random numbers whichever statistics it needs, and builds
#' two group summaries per shuffle: plain Hellinger distances for \eqn{R} and
#' squared distances for \eqn{E}. The permutation test keeps the ratio-only
#' helper, which costs half as much per shuffle.
#'
#' @param contrast_group Named factor aligned to `dist_matrix` rownames.
#' @param dist_matrix Symmetric subject Hellinger distance matrix.
#' @param n_shuffles Number of label shuffles.
#' @return Matrix with `n_shuffles` rows and columns `ratio` and `energy`.
#' @keywords internal
compute_shuffled_statistics <- function(contrast_group, dist_matrix, n_shuffles) {
    squared_dist_matrix <- dist_matrix^2
    shuffled <- vapply(seq_len(n_shuffles), function(shuffle_idx) {
        shuffled_group <- sample(contrast_group)
        names(shuffled_group) <- names(contrast_group)
        c(
            ratio = compute_between_within_ratio(
                build_group_distance_summary_matrix(shuffled_group, dist_matrix)
            ),
            energy = compute_energy_distance(
                build_group_distance_summary_matrix(shuffled_group, squared_dist_matrix)
            )
        )
    }, FUN.VALUE = c(ratio = 0, energy = 0))
    t(shuffled)
}

#' Check whether a contrast has enough distinct label shuffles for a NOPE.
#'
#' Each NOPE edge is the `nope_quantile` of its statistic (\eqn{R} or \eqn{E})
#' across label shuffles, so it only works as a bar when the true grouping is a
#' small fraction of all shuffles. With few subjects there are few distinct
#' ways to split them; the true split then sets the quantile itself and neither
#' statistic can clear it, the
#' same limit that keeps a 3 versus 3 permutation p-value at or above 0.1.
#'
#' \eqn{R} and \eqn{E} depend only on which subjects share a group, so relabelings that swap
#' two equal-sized groups count as one split. The number of distinct splits is
#' `n! / (prod(n_g!) * prod(m_k!))`, where `n_g` are group sizes and `m_k`
#' counts groups of each size (10 for 3 versus 3, 35 for 4 versus 4). The
#' contrast passes when the true split is under `1 - nope_quantile` of all
#' shuffles, i.e. more than 20 distinct splits at the default 0.95. Depends
#' only on group sizes, so it uses no random numbers and runs once per contrast.
#'
#' @param contrast_group Factor of group labels for one contrast.
#' @param nope_quantile Quantile of each shuffled statistic used as its NOPE
#'   edge.
#' @return `TRUE` when \eqn{R} and \eqn{E} can clear their NOPEs, `FALSE`
#'   otherwise.
#' @keywords internal
has_enough_shuffles_for_nope <- function(contrast_group, nope_quantile = 0.95) {
    group_sizes <- as.vector(table(droplevels(contrast_group)))
    n_groups_per_size <- as.vector(table(group_sizes))

    # Log scale because n! overflows past about 170 subjects.
    log_n_distinct_splits <- lfactorial(sum(group_sizes)) -
        sum(lfactorial(group_sizes)) -
        sum(lfactorial(n_groups_per_size))
    min_distinct_splits <- 1 / (1 - nope_quantile)

    # The tolerance keeps exact boundary cases (20 splits at 0.95) from passing
    # on floating-point error in 1 - nope_quantile.
    log_n_distinct_splits > log(min_distinct_splits) + 1e-8
}

#' Run permutation null for one contrast.
#' @keywords internal
run_contrast_permutation <- function(counts, group, spec, n_permutations, prior_pseudocounts) {
    observed <- compute_observed_contrast(counts, group, spec, prior_pseudocounts)
    subset_data <- subset_for_contrast(counts, group, spec)
    softened <- apply_prior_pseudocounts(subset_data$counts, prior_pseudocounts)
    dist_matrix <- build_hellinger_distance_matrix(softened)
    contrast_group <- subset_data$group

    permuted_ratios <- compute_shuffled_ratios(contrast_group, dist_matrix, n_permutations)

    p_value <- permutation_p_value(observed$effectSize, permuted_ratios)

    list(
        effectSize = observed$effectSize,
        energyDistance = observed$energyDistance,
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
            energyDistance = result$energyDistance,
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
            energyDistance = omnibus_row$energyDistance[1],
            pValue = omnibus_row$pValue[1]
        ),
        contrasts = contrast_only,
        draws = NULL,
        diagnostics = NULL
    )
}

#' One Bayes draw: sample compositions, rebuild \eqn{R} and \eqn{E}, and set this
#' draw's NOPE upper edges.
#'
#' When `nope` or `energy_nope` is `"shuffle"`, the function shuffles labels on this
#' draw's own distance matrix, so each edge carries the same count noise as
#' this draw's statistic rather than the noise of one fixed point-estimate
#' matrix. One set of shuffles serves both statistics. Each NOPE edge is the
#' `nope_quantile` of its shuffled statistic, and the function floors it at the
#' no-difference value (1 for \eqn{R}, 0 for \eqn{E}): a quantile below that value
#' would invert the interval and let a draw below no difference count as clearing
#' it. A numeric setting is a fixed edge and needs no shuffles.
#' @keywords internal
run_bayes_draw <- function(counts,
                           group,
                           spec,
                           prior_pseudocounts,
                           nope = "shuffle",
                           energy_nope = "shuffle",
                           n_nope_shuffles = 200,
                           nope_quantile = 0.95) {
    subset_data <- subset_for_contrast(counts, group, spec)
    draw_matrix <- sample_dirichlet_compositions(subset_data$counts, prior_pseudocounts)
    dist_matrix <- build_distance_matrix_from_draw(draw_matrix)
    contrast_group <- subset_data$group

    summary_matrix <- build_group_distance_summary_matrix(contrast_group, dist_matrix)
    observed_ratio <- compute_between_within_ratio(summary_matrix)
    observed_energy <- compute_energy_distance(
        build_group_distance_summary_matrix(contrast_group, dist_matrix^2)
    )

    shuffled_statistics <- NULL
    if (identical(nope, "shuffle") || identical(energy_nope, "shuffle")) {
        shuffled_statistics <- compute_shuffled_statistics(contrast_group, dist_matrix, n_nope_shuffles)
    }

    list(
        effectSize = observed_ratio,
        energyDistance = observed_energy,
        ratioNopeHigh = resolve_nope_edge(nope, shuffled_statistics, "ratio", 1, nope_quantile),
        energyNopeHigh = resolve_nope_edge(energy_nope, shuffled_statistics, "energy", 0, nope_quantile)
    )
}

#' Upper NOPE edge for one statistic in one posterior draw.
#'
#' @param nope_setting `"shuffle"`, `"none"`, or a fixed numeric edge.
#' @param shuffled_statistics Matrix from compute_shuffled_statistics, or
#'   `NULL` when no shuffles ran.
#' @param statistic_column `"ratio"` or `"energy"`.
#' @param floor_value The statistic's no-difference value (1 for \eqn{R}, 0 for
#'   \eqn{E}); the NOPE edge never falls below it.
#' @param nope_quantile Quantile of the shuffled statistic used as the edge.
#' @return Scalar edge, or `NA` when `nope_setting` is `"none"`.
#' @keywords internal
resolve_nope_edge <- function(nope_setting, shuffled_statistics, statistic_column, floor_value, nope_quantile) {
    if (identical(nope_setting, "shuffle")) {
        shuffled_values <- shuffled_statistics[, statistic_column]
        return(max(floor_value, stats::quantile(shuffled_values, probs = nope_quantile, names = FALSE)))
    }
    if (is.numeric(nope_setting)) {
        return(nope_setting)
    }
    NA_real_
}

#' Execute Bayes enrichment from Dirichlet posterior draws of R and E.
#'
#' For each contrast, returns both a one-row posterior summary and the long
#' draw table used by half-eye plots. Draw storage scales as
#' `n_posterior * n_contrasts` effect sizes, including omnibus, so the
#' posterior of the global ratio can be inspected alongside pairwise draws.
#'
#' `R_PPGTnope` compares each draw's \eqn{R} with that same draw's \eqn{R} edge,
#' and `E_PPGTnope` does the same for \eqn{E}; `R_nopeHigh` and `E_nopeHigh` in
#' the summary are the medians of the per-draw edges. Shuffles consume the
#' random stream between Dirichlet draws, so a fixed seed gives different draws
#' when `nope` or `energy_nope` is `"shuffle"` than when neither is.
#' @keywords internal
run_bayes_enrichment <- function(composition,
                                 contrast_specs,
                                 n_posterior,
                                 n_cores,
                                 prior_pseudocounts,
                                 credible_interval = "hpdi",
                                 nope = "shuffle",
                                 energy_nope = "shuffle",
                                 n_nope_shuffles = 200,
                                 nope_quantile = 0.95) {
    counts <- composition$counts
    group <- composition$group
    contrast_names <- names(contrast_specs)
    use_ratio_nope <- !identical(nope, "none")
    use_energy_nope <- !identical(energy_nope, "none")
    any_shuffle_nope <- identical(nope, "shuffle") || identical(energy_nope, "shuffle")

    # Group sizes stay fixed across draws, so the main process decides once
    # here whether each contrast has enough shuffles for a NOPE, where
    # forked workers cannot swallow warnings. R and E depend only on which
    # subjects share a group, so one check covers both. A fixed edge needs no
    # shuffles.
    enough_shuffles <- stats::setNames(rep(TRUE, length(contrast_names)), contrast_names)
    if (any_shuffle_nope) {
        for (contrast_name in contrast_names) {
            contrast_group <- prepare_contrast_group(group, contrast_specs[[contrast_name]])
            enough_shuffles[[contrast_name]] <- has_enough_shuffles_for_nope(contrast_group, nope_quantile)
        }
        if (!all(enough_shuffles)) {
            warning(
                sprintf(
                    paste0(
                        "groups too small for a NOPE in %s: too few distinct label shuffles, ",
                        "so the NOPE probabilities are NA and enoughShufflesForNope is FALSE ",
                        "(the result still includes PPGT1); a fixed edge (nope = <number above 1>, ",
                        "energyNope = <number above 0>) works at any group size"
                    ),
                    paste(contrast_names[!enough_shuffles], collapse = ", ")
                ),
                call. = FALSE
            )
        }
    }

    # Share of draws whose statistic clears that draw's edge. A NOPE on
    # groups too small to shuffle carries no information, so it reports NA.
    summarize_nope_probability <- function(statistic_values, nope_highs, nope_setting, contrast_name) {
        if (identical(nope_setting, "shuffle") && !enough_shuffles[[contrast_name]]) {
            return(NA_real_)
        }
        mean(statistic_values > nope_highs)
    }

    run_contrast_bayes <- function(contrast_name) {
        spec <- contrast_specs[[contrast_name]]
        draw_indices <- seq_len(n_posterior)
        compute_draw <- function(draw_idx) {
            run_bayes_draw(
                counts,
                group,
                spec,
                prior_pseudocounts,
                nope = nope,
                energy_nope = energy_nope,
                n_nope_shuffles = n_nope_shuffles,
                nope_quantile = nope_quantile
            )
        }

        draw_results <- if (n_cores > 1) {
            parallel::mclapply(draw_indices, compute_draw, mc.cores = min(n_cores, n_posterior))
        } else {
            lapply(draw_indices, compute_draw)
        }

        effect_sizes <- vapply(draw_results, function(x) x$effectSize, numeric(1))
        energy_distances <- vapply(draw_results, function(x) x$energyDistance, numeric(1))
        ci_bounds <- compute_posterior_interval(
            effect_sizes,
            prob = 0.95,
            type = credible_interval
        )
        energy_ci_bounds <- compute_posterior_interval(
            energy_distances,
            prob = 0.95,
            type = credible_interval
        )

        summary_cols <- list(
            contrastId = spec$contrastId,
            effectSize = mean(effect_sizes),
            PPGT1 = mean(effect_sizes > 1)
        )
        draw_cols <- list(
            contrastId = spec$contrastId,
            draw = draw_indices,
            effectSize = effect_sizes,
            energyDistance = energy_distances
        )
        if (use_ratio_nope) {
            ratio_nope_highs <- vapply(draw_results, function(x) x$ratioNopeHigh, numeric(1))
            summary_cols$R_PPGTnope <- summarize_nope_probability(effect_sizes, ratio_nope_highs, nope, contrast_name)
            summary_cols$R_nopeHigh <- stats::median(ratio_nope_highs)
            draw_cols$R_nopeHigh <- ratio_nope_highs
        }
        if (use_energy_nope) {
            energy_nope_highs <- vapply(draw_results, function(x) x$energyNopeHigh, numeric(1))
            summary_cols$E_PPGTnope <- summarize_nope_probability(
                energy_distances,
                energy_nope_highs,
                energy_nope,
                contrast_name
            )
            summary_cols$E_nopeHigh <- stats::median(energy_nope_highs)
            draw_cols$E_nopeHigh <- energy_nope_highs
        }
        if (any_shuffle_nope) {
            summary_cols$enoughShufflesForNope <- enough_shuffles[[contrast_name]]
        }
        summary_cols$effectCiLow <- ci_bounds[1]
        summary_cols$effectCiHigh <- ci_bounds[2]
        summary_cols$energyDistance <- mean(energy_distances)
        summary_cols$energyCiLow <- energy_ci_bounds[1]
        summary_cols$energyCiHigh <- energy_ci_bounds[2]

        summary_row <- as.data.frame(summary_cols, stringsAsFactors = FALSE)

        # Keep per-draw effect sizes for PlotCompositionContrasts half-eyes;
        # summarizing to CI alone would force geom_interval and lose the slab.
        draws_df <- as.data.frame(draw_cols, stringsAsFactors = FALSE)

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

    omnibus <- list(
        effectSize = omnibus_row$effectSize[1],
        PPGT1 = omnibus_row$PPGT1[1]
    )
    nope_fields <- c("R_PPGTnope", "R_nopeHigh", "E_PPGTnope", "E_nopeHigh", "enoughShufflesForNope")
    for (nope_field in intersect(nope_fields, colnames(omnibus_row))) {
        omnibus[[nope_field]] <- omnibus_row[[nope_field]][1]
    }
    omnibus$effectCiLow <- omnibus_row$effectCiLow[1]
    omnibus$effectCiHigh <- omnibus_row$effectCiHigh[1]
    omnibus$energyDistance <- omnibus_row$energyDistance[1]
    omnibus$energyCiLow <- omnibus_row$energyCiLow[1]
    omnibus$energyCiHigh <- omnibus_row$energyCiHigh[1]
    omnibus$draws <- omnibus_draws

    list(
        omnibus = omnibus,
        contrasts = contrast_only,
        draws = draws_df,
        diagnostics = list(nPosterior = n_posterior)
    )
}
