#' Wilson score interval for a binomial proportion.
#'
#' Centered on the score-test mean so bounds stay in (0, 1) when the observed
#' proportion is 0 or 1, then clipped to the unit interval.
#'
#' @param n_successes Number of successes.
#' @param n_trials Number of trials.
#' @param confidence Confidence level (default 0.95).
#' @return data.frame with `intervalLow` and `intervalHigh`.
#' @keywords internal
wilson_proportion_interval <- function(n_successes, n_trials, confidence = 0.95) {
    if (n_trials < 1) {
        stop("n_trials must be at least 1")
    }
    if (n_successes < 0 || n_successes > n_trials) {
        stop("n_successes must lie in [0, n_trials]")
    }
    if (confidence <= 0 || confidence >= 1) {
        stop("confidence must lie in (0, 1)")
    }

    z <- stats::qnorm((1 + confidence) / 2)
    p_hat <- n_successes / n_trials
    z2 <- z * z
    denom <- 1 + z2 / n_trials
    center <- (p_hat + z2 / (2 * n_trials)) / denom
    half_width <- z * sqrt(
        p_hat * (1 - p_hat) / n_trials + z2 / (4 * n_trials * n_trials)
    ) / denom

    data.frame(
        intervalLow = max(0, center - half_width),
        intervalHigh = min(1, center + half_width)
    )
}

#' Monte Carlo error for a permutation p-value or the posterior probability that R exceeds 1.
#'
#' A reported permutation p-value or posterior probability that the
#' Hellinger ratio \(R\) exceeds 1 is a binomial Monte Carlo proportion.
#' This helper does not draw additional permutations or posterior
#' samples; it diagnoses the result you already have.
#'
#' For permutation, the observed-inclusive p-value is a proportion of
#' `nPermutations + 1` trials, and the smallest reportable p is
#' \(1/(B+1)\). For Bayes, the posterior probability that \(R > 1\)
#' (stored as `PPGT1`) is a proportion of `nPosterior` Dirichlet draws
#' with \(R > 1\).
#'
#' Two flags answer different questions. `decisionStable` is TRUE when the
#' method-specific cutoff (0.05 for permutation p, 0.95 for the posterior
#' probability that \(R > 1\) by default) lies outside the Wilson interval.
#' `samplingSufficient` is TRUE when the Wilson half-width is at most
#' `maxHalfWidth`. A p-value sitting on the resolution floor at small \(B\)
#' can look significant while the interval still contains 0.05.
#'
#' Permutation uses the unadjusted p-value. Holm-adjusted `pAdj` is not a
#' binomial proportion of \(B+1\) trials.
#'
#' @param result A HellingerEnrichmentResult from CompareGroupCompositions.
#' @param alpha Decision cutoff for permutation p-values (default 0.05).
#' @param ppgt1Cutoff Decision cutoff for the posterior probability that
#'   \(R > 1\) (default 0.95).
#' @param maxHalfWidth Maximum Wilson half-width treated as sufficient
#'   (default 0.02).
#' @param contrastId Optional contrast identifier. Default is the omnibus
#'   result. Pass a `contrasts$contrastId` value to diagnose that row.
#' @return A one-row data.frame with the estimate, trial counts, resolution
#'   floor (NA for Bayes), Wilson interval, cutoff, `decisionStable`, and
#'   `samplingSufficient`.
#' @export
#' @examples
#' long_table <- data.frame(
#'     subjectId = rep(c("S1", "S2", "S3", "S4"), each = 2),
#'     category = rep(c("A", "B"), 4),
#'     group = rep(c("Control", "Control", "Treatment", "Treatment"), each = 2),
#'     n = c(20, 10, 18, 12, 8, 22, 9, 21),
#'     stringsAsFactors = FALSE
#' )
#' result <- CompareGroupCompositions(long_table, nPermutations = 50, seed = 1)
#' SummarizeSamplingError(result)
SummarizeSamplingError <- function(result,
                                   alpha = 0.05,
                                   ppgt1Cutoff = 0.95,
                                   maxHalfWidth = 0.02,
                                   contrastId = NULL) {
    if (!inherits(result, "HellingerEnrichmentResult")) {
        stop("result must be a HellingerEnrichmentResult")
    }
    if (!result$method %in% c("permutation", "bayes")) {
        stop("result$method must be \"permutation\" or \"bayes\"")
    }
    if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
        stop("alpha must be a single number in (0, 1)")
    }
    if (!is.numeric(ppgt1Cutoff) || length(ppgt1Cutoff) != 1 ||
        ppgt1Cutoff <= 0 || ppgt1Cutoff >= 1) {
        stop("ppgt1Cutoff must be a single number in (0, 1)")
    }
    if (!is.numeric(maxHalfWidth) || length(maxHalfWidth) != 1 ||
        maxHalfWidth <= 0) {
        stop("maxHalfWidth must be a single positive number")
    }

    use_omnibus <- is.null(contrastId) || identical(contrastId, "omnibus")
    if (use_omnibus) {
        selected_contrast <- "omnibus"
        if (identical(result$method, "permutation")) {
            estimate <- result$omnibus$pValue
        } else {
            estimate <- result$omnibus$PPGT1
        }
    } else {
        if (length(contrastId) != 1 || !is.character(contrastId)) {
            stop("contrastId must be a single character string")
        }
        contrast_row <- result$contrasts[
            result$contrasts$contrastId == contrastId,
            ,
            drop = FALSE
        ]
        if (nrow(contrast_row) != 1) {
            stop(sprintf(
                "contrastId '%s' not found; available: %s",
                contrastId,
                paste(c("omnibus", result$contrasts$contrastId), collapse = ", ")
            ))
        }
        selected_contrast <- contrastId
        if (identical(result$method, "permutation")) {
            estimate <- contrast_row$pValue[1]
        } else {
            estimate <- contrast_row$PPGT1[1]
        }
    }

    if (!is.numeric(estimate) || length(estimate) != 1 || is.na(estimate)) {
        stop("selected estimate is missing")
    }

    if (identical(result$method, "permutation")) {
        n_permutations <- result$settings$nPermutations
        n_trials <- as.integer(n_permutations) + 1L
        resolution_floor <- 1 / n_trials
        on_floor <- abs(estimate - resolution_floor) < 1e-12
        cutoff <- alpha
    } else {
        n_trials <- as.integer(result$settings$nPosterior)
        resolution_floor <- NA_real_
        on_floor <- NA
        cutoff <- ppgt1Cutoff
    }

    n_successes <- as.integer(round(estimate * n_trials))
    if (n_successes < 0L || n_successes > n_trials) {
        stop("could not recover binomial counts from the reported estimate")
    }

    interval <- wilson_proportion_interval(n_successes, n_trials)
    half_width <- (interval$intervalHigh - interval$intervalLow) / 2
    decision_stable <- cutoff < interval$intervalLow || cutoff > interval$intervalHigh
    sampling_sufficient <- half_width <= maxHalfWidth

    data.frame(
        method = result$method,
        contrastId = selected_contrast,
        estimate = estimate,
        nTrials = n_trials,
        nSuccesses = n_successes,
        floor = resolution_floor,
        onFloor = on_floor,
        intervalLow = interval$intervalLow,
        intervalHigh = interval$intervalHigh,
        halfWidth = half_width,
        cutoff = cutoff,
        decisionStable = decision_stable,
        samplingSufficient = sampling_sufficient,
        stringsAsFactors = FALSE
    )
}
