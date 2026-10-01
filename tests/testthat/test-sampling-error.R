library(testthat)
library(HellingerDistanceEnrichment)

fake_permutation_result <- function(p_value,
                                    n_permutations,
                                    contrast_p = 0.4) {
    structure(
        list(
            omnibus = list(effectSize = 1.2, pValue = p_value),
            contrasts = data.frame(
                contrastId = "Control_vs_Treatment",
                effectSize = 1.1,
                pValue = contrast_p,
                pAdj = contrast_p,
                stringsAsFactors = FALSE
            ),
            method = "permutation",
            settings = list(
                nPermutations = n_permutations,
                nPosterior = 100
            )
        ),
        class = "HellingerEnrichmentResult"
    )
}

fake_bayes_result <- function(ppgt1, n_posterior, contrast_ppgt1 = 0.6) {
    structure(
        list(
            omnibus = list(
                effectSize = 1.4,
                PPGT1 = ppgt1,
                effectCiLow = 1.1,
                effectCiHigh = 1.7
            ),
            contrasts = data.frame(
                contrastId = "Control_vs_Treatment",
                effectSize = 1.3,
                PPGT1 = contrast_ppgt1,
                effectCiLow = 0.9,
                effectCiHigh = 1.6,
                stringsAsFactors = FALSE
            ),
            method = "bayes",
            settings = list(
                nPermutations = 1000,
                nPosterior = n_posterior
            )
        ),
        class = "HellingerEnrichmentResult"
    )
}

test_that("permutation trial count is B + 1 and floor is 1 / (B + 1)", {
    n_permutations <- 10L
    p_floor <- 1 / (n_permutations + 1)
    result <- fake_permutation_result(p_floor, n_permutations)
    summary_row <- SummarizeSamplingError(result)

    expect_equal(summary_row$nTrials, n_permutations + 1L)
    expect_equal(summary_row$floor, p_floor)
    expect_true(summary_row$onFloor)
    expect_equal(summary_row$nSuccesses, 1L)
    expect_equal(summary_row$cutoff, 0.05)
    expect_equal(summary_row$contrastId, "omnibus")
    expect_equal(summary_row$method, "permutation")
})

test_that("Bayes trial count is nPosterior and floor is NA", {
    result <- fake_bayes_result(ppgt1 = 1, n_posterior = 10L)
    summary_row <- SummarizeSamplingError(result)

    expect_equal(summary_row$nTrials, 10L)
    expect_equal(summary_row$nSuccesses, 10L)
    expect_true(is.na(summary_row$floor))
    expect_true(is.na(summary_row$onFloor))
    expect_equal(summary_row$cutoff, 0.95)
    expect_equal(summary_row$method, "bayes")
})

test_that("decisionStable is FALSE when the Wilson interval contains the cutoff", {
    perm_floor <- SummarizeSamplingError(
        fake_permutation_result(1 / 11, 10L)
    )
    expect_false(perm_floor$decisionStable)
    expect_true(perm_floor$intervalLow <= 0.05)
    expect_true(perm_floor$intervalHigh >= 0.05)
    expect_false(perm_floor$samplingSufficient)

    bayes_small <- SummarizeSamplingError(
        fake_bayes_result(ppgt1 = 1, n_posterior = 10L)
    )
    expect_false(bayes_small$decisionStable)
    expect_true(bayes_small$intervalLow <= 0.95)
    expect_true(bayes_small$intervalHigh >= 0.95)
    expect_false(bayes_small$samplingSufficient)
})

test_that("samplingSufficient follows Wilson half-width", {
    precise <- SummarizeSamplingError(
        fake_bayes_result(ppgt1 = 1, n_posterior = 5000L)
    )
    expect_true(precise$decisionStable)
    expect_true(precise$samplingSufficient)
    expect_lte(precise$halfWidth, 0.02)
    expect_gt(precise$intervalLow, 0.95)

    coarse <- SummarizeSamplingError(
        fake_permutation_result(p_value = 0.5, n_permutations = 10L)
    )
    expect_false(coarse$samplingSufficient)
    expect_gt(coarse$halfWidth, 0.02)
})

test_that("contrastId selects the requested unadjusted p-value", {
    result <- fake_permutation_result(
        p_value = 0.2,
        n_permutations = 99L,
        contrast_p = 0.4
    )
    omnibus_row <- SummarizeSamplingError(result)
    contrast_row <- SummarizeSamplingError(
        result,
        contrastId = "Control_vs_Treatment"
    )

    expect_equal(omnibus_row$estimate, 0.2)
    expect_equal(omnibus_row$contrastId, "omnibus")
    expect_equal(contrast_row$estimate, 0.4)
    expect_equal(contrast_row$contrastId, "Control_vs_Treatment")
    expect_equal(contrast_row$nTrials, 100L)
})

test_that("SummarizeSamplingError rejects bad input", {
    expect_error(SummarizeSamplingError(list()), "HellingerEnrichmentResult")
    expect_error(
        SummarizeSamplingError(
            fake_permutation_result(0.2, 20L),
            contrastId = "missing_contrast"
        ),
        "not found"
    )
    expect_error(
        SummarizeSamplingError(
            fake_permutation_result(0.2, 20L),
            alpha = 0
        ),
        "alpha"
    )
    expect_error(
        SummarizeSamplingError(
            fake_bayes_result(1, 20L),
            ppgt1Cutoff = 1
        ),
        "ppgt1Cutoff"
    )
})

test_that("Wilson interval stays in [0, 1] at the boundaries", {
    zero <- HellingerDistanceEnrichment:::wilson_proportion_interval(0, 20)
    expect_gte(zero$intervalLow, 0)
    expect_lte(zero$intervalHigh, 1)
    expect_lt(zero$intervalLow, zero$intervalHigh)

    one <- HellingerDistanceEnrichment:::wilson_proportion_interval(20, 20)
    expect_gte(one$intervalLow, 0)
    expect_lte(one$intervalHigh, 1)
    expect_gt(one$intervalLow, 0.8)
})

test_that("real CompareGroupCompositions results have matching trial counts", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 4,
        n_categories = 3,
        seed = 21
    )
    perm_result <- CompareGroupCompositions(
        long_table,
        method = "permutation",
        nPermutations = 40,
        seed = 22
    )
    perm_summary <- SummarizeSamplingError(perm_result)
    expect_equal(perm_summary$nTrials, 41L)
    expect_equal(perm_summary$estimate, perm_result$omnibus$pValue)
    expect_gte(perm_summary$floor, 1 / 41)

    bayes_result <- CompareGroupCompositions(
        long_table,
        method = "bayes",
        nPosterior = 15,
        seed = 23
    )
    bayes_summary <- SummarizeSamplingError(bayes_result)
    expect_equal(bayes_summary$nTrials, 15L)
    expect_equal(bayes_summary$estimate, bayes_result$omnibus$PPGT1)
    expect_equal(bayes_summary$cutoff, 0.95)
})
