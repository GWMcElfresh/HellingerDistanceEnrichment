library(testthat)
library(HellingerDistanceEnrichment)

test_that("long table coerces to CategoryComposition and round-trips", {
    fixture_path <- testthat::test_path("fixtures", "synthetic_long_counts.csv")
    long_table <- read.csv(fixture_path, stringsAsFactors = FALSE)

    composition <- CompareGroupCompositions(long_table, nPermutations = 10)
    expect_s3_class(composition, "HellingerEnrichmentResult")

    comp <- HellingerDistanceEnrichment:::coerce_long_table_to_composition(long_table)
    expect_s3_class(comp, "CategoryComposition")
    expect_equal(nrow(comp$counts), 4)
    expect_equal(ncol(comp$counts), 4)
    expect_equal(comp$counts["Subj01", "Category0"], 25)
    expect_equal(comp$counts["Subj01", "Category2"], 20)
})

test_that("missing subject-category cells are zero-filled", {
    long_table <- data.frame(
        subjectId = c("A", "A"),
        category = c("X", "Y"),
        group = c("G1", "G1"),
        n = c(5, 3),
        stringsAsFactors = FALSE
    )
    comp <- HellingerDistanceEnrichment:::coerce_long_table_to_composition(long_table)
    expect_equal(comp$counts["A", "X"], 5)
    expect_equal(comp$counts["A", "Y"], 3)
})

test_that("Hellinger distance matches known toy vectors", {
    freqs1 <- c(0.5, 0.5)
    freqs2 <- c(1, 0)
    distance <- HellingerDistanceEnrichment:::compute_hellinger_distance_freqs(freqs1, freqs2)
    diff_sqrt <- sqrt(freqs1) - sqrt(freqs2)
    expected <- sqrt(sum(diff_sqrt^2)) / sqrt(2)
    expect_equal(distance, expected, tolerance = 1e-10)

    identical_freqs <- HellingerDistanceEnrichment:::compute_hellinger_distance_freqs(freqs1, freqs1)
    expect_equal(identical_freqs, 0)
})

test_that("permutation null does not systematically yield tiny p under exchangeability", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 4,
        n_categories = 3,
        groups = c("A", "B"),
        seed = 42
    )

    set.seed(99)
    p_values <- replicate(20, {
        result <- CompareGroupCompositions(
            long_table,
            method = "permutation",
            nPermutations = 50,
            seed = sample.int(1e6, 1)
        )
        result$omnibus$pValue
    })

    expect_gt(median(p_values), 0.05)
})

test_that("planted structure yields enrichment signal", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 8,
        n_categories = 4,
        seed = 7
    )
    structured <- HellingerDistanceEnrichment:::plant_group_structure(
        long_table,
        target_group = "Treatment",
        target_category = "Category0",
        boost = 12
    )

    result <- CompareGroupCompositions(
        structured,
        method = "permutation",
        nPermutations = 300,
        seed = 123
    )

    expect_lt(result$omnibus$pValue, 0.15)
    expect_gt(result$omnibus$effectSize, 1)
})

test_that("three-group design plants one arm without disturbing the other", {
    groups <- c("Control", "Condition1", "Condition2")
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 8,
        n_categories = 4,
        groups = groups,
        seed = 21
    )
    structured <- HellingerDistanceEnrichment:::plant_group_structure(
        long_table,
        target_group = "Condition1",
        target_category = "Category0",
        boost = 12
    )

    result <- CompareGroupCompositions(
        structured,
        method = "permutation",
        nPermutations = 300,
        seed = 88
    )

    contrasts <- result$contrasts
    planted <- contrasts$contrastId == "Condition1_vs_Control"
    null_arm <- contrasts$contrastId == "Condition2_vs_Control"

    expect_lt(contrasts$pValue[planted], 0.15)
    expect_gt(contrasts$pValue[null_arm], 0.2)
    expect_lt(result$omnibus$pValue, 0.2)
})

test_that("plant_group_structure rejects unknown group or category", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(seed = 1)

    expect_error(
        HellingerDistanceEnrichment:::plant_group_structure(
            long_table,
            target_group = "MissingGroup"
        ),
        "target_group"
    )
    expect_error(
        HellingerDistanceEnrichment:::plant_group_structure(
            long_table,
            target_category = "MissingCategory"
        ),
        "target_category"
    )
    expect_error(
        HellingerDistanceEnrichment:::plant_group_structure(long_table, boost = 0),
        "boost must be positive"
    )
})

test_that("Bayes smoke test returns finite PPGT1", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 5,
        n_categories = 4,
        seed = 11
    )

    result <- CompareGroupCompositions(
        long_table,
        method = "bayes",
        nPosterior = 5,
        nCores = 1,
        seed = 55
    )

    expect_s3_class(result, "HellingerEnrichmentResult")
    expect_true(is.finite(result$omnibus$PPGT1))
    expect_true(result$omnibus$PPGT1 >= 0 && result$omnibus$PPGT1 <= 1)
    expect_true(all(is.finite(result$contrasts$PPGT1)))
    expect_true(all(result$contrasts$PPGT1 >= 0 & result$contrasts$PPGT1 <= 1))
    expect_false("pValue" %in% names(result$contrasts))
    expect_false("pAdj" %in% names(result$contrasts))
    expect_null(result$omnibus$pValue)
    expect_true(all(is.finite(result$contrasts$effectCiLow)))
    expect_true(!is.null(result$draws))
    expect_gt(nrow(result$draws), 0)
    expect_true(all(c("contrastId", "draw", "effectSize") %in% colnames(result$draws)))
    expect_true("omnibus" %in% result$draws$contrastId)
    expect_equal(result$settings$credibleInterval, "hpdi")
    expect_equal(length(result$omnibus$draws), 5)
    expect_equal(
        result$omnibus$draws,
        result$draws$effectSize[result$draws$contrastId == "omnibus"]
    )
})

test_that("compute_posterior_interval distinguishes HPDI from CPI on an asymmetric sample", {
    x <- c(rep(0, 90), 10:19)
    cpi <- HellingerDistanceEnrichment:::compute_posterior_interval(x, prob = 0.95, type = "cpi")
    hpdi <- HellingerDistanceEnrichment:::compute_posterior_interval(x, prob = 0.95, type = "hpdi")

    expect_equal(
        cpi,
        as.numeric(stats::quantile(x, probs = c(0.025, 0.975), names = FALSE))
    )
    expect_false(isTRUE(all.equal(hpdi, cpi)))
    expect_lt(diff(hpdi), diff(cpi))
})

test_that("Bayes CPI option is recorded and returns finite intervals", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 4,
        n_categories = 3,
        seed = 13
    )

    result <- CompareGroupCompositions(
        long_table,
        method = "bayes",
        nPosterior = 8,
        nCores = 1,
        seed = 2,
        credibleInterval = "cpi"
    )

    expect_equal(result$settings$credibleInterval, "cpi")
    expect_true(all(is.finite(result$contrasts$effectCiLow)))
    expect_true(all(is.finite(result$contrasts$effectCiHigh)))
    expect_true(is.finite(result$omnibus$effectCiLow))
    expect_true(is.finite(result$omnibus$effectCiHigh))
})

test_that("Bayes warns when unused permutation arguments are non-default", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 4,
        n_categories = 3,
        seed = 12
    )

    expect_warning(
        CompareGroupCompositions(
            long_table,
            method = "bayes",
            nPosterior = 3,
            nPermutations = 20,
            nCores = 1,
            seed = 1
        ),
        "nPermutations"
    )
})

test_that("Bayes NOPE adds R_PPGTnope bounded above by PPGT1", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 5,
        n_categories = 4,
        seed = 11
    )

    result <- CompareGroupCompositions(
        long_table,
        method = "bayes",
        nPosterior = 10,
        nNopeShuffles = 30,
        seed = 55
    )

    expect_equal(result$settings$nope, "shuffle")
    expect_true(all(c("R_PPGTnope", "R_nopeHigh", "R_hdiAboveNope", "enoughShufflesForNope") %in% colnames(result$contrasts)))
    expect_true(all(result$contrasts$enoughShufflesForNope))
    expect_true(result$omnibus$enoughShufflesForNope)
    expect_true("R_nopeHigh" %in% colnames(result$draws))
    expect_true(all(result$draws$R_nopeHigh >= 1))
    expect_true(all(result$contrasts$R_PPGTnope <= result$contrasts$PPGT1))
    expect_lte(result$omnibus$R_PPGTnope, result$omnibus$PPGT1)
    expect_gte(result$omnibus$R_nopeHigh, 1)
    expect_type(result$omnibus$R_hdiAboveNope, "logical")
    expect_equal(result$omnibus$R_hdiAboveNope, result$omnibus$effectCiLow > result$omnibus$R_nopeHigh)
    expect_equal(result$contrasts$R_hdiAboveNope, result$contrasts$effectCiLow > result$contrasts$R_nopeHigh)
})

test_that("nope and energyNope = 'none' drop their NOPE columns", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 4,
        n_categories = 3,
        seed = 13
    )

    ratio_none <- CompareGroupCompositions(long_table, method = "bayes", nPosterior = 5, nope = "none", seed = 3)
    expect_false(any(c("R_PPGTnope", "R_nopeHigh", "R_hdiAboveNope") %in% colnames(ratio_none$contrasts)))
    expect_true(all(c("E_PPGTnope", "E_nopeHigh", "E_hdiAboveNope", "enoughShufflesForNope") %in% colnames(ratio_none$contrasts)))

    both_none <- CompareGroupCompositions(
        long_table,
        method = "bayes",
        nPosterior = 5,
        nope = "none",
        energyNope = "none",
        seed = 3
    )
    nope_columns <- c(
        "R_PPGTnope", "R_nopeHigh", "R_hdiAboveNope",
        "E_PPGTnope", "E_nopeHigh", "E_hdiAboveNope",
        "enoughShufflesForNope"
    )
    expect_true("PPGT1" %in% colnames(both_none$contrasts))
    expect_false(any(nope_columns %in% colnames(both_none$contrasts)))
    expect_false(any(c("R_nopeHigh", "E_nopeHigh") %in% colnames(both_none$draws)))
    expect_null(both_none$omnibus$R_PPGTnope)
    expect_null(both_none$omnibus$E_PPGTnope)
    expect_output(print(both_none), "PPGT1")
})

test_that("the E NOPE reports E_PPGTnope bounded by the share of draws with E above 0", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 5,
        n_categories = 4,
        groups = c("A", "B", "C"),
        seed = 11
    )

    result <- CompareGroupCompositions(long_table, method = "bayes", nPosterior = 12, nNopeShuffles = 30, seed = 55)

    expect_true(all(c("E_PPGTnope", "E_nopeHigh", "E_hdiAboveNope") %in% colnames(result$contrasts)))
    expect_true(all(result$draws$E_nopeHigh >= 0))
    expect_equal(result$contrasts$E_hdiAboveNope, result$contrasts$energyCiLow > result$contrasts$E_nopeHigh)
    for (contrast_id in result$contrasts$contrastId) {
        contrast_draws <- result$draws[result$draws$contrastId == contrast_id, ]
        contrast_row <- result$contrasts[result$contrasts$contrastId == contrast_id, ]
        expect_lte(contrast_row$E_PPGTnope, mean(contrast_draws$energyDistance > 0))
        expect_equal(contrast_row$E_PPGTnope, mean(contrast_draws$energyDistance > contrast_draws$E_nopeHigh))
    }
    expect_output(print(result), "E_PPGTnope")
})

test_that("a fixed energyNope edge applies to every row and draw, alone or mixed", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 3,
        n_categories = 4,
        seed = 7
    )

    # 3 versus 3 is too small for a NOPE, but fixed edges on both scales
    # still apply and draw no random numbers.
    fixed_both <- CompareGroupCompositions(
        long_table,
        method = "bayes",
        nPosterior = 12,
        nope = 1.2,
        energyNope = 0.02,
        seed = 9
    )
    none_both <- CompareGroupCompositions(
        long_table,
        method = "bayes",
        nPosterior = 12,
        nope = "none",
        energyNope = "none",
        seed = 9
    )
    expect_true(all(fixed_both$contrasts$E_nopeHigh == 0.02))
    expect_true(all(fixed_both$draws$E_nopeHigh == 0.02))
    expect_true(all(is.finite(fixed_both$contrasts$E_PPGTnope)))
    expect_false("enoughShufflesForNope" %in% colnames(fixed_both$contrasts))
    expect_identical(fixed_both$draws$energyDistance, none_both$draws$energyDistance)

    # Mixed: a fixed R edge with an E NOPE still shuffles, so at 3 versus 3
    # only the NOPE side reports NA.
    expect_warning(
        mixed <- CompareGroupCompositions(
            long_table,
            method = "bayes",
            nPosterior = 12,
            nope = 1.2,
            energyNope = "shuffle",
            nNopeShuffles = 20,
            seed = 9
        ),
        "too small for a NOPE"
    )
    expect_true(all(is.finite(mixed$contrasts$R_PPGTnope)))
    expect_true(all(is.na(mixed$contrasts$E_PPGTnope)))
    expect_equal(mixed$contrasts$R_hdiAboveNope, mixed$contrasts$effectCiLow > 1.2)
    expect_true(all(is.na(mixed$contrasts$E_hdiAboveNope)))
    expect_false(any(mixed$contrasts$enoughShufflesForNope))
})

test_that("deeply sampled null cohorts do not clear the NOPE", {
    # At 2000 cells per subject the Dirichlet posterior collapses onto the
    # sample R, so PPGT1 tracks the sign of chance noise around 1. The NOPE
    # edge from label shuffling absorbs that noise.
    simulate_null_cohort <- function(n_per_group, n_categories, cells) {
        rows <- list()
        for (group_label in c("A", "B")) {
            for (subject_idx in seq_len(n_per_group)) {
                gamma_draws <- stats::rgamma(n_categories, shape = 5)
                subject_counts <- stats::rmultinom(1, cells, gamma_draws / sum(gamma_draws))
                rows[[length(rows) + 1]] <- data.frame(
                    subjectId = paste0(group_label, subject_idx),
                    category = paste0("Category", seq_len(n_categories)),
                    group = group_label,
                    n = as.vector(subject_counts),
                    stringsAsFactors = FALSE
                )
            }
        }
        do.call(rbind, rows)
    }

    set.seed(2026)
    omnibus_summaries <- t(replicate(4, {
        null_table <- simulate_null_cohort(n_per_group = 8, n_categories = 10, cells = 2000)
        result <- CompareGroupCompositions(
            null_table,
            method = "bayes",
            nPosterior = 20,
            nNopeShuffles = 50,
            contrasts = c("A", "B"),
            seed = sample.int(1e6, 1)
        )
        c(PPGT1 = result$omnibus$PPGT1, R_PPGTnope = result$omnibus$R_PPGTnope)
    }))

    expect_true(all(omnibus_summaries[, "R_PPGTnope"] < 0.95))
    expect_true(all(omnibus_summaries[, "R_PPGTnope"] <= omnibus_summaries[, "PPGT1"]))
})

test_that("planted enrichment still clears the NOPE", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 8,
        n_categories = 4,
        seed = 7
    )
    structured <- HellingerDistanceEnrichment:::plant_group_structure(
        long_table,
        target_group = "Treatment",
        target_category = "Category0",
        boost = 12
    )

    result <- CompareGroupCompositions(
        structured,
        method = "bayes",
        nPosterior = 20,
        nNopeShuffles = 50,
        seed = 123
    )

    expect_gt(result$omnibus$R_PPGTnope, 0.9)
})

test_that("energy distance matches hand calculation and averages over group pairs", {
    two_group <- matrix(
        c(0.10, 0.30,
          0.30, 0.20),
        nrow = 2,
        dimnames = list(c("A", "B"), c("A", "B"))
    )
    expect_equal(HellingerDistanceEnrichment:::compute_energy_distance(two_group), 2 * 0.30 - 0.10 - 0.20)

    three_group <- matrix(
        c(0.10, 0.30, 0.25,
          0.30, 0.20, 0.40,
          0.25, 0.40, 0.15),
        nrow = 3,
        dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
    )
    pair_energies <- c(
        2 * 0.30 - 0.10 - 0.20,
        2 * 0.25 - 0.10 - 0.15,
        2 * 0.40 - 0.20 - 0.15
    )
    expect_equal(HellingerDistanceEnrichment:::compute_energy_distance(three_group), mean(pair_energies))
})

test_that("energy distance uses squared Hellinger distances", {
    omnibus_spec <- list(contrastId = "omnibus", groups = c("A", "B"), collapse = NULL)

    # Subjects identical within each group: E reduces to 2 * (1 - BC) between
    # the two group compositions, which plain Hellinger distances would miss.
    p_a <- c(40, 30, 20, 10)
    p_b <- c(10, 20, 30, 40)
    identical_counts <- rbind(p_a, p_a, p_a, p_b, p_b, p_b)
    rownames(identical_counts) <- paste0("S", 1:6)
    identical_group <- factor(rep(c("A", "B"), each = 3))
    names(identical_group) <- rownames(identical_counts)

    identical_result <- HellingerDistanceEnrichment:::compute_observed_contrast(
        identical_counts, identical_group, omnibus_spec, prior_pseudocounts = 0
    )
    bhattacharyya <- sum(sqrt((p_a / sum(p_a)) * (p_b / sum(p_b))))
    expect_equal(identical_result$energyDistance, 2 * (1 - bhattacharyya))

    # Random compositions: E equals the squared distance between square-root
    # centroids minus each group's summed unbiased variance over its size.
    set.seed(31)
    random_counts <- matrix(stats::rpois(7 * 5, lambda = 30) + 1, nrow = 7)
    rownames(random_counts) <- paste0("S", 1:7)
    random_group <- factor(c("A", "A", "A", "B", "B", "B", "B"))
    names(random_group) <- rownames(random_counts)

    sqrt_freqs <- sqrt(random_counts / rowSums(random_counts))
    sqrt_a <- sqrt_freqs[random_group == "A", , drop = FALSE]
    sqrt_b <- sqrt_freqs[random_group == "B", , drop = FALSE]
    centroid_term <- sum((colMeans(sqrt_a) - colMeans(sqrt_b))^2)
    spread_term <- sum(apply(sqrt_a, 2, stats::var)) / nrow(sqrt_a) +
        sum(apply(sqrt_b, 2, stats::var)) / nrow(sqrt_b)

    random_result <- HellingerDistanceEnrichment:::compute_observed_contrast(
        random_counts, random_group, omnibus_spec, prior_pseudocounts = 0
    )
    expect_equal(random_result$energyDistance, centroid_term - spread_term)
})

test_that("both methods report energy distance beside the ratio", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 5,
        n_categories = 4,
        seed = 11
    )

    perm_result <- CompareGroupCompositions(long_table, method = "permutation", nPermutations = 20, seed = 5)
    expect_true(is.finite(perm_result$omnibus$energyDistance))
    expect_true(all(is.finite(perm_result$contrasts$energyDistance)))
    expect_output(print(perm_result), "energyDistance")

    bayes_result <- CompareGroupCompositions(
        long_table,
        method = "bayes",
        nPosterior = 12,
        nNopeShuffles = 20,
        seed = 5
    )
    expect_true(all(c("energyDistance", "energyCiLow", "energyCiHigh") %in% colnames(bayes_result$contrasts)))
    expect_true(all(is.finite(bayes_result$contrasts$energyCiLow)))
    expect_true(is.finite(bayes_result$omnibus$energyDistance))
    expect_true("energyDistance" %in% colnames(bayes_result$draws))
})

test_that("has_enough_shuffles_for_nope needs more than 20 distinct splits at 0.95", {
    enough <- function(sizes, nope_quantile = 0.95) {
        group <- factor(rep(paste0("G", seq_along(sizes)), sizes))
        HellingerDistanceEnrichment:::has_enough_shuffles_for_nope(group, nope_quantile)
    }

    # Distinct splits: 2v2 = 3, 3v3 = 10, 2v4 = 15, 2v5 = 21, 4v4 = 35, 3v3v3 = 280.
    expect_false(enough(c(2, 2)))
    expect_false(enough(c(3, 3)))
    expect_false(enough(c(2, 4)))
    expect_true(enough(c(2, 5)))
    expect_true(enough(c(4, 4)))
    expect_true(enough(c(3, 3, 3)))

    # A split count exactly at the threshold does not pass.
    expect_false(enough(c(4, 4), nope_quantile = 1 - 1 / 35))
})

test_that("groups too small for the NOPE report NA with a warning", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 3,
        n_categories = 4,
        seed = 7
    )
    structured <- HellingerDistanceEnrichment:::plant_group_structure(
        long_table,
        target_group = "Treatment",
        target_category = "Category0",
        boost = 12
    )

    expect_warning(
        result <- CompareGroupCompositions(
            structured,
            method = "bayes",
            nPosterior = 10,
            nNopeShuffles = 30,
            seed = 4
        ),
        "too small for a NOPE"
    )

    expect_gt(result$omnibus$PPGT1, 0.9)
    expect_true(is.na(result$omnibus$R_PPGTnope))
    expect_true(is.na(result$omnibus$R_hdiAboveNope))
    expect_true(all(is.na(result$contrasts$R_PPGTnope)))
    expect_true(all(is.na(result$contrasts$R_hdiAboveNope)))
    expect_true(all(is.finite(result$contrasts$R_nopeHigh)))
    expect_false(result$omnibus$enoughShufflesForNope)
    expect_false(any(result$contrasts$enoughShufflesForNope))
    expect_output(print(result), "too few distinct label shuffles")
})

test_that("a numeric nope sets a fixed edge without consuming random numbers", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 3,
        n_categories = 4,
        seed = 7
    )

    fixed_result <- CompareGroupCompositions(
        long_table, method = "bayes", nPosterior = 12, nope = 1.2, energyNope = "none", seed = 9
    )
    none_result <- CompareGroupCompositions(
        long_table, method = "bayes", nPosterior = 12, nope = "none", energyNope = "none", seed = 9
    )

    expect_equal(fixed_result$settings$nope, 1.2)
    expect_true(all(fixed_result$contrasts$R_nopeHigh == 1.2))
    expect_true(all(fixed_result$draws$R_nopeHigh == 1.2))
    expect_false("enoughShufflesForNope" %in% colnames(fixed_result$contrasts))
    # 3 versus 3 is too small for the NOPE, but a fixed edge still applies.
    expect_true(all(is.finite(fixed_result$contrasts$R_PPGTnope)))
    expect_true(all(fixed_result$contrasts$R_PPGTnope <= fixed_result$contrasts$PPGT1))
    expect_equal(fixed_result$contrasts$R_hdiAboveNope, fixed_result$contrasts$effectCiLow > 1.2)
    expect_output(print(fixed_result), "R_hdiAboveNope")
    expect_identical(fixed_result$draws$effectSize, none_result$draws$effectSize)
    expect_identical(fixed_result$contrasts$PPGT1, none_result$contrasts$PPGT1)
    expect_output(print(fixed_result), "fixed NOPE bound")
})

test_that("invalid nope values stop", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 4,
        n_categories = 3,
        seed = 12
    )

    for (bad_nope in list(0.9, 1, c(1, 2), NA_real_, "foo", Inf)) {
        expect_error(
            CompareGroupCompositions(long_table, method = "bayes", nPosterior = 3, nope = bad_nope, seed = 1),
            "nope must be"
        )
    }
    for (bad_energy_nope in list(0, -0.1, c(0.01, 0.02), NA_real_, "foo", Inf)) {
        expect_error(
            CompareGroupCompositions(
                long_table, method = "bayes", nPosterior = 3, energyNope = bad_energy_nope, seed = 1
            ),
            "energyNope must be"
        )
    }
    expect_warning(
        CompareGroupCompositions(long_table, method = "permutation", nPermutations = 10, energyNope = 0.02, seed = 1),
        "energyNope"
    )
    expect_warning(
        CompareGroupCompositions(long_table, method = "permutation", nPermutations = 10, nope = 1.2, seed = 1),
        "nope unused"
    )
})

test_that("PlotCompositionContrasts shades the NOPE behind Bayes half-eyes", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 5,
        n_categories = 4,
        groups = c("A", "B", "C"),
        seed = 3
    )
    rect_layer_count <- function(plot) {
        sum(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomRect"), logical(1)))
    }

    shuffle_result <- CompareGroupCompositions(long_table, method = "bayes", nPosterior = 10, nNopeShuffles = 20, seed = 4)
    shuffle_plot <- PlotCompositionContrasts(shuffle_result)
    expect_silent(ggplot2::ggplot_build(shuffle_plot))
    expect_equal(rect_layer_count(shuffle_plot), 1)

    fixed_plot <- PlotCompositionContrasts(
        CompareGroupCompositions(long_table, method = "bayes", nPosterior = 10, nope = 1.25, seed = 4)
    )
    expect_silent(ggplot2::ggplot_build(fixed_plot))
    expect_equal(fixed_plot$layers[[1]]$data$xmax[1], 1.25)

    none_plot <- PlotCompositionContrasts(
        CompareGroupCompositions(long_table, method = "bayes", nPosterior = 10, nope = "none", energyNope = "none", seed = 4)
    )
    expect_equal(rect_layer_count(none_plot), 0)

    energy_only_plot <- PlotCompositionContrasts(
        CompareGroupCompositions(long_table, method = "bayes", nPosterior = 10, nope = "none", energyNope = 0.01, seed = 4)
    )
    expect_silent(ggplot2::ggplot_build(energy_only_plot))
    energy_band <- Filter(function(layer) inherits(layer$geom, "GeomRect"), energy_only_plot$layers)[[1]]$data
    expect_true(all(energy_band$xmin == 0 & energy_band$xmax == 0.01))

    small_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 3,
        n_categories = 4,
        seed = 3
    )
    small_result <- suppressWarnings(
        CompareGroupCompositions(small_table, method = "bayes", nPosterior = 10, nNopeShuffles = 20, seed = 4)
    )
    small_plot <- PlotCompositionContrasts(small_result)
    expect_equal(rect_layer_count(small_plot), 0)
    expect_silent(ggplot2::ggplot_build(small_plot))
})

test_that("PlotCompositionContrasts stacks R and E panels unless asked not to", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 5,
        n_categories = 4,
        groups = c("A", "B", "C"),
        seed = 3
    )
    panel_count <- function(plot) nrow(ggplot2::ggplot_build(plot)$layout$layout)

    bayes_result <- CompareGroupCompositions(long_table, method = "bayes", nPosterior = 10, nNopeShuffles = 20, seed = 4)
    stacked_plot <- PlotCompositionContrasts(bayes_result)
    expect_equal(panel_count(stacked_plot), 2)
    expect_equal(ggplot2::ggplot_build(stacked_plot)$layout$layout$ROW, c(1, 2))
    halfeye_layer <- Filter(function(layer) inherits(layer$stat, "StatSlabinterval"), stacked_plot$layers)[[1]]
    expect_equal(halfeye_layer$geom_params$normalize, "panels")

    expect_equal(panel_count(PlotCompositionContrasts(bayes_result, showEnergyDistance = FALSE)), 1)

    perm_result <- CompareGroupCompositions(long_table, method = "permutation", nPermutations = 20, seed = 4)
    expect_equal(panel_count(PlotCompositionContrasts(perm_result)), 2)
})

test_that("ComputeEnergyNopeBound turns a composition shift into an E bound", {
    reference <- c(TypeA = 0.5, TypeB = 0.3, TypeC = 0.2)
    shifted <- c(TypeA = 0.45, TypeB = 0.35, TypeC = 0.2)
    expected <- 2 * (1 - sum(sqrt(reference * shifted)))

    expect_equal(ComputeEnergyNopeBound(reference, shifted), expected)
    # Counts normalize to the same frequencies; named categories match by name.
    expect_equal(ComputeEnergyNopeBound(reference * 1000, shifted * 40), expected)
    expect_equal(ComputeEnergyNopeBound(reference, rev(shifted)), expected)

    expect_error(ComputeEnergyNopeBound(reference, reference), "same composition")
    expect_error(ComputeEnergyNopeBound(reference, c(TypeA = 0.5, TypeB = 0.5, TypeD = 0)), "same categories")
    expect_error(ComputeEnergyNopeBound(c(0.5, 0.5), c(0.2, 0.3, 0.5)), "same number of categories")
    expect_error(ComputeEnergyNopeBound(c(0.5, -0.1, 0.6), c(0.2, 0.3, 0.5)), "non-negative")

    # The bound goes straight into energyNope.
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 4,
        n_categories = 3,
        seed = 12
    )
    energy_bound <- ComputeEnergyNopeBound(c(1, 1, 1), c(0.9, 1.1, 1))
    result <- CompareGroupCompositions(
        long_table,
        method = "bayes",
        nPosterior = 5,
        nope = "none",
        energyNope = energy_bound,
        seed = 1
    )
    expect_true(all(result$contrasts$E_nopeHigh == energy_bound))
})

test_that("NOPE arguments warn when they cannot apply", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 4,
        n_categories = 3,
        seed = 12
    )

    expect_warning(
        CompareGroupCompositions(
            long_table,
            method = "permutation",
            nPermutations = 10,
            nope = "none",
            seed = 1
        ),
        "nope unused"
    )
    expect_warning(
        CompareGroupCompositions(
            long_table,
            method = "bayes",
            nPosterior = 3,
            nope = "none",
            energyNope = "none",
            nNopeShuffles = 20,
            seed = 1
        ),
        "nNopeShuffles"
    )
    expect_error(
        CompareGroupCompositions(
            long_table,
            method = "bayes",
            nPosterior = 3,
            nNopeShuffles = 0,
            seed = 1
        ),
        "nNopeShuffles"
    )
})

test_that("Holm adjustment is monotonic on contrasts", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 5,
        n_categories = 4,
        seed = 3
    )
    structured <- HellingerDistanceEnrichment:::plant_group_structure(long_table, boost = 6)

    result <- CompareGroupCompositions(
        structured,
        method = "permutation",
        nPermutations = 100,
        pAdjustMethod = "holm",
        seed = 1
    )

    expect_true(all(result$contrasts$pAdj >= result$contrasts$pValue))
    expect_equal(
        result$contrasts$pAdj,
        stats::p.adjust(result$contrasts$pValue, method = "holm")
    )
})

test_that("study raw data is excluded from package and ignore rules", {
    expect_false(file.exists(system.file("raw_data.csv", package = "HellingerDistanceEnrichment")))

    pkg_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = FALSE)
    buildignore_path <- file.path(pkg_root, ".Rbuildignore")

    if (file.exists(buildignore_path)) {
        buildignore_lines <- readLines(buildignore_path, warn = FALSE)
        expect_true(any(grepl("raw_data", buildignore_lines)))
    }
})

test_that("PlotCompositionContrasts returns ggplot for permutation result", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(seed = 2)
    result <- CompareGroupCompositions(
        long_table,
        method = "permutation",
        nPermutations = 30,
        seed = 4
    )
    p <- PlotCompositionContrasts(result)
    expect_s3_class(p, "ggplot")
})

test_that("PlotCompositionContrasts uses halfeye for Bayes result", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 4,
        n_categories = 3,
        seed = 7
    )
    result <- CompareGroupCompositions(
        long_table,
        method = "bayes",
        nPosterior = 5,
        nCores = 1,
        seed = 8
    )
    p <- PlotCompositionContrasts(result)
    expect_s3_class(p, "ggplot")
    layer_classes <- vapply(p$layers, function(layer) class(layer$stat)[1], character(1))
    expect_true(any(grepl("halfeye|slabinterval", layer_classes, ignore.case = TRUE)))
})

test_that("custom contrast subset runs", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(seed = 5)
    result <- CompareGroupCompositions(
        long_table,
        method = "permutation",
        nPermutations = 30,
        contrasts = list(
            pair_only = c("Control", "Treatment")
        ),
        seed = 6
    )
    expect_true("pair_only" %in% result$contrasts$contrastId)
    expect_equal(result$contrasts$contrastId, "pair_only")
    expect_true(!is.null(result$omnibus$effectSize))
})

test_that("a length-2 contrasts vector is a planned pairwise family", {
    groups <- c("Control", "Condition1", "Condition2")
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 4,
        n_categories = 3,
        groups = groups,
        seed = 14
    )

    result <- CompareGroupCompositions(
        long_table,
        method = "permutation",
        nPermutations = 20,
        contrasts = c("Control", "Condition1"),
        seed = 3
    )

    expect_equal(result$contrasts$contrastId, "Control_vs_Condition1")
    expect_false("Condition1_vs_Condition2" %in% result$contrasts$contrastId)
    expect_false("Condition2_vs_Control" %in% result$contrasts$contrastId)
    expect_true(is.finite(result$omnibus$effectSize))
})

test_that("unnamed contrast list keeps only the supplied pairs", {
    groups <- c("Control", "Condition1", "Condition2")
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 4,
        n_categories = 3,
        groups = groups,
        seed = 15
    )

    result <- CompareGroupCompositions(
        long_table,
        method = "permutation",
        nPermutations = 20,
        contrasts = list(
            c("Control", "Condition1"),
            c("Control", "Condition2")
        ),
        seed = 4
    )

    expect_setequal(
        result$contrasts$contrastId,
        c("Control_vs_Condition1", "Control_vs_Condition2")
    )
    expect_false("Condition1_vs_Condition2" %in% result$contrasts$contrastId)
})

test_that("unknown contrast group names stop", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(seed = 16)

    expect_error(
        CompareGroupCompositions(
            long_table,
            method = "permutation",
            nPermutations = 10,
            contrasts = c("Control", "MissingGroup"),
            seed = 5
        ),
        "MissingGroup"
    )
})

test_that("includePairwise TRUE restores additive pairwise overlay", {
    groups <- c("Control", "Condition1", "Condition2")
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 4,
        n_categories = 3,
        groups = groups,
        seed = 17
    )

    result <- CompareGroupCompositions(
        long_table,
        method = "permutation",
        nPermutations = 20,
        contrasts = list(primary = c("Control", "Condition1")),
        includePairwise = TRUE,
        seed = 6
    )

    expect_true("primary" %in% result$contrasts$contrastId)
    expect_true("Condition1_vs_Condition2" %in% result$contrasts$contrastId)
    expect_true("Condition2_vs_Control" %in% result$contrasts$contrastId)
})

test_that("includePairwise FALSE without contrasts is omnibus only", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(seed = 18)

    expect_warning(
        result <- CompareGroupCompositions(
            long_table,
            method = "permutation",
            nPermutations = 20,
            includePairwise = FALSE,
            seed = 7
        ),
        "only the omnibus test will run"
    )

    expect_equal(nrow(result$contrasts), 0)
    expect_true(is.finite(result$omnibus$effectSize))
    expect_error(
        PlotCompositionContrasts(result),
        "omnibus-only"
    )
})

test_that("packaged h5ad extracts to CategoryComposition", {
    skip_if_not_installed("anndata")
    h5ad_path <- system.file(
        "extdata", "tiny_obs.h5ad",
        package = "HellingerDistanceEnrichment"
    )
    skip_if_not(nzchar(h5ad_path) && file.exists(h5ad_path))
    adata <- try(anndata::read_h5ad(h5ad_path), silent = TRUE)
    skip_if(inherits(adata, "try-error"), "Python anndata/read_h5ad unavailable")

    composition <- ExtractClusterComposition(
        adata,
        subjectCol = "subjectId",
        categoryCol = "category",
        groupCol = "group"
    )
    expect_s3_class(composition, "CategoryComposition")
    expect_equal(length(composition$subjectIds), 4)
    expect_setequal(as.character(composition$group), c("A", "B"))
})

test_that("subject with multiple groups errors", {
    long_table <- data.frame(
        subjectId = c("S1", "S1"),
        category = c("A", "A"),
        group = c("G1", "G2"),
        n = c(1, 2),
        stringsAsFactors = FALSE
    )
    expect_error(
        HellingerDistanceEnrichment:::coerce_long_table_to_composition(long_table),
        "multiple groups"
    )
})
