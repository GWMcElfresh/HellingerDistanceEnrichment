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

test_that("Bayes nesting smoke test returns finite posterior-mean p", {
    long_table <- HellingerDistanceEnrichment:::build_synthetic_long_table(
        n_subjects_per_group = 5,
        n_categories = 4,
        seed = 11
    )

    result <- CompareGroupCompositions(
        long_table,
        method = "bayes",
        nPosterior = 5,
        nPermutations = 20,
        nCores = 1,
        seed = 55
    )

    expect_s3_class(result, "HellingerEnrichmentResult")
    expect_true(is.finite(result$omnibus$pValue))
    expect_true(all(is.finite(result$contrasts$pValue)))
    expect_true(all(is.finite(result$contrasts$effectCiLow)))
    expect_true(!is.null(result$draws))
    expect_gt(nrow(result$draws), 0)
    expect_true(all(c("contrastId", "draw", "effectSize") %in% colnames(result$draws)))
    expect_false("omnibus" %in% result$draws$contrastId)
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
        nPermutations = 10,
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
