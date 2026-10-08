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
