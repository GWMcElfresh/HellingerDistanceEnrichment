# Regenerate inst/extdata/sampling_study.rds for the sampling-depth vignette.
#
# Run from the package root:
#   Rscript data-raw/simulate_sampling_study.R
#
# Every cohort takes its seed from a fixed formula, 1e6 * scenario +
# 1e3 * cell + replicate, set before the first run; we did not tune any seed
# to favor an outcome. The full grid runs about 800 Bayes cohorts at the package
# defaults (nPosterior = 100, nNopeShuffles = 200), roughly 10 s each with
# both NOPEs, so the script spreads cohorts over a PSOCK cluster (mclapply
# cannot fork on Windows). Expect about 35 minutes on 5 workers.

library(parallel)

n_workers <- max(1, min(5, detectCores(logical = FALSE) - 1))
n_cohorts_nope <- 100
n_cohorts_depth <- 200
n_cohorts_heterogeneity <- 100
n_cohorts_ranking <- 200
n_cohorts_spread <- 100
n_permutations <- 999

package_root <- normalizePath(".")
output_path <- file.path(package_root, "inst", "extdata", "sampling_study.rds")

cluster <- makeCluster(n_workers)
on.exit(stopCluster(cluster), add = TRUE)
clusterExport(cluster, "package_root")
invisible(clusterEvalQ(cluster, {
    suppressMessages(pkgload::load_all(package_root, quiet = TRUE))
    NULL
}))

worker_functions <- function() {
    omnibus_spec <- list(contrastId = "omnibus", groups = c("Control", "Treatment"), collapse = NULL)

    # Omnibus-only inference through the package internals; with two groups the
    # pairwise contrast would repeat the omnibus and double the cost.
    run_omnibus_inference <- function(long_table, method, seed) {
        composition <- coerce_long_table_to_composition(long_table)
        contrast_specs <- build_contrast_specifications(composition$group, include_pairwise = FALSE)
        set.seed(seed)
        if (method == "bayes") {
            run_bayes_enrichment(
                composition = composition,
                contrast_specs = contrast_specs,
                n_posterior = 100,
                n_cores = 1,
                prior_pseudocounts = 1 / 2,
                nope = "shuffle",
                energy_nope = "shuffle",
                n_nope_shuffles = 200
            )$omnibus
        } else {
            run_permutation_enrichment(
                composition = composition,
                contrast_specs = contrast_specs,
                n_permutations = n_permutations,
                n_cores = 1,
                prior_pseudocounts = 1 / 2
            )$omnibus
        }
    }

    # Noise-free reference: R and E on the latent compositions themselves.
    compute_reference <- function(cohort) {
        compute_observed_contrast(cohort$latent_counts, cohort$group, omnibus_spec, prior_pseudocounts = 0)
    }

    resample_counts_at_depth <- function(cohort, reads_per_subject) {
        latent_freqs <- cohort$latent_counts / rowSums(cohort$latent_counts)
        counts <- t(apply(latent_freqs, 1, function(freqs) stats::rmultinom(1, reads_per_subject, freqs)))
        dimnames(counts) <- dimnames(cohort$latent_counts)
        counts
    }

    list(
        omnibus_spec = omnibus_spec,
        run_omnibus_inference = run_omnibus_inference,
        compute_reference = compute_reference,
        resample_counts_at_depth = resample_counts_at_depth
    )
}
clusterExport(cluster, c("worker_functions", "n_permutations"))
invisible(clusterEvalQ(cluster, {
    helpers <- worker_functions()
    NULL
}))

run_grid <- function(task_table, task_function) {
    task_rows <- split(task_table, seq_len(nrow(task_table)))
    result_rows <- parLapplyLB(cluster, task_rows, task_function)
    out <- do.call(rbind, result_rows)
    rownames(out) <- NULL
    out
}

# Scenario 1: NOPE calibration and power across depth. Null cohorts (no
# shift) and 5-point-shift cohorts at three depths, 8 versus 8 subjects.
nope_tasks <- expand.grid(
    readsPerSubject = c(50, 500, 2000),
    shift = c(0, 0.05),
    replicate = seq_len(n_cohorts_nope)
)
nope_tasks$cell <- as.integer(factor(paste(nope_tasks$readsPerSubject, nope_tasks$shift)))
nope_tasks$seed <- 1e6 * 1 + 1e3 * nope_tasks$cell + nope_tasks$replicate

nope_calibration <- run_grid(nope_tasks, function(task) {
    cohort <- simulate_composition_cohort(
        reads_per_subject = task$readsPerSubject,
        concentration = 50,
        shift = task$shift,
        seed = task$seed
    )
    bayes <- helpers$run_omnibus_inference(cohort$long_table, "bayes", task$seed + 1)
    permutation <- helpers$run_omnibus_inference(cohort$long_table, "permutation", task$seed + 2)
    data.frame(
        readsPerSubject = task$readsPerSubject,
        shift = task$shift,
        replicate = task$replicate,
        PPGT1 = bayes$PPGT1,
        R_PPGTnope = bayes$R_PPGTnope,
        R_nopeHigh = bayes$R_nopeHigh,
        E_PPGTnope = bayes$E_PPGTnope,
        E_nopeHigh = bayes$E_nopeHigh,
        effectSize = bayes$effectSize,
        energyDistance = bayes$energyDistance,
        energyCiLow = bayes$energyCiLow,
        pValue = permutation$pValue
    )
})

# Scenario 2: drift across depth for one fixed cohort. Each replicate draws
# latent compositions once, then resamples counts at every depth, so the only
# thing that changes across depth is count noise.
depth_tasks <- data.frame(replicate = seq_len(n_cohorts_depth))
depth_tasks$seed <- 1e6 * 2 + depth_tasks$replicate

depth_drift <- run_grid(depth_tasks, function(task) {
    cohort <- simulate_composition_cohort(concentration = 50, shift = 0.05, seed = task$seed)
    reference <- helpers$compute_reference(cohort)
    depths <- c(25, 50, 100, 250, 500, 1000, 2000, 5000)
    do.call(rbind, lapply(depths, function(reads_per_subject) {
        counts <- helpers$resample_counts_at_depth(cohort, reads_per_subject)
        observed <- compute_observed_contrast(counts, cohort$group, helpers$omnibus_spec, prior_pseudocounts = 1 / 2)
        data.frame(
            replicate = task$replicate,
            readsPerSubject = reads_per_subject,
            effectSize = observed$effectSize,
            energyDistance = observed$energyDistance,
            referenceEffectSize = reference$effectSize,
            referenceEnergyDistance = reference$energyDistance
        )
    }))
})

# Scenario 3: heterogeneity sweep. One fixed 5-point shift at 2000 reads while
# each group's concentration varies independently; the diagonal holds equal
# spread, off-diagonal cells hold one group more variable than the other.
heterogeneity_tasks <- expand.grid(
    controlConcentration = c(10, 30, 100, 300, 1000),
    treatmentConcentration = c(10, 30, 100, 300, 1000),
    replicate = seq_len(n_cohorts_heterogeneity)
)
heterogeneity_tasks$cell <- as.integer(factor(paste(
    heterogeneity_tasks$controlConcentration,
    heterogeneity_tasks$treatmentConcentration
)))
heterogeneity_tasks$seed <- 1e6 * 3 + 1e3 * heterogeneity_tasks$cell + heterogeneity_tasks$replicate

heterogeneity <- run_grid(heterogeneity_tasks, function(task) {
    cohort <- simulate_composition_cohort(
        concentration = c(task$controlConcentration, task$treatmentConcentration),
        shift = 0.05,
        seed = task$seed
    )
    composition <- coerce_long_table_to_composition(cohort$long_table)
    observed <- compute_observed_contrast(
        composition$counts,
        composition$group,
        helpers$omnibus_spec,
        prior_pseudocounts = 1 / 2
    )
    permutation <- helpers$run_omnibus_inference(cohort$long_table, "permutation", task$seed + 1)
    data.frame(
        controlConcentration = task$controlConcentration,
        treatmentConcentration = task$treatmentConcentration,
        replicate = task$replicate,
        effectSize = observed$effectSize,
        energyDistance = observed$energyDistance,
        pValue = permutation$pValue
    )
})

# Scenario 4: ranking three contrasts. The largest shift sits in the most
# heterogeneous groups and the smallest in the most homogeneous, so ordering by
# shift and ordering by separation relative to spread disagree.
ranking_design <- data.frame(
    contrastLabel = c("small shift, homogeneous", "medium shift, moderate", "large shift, heterogeneous"),
    shift = c(0.02, 0.05, 0.08),
    concentration = c(1000, 100, 20)
)
ranking_tasks <- data.frame(replicate = seq_len(n_cohorts_ranking))
ranking_tasks$seed <- 1e6 * 4 + 1e3 * ranking_tasks$replicate
clusterExport(cluster, "ranking_design")

ranking <- run_grid(ranking_tasks, function(task) {
    do.call(rbind, lapply(seq_len(nrow(ranking_design)), function(design_idx) {
        cohort <- simulate_composition_cohort(
            concentration = ranking_design$concentration[design_idx],
            shift = ranking_design$shift[design_idx],
            seed = task$seed + design_idx
        )
        composition <- coerce_long_table_to_composition(cohort$long_table)
        observed <- compute_observed_contrast(
            composition$counts,
            composition$group,
            helpers$omnibus_spec,
            prior_pseudocounts = 1 / 2
        )
        data.frame(
            replicate = task$replicate,
            contrastLabel = ranking_design$contrastLabel[design_idx],
            shift = ranking_design$shift[design_idx],
            concentration = ranking_design$concentration[design_idx],
            effectSize = observed$effectSize,
            energyDistance = observed$energyDistance
        )
    }))
})

# Scenario 5: spread difference with no shift. Both groups share one mean
# composition; the treatment group is far more variable than control. A
# matched arm with equal spread anchors the comparison.
spread_tasks <- expand.grid(
    treatmentConcentration = c(200, 15),
    replicate = seq_len(n_cohorts_spread)
)
spread_tasks$cell <- as.integer(factor(spread_tasks$treatmentConcentration))
spread_tasks$seed <- 1e6 * 5 + 1e3 * spread_tasks$cell + spread_tasks$replicate

spread_only <- run_grid(spread_tasks, function(task) {
    cohort <- simulate_composition_cohort(
        concentration = c(200, task$treatmentConcentration),
        shift = 0,
        seed = task$seed
    )
    bayes <- helpers$run_omnibus_inference(cohort$long_table, "bayes", task$seed + 1)
    permutation <- helpers$run_omnibus_inference(cohort$long_table, "permutation", task$seed + 2)
    data.frame(
        controlConcentration = 200,
        treatmentConcentration = task$treatmentConcentration,
        replicate = task$replicate,
        effectSize = permutation$effectSize,
        energyDistance = permutation$energyDistance,
        pValue = permutation$pValue,
        PPGT1 = bayes$PPGT1,
        R_PPGTnope = bayes$R_PPGTnope,
        E_PPGTnope = bayes$E_PPGTnope,
        energyCiLow = bayes$energyCiLow
    )
})

git_describe <- tryCatch(
    system2("git", c("describe", "--always", "--dirty"), stdout = TRUE),
    error = function(e) NA_character_
)

sampling_study <- list(
    metadata = list(
        packageVersion = unname(read.dcf(file.path(package_root, "DESCRIPTION"), fields = "Version")[1, 1]),
        gitDescribe = git_describe,
        generatedOn = as.character(Sys.Date()),
        settings = list(
            subjectsPerGroup = 8,
            nCategories = 10,
            nPosterior = 100,
            nNopeShuffles = 200,
            nPermutations = n_permutations,
            priorPseudocounts = 1 / 2,
            seedFormula = "1e6 * scenario + 1e3 * cell + replicate"
        )
    ),
    nopeCalibration = nope_calibration,
    depthDrift = depth_drift,
    heterogeneity = heterogeneity,
    ranking = ranking,
    spreadOnly = spread_only
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(sampling_study, output_path, compress = "xz")
cat(sprintf("Wrote %s\n", output_path))
