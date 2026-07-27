#!/usr/bin/env Rscript

# CLI entrypoint for Hellinger distance enrichment on bind-mounted inputs.
# Study data must never be baked into the image; pass paths at runtime.

suppressPackageStartupMessages({
    library(HellingerDistanceEnrichment)
})

parse_args <- function(args) {
    if (length(args) == 0 || args[1] %in% c("-h", "--help")) {
        cat(paste0(
            "Usage: hellinger-enrichment --input <path> --format <long|seurat|h5ad> ",
            "[options]\n\n",
            "Options:\n",
            "  --input PATH          Input file (CSV long table, Seurat RDS, or h5ad)\n",
            "  --format FORMAT       long, seurat, or h5ad (default: long)\n",
            "  --subject-col NAME    Subject column (default: subjectId)\n",
            "  --category-col NAME   Category column (default: category)\n",
            "  --group-col NAME      Group column (default: group)\n",
            "  --count-col NAME      Count column for long format (default: n)\n",
            "  --method METHOD       permutation or bayes (default: permutation)\n",
            "  --n-perm N            Permutations (default: 1000)\n",
            "  --n-posterior N       Posterior draws for bayes (default: 100)\n",
            "  --n-cores N           Parallel workers (default: 1)\n",
            "  --seed N              Random seed\n",
            "  --output-dir PATH     Output directory (default: ./hellinger_output)\n",
            "  --plot                Write contrast PDF via ggplot2\n"
        ))
        quit(status = 0)
    }

    defaults <- list(
        format = "long",
        subject_col = "subjectId",
        category_col = "category",
        group_col = "group",
        count_col = "n",
        method = "permutation",
        n_perm = 1000,
        n_posterior = 100,
        n_cores = 1,
        seed = NA,
        output_dir = "./hellinger_output",
        plot = FALSE
    )

    i <- 1
    while (i <= length(args)) {
        key <- args[i]
        if (key == "--input") {
            defaults$input <- args[i + 1]; i <- i + 2; next
        }
        if (key == "--format") {
            defaults$format <- args[i + 1]; i <- i + 2; next
        }
        if (key == "--subject-col") {
            defaults$subject_col <- args[i + 1]; i <- i + 2; next
        }
        if (key == "--category-col") {
            defaults$category_col <- args[i + 1]; i <- i + 2; next
        }
        if (key == "--group-col") {
            defaults$group_col <- args[i + 1]; i <- i + 2; next
        }
        if (key == "--count-col") {
            defaults$count_col <- args[i + 1]; i <- i + 2; next
        }
        if (key == "--method") {
            defaults$method <- args[i + 1]; i <- i + 2; next
        }
        if (key == "--n-perm") {
            defaults$n_perm <- as.integer(args[i + 1]); i <- i + 2; next
        }
        if (key == "--n-posterior") {
            defaults$n_posterior <- as.integer(args[i + 1]); i <- i + 2; next
        }
        if (key == "--n-cores") {
            defaults$n_cores <- as.integer(args[i + 1]); i <- i + 2; next
        }
        if (key == "--seed") {
            defaults$seed <- as.integer(args[i + 1]); i <- i + 2; next
        }
        if (key == "--output-dir") {
            defaults$output_dir <- args[i + 1]; i <- i + 2; next
        }
        if (key == "--plot") {
            defaults$plot <- TRUE; i <- i + 1; next
        }
        stop(sprintf("unknown argument: %s", key))
    }

    if (is.null(defaults$input)) {
        stop("--input is required")
    }

    defaults
}

load_input <- function(path, format, subject_col, category_col, group_col, count_col) {
    if (format == "long") {
        table <- read.csv(path, stringsAsFactors = FALSE)
        return(ExtractClusterComposition(
            table,
            subjectCol = subject_col,
            categoryCol = category_col,
            groupCol = group_col,
            countCol = count_col
        ))
    }

    if (format == "seurat") {
        object <- readRDS(path)
        return(ExtractClusterComposition(
            object,
            subjectCol = subject_col,
            categoryCol = category_col,
            groupCol = group_col
        ))
    }

    if (format == "h5ad") {
        if (!requireNamespace("anndata", quietly = TRUE)) {
            stop("anndata package required for h5ad input")
        }
        object <- anndata::read_h5ad(path)
        return(ExtractClusterComposition(
            object,
            subjectCol = subject_col,
            categoryCol = category_col,
            groupCol = group_col
        ))
    }

    stop(sprintf("unsupported format: %s", format))
}

main <- function() {
    opts <- parse_args(commandArgs(trailingOnly = TRUE))

    dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)

    composition <- load_input(
        opts$input,
        opts$format,
        opts$subject_col,
        opts$category_col,
        opts$group_col,
        opts$count_col
    )

    result <- CompareGroupCompositions(
        composition,
        method = opts$method,
        nPermutations = opts$n_perm,
        nPosterior = opts$n_posterior,
        nCores = opts$n_cores,
        seed = if (is.na(opts$seed)) NULL else opts$seed
    )

    saveRDS(result, file.path(opts$output_dir, "hellinger_result.rds"))
    write.csv(result$contrasts, file.path(opts$output_dir, "contrasts.csv"), row.names = FALSE)

    omnibus_df <- data.frame(
        effectSize = result$omnibus$effectSize,
        pValue = result$omnibus$pValue
    )
    write.csv(omnibus_df, file.path(opts$output_dir, "omnibus.csv"), row.names = FALSE)

    if (isTRUE(opts$plot)) {
        p <- PlotCompositionContrasts(result)
        ggplot2::ggsave(
            filename = file.path(opts$output_dir, "contrasts.pdf"),
            plot = p,
            width = 7,
            height = 5
        )
    }

    cat(sprintf("Wrote results to %s\n", opts$output_dir))
}

main()
