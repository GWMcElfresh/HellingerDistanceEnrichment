#' Plot Hellinger enrichment contrasts with uncertainty.
#'
#' Takes a HellingerEnrichmentResult, maps each pairwise contrast onto the
#' between/within Hellinger ratio axis, and returns a ggplot. Bayes results
#' are drawn from the retained posterior effect-size draws
#' (`nPosterior` rows per contrast) via `ggdist::stat_halfeye`, because a
#' density slab needs samples rather than summarized `effectCiLow`/`effectCiHigh`
#' intervals. Half-eye nested intervals follow `result$settings$credibleInterval`
#' (HPDI via `ggdist::mean_hdi`, CPI via `ggdist::mean_qi`) unless `...`
#' overrides `point_interval`. Permutation results have no draws, so they
#' render as points colored by adjusted significance. When `showOmnibus` is
#' TRUE, the omnibus ratio is overlaid as a dashed reference with its p-value
#' (permutation) or PPGT1 (Bayes) labeled beside the line. Omnibus draws stay
#' in `result$omnibus$draws` rather than as a plot row.
#'
#' @param result A HellingerEnrichmentResult from CompareGroupCompositions.
#' @param theme ggplot2 theme (default egg::theme_article()).
#' @param showOmnibus If TRUE, annotate the omnibus effect size as a dashed line.
#' @param ... Passed to ggdist::stat_halfeye for Bayes plots.
#' @return A ggplot object.
#' @export
#' @importFrom ggdist stat_halfeye
PlotCompositionContrasts <- function(result,
                                     theme = egg::theme_article(),
                                     showOmnibus = TRUE,
                                     ...) {
    if (!inherits(result, "HellingerEnrichmentResult")) {
        stop("result must be a HellingerEnrichmentResult")
    }

    plot_data <- result$contrasts
    if (nrow(plot_data) == 0) {
        stop("result has no contrasts to plot; omnibus-only results are in result$omnibus")
    }

    contrast_levels <- rev(plot_data$contrastId)
    plot_data$contrastId <- factor(plot_data$contrastId, levels = contrast_levels)

    if (result$method == "bayes") {
        # Half-eye needs one effectSize per posterior draw; summarized CI
        # columns alone cannot rebuild the density slab.
        if (is.null(result$draws) || nrow(result$draws) == 0) {
            stop("Bayes result has no posterior draws to plot; re-run CompareGroupCompositions")
        }

        # Omnibus draws are retained for extraction but are not a contrast row.
        draws_data <- result$draws[result$draws$contrastId %in% plot_data$contrastId, , drop = FALSE]
        draws_data$contrastId <- factor(draws_data$contrastId, levels = contrast_levels)

        halfeye_args <- list(...)
        if (is.null(halfeye_args$point_interval)) {
            interval_type <- result$settings$credibleInterval
            halfeye_args$point_interval <- if (identical(interval_type, "cpi")) {
                ggdist::mean_qi
            } else {
                ggdist::mean_hdi
            }
        }

        p <- ggplot2::ggplot(
            draws_data,
            ggplot2::aes(x = .data$effectSize, y = .data$contrastId)
        ) +
            do.call(ggdist::stat_halfeye, halfeye_args)
    } else {
        plot_data$significant <- plot_data$pAdj < 0.05
        p <- ggplot2::ggplot(
            plot_data,
            ggplot2::aes(x = .data$effectSize, y = .data$contrastId)
        ) +
            ggplot2::geom_point(
                ggplot2::aes(color = .data$significant),
                size = 3
            ) +
            ggplot2::scale_color_manual(
                values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2"),
                name = "Adj. p < 0.05"
            )
    }

    p <- p +
        ggplot2::geom_vline(xintercept = 1, linetype = "dotted", color = "gray50") +
        # Extra right-side expand so the omnibus label (hjust slightly past the
        # vline) is not clipped when the omnibus ratio sits near the data max.
        ggplot2::scale_x_continuous(
            expand = ggplot2::expansion(mult = c(0.05, 0.15))
        ) +
        ggplot2::labs(
            x = "Between / within Hellinger ratio",
            y = NULL,
            title = sprintf("Composition contrasts (%s)", result$method)
        ) +
        theme

    if (showOmnibus) {
        p <- p + ggplot2::geom_vline(
            xintercept = result$omnibus$effectSize,
            linetype = "dashed",
            color = "#009E73",
            linewidth = 0.8
        ) +
            ggplot2::annotate(
                "text",
                x = result$omnibus$effectSize,
                # Discrete contrast levels start at 1; 0.5 sits just below the
                # lowest row without overlapping contrast labels.
                y = 0.5,
                label = if (identical(result$method, "bayes")) {
                    sprintf("omnibus (PPGT1=%.3f)", result$omnibus$PPGT1)
                } else {
                    sprintf("omnibus (p=%.3f)", result$omnibus$pValue)
                },
                hjust = -0.05,
                size = 3,
                color = "#009E73"
            )
    }

    p
}
