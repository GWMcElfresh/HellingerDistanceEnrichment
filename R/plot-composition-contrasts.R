#' Plot Hellinger enrichment contrasts with uncertainty.
#'
#' Displays contrast-level effect sizes (between/within Hellinger ratio) with
#' ggdist half-eye densities for Bayes results or point estimates for
#' permutation results. The omnibus statistic is shown as a reference line.
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
        stop("result has no contrasts to plot")
    }

    contrast_levels <- rev(plot_data$contrastId)
    plot_data$contrastId <- factor(plot_data$contrastId, levels = contrast_levels)
    plot_data$significant <- plot_data$pAdj < 0.05

    if (result$method == "bayes") {
        if (is.null(result$draws) || nrow(result$draws) == 0) {
            stop("Bayes result has no posterior draws to plot; re-run CompareGroupCompositions")
        }

        draws_data <- result$draws
        draws_data$contrastId <- factor(draws_data$contrastId, levels = contrast_levels)

        p <- ggplot2::ggplot(
            draws_data,
            ggplot2::aes(x = .data$effectSize, y = .data$contrastId)
        ) +
            ggdist::stat_halfeye(...)
    } else {
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
                y = 0.5,
                label = sprintf("omnibus (p=%.3f)", result$omnibus$pValue),
                hjust = -0.05,
                size = 3,
                color = "#009E73"
            )
    }

    p
}
