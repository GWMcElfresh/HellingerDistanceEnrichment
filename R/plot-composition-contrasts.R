#' Plot Hellinger enrichment contrasts with uncertainty.
#'
#' Takes a HellingerEnrichmentResult and returns a ggplot with one row per
#' pairwise contrast. By default the plot stacks two panels that share those
#' rows: the between/within Hellinger ratio \eqn{R} on top and the energy
#' distance \eqn{E} below, each with its own x axis, so both read left to right. \eqn{R}
#' carries the p-values, `PPGT1`, and `R_PPGTnope`; \eqn{E}
#' shows the size of the shift in squared Hellinger units. A dotted line marks
#' each statistic's no-difference value (\eqn{R = 1}, \eqn{E = 0}). Set
#' `showEnergyDistance = FALSE` for the ratio alone.
#'
#' Bayes results draw from the retained posterior draws (`nPosterior` rows per
#' contrast) via `ggdist::stat_halfeye`, because a density slab needs samples
#' rather than summarized intervals. Half-eye nested intervals follow
#' `result$settings$credibleInterval` (HPDI via `ggdist::mean_hdi`, CPI via
#' `ggdist::mean_qi`) unless `...` overrides `point_interval`. When a Bayes
#' result carries a NOPE, permutation-calibrated or fixed, a gray band behind each half-eye
#' spans that statistic's interval: \eqn{R} from 1 to `R_nopeHigh` in the top panel, \eqn{E} from
#' 0 to `E_nopeHigh` in the bottom panel. You can then see how much of each
#' posterior clears its edge. The legend names each band: "NOPE (permutation)"
#' for the default permutation-calibrated bound, whose band shows the median of
#' the per-draw bounds, or "NOPE (fixed)" for a bound you supplied. Permutation
#' results have no draws, so they render as points colored by adjusted
#' significance. When `showOmnibus` is TRUE, a dashed line in each panel marks
#' the omnibus value; the \eqn{R} panel labels it with its p-value (permutation)
#' or `PPGT1` and `R_PPGTnope` (Bayes), and the \eqn{E} panel labels it with \eqn{E}
#' and, under Bayes, `E_PPGTnope`.
#' Omnibus draws stay in `result$omnibus$draws` rather than as a plot row.
#'
#' @param result A HellingerEnrichmentResult from CompareGroupCompositions.
#' @param theme ggplot2 theme (default egg::theme_article()). With the
#'   default, the legend sits above the panels; a theme you supply keeps its
#'   own legend placement.
#' @param showOmnibus If TRUE, mark the omnibus value in each panel with a
#'   dashed line.
#' @param showEnergyDistance If TRUE (default), stack an energy distance
#'   panel below the ratio panel.
#' @param ... Passed to ggdist::stat_halfeye for Bayes plots.
#' @return A ggplot object.
#' @export
#' @importFrom ggdist stat_halfeye
PlotCompositionContrasts <- function(result,
                                     theme = egg::theme_article(),
                                     showOmnibus = TRUE,
                                     showEnergyDistance = TRUE,
                                     ...) {
    use_default_theme <- missing(theme)
    if (!inherits(result, "HellingerEnrichmentResult")) {
        stop("result must be a HellingerEnrichmentResult")
    }

    plot_data <- result$contrasts
    if (nrow(plot_data) == 0) {
        stop("result has no contrasts to plot; omnibus-only results are in result$omnibus")
    }

    contrast_levels <- rev(plot_data$contrastId)
    plot_data$contrastId <- factor(plot_data$contrastId, levels = contrast_levels)
    is_bayes <- identical(result$method, "bayes")

    if (is_bayes) {
        # Half-eye needs one value per posterior draw; summarized CI columns
        # alone cannot rebuild the density slab.
        if (is.null(result$draws) || nrow(result$draws) == 0) {
            stop("Bayes result has no posterior draws to plot; re-run CompareGroupCompositions")
        }
        # The result keeps omnibus draws for extraction, but they are not a contrast row.
        source_data <- result$draws[result$draws$contrastId %in% plot_data$contrastId, , drop = FALSE]
        source_data$contrastId <- factor(source_data$contrastId, levels = contrast_levels)
    } else {
        source_data <- plot_data
        source_data$significant <- source_data$pAdj < 0.05
    }

    statistic_labels <- c(ratio = "Between / within ratio (R)", energy = "Energy distance (E)")
    show_energy <- isTRUE(showEnergyDistance)
    long_data <- stack_plot_statistics(source_data, show_energy, statistic_labels)
    shown_labels <- levels(long_data$statistic)

    p <- ggplot2::ggplot(
        long_data,
        ggplot2::aes(x = .data$value, y = .data$contrastId)
    ) +
        # An explicit discrete scale lets numeric positions (NOPE bands, the
        # omnibus labels) share the y axis with the factor-valued rows.
        ggplot2::scale_y_discrete(limits = contrast_levels)

    if (is_bayes) {
        halfeye_args <- list(...)
        # Scale slab heights within each panel; across panels, the narrow E
        # densities would flatten the R slabs.
        if (is.null(halfeye_args$normalize)) {
            halfeye_args$normalize <- "panels"
        }
        if (is.null(halfeye_args$point_interval)) {
            halfeye_args$point_interval <- if (identical(result$settings$credibleInterval, "cpi")) {
                ggdist::mean_qi
            } else {
                ggdist::mean_hdi
            }
        }

        nope_band <- rbind(
            build_nope_band(plot_data, contrast_levels, "R", result$settings$nope, halfeye_args$scale),
            if (show_energy) {
                build_nope_band(plot_data, contrast_levels, "E", result$settings$energyNope, halfeye_args$scale)
            }
        )
        if (!is.null(nope_band)) {
            nope_band$statistic <- factor(
                statistic_labels[ifelse(nope_band$statistic == "R", "ratio", "energy")],
                levels = shown_labels
            )
            nope_fill <- c("NOPE (permutation)" = "gray70", "NOPE (fixed)" = "gray40")[unique(nope_band$nopeLabel)]
            p <- p +
                ggplot2::geom_rect(
                    data = nope_band,
                    ggplot2::aes(
                        xmin = .data$xmin,
                        xmax = .data$xmax,
                        ymin = .data$ymin,
                        ymax = .data$ymax,
                        fill = .data$nopeLabel
                    ),
                    inherit.aes = FALSE,
                    alpha = 0.35
                ) +
                # A darker gray marks a fixed edge, so a mixed setting (e.g. a
                # fixed R edge with an E NOPE) shows two distinct swatches.
                ggplot2::scale_fill_manual(values = nope_fill, name = NULL)
        }

        p <- p + do.call(ggdist::stat_halfeye, halfeye_args)
    } else {
        p <- p +
            ggplot2::geom_point(
                ggplot2::aes(color = .data$significant),
                size = 3
            ) +
            ggplot2::scale_color_manual(
                values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2"),
                name = "Adj. p < 0.05"
            )
    }

    reference_lines <- data.frame(
        statistic = factor(unname(statistic_labels), levels = shown_labels),
        xintercept = c(1, 0)
    )
    reference_lines <- reference_lines[!is.na(reference_lines$statistic), , drop = FALSE]

    p <- p +
        ggplot2::geom_vline(
            data = reference_lines,
            ggplot2::aes(xintercept = .data$xintercept),
            linetype = "dotted",
            color = "gray50"
        ) +
        # Extra expansion on both sides leaves room for the omnibus labels.
        ggplot2::scale_x_continuous(
            expand = ggplot2::expansion(mult = c(0.08, 0.15))
        ) +
        ggplot2::labs(
            x = if (show_energy) NULL else "Between / within Hellinger ratio",
            y = NULL,
            title = sprintf("Composition contrasts (%s)", result$method)
        )
    if (show_energy) {
        p <- p + ggplot2::facet_wrap(~ statistic, ncol = 1, scales = "free_x")
    }
    p <- p + theme

    # Only with the default theme: a theme the caller supplies keeps its own
    # legend placement.
    has_nope_columns <- any(c("R_nopeHigh", "E_nopeHigh") %in% colnames(plot_data))
    if (use_default_theme && (!is_bayes || has_nope_columns)) {
        p <- p + ggplot2::theme(legend.position = "top")
    }

    if (showOmnibus) {
        omnibus_marks <- build_omnibus_marks(result, long_data, statistic_labels, shown_labels)
        p <- p +
            ggplot2::geom_vline(
                data = omnibus_marks,
                ggplot2::aes(xintercept = .data$value),
                linetype = "dashed",
                color = "#009E73",
                linewidth = 0.8
            ) +
            ggplot2::geom_text(
                data = omnibus_marks,
                ggplot2::aes(x = .data$value, y = 0.5, label = .data$label, hjust = .data$hjust),
                inherit.aes = FALSE,
                size = 3,
                color = "#009E73"
            )
    }

    p
}

#' Stack the ratio and energy distance into one long table for faceting.
#'
#' Each input row (a posterior draw or a contrast summary) appears once per
#' shown statistic, with the statistic's value in `value` and its panel label
#' in `statistic`. The ratio panel always comes first.
#'
#' @param source_data Draws (Bayes) or contrast rows (permutation) carrying
#'   `effectSize` and, when `show_energy` is TRUE, `energyDistance`.
#' @param show_energy Whether to include the energy distance panel.
#' @param statistic_labels Named panel labels (`ratio`, `energy`).
#' @return data.frame with `value` and a factor `statistic`.
#' @keywords internal
stack_plot_statistics <- function(source_data, show_energy, statistic_labels) {
    ratio_rows <- source_data
    ratio_rows$value <- ratio_rows$effectSize
    ratio_rows$statistic <- statistic_labels[["ratio"]]
    stacked <- ratio_rows

    if (show_energy) {
        energy_rows <- source_data
        energy_rows$value <- energy_rows$energyDistance
        energy_rows$statistic <- statistic_labels[["energy"]]
        stacked <- rbind(ratio_rows, energy_rows)
    }

    shown_labels <- if (show_energy) unname(statistic_labels) else statistic_labels[["ratio"]]
    stacked$statistic <- factor(stacked$statistic, levels = shown_labels)
    stacked
}

#' Omnibus line positions and labels for each shown panel.
#'
#' The ratio panel's label carries the decision statistics (the permutation
#' p-value, or `PPGT1` and `R_PPGTnope` under Bayes); the energy panel's label
#' carries \eqn{E} and, under Bayes, `E_PPGTnope`. Each label sits on whichever side of its line has
#' more room within the panel, counting the reference line, so the panel edge
#' does not clip it.
#'
#' @param result A HellingerEnrichmentResult.
#' @param long_data Output of stack_plot_statistics.
#' @param statistic_labels Named panel labels (`ratio`, `energy`).
#' @param shown_labels Panel labels present in the plot.
#' @return data.frame with `statistic`, `value`, `label`, and `hjust`.
#' @keywords internal
build_omnibus_marks <- function(result, long_data, statistic_labels, shown_labels) {
    ratio_label <- if (identical(result$method, "bayes") && !is.null(result$omnibus$R_PPGTnope)) {
        sprintf("omnibus (PPGT1=%.3f, R_PPGTnope=%.3f)", result$omnibus$PPGT1, result$omnibus$R_PPGTnope)
    } else if (identical(result$method, "bayes")) {
        sprintf("omnibus (PPGT1=%.3f)", result$omnibus$PPGT1)
    } else {
        sprintf("omnibus (p=%.3f)", result$omnibus$pValue)
    }
    marks <- data.frame(
        statistic = statistic_labels[["ratio"]],
        value = result$omnibus$effectSize,
        label = ratio_label,
        stringsAsFactors = FALSE
    )

    if (statistic_labels[["energy"]] %in% shown_labels && !is.null(result$omnibus$energyDistance)) {
        marks <- rbind(marks, data.frame(
            statistic = statistic_labels[["energy"]],
            value = result$omnibus$energyDistance,
            label = if (!is.null(result$omnibus$E_PPGTnope)) {
                sprintf("omnibus (E=%.3g, E_PPGTnope=%.3f)", result$omnibus$energyDistance, result$omnibus$E_PPGTnope)
            } else {
                sprintf("omnibus (E=%.3g)", result$omnibus$energyDistance)
            },
            stringsAsFactors = FALSE
        ))
    }

    # The reference line (R = 1, E = 0) widens each panel's range, so count it
    # when measuring the room on each side of the omnibus line.
    reference_values <- c(1, 0)
    names(reference_values) <- unname(statistic_labels)
    marks$hjust <- vapply(seq_len(nrow(marks)), function(mark_idx) {
        panel_values <- long_data$value[long_data$statistic == marks$statistic[mark_idx]]
        panel_range <- range(c(panel_values, reference_values[[marks$statistic[mark_idx]]]), na.rm = TRUE)
        room_left <- marks$value[mark_idx] - panel_range[1]
        room_right <- panel_range[2] - marks$value[mark_idx]
        if (room_left > room_right) 1.05 else -0.05
    }, numeric(1))
    marks$statistic <- factor(marks$statistic, levels = shown_labels)
    marks
}

#' Build the shaded NOPE band for one statistic in each Bayes contrast row.
#'
#' Each band runs from the statistic's no-difference value (1 for \eqn{R}, 0 for
#' \eqn{E}) to the contrast's `R_nopeHigh` or `E_nopeHigh` and covers the height
#' of that row's half-eye slab, which `ggdist::stat_halfeye` draws upward from
#' the row by `slab_scale` (0.9 by default). For the NOPE, the edge is the
#' median of the per-draw edges. Contrasts whose NOPE probability is `NA`
#' (groups too small for the NOPE) get no band, because the edge there
#' carries no information.
#'
#' @param plot_data Contrast summary rows with `contrastId` as a factor.
#' @param contrast_levels Factor levels in plotting order.
#' @param statistic_prefix `"R"` or `"E"`, matching the column prefixes.
#' @param nope_setting `result$settings$nope` or `result$settings$energyNope`.
#' @param slab_scale Half-eye `scale` argument, or `NULL` for the default.
#' @return data.frame of rectangles with a `statistic` column, or `NULL` when
#'   the result has no NOPE or fixed edge on that statistic.
#' @keywords internal
build_nope_band <- function(plot_data, contrast_levels, statistic_prefix, nope_setting, slab_scale = NULL) {
    edge_column <- paste0(statistic_prefix, "_nopeHigh")
    probability_column <- paste0(statistic_prefix, "_PPGTnope")
    if (!edge_column %in% colnames(plot_data)) {
        return(NULL)
    }
    has_edge <- is.finite(plot_data[[edge_column]]) & !is.na(plot_data[[probability_column]])
    if (!any(has_edge)) {
        return(NULL)
    }
    if (is.null(slab_scale)) {
        slab_scale <- 0.9
    }

    nope_label <- if (identical(nope_setting, "shuffle")) {
        "NOPE (permutation)"
    } else {
        "NOPE (fixed)"
    }
    row_positions <- match(as.character(plot_data$contrastId[has_edge]), contrast_levels)
    data.frame(
        statistic = statistic_prefix,
        xmin = if (identical(statistic_prefix, "R")) 1 else 0,
        xmax = plot_data[[edge_column]][has_edge],
        ymin = row_positions,
        ymax = row_positions + slab_scale,
        nopeLabel = nope_label,
        stringsAsFactors = FALSE
    )
}
