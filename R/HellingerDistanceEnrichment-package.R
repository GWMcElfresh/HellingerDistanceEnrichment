#' @keywords internal
"_PACKAGE"

#' Hellinger distance enrichment for categorical composition data.
#'
#' @description
#' Compare subject-level categorical compositions (cluster counts, cell-type
#' proportions) across experimental groups using Hellinger-distance enrichment
#' ratios with permutation or conjugate Bayesian posterior-predictive inference.
#'
#' @import ggplot2 ggdist egg parallel reticulate Seurat anndata
#' @name HellingerDistanceEnrichment
NULL

utils::globalVariables(c(
    "contrastId",
    "effectCiHigh",
    "effectCiLow",
    "effectSize",
    "significant"
))
