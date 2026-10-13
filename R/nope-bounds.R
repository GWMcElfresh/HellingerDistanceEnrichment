#' Energy distance NOPE bound for a composition shift.
#'
#' Translates the smallest shift in mean composition that you would call
#' meaningful into a fixed NOPE bound on the energy distance \eqn{E}, ready to
#' pass as `energyNope` to CompareGroupCompositions. You give the reference
#' composition (e.g. control) and the shifted composition you have in mind;
#' the function normalizes both to frequencies and returns
#' \eqn{2(1 - \mathrm{BC})}{2(1 - BC)}, where \eqn{\mathrm{BC} = \sum_k \sqrt{p_k q_k}}{BC = sum_k sqrt(p_k q_k)} is the
#' Bhattacharyya coefficient between them, in squared Hellinger units.
#'
#' The bound equals \eqn{E} exactly when every subject sits at its group's mean
#' composition. With within-group spread, \eqn{E} is the squared distance between
#' group centroids in square-root composition space, so treat the bound as a
#' close guide rather than an exact translation. No such helper exists for the
#' ratio \eqn{R}, because \eqn{R} also divides by within-group spread; a fixed
#' `nope` bound on \eqn{R} needs a simulation at your own depths and group sizes.
#'
#' @param referenceComposition Non-negative numeric vector of category
#'   frequencies or counts for the reference group. The function normalizes it
#'   to sum to 1.
#' @param shiftedComposition Non-negative numeric vector for the shifted group,
#'   over the same categories. When both vectors are named, the function
#'   matches categories by name; otherwise it matches them by position.
#' @return Scalar energy distance in squared Hellinger units, above 0.
#' @examples
#' reference <- c(TypeA = 0.50, TypeB = 0.30, TypeC = 0.20)
#' shifted <- c(TypeA = 0.45, TypeB = 0.35, TypeC = 0.20)
#' ComputeEnergyNopeBound(reference, shifted)
#'
#' # Counts work as well; the function normalizes them.
#' ComputeEnergyNopeBound(c(500, 300, 200), c(450, 350, 200))
#' @seealso [CompareGroupCompositions()], whose `energyNope` argument takes
#'   the bound.
#' @export
ComputeEnergyNopeBound <- function(referenceComposition, shiftedComposition) {
    reference_freqs <- normalize_composition_vector(referenceComposition, "referenceComposition")
    shifted_freqs <- normalize_composition_vector(shiftedComposition, "shiftedComposition")

    reference_names <- names(reference_freqs)
    shifted_names <- names(shifted_freqs)
    if (!is.null(reference_names) && !is.null(shifted_names)) {
        if (!setequal(reference_names, shifted_names) || anyDuplicated(reference_names) || anyDuplicated(shifted_names)) {
            stop("referenceComposition and shiftedComposition must name the same categories, each once")
        }
        shifted_freqs <- shifted_freqs[reference_names]
    } else if (length(reference_freqs) != length(shifted_freqs)) {
        stop("referenceComposition and shiftedComposition must have the same number of categories")
    }

    bhattacharyya_coefficient <- sum(sqrt(reference_freqs * shifted_freqs))
    # Floating-point rounding can push BC a hair above 1 for identical inputs.
    energy_bound <- max(0, 2 * (1 - bhattacharyya_coefficient))
    if (energy_bound <= .Machine$double.eps * 8) {
        stop("referenceComposition and shiftedComposition describe the same composition, so they give no bound; energyNope needs a number above 0")
    }
    energy_bound
}

#' Normalize one composition vector to frequencies, with input checks.
#'
#' @param composition Non-negative numeric vector of frequencies or counts.
#' @param argument_name Argument name for error messages.
#' @return Numeric vector summing to 1, names kept.
#' @keywords internal
normalize_composition_vector <- function(composition, argument_name) {
    if (!is.numeric(composition) || length(composition) < 2 || anyNA(composition) ||
        any(!is.finite(composition)) || any(composition < 0)) {
        stop(sprintf("%s must be a numeric vector of at least 2 non-negative, finite values", argument_name))
    }
    total <- sum(composition)
    if (total <= 0) {
        stop(sprintf("%s must contain at least one positive value", argument_name))
    }
    composition / total
}
