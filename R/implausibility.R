#' nth Maximum Implausibility
#'
#' For a collection of emulators, it can be helpful to combine the implausibility
#' measures for a given set of observations. The maximum implausibility of a point,
#' given a set of univariate emulators and an associated collection of target values,
#' is the largest implausibility of the collected set of implasusibilities. The 2nd
#' maximum is the maximum of the set without the largest value, and so on.
#'
#' @param ems A set of \code{\link{Emulator}} objects.
#' @param x An input point, or \code{data.frame} of points.
#' @param z The target values.
#' @param n The implausibility level to return.
#' @param max_imp A maximum implausibility to return (often used with plotting)
#' @param cutoff A numeric value, or vector of such, representing allowed implausibility
#'
#' @return Either the nth maximum implausibilities, or booleans (if cutoff is given).
#' @export
#'
#' @examples
#' # A single point
#' nth_implausible(sample_emulators$ems, data.frame(aSI = 0.4, aIR = 0.25, aSR = 0.025),
#'  sample_emulators$targets)
#' # A data.frame of points
#' grid <- expand.grid(
#'  aSI = seq(0.1, 0.8, length.out = 4),
#'  aIR = seq(0, 0.5, length.out = 4),
#'  aSR = seq(0, 0.05, length.out = 4)
#' )
#' # Vector of numerics
#' i1 <- nth_implausible(sample_emulators$ems, grid, sample_emulators$targets)
#' # Vector of booleans (same as i1 <= 3)
#' i2 <- nth_implausible(sample_emulators$ems, grid, sample_emulators$targets, cutoff = 3)
#' # Throws a warning as n > no. of targets
#' i3 <- nth_implausible(sample_emulators$ems, grid, sample_emulators$targets, n = 4)
#' # Vector of booleans (note different output to i2)
#' i4 <- nth_implausible(sample_emulators$ems, grid, sample_emulators$targets,
#'  cutoff = c(4, 2.5, 1))
nth_implausible <- function(ems, x, z, n = 1, max_imp = 20, cutoff = NULL) {
  if ("Emulator" %in% class(ems)) {
    if (!ems$output_name %in% names(z)) stop("Target not found corresponding to named emulator.")
    else return(ems$implausibility(x, z[[ems$output_name]], cutoff))
  }
  for (i in 1:length(z)) {
    if (length(z[[i]]) == 1) {
      warning(paste("Target", names(z)[i], "is a single value. Assuming it's a value with sigma = 0.01."))
      z[[i]] <- list(val = z[[i]], sigma = 0.01)
    }
  }
  if (n > length(unique(purrr::map_chr(ems, ~.$output_name)))) {
    warning("n cannot be greater than the number of targets to match to. Switching to minimum implausibility.")
    n <- length(unique(purrr::map_chr(ems, ~.$output_name)))
  }
  if (length(cutoff) == 1) cutoff <- rep(cutoff, length(ems))
  implausibles <- do.call('cbind', purrr::map(seq_along(ems), ~ems[[.]]$implausibility(x, z[[ems[[.]]$output_name]], cutoff[[.]])))
  d_implausibles <- setNames(data.frame(implausibles), purrr::map_chr(ems, ~.$output_name))
  if (!is.null(cutoff)) {
    implausibles <- t(apply(d_implausibles, 1, function(x) purrr::map_lgl(unique(names(x)), ~all(x[names(x) == .]))))
    if (length(ems) == 1) return(t(implausibles))
    return(apply(implausibles, 1, function(x) sum(x) > length(unique(purrr::map_chr(ems, ~.$output_name)))-n))
  }
  if (length(ems) == 1) return(t(implausibles))
  implausibles <- t(apply(d_implausibles, 1, function(x) purrr::map_dbl(unique(names(x)), ~max(x[names(x) == .]))))
  if (nrow(implausibles) == 1 && nrow(x) != 1) implausibles <- t(implausibles)
  if (n == 1) imps <- apply(implausibles, 1, max)
  else imps <- apply(implausibles, 1, function(x) -sort(-x, partial = 1:n)[n])
  return(purrr::map_dbl(imps, ~min(., max_imp)))
}

# nth_implausible <- function(ems, x, z, n = 1, max_imp = 20, cutoff = NULL) {
#   if ("Emulator" %in% class(ems)) {
#     if (is.null(z$val)) z <- z[[ems$output_name]]
#     return(ems$implausibility(x, z, cutoff))
#   }
#   if (n > length(ems)) {
#     warning("n cannot be greater than the number of targets to match to. Switching to minimum implausibility.")
#     n <- length(ems)
#   }
#   outputs <- if (is.numeric(z)) purrr::map(z, ~list(val = ., sigma = 0.01)) else z
#   if (length(cutoff) == 1) cutoff <- rep(cutoff, length(ems))
#   implausibles <- do.call('cbind', purrr::map(seq_along(ems), ~ems[[.]]$implausibility(x, outputs[[ems[[.]]$output_name]], cutoff[[.]])))
#   d_implausibles <- setNames(data.frame(implausibles), purrr::map_chr(ems, ~.$output_name))
#   if (!is.null(cutoff))
#     implausibles <- t(apply(d_implausibles, 1, function(x) purrr::map_dbl(unique(names(x)), ~all(x[names(x) == .]))))
#   else implausibles <- t(apply(d_implausibles, 1, function(x) purrr::map_dbl(unique(names(x)), ~max(x[names(x) == .]))))
#   if (!is.null(cutoff))
#     return(apply(implausibles, 1, function(x) sum(x) > length(unique(purrr::map_chr(ems, ~.$output_name)))-n))
#   if (n == 1)
#     imps <- apply(implausibles, 1, max)
#   else
#     imps <- apply(implausibles, 1, function(x) -sort(-x, partial = 1:n)[n])
#   return(purrr::map_dbl(imps, ~min(., max_imp)))
# }

