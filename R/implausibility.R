sequential_imp <- function(ems, x, z, n = 1, cutoff = 3) {
  t_ems <- ems
  results <- purrr::map_lgl(1:nrow(x), function(a) {
    current_fails <- 0
    i <- 1
    while (i <= length(t_ems)) {
      if (!t_ems[[i]]$implausibility(x[a,], z[[t_ems[[i]]$output_name]], cutoff)) {
        current_fails <- current_fails + 1
        this_output <- which(purrr::map_chr(t_ems, ~.$output_name) == t_ems[[i]]$output_name)
        t_ems <- t_ems[-this_output]
        how_many_before <- sum(this_output <= i)
        i <- i - how_many_before
      }
      if (current_fails >= n) break
      i <- i + 1
    }
    return(current_fails < n)
  })
  return(results)
}

#' nth Maximum Implausibility
#'
#' For a collection of emulators, it can be helpful to combine the implausibility
#' measures for a given set of observations. The maximum implausibility of a point,
#' given a set of univariate emulators and an associated collection of target values,
#' is the largest implausibility of the collected set of implasusibilities. The 2nd
#' maximum is the maximum of the set without the largest value, and so on.
#'
#' If \code{sequential = TRUE} and a specific \code{cutoff} has been provided, then the
#' emulators' implausibility will be evaluated one emulator at a time. If a point
#' is judged non-implausible by more than \code{n} emulators, \code{FALSE} is
#' returned without evaluating any more. Due to R efficiencies, this is more efficient
#' than the 'evaluate all' method once more than around 10 emulators are considered.
#'
#' @param ems A set of \code{\link{Emulator}} objects.
#' @param x An input point, or \code{data.frame} of points.
#' @param z The target values.
#' @param n The implausibility level to return.
#' @param max_imp A maximum implausibility to return (often used with plotting)
#' @param cutoff A numeric value, or vector of such, representing allowed implausibility
#' @param sequential Should the emulators be evaluated sequentially?
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
nth_implausible <- function(ems, x, z, n = 1, max_imp = 20, cutoff = NULL, sequential = FALSE) {
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
  if (!is.null(cutoff) && (length(ems) > 10 || sequential))
    return(sequential_imp(ems, x, z, n, cutoff[1]))
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

