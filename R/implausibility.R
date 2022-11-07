# Code to check emulator implausibilities sequentially
sequential_imp <- function(ems, x, z, n = 1, cutoff = 3) {
  t_ems <- ems
  results <- purrr::map_lgl(seq_len(nrow(x)), function(a) {
    current_fails <- 0
    i <- 1
    while (i <= length(t_ems)) {
      if (!t_ems[[i]]$implausibility(x[a,],
                                     z[[t_ems[[i]]$output_name]], cutoff)) {
        current_fails <- current_fails + 1
        this_output <- which(
          purrr::map_chr(t_ems, ~.$output_name) == t_ems[[i]]$output_name)
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
#' Computes the nth-maximum implausibility of points relative to a set of emulators.
#'
#' For a collection of emulators, we often combine the implausibility
#' measures for a given set of observations. The maximum implausibility of a point,
#' given a set of univariate emulators and an associated collection of target values,
#' is the largest implausibility of the collected set of implausibilities. The 2nd
#' maximum is the maximum of the set without the largest value, and so on. By default,
#' maximum implausibility will be considered when there are fewer than 10 targets to
#' match to; otherwise second-maximum implausibility is considered.
#'
#' If \code{sequential = TRUE} and a specific \code{cutoff} has been provided, then the
#' emulators' implausibility will be evaluated one emulator at a time. If a point
#' is judged non-implausible by more than \code{n} emulators, \code{FALSE} is
#' returned without evaluating any more. Due to R efficiencies, this is more efficient
#' than the 'evaluate all' method once more than around 10 emulators are considered.
#'
#' This function also deals with variance emulators and bimodal emulators, working in a nested
#' fashion. If targets are provided for both the expectation and variance as a list, then
#' given \code{ems = list(expectation = ..., variance = ...)} the implausibility is calculated
#' with respect to both sets of emulators, maximising as relevant. If targets are provided in
#' the 'normal' fashion, then only the mean emulators are used. The bimodal case is similar;
#' given a set of emulators \code{list(mode1 = list(expectation = ..., variance = ...), ...)}
#' then each mode has implausibility evaluated separately. The results from the two modes are
#' combined via piecewise minimisation.
#'
#' @param ems A set of \code{\link{Emulator}} objects or nested sets thereof (see description)
#' @param x An input point, or \code{data.frame} of points.
#' @param z The target values, in the usual form or nested thereof.
#' @param n The implausibility level to return.
#' @param max_imp A maximum implausibility to return (often used with plotting)
#' @param cutoff A numeric value, or vector of such, representing allowed implausibility
#' @param sequential Should the emulators be evaluated sequentially?
#' @param get_raw Boolean - determines whether nth-implausibility should be applied.
#' @param ordered If FALSE, emulators are ordered according to restrictiveness.
#'
#' @return Either the nth maximum implausibilities, or booleans (if cutoff is given).
#' @export
#'
#' @examples
#' # A single point
#' nth_implausible(SIREmulators$ems, data.frame(aSI = 0.4, aIR = 0.25, aSR = 0.025),
#'  SIREmulators$targets)
#' # A data.frame of points
#' grid <- expand.grid(
#'  aSI = seq(0.1, 0.8, length.out = 4),
#'  aIR = seq(0, 0.5, length.out = 4),
#'  aSR = seq(0, 0.05, length.out = 4)
#' )
#' # Vector of numerics
#' i1 <- nth_implausible(SIREmulators$ems, grid, SIREmulators$targets)
#' # Vector of booleans (same as i1 <= 3)
#' i2 <- nth_implausible(SIREmulators$ems, grid, SIREmulators$targets, cutoff = 3)
#' # Throws a warning as n > no. of targets
#' i3 <- nth_implausible(SIREmulators$ems, grid, SIREmulators$targets, n = 4)
#' # Vector of booleans (note different output to i2)
#' i4 <- nth_implausible(SIREmulators$ems, grid, SIREmulators$targets,
#'  cutoff = c(4, 2.5, 2))
#'
#' # Variance Emulators
#' v_ems <- variance_emulator_from_data(BirthDeath$training, c('Y'),
#'  list(lambda = c(0, 0.08), mu = c(0.04, 0.13)))
#' v_targs = list(expectation = list(Y = c(90, 110)), variance = list(Y = c(55, 95)))
#' nth_implausible(v_ems, unique(BirthDeath$validation[,1:2]), v_targs)
#' ## If there is a mismatch between emulators and targets, expectation is assumed
#' nth_implausible(v_ems$expectation, unique(BirthDeath$validation[,1:2]), v_targs)
#' nth_implausible(v_ems, unique(BirthDeath$validation[,1:2]), v_targs$expectation)
#'
nth_implausible <- function(ems, x, z, n = NULL,
                            max_imp = Inf, cutoff = NULL,
                            sequential = FALSE, get_raw = FALSE,
                            ordered = FALSE) {
  if (!"data.frame" %in% class(x))  {
    if (!is.null(dim(x)) && !is.null(colnames(x)))
      x <- data.frame(x)
    else stop("Named array or data.frame of points required.")
  }
  if (!ordered) ems <- collect_emulators(ems, z)
  ## Preprocessing for variance emulation
  if (!is.null(ems$expectation) && !is.null(ems$variance)) {
    if (is.null(n))
      n <- ifelse(length(unique(purrr::map_chr(
        ems$expectation, ~.$output_name))) > 10, 2, 1)
    if (!is.null(z$expectation) && !is.null(z$variance)) {
      imps_list <- list(
        expectation = nth_implausible(ems$expectation, x,
                                      z$expectation, n, max_imp,
                                      cutoff, FALSE, TRUE),
        variance = nth_implausible(ems$variance, x,
                                   z$variance, n, max_imp,
                                   cutoff, FALSE, TRUE))
      imp_mat <- cbind.data.frame(imps_list$expectation, imps_list$variance)
      if (get_raw) {
        return(setNames(
          imp_mat,
          c(paste0(purrr::map_chr(ems$expectation, ~.$output_name), "Exp"),
            paste0(purrr::map_chr(ems$variance, ~.$output_name), "Var"))))
      }
      else {
        if (n == 1) imps <- apply(imp_mat, 1, max)
        else imps <- apply(imp_mat, 1, function(x) -sort(-x, partial = 1:n)[n])
        return(purrr::map_dbl(imps, ~min(., max_imp)))
      }
    }
    else {
      return(nth_implausible(ems$expectation, x, z, n,
                             max_imp, cutoff, sequential, get_raw))
    }
  }
  else if (!is.null(z$expectation) && !is.null(z$variance)) {
    return(nth_implausible(ems, x, z$expectation, n,
                           max_imp, cutoff, sequential, get_raw))
  }
  else if (!is.null(ems$mode1) && !is.null(ems$mode2)) {
    if (is.null(n))
      n <- ifelse(length(unique(purrr::map_chr(
        ems$mode1$expectation, ~.$output_name))) > 10, 2, 1)
    imps1 <- nth_implausible(ems$mode1, x, z, n, max_imp, cutoff, FALSE, TRUE)
    imps2 <- nth_implausible(ems$mode2, x, z, n, max_imp, cutoff, FALSE, TRUE)
    imps2 <- imps2[,names(imps1)]
    get_min_concrete <- function(v1, v2) {
      if (all(is.numeric(v1))) {
        output <- purrr::map_dbl(seq_along(v1), function(i) {
          if (is.nan(v1[i]) && is.nan(v2[i])) return(0)
          if (is.nan(v1[i])) return(v2[i])
          if (is.nan(v2[i])) return(v1[i])
          return(min(v1[i], v2[i]))
        })
      }
      if (all(is.logical(v1))) {
        output <- v1 | v2
      }
      return(output)
    }
    imp_mat <- as.data.frame(Map(get_min_concrete, imps1, imps2))
    if (get_raw) return(imp_mat)
  }
  else {
    if (is.null(n))
      n <- ifelse(length(unique(purrr::map_chr(ems, ~.$output_name))) > 10, 2, 1)
    for (i in seq_along(z)) {
      if (length(z[[i]]) == 1) {
        warning(paste("Target",
                      names(z)[i],
                      "is a single value. Assuming it's a val with sigma = 5%."))
        z[[i]] <- list(val = z[[i]], sigma = 0.05*z[[i]])
      }
    }
    if ("Emulator" %in% class(ems)) {
      if (!ems$output_name %in% names(z))
        stop("Target not found corresponding to named emulator.")
      else return(ems$implausibility(x, z[[ems$output_name]], cutoff))
    }
    if (n > length(unique(purrr::map_chr(ems, ~.$output_name)))) {
      warning(paste("n cannot be greater than the number of targets to match to.",
                    "Switching to minimum implausibility."))
      n <- length(unique(purrr::map_chr(ems, ~.$output_name)))
    }
    if (length(cutoff) == 1) cutoff <- rep(cutoff, length(ems))
    if (!is.null(cutoff) && (length(ems) > 10 || sequential)) {
      return(sequential_imp(ems, x, z, n, cutoff[1]))
    }
    implausibles <- do.call(
      'cbind', purrr::map(seq_along(ems),
                          ~ems[[.]]$implausibility(x,
                                                   z[[ems[[.]]$output_name]],
                                                   cutoff[[.]])))
    d_implausibles <- setNames(
      data.frame(implausibles),
      purrr::map_chr(ems, ~.$output_name))
    if (!is.null(cutoff))
      imp_mat <- t(apply(
        d_implausibles, 1,
        function(x) purrr::map_lgl(unique(names(x)), ~all(x[names(x) == .]))))
    else
      imp_mat <- t(apply(
        d_implausibles, 1,
        function(x) purrr::map_dbl(unique(names(x)), ~max(x[names(x) == .]))))
  }
  if (length(ems) == 1 ||
      (!is.null(ems$expectation) &&
       length(ems$expectation) == 1)) return(t(imp_mat))
  if (nrow(imp_mat) == 1 && nrow(x) != 1) imp_mat <- t(imp_mat)
  if (get_raw)
    return(setNames(data.frame(imp_mat),
                    unique(purrr::map_chr(ems, ~.$output_name))))
  if (!is.null(cutoff)) {
    return(apply(imp_mat, 1, function(x) sum(x) > length(x)-n))
  }
  if (n == 1) imps <- apply(imp_mat, 1, max)
  else imps <- apply(imp_mat, 1, function(x) -sort(-x, partial = 1:n)[n])
  return(purrr::map_dbl(imps, ~min(., max_imp)))
}

