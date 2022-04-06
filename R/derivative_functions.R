get_deriv_info <- function(em, x, var = FALSE, ...) {
  deriv_exp <- purrr::map_dbl(names(em$ranges),
                              function(y) em$get_exp_d(x, p = y, ...))
  if (var) {
    deriv_var <- matrix(
      unlist
      (purrr::map(names(em$ranges),
                  function(a)
                    purrr::map_dbl(names(em$ranges),
                                   function(b)
                                     em$get_cov_d(x, p1 = a, p2 = b, ...)))),
      nrow = length(em$ranges), byrow = TRUE)
    return(list(exp = deriv_exp, var = deriv_var))
  }
  return(list(exp = deriv_exp))
}

#' Derivative inner product
#'
#' Find the (uncertainty modified) inner product between the derivative at a point \code{x}
#' and a proposed direction \code{v}.
#'
#' Given a point \code{x} and a direction \code{v}, we find the overlap between E[f'(x)] and
#' \code{v}. The emulated derivative has uncertainty associated with it: the variance is taken
#' into account using \eqn{v^{T} Var[f'(x)] v}.
#'
#' If \code{sd == NULL}, then only the (normed) overlap between the derivative and the direction
#' vector is returned. Otherwise a pair of values are returned: these are the normed overlap plus
#' or minus \code{sd} times the uncertainty.
#'
#' This function is concerned with ascertaining whether a direction is oriented in the direction
#' of the emulator gradient, subject to the uncertainty around the estimate of the derivative.
#' It allows for a consideration of "emulated gradient descent".
#'
#' @param em The emulator in question
#' @param x The point in input space to evaluate at
#' @param v The direction to assess
#' @param sd How many standard deviations to consider.
#' @param ... Additional arguments to pass through (eg local.var to the emulator functions)
#'
#' @return Either a single numeric or a pair of numerics (see description)
#'
#' @export
#'
#' @examples
#'  directional_deriv(SIREmulators$ems[[1]], SIRSample$validation[1,], c(1,1,1))
#'
directional_deriv <- function(em, x, v, sd = NULL, ...) {
  normed_v <- v/sqrt(sum(v^2))
  if (!is.null(sd))
    deriv_info <- get_deriv_info(em, x, var = TRUE)
  else
    deriv_info <- get_deriv_info(em, x, ...)
  normed_info <- deriv_info$exp/sqrt(sum(deriv_info$exp^2))
  suitability <- normed_info %*% normed_v
  if (is.null(sd)) return(suitability)
  else {
    uncertainty <- sqrt(normed_v %*% deriv_info$var %*% normed_v)
    return(c(suitability - sd * uncertainty, suitability + sd * uncertainty))
  }
}

#' Emulated Derivative Point Proposal
#'
#' Proposes a new point by applying `emulated gradient descent' on an existing point.
#'
#' Given a point (preferably close to the implausibility boundary) \code{x}, we can calculate
#' the emulated gradient at this point for each emulator. If the estimate of the expectation
#' at this point for a given emulator is larger than the target value, then we would like to
#' move in the direction of greatest decrease for this emulator, and conversely for an estimate
#' of the expectation that's smaller than the target value. The combination of this information
#' for every emulator under consideration defines a preferred set of directions of travel from
#' this point.
#'
#' We may try to find a shared direction which improves (or at least does not worsen) all
#' emulator evaluations. If a point is already well inside the implausibility boundary for a given
#' output (where `well inside' is defined by the value of \code{accept}), we may allow this
#' output to worsen in order to improve the others.
#'
#' Provided a shared direction, v, can be identified, we iteratively move in this direction. Define
#' the new proposed point x' = x + h*v, where h is a step-size given by \code{hstart}. Compare
#' the summary statistic (either expectational difference or implausibility) to that provided by
#' the original point; if the new point gives improvement, then continue to move in this direction
#' until no further improvement is possible for this step-size. The step-size is reduced (up to
#' a minimum of \code{hcutoff}) and the process is repeated. Only finitely many iteration steps
#' are permitted; this can be tuned by supplying a value of \code{iteration.steps}.
#'
#' @param ems The emulators to evaluate with respect to.
#' @param x The original point.
#' @param targets The list of emulator targets.
#' @param accept The implausibility below which we allow an output to worsen.
#' @param hstart The initial step size.
#' @param hcutoff The minimum allowed step size.
#' @param iteration.measure Either `exp' for expectation or `imp' for implausibility.
#' @param iteration.steps The number of allowed iterations.
#' @param nv The number of directions on the n-sphere to try.
#'
#' @return Either a new proposal point, or the original point if an improvement could not be found.
#' @export
#'
#' @examples
#'  # Take a point from the SIR system at later waves with low (but >3) implausibility
#'  start_point <- SIRMultiWaveData[[2]][90,1:3]
#'  ems <- SIRMultiWaveEmulators[[3]]
#'  targs <- SIREmulators$targets
#'  # Using expected error as measure
#'  new_point1 <- directional_proposal(ems, start_point, targs)
#'  # Using implausibility as measure
#'  new_point2 <- directional_proposal(ems, start_point, targs, iteration.measure = 'imp')
#'  all_points <- do.call('rbind.data.frame', list(start_point, new_point1, new_point2))
#'  nth_implausible(ems, all_points, targs)
#'
directional_proposal <- function(ems, x, targets, accept = 2, hstart = 1e-04,
                                 hcutoff = 1e-09, iteration.measure = 'exp',
                                 iteration.steps = 100, nv = 500) {
  if (length(x) > length(ems[[1]]$ranges))
    x <- x[,names(ems[[1]]$ranges)]
  if (any(!names(targets) %in% purrr::map_chr(ems, ~.$output_name))) {
    warning("Not all targets have a corresponding emulator.
            Restricting to only emulated outputs.")
    targets <- targets[names(targets) %in% purrr::map_chr(ems, ~.$output_name)]
  }
  point_implaus <- purrr::map_dbl(seq_along(ems),
                                  ~ems[[.]]$implausibility(x, z = targets[[.]]))
  x_predict <- purrr::map_dbl(ems, ~.$get_exp(x))
  is_bigger <- purrr::map_lgl(seq_along(targets), function(y) {
    if (!is.numeric(targets[[y]]))
      comparative <- targets[[y]]$val < x_predict[[y]]
    else comparative <- targets[[y]][2] < x_predict[[y]]
  })
  x_diffs <- do.call('rbind', purrr::map(ems, ~get_deriv_info(., x)$exp))
  x_dir <- do.call('rbind', purrr::map(seq_along(is_bigger), function(i) {
                                       if(is_bigger[[i]])
                                         return(-1*x_diffs[i,])
                                       else
                                         return(x_diffs[i,])}))
  x_norms <- apply(x_dir, 1, function(y) sqrt(sum(y^2)))
  x_dir <- sweep(x_dir, 1, x_norms, "/")
  test_dirs <- runifs(nv*length(x), length(x))
  suits <- apply(x_dir, 1, function(y) apply(test_dirs, 1, function(z) z %*% y))
  suit_means <- apply(suits, 1, mean)
  order_dirs <- test_dirs[order(suit_means, decreasing = TRUE),]
  order_suits <- suits[order(suit_means, decreasing = TRUE),]
  restrict_dirs <- order_dirs[apply(order_suits, 1,
                                    function(y)
                                      all(y >= 0 | point_implaus < accept)),]
  nth_discrepancy <- function(ems, x, targets, n = 1) {
    discs <- purrr::map_dbl(seq_along(ems), function(y) {
      if (!is.numeric(targets[[y]]))
        return(abs((ems[[y]]$get_exp(x) - targets[[y]]$val)/targets[[y]]$val))
      else {
        emval <- ems[[y]]$get_exp(x)
        if (emval < targets[[y]][1])
          return((targets[[y]][1]-emval)/targets[[y]][1])
        if (emval > targets[[y]][2])
          return((emval - targets[[y]][2])/targets[[y]][2])
        return(0)
      }
    })
    if (n==1) return(max(discs))
    return(order(discs, decreasing = TRUE)[[n]])
  }
  if (nrow(restrict_dirs) == 0) {
    warning("No direction improves the performance on the relevant targets.")
    return(x)
  }
  range_dists <- purrr::map_dbl(ems[[1]]$ranges, ~diff(.)/2)
  best_dir <- restrict_dirs[1,] * range_dists
  track_measure <- if (iteration.measure == "exp")
    nth_discrepancy(ems, x, targets) else max(point_implaus)
  better_point <- NULL
  attempts <- 0
  gap <- 1e-04
  index <- 1
  while (attempts < iteration.steps) {
    new_point <- x + gap * index * best_dir
    new_measure <- if (iteration.measure == "exp")
      nth_discrepancy(ems, new_point, targets)
    else
      nth_implausible(ems, new_point, targets)
    if (new_measure < track_measure) {
      better_point <- new_point
      track_measure <- new_measure
      index <- index + 1
    }
    else {
      if (is.null(better_point) && gap > hcutoff) gap <- gap * 0.1
      else break
    }
    attempts <- attempts + 1
  }
  if (is.null(better_point)) return(x)
  return(better_point)
}
