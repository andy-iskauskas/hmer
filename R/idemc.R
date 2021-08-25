#' Crossover step
#'
#' Takes two points, chooses a pivot along the points, and exchanges
#' the values after the pivot. Concretely, given x = (x1, x2, ..., xd) and
#' y = (y1, y2, ..., yd), a value k in 1:d is picked, and the resulting points
#' are (x1, x2, ..., xk, y(k+1), ... yd) and (y1, y2, ..., yk, x(k+1), ... xd).
#' If the resulting points are non-implausible with respect to their ladder rungs,
#' they are kept.
#'
#' @param x The current set of points, as a list of data.frames
#' @param imp_func An implausibility function (defined in \code{IDEMC})
#' @param imps The list of implausibilities for the current ladder
#'
#' @keywords internal
#' @noRd
#'
#' @return A list of data.frames containing the resulting points.
crossover <- function(x, imp_func, imps) {
  n <- length(x)
  indices <- sample(n, 2)
  xi <- x[[indices[1]]][sample(nrow(x[[indices[1]]]), 1),]
  xj <- x[[indices[2]]][sample(nrow(x[[indices[2]]]), 1),]
  d <- length(names(xi))
  c_point <- sample(d-1, 1)
  yi <- unlist(c(xi[,1:c_point], xj[,(c_point+1):d]), use.names = FALSE)
  yj <- unlist(c(xj[,1:c_point], xi[,(c_point+1):d]), use.names = FALSE)
  if (imp_func(setNames(data.frame(matrix(yi, nrow = 1)), names(xi)), imps[indices[1]]) && imp_func(setNames(data.frame(matrix(yj, nrow = 1)), names(xj)), imps[indices[2]])) {
    x[[indices[1]]] <- rbind(x[[indices[1]]], yi)
    x[[indices[2]]] <- rbind(x[[indices[2]]], yj)
  }
  return(x)
}

#' Exchange step
#'
#' Takes two points from neighbouring rungs of the implausibility
#' ladder and swaps them. The exchange is kept if the point moved into
#' a lower implausibility ladder satisfies the new, more stringent, implausibility
#' constraint.
#'
#' @param x The current set of points, as a list of data.frames
#' @param imp_func An implausibility function (defined in \code{IDEMC})
#' @param imps The list of implausibilities for the ladder
#'
#' @keywords internal
#' @noRd
#'
#' @return A list of data.frames containing the resulting points
exchange <- function(x, imp_func, imps) {
  n <- length(x)
  a <- sample(n, 1)
  if (a != n && a != 1) {
    if (runif(1) < 0.5) {
      b <- a
      a <- b - 1
    }
    else b <- a + 1
  }
  else if (a == n) {
    b <- n
    a <- n - 1
  }
  else b <- 2
  samp1 <- sample(nrow(x[[a]]),1)
  samp2 <- sample(nrow(x[[b]]),1)
  xa <- x[[a]][samp1,]
  xb <- x[[b]][samp2,]
  if (imp_func(xa, imps[b])) {
    x[[a]][samp1,] <- xb
    x[[b]][samp2,] <- xa
  }
  return(x)
}

##### Need to think about a more robust mutation method. Maybe take lessons from importance_sample?
#' Mutation step
#'
#' Takes a point and mutates it using Metropolis-Hastings. The distribution with which
#' the point is mutated is a weighted uniform spherical distribution.
#'
#' @param ems The emulators in question
#' @param x The current points, as a list of data.frames
#' @param imp_func An implausibility function (defined in \code{IDEMC})
#' @param imps The list of implausibilities for the ladder
#' @param index The index of the ladder rung to mutate from. Default: \code{NULL}
#'
#' @keywords internal
#' @noRd
#'
#' @return A list of data.frames with the resulting points.
mutate <- function(ems, x, imp_func, imps, index) {
  pts <- x[[index]]
  ranges <- ems[[1]]$ranges
  if (index == 1) {
    y <- purrr::map_dbl(ranges, ~runif(1, .[1], .[2]))
    x[[index]] <- rbind(x[[index]], y)
    return(x)
  }
  sampled_point <- sample(nrow(pts), 1)
  mutate_point <- unlist(pts[sampled_point,], use.names = FALSE)
  imp <- imps[index]
  y <- runifs(1, length(mutate_point), mutate_point, r = purrr::map_dbl(ranges, ~.[[2]]-.[[1]])/4)
  if (imp_func(setNames(data.frame(matrix(y, nrow = 1)), names(pts)), imps[index])) {
    p_weight <- sum(apply(pts, 1, function(a) punifs(y, a, purrr::map_dbl(ranges, ~.[[2]]-.[[1]])/4)))
    inv_pts <- rbind(pts, y)
    inv_p_weight <- sum(apply(inv_pts, 1, function(a) punifs(mutate_point, a, purrr::map_dbl(ranges, ~.[[2]]-.[[1]])/8)))-1
    if (runif(1) < (inv_p_weight/p_weight)) {
      x[[index]] <- rbind(x[[index]], y)
    }
  }
  return(x)
}


#' IDEMC step
#'
#' Performs a single IDEMC step: given a set of points in an implausibility ladder,
#' it performs some number of mutations, crossovers, and exchanges. Continues until
#' all rungs of the ladder have at least \code{s} points in them.
#'
#' @param ems The emulators in question
#' @param points The current set of points
#' @param imp_func An implausibility function (defined in \code{IDEMC})
#' @param imp_levels The implausibility ladder levels
#' @param s The desired number of points per rung
#' @param pm The probability of mutation rather than crossover. Default: 0.9
#' @param M The number of mutations to perform each time. Default: 10
#'
#' @keywords internal
#' @noRd
#'
#' @return A list of data.frames corresponding to the ladder points.
IDEMC_step <- function(ems, points, imp_func, imp_levels, s, pm = 0.9, M = 10) {
  while(any(purrr::map(points, nrow) < s)) {
    r <- runif(1)
    if (r < pm) {
      for (i in 1:M) {
        for (j in 1:length(points)) points <- mutate(ems, points, imp_func, imp_levels, index = j)
      }
    }
    else {
      for (j in 1:floor((length(points)+1)/2)) points <- crossover(points, imp_func, imp_levels)
    }
    for (i in 1:(length(points)+1)) points <- exchange(points, imp_func, imp_levels)
  }
  points <- purrr::map(points, ~.[sample(nrow(.), s),])
}

#' IDEMC Point Generation
#'
#' Performs Implausibility-Driven Evolutionary Monte Carlo.
#'
#' Given a set of initial points (preferably space-filling across the space in question),
#' the implausibility ladder is set up via a burn-in phase. Once the ladder of implausibilities
#' has been determined, these are used to generate the full set of points.
#'
#' This is a very computationally expensive procedure for generating points, and should only be
#' used when it is strongly suspected that the target region is extremely small or has an
#' interesting disconnected structure. For more mundane point generation, a more suitable
#' approach is the functionality in \code{\link{generate_new_runs}}.
#'
#' The burn-in phase starts with a rung defined over the full space (i.e. the implausibility
#' for this rung is simply the maximum implausibility over the space). The implausibility of
#' the next rung is chosen to be such that 30% of the previous points are below this threshold.
#' A full complement of points are generated at this new level, and the process is repeated
#' using these new points to find the next rung of the ladder. This iterates until the desired
#' final implausibility, \code{imp} has been reached as a ladder rung.
#'
#' Once the burn in has been performed, a full set of points are generated at each rung of the
#' ladder (where it is assumed that the burn-in generated fewer points than required, so that
#' \code{sn}<\code{s}).
#'
#' @param xsamp The initial sample of points, as a data.frame
#' @param ems The emulators with which to evaluate implausibility
#' @param targets The corresponding output targets
#' @param s The number of points to generate in the burn-in phase
#' @param sn The final number of points to generate
#' @param p The proportion of points to keep in each new ladder rung
#' @param imp The value of implausibility that is ultimately desired
#' @param verbose Should logging messages be outputted? Default: F
#' @param ... Any additional parameters to pass to \code{IDEMC_step}
#'
#' @return A list of data.frames, corresponding to the points generated.
#'
#' @seealso \code{\link{generate_new_runs}} for other point generation methods.
#'
#' @references
#' Vernon, I. & Williamson, D. (2013) Efficient uniform designs for multi-wave computer experiments. arXiv:1309.3520
#'
#' @examples
#' \dontrun{
#' ranges <- list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05))
#' out_vars <- c('nS', 'nI', 'nR')
#' ems <- emulator_from_data(GillespieSIR, out_vars, ranges)
#' z <- list(
#'  nS = list(val = 281, sigma = 10.43),
#'  nI = list(val = 30, sigma = 11.16),
#'  nR = list(val = 689, sigma = 14.32)
#' )
#' start_pts <- data.frame(
#'   aSI = runif(500, ranges$aSI[1], ranges$aSI[2]),
#'   aIR = runif(500, ranges$aIR[1], ranges$aIR[2]),
#'   aSR = runif(500, ranges$aSR[1], ranges$aSR[2])
#' )
#' result <- IDEMC(start_pts, ems, z, 50, 100, 0.3, imp = 2)
#' }
#'
#' @export
#'
IDEMC <- function(xsamp, ems, targets, s, sn, p, imp = 3, verbose = F, ...) {
  xsamp <- xsamp[,names(ems[[1]]$ranges)]
  sample_imps <- nth_implausible(ems, xsamp, targets, max_imp = Inf)
  test_i <- max(sample_imps)
  range_func <- function(x, ranges) {
    all(purrr::map_lgl(seq_along(ranges), ~x[.]>=ranges[[.]][1] && x[.] <= ranges[[.]][2]))
  }
  check_imp <- function(x, imp) {
    if (!range_func(x, ems[[1]]$ranges)) return(FALSE)
    for (i in 1:length(ems)) if (!ems[[i]]$implausibility(x, targets[[ems[[i]]$output_name]], imp)) return(FALSE)
    return(TRUE)
  }
  imps <- c(test_i)
  if (verbose) print(test_i)
  samp_with_imps <- setNames(cbind(xsamp, sample_imps), c(names(xsamp), "I"))
  samp_with_imps <- samp_with_imps[order(samp_with_imps$I),]
  cutoff_index <- floor(p*nrow(samp_with_imps)) + 1
  test_i <- samp_with_imps[cutoff_index, "I"]
  xsub <- samp_with_imps[1:cutoff_index, names(xsamp)]
  if (nrow(xsub) > s) xsub <- xsub[sample(nrow(xsub), s),]
  chroms <- list(xsamp, xsub)
  imps <- c(imps, test_i)
  while(test_i > imp) {
    if (verbose) print(test_i)
    new_sample <- IDEMC_step(ems, chroms, check_imp, imps, s, ...)
    new_sample <- new_sample[[length(new_sample)]]
    sample_imps <- nth_implausible(ems, new_sample, targets, max_imp = Inf)
    samp_with_imps <- setNames(cbind(new_sample, sample_imps), c(names(new_sample), "I"))
    samp_with_imps <- samp_with_imps[order(samp_with_imps$I),]
    cutoff_index <- floor(p*nrow(samp_with_imps)) + 1
    test_i <- samp_with_imps[cutoff_index, "I"]
    if (test_i < imp) test_i <- imp
    xsub <- samp_with_imps[1:cutoff_index, names(new_sample)]
    chroms[[length(chroms)+1]] <- xsub
    imps <- c(imps, test_i)
  }
  cat("Finished burn-in. Implausibility ladder:\n")
  cat(imps)
  cat("\n")
  final_out <- IDEMC_step(ems, chroms, check_imp, imps, sn, ...)
  for (i in 1:length(final_out)) row.names(final_out[[i]]) <- NULL
  return(final_out)
}
