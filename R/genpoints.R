## Helper functions for uniform sphere sampling
#'
#' @importFrom stats rexp rnorm
runifs <- function(n, d, c = rep(0, d), r = 1) {
  init_points <- matrix(rnorm(n*d, 0, 1), nrow = n, byrow = T)
  exps <- rexp(n, 0.5)
  denoms <- (exps + apply(init_points, 1, function(x) sum(x^2)))^(0.5)
  swept <- sweep(init_points, 1, denoms, "/")
  return(t(apply(swept, 1, function(x) x*r + c)))
}
punifs <- function(x, c = rep(0, length(x)), r = 1) {
  return(ifelse(sum((x-c)^2 / r^2) <= 1, 1, 0))
}
pca_transform <- function(x, s_points, forward = TRUE) {
  if ("data.frame" %in% class(s_points)) s_points <- data.matrix(s_points)
  if ("data.frame" %in% class(x)) x <- data.matrix(x)
  s_trafo <- sweep(
    sweep(
      s_points, 2,
      apply(s_points, 2, mean), "-"), 2, apply(s_points, 2, sd), "/")

  s_estruct <- eigen(cov(s_trafo))
  s_estruct$values <- purrr::map_dbl(
    s_estruct$values, function(x) {if(x < 1e-10) 1e-10 else x})
  if (forward) x <- sweep(
    sweep(
      x, 2,
      apply(s_points, 2, mean), "-"), 2, apply(s_points, 2, sd), "/")
  if (forward) return (x %*% s_estruct$vectors %*% diag(1/sqrt(s_estruct$values)))
  pre_traf <- x %*% diag(sqrt(s_estruct$values)) %*% t(s_estruct$vectors)
  return(sweep(sweep(pre_traf, 2, apply(s_points, 2, sd), "*"), 2, apply(s_points, 2, mean), "+"))
}
maximin_sample <- function(points, n, reps = 1000, nms) {
  c_measure <- op <- NULL
  options <- purrr::map(1:reps, function(rep) {
    tp <- points[sample(nrow(points), n),]
    if (!"data.frame" %in% class(tp)) tp <- setNames(tp, nms)
    measure <- min(dist(tp))
    return(list(points = tp, value = measure))
  })
  return(options[[which.max(purrr::map_dbl(options, "value"))]]$points)
  for (i in 1:reps) {
    tp <- points[sample(nrow(points), n),]
    if (!"data.frame" %in% class(tp)) tp <- setNames(tp, nms)
    measure <- min(dist(tp))
    if (is.null(c_measure) || measure > c_measure) {
      op <- tp
      c_measure <- measure
    }
  }
  return(op)
}

#' Generate Proposal Points
#'
#' Given a set of trained emulators, this finds the next set of points that will be
#' informative for a subsequent wave of emulation.
#'
#' If the method is \code{'lhs'}, a Latin hypercube is generated and non-implausible
#' points from this design are retained. If enough points are accepted, the points
#' outputted are chosen using either a maximin or V-optimality criteria (chosen by
#' \code{measure.method}).
#'
#' The methods \code{'line'} and \code{'importance'} both require a predetermined set
#' of non-implausible points \code{s_points}; if they are not provided then lhs sampling
#' is performed first.
#'
#' The method \code{'line'} performs line sampling for boundary detection. Given a set of
#' non-implausible points, rays are drawn between pairs of points (selected so as to
#' maximise the distance between them), and more points are sampled along these lines.
#' Points are kept if they lie near a boundary of the non-implausible space.
#'
#' The method \code{'importance'} performs importance sampling, using a mixture distribution
#' of multivariate normal or uniform spherical proposals around the current non-implausible
#' points (determined by \code{distro}). The optimal standard deviation or radius is found
#' using a burn-in phase, before the full set of points is generated.
#'
#' The method \code{'slice'} performs slice sampling. Given one known non-implausible
#' point, it attempts to find a minimum enclosing hyperrectangle for the non-implausible
#' region around this point and samples uniformly from this, shrinking the hyperrectangle
#' as appropriate. This method is also called if LH sampling has only generated one point
#' (since any later methods require at least two points to be useful).
#'
#' The method \code{'optical'} uses optical depth sampling: given a set of known
#' non-implausible points, the an approximation of the one-dimensional marginal distributions
#' in each parameter direction can be determined. From these derived marginals, points are
#' sampled.
#'
#' For any sampling strategy, the parameters \code{ems} and \code{z} must be provided.
#' All of these methods depend on a means of assessing point suitability (henceforth
#' referred to as 'implausibility'). By default, this is uses nth-maximum implausibility
#' as provided by \code{\link{nth_implausible}}: a user-defined method can be substituted
#' for this by supplying the function call to \code{imp_func}. As a minimum, the function
#' should take four arguments: the emulators, the points, and the targets, and a \code{...}
#' argument to allow compatibility with the standard behaviour.
#'
#' The option \code{seek} determines how many points should be chosen that have a higher
#' probability of matching targets, as opposed to not missing targets. Due to the danger of
#' such an approach in terms of obtaining a representative space-filling design over the
#' space, this value should not be too high: a rough guide is that it should be no larger
#' than 10\% of the desired number of points. The default is \code{seek = 0}. If \code{seek}
#' is in the range [0, 1], it is assumed to represent a proportion of the total number of
#' points that should be sought using this method; otherwise, it is assumed to be the number
#' of points that are desired.
#'
#' The default behaviour is as follows. A set of initial points are generated from an LHD;
#' line sampling is performed to find the boundaries; and finally this collection of points
#' is augmented to the desired number of points by importance sampling using uniform
#' spherical proposals.
#'
#' In regions where the non-implausible space (at the given cutoff) is very hard to find,
#' the function will start at a higher implausibility where it can find a space-filling
#' design; using this as a starting point any other methods are performed. From this new
#' proposal, a subset of lower-implausibility points are selected. This process iterates
#' until either the desired implausibility has been reached or the process has reached a
#' barrier to further reductions in implausibility. The argument \code{c_tol} is used to
#' determine if the improvement in implausibility is small enough to justify stopping the
#' process; by default this is 0.1. The process will also stop if it has produced points
#' `close to' the desired implausibility: the level of closeness is defined using the
#' \code{i_tol} argument.
#'
#' These methods may not work, or may work slowly, if the target space is very small
#' compared to the current not-yet-ruled-out space, or it may miss small disconnected
#' regions of parameter space.
#'
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom stats setNames runif dist cov
#' @importFrom utils write.csv
#' @importFrom tidyr pivot_longer
#'
#' @param ems A list of \code{\link{Emulator}} objects, trained on previous design points.
#' @param n_points The desired number of points to propose for the next wave.
#' @param z The targets to match to.
#' @param method The method(s) to use.
#' @param cutoff The implausibility cutoff(s) to compare outputs to.
#' @param nth A parameter to be passed to the \code{n} argument of \code{\link{nth_implausible}}.
#' @param plausible_set An optional set of known non-implausible points (for eg line sampling).
#' @param verbose Should progress statements be printed to the console?
#' @param cluster Should emulator clustering be considered in the LHS generation?
#' @param resample Number of times to resample using line and/or importance sampling.
#' @param seek How many `good' points to search for
#' @param c_tol The tolerance with which to determine that best implausibility has been reached.
#' @param i_tol The tolerance on final desired implausibility
#' @param to_file The filename to write to sequentially during proposal. Default is NULL (no writing)
#' @param imp_func The implausibility measure to use to determine point acceptance.
#' @param ... Any parameters to pass to individual sampling functions, eg \code{distro} for importance sampling.
#'
#' @return A data.frame containing the set of new points to run the model at.
#'
#' @export
#'
#' @examples
#' \donttest{
#'  # A simple example that uses a number of the native and ... parameter options
#'  pts <- generate_new_runs(SIREmulators$ems, 100, SIREmulators$targets,
#'  measure.method = 'maximin', distro = 'sphere', resample = 0)
#'  pts_optical <- generate_new_runs(SIREmulators$ems, 100, SIREmulators$targets,
#'   method = c('optical'))
#'  pts_slice <- generate_new_runs(SIREmulators$ems, 100, SIREmulators$targets,
#'   method = c('slice'))
#'  pts_no_importance <- generate_new_runs(SIREmulators$ems, 100, SIREmulators$targets,
#'   method = c('line'))
#' }
generate_new_runs <- function(ems, n_points, z,
                              method = 'default',
                              cutoff = 3,
                              nth = NULL,
                              plausible_set, verbose = interactive(),
                              cluster = FALSE, resample = 1, seek = 0,
                              c_tol = 0.5, i_tol = 0.01, to_file = NULL,
                              imp_func = function(ems, x, z, ...) nth_implausible(ems, x, z, ...),
                              ...) {
  lhd_pca <- FALSE
  if (!is.null(to_file)) { #nocov start
    tryCatch(
      write.csv(data.frame(), file = to_file, row.names = FALSE),
      error = function(e) {
        warning("Cannot open directory provided in to_file; output will not be saved.")
        to_file <<- NULL
      }
    )
  } #nocov end
  ems <- collect_emulators(ems, z)
  ranges <- getRanges(ems)
  if (is.null(nth)) {
    if (!is.null(ems$expectation))
      nems <- length(unique(purrr::map_chr(
        ems$expectation, ~.$output_name
      )))
    else if (!is.null(ems$mode1))
      nems <- length(unique(purrr::map_chr(
        ems$mode1$expectation, ~.$output_name
      )))
    else
      nems <- length(unique(purrr::map_chr(
        ems, ~.$output_name
      )))
    nth <- ifelse(nems > 10, 2, 1)
  }
  if (length(method) == 1 && method == "default") {
    method <- c('lhs', 'line', 'importance')
    user_select <- FALSE
  }
  else user_select <- TRUE
  possible_methods <- c('lhs', 'line', 'importance', 'slice', 'optical')
  which_methods <- possible_methods[possible_methods %in% method]
  n_current <- 0
  if (any(!method %in% possible_methods))
    warning(paste("Unrecognised method name(s)",
                  method[!method %in% possible_methods], "ignored."))
  if (missing(plausible_set) || 'lhs' %in% which_methods) {
    if (verbose) cat("Proposing from LHS...\n") #nocov
    if (!cluster) {
      lh_gen <- lhs_gen(ems, ranges, max(n_points, 10*length(ranges)),
                        z, cutoff, nth, use_pca = lhd_pca, imp_func = imp_func, ...)
      points <- lh_gen$points
      this_cutoff <- lh_gen$cutoff
    }
    else {
      recent_ems <- ems[!duplicated(purrr::map_chr(ems, ~.$output_name))]
      cluster_gen <- lhs_gen_cluster(recent_ems, ranges,
                                     max(n_points, 10*length(ranges)),
                                     z, cutoff, nth, verbose,
                                     c_tol = c_tol, imp_func = imp_func, ...)
      if (length(recent_ems) != length(ems)) {
        leftover_imps <- imp_func(
          ems[duplicated(purrr::map_chr(ems, ~.$output_name))],
          cluster_gen$points, z, n = nth, ordered = TRUE)
        this_cutoff <- max(cluster_gen$cutoff,
                           sort(leftover_imps)[5*length(ranges)])
        points <- cluster_gen$points[leftover_imps <= this_cutoff,]
      }
      else {
        points <- cluster_gen$points
        this_cutoff <- cluster_gen$cutoff
      }
    }
    maybe_verbose <- TRUE
    if (this_cutoff == cutoff && nrow(points) >= 0.25 * n_points && !user_select) {
      if (verbose) cat("LHS has high yield - no other methods required.\n") #nocov
      maybe_verbose <- FALSE
      all_points <- points
      while(nrow(all_points) < n_points) {
        if (verbose) cat("Proposing from LHS...\n") #nocov
        if (!cluster) {
          lh_gen <- lhs_gen(ems, ranges, max(n_points, 10*length(ranges)),
                            z, cutoff, nth, use_pca = lhd_pca, imp_func = imp_func, ...)
          points <- lh_gen$points
          this_cutoff <- lh_gen$cutoff
        }
        else {
          recent_ems <- ems[!duplicated(purrr::map_chr(ems, ~.$output_name))]
          cluster_gen <- lhs_gen_cluster(recent_ems, ranges,
                                         max(n_points, 10*length(ranges)),
                                         z, cutoff, nth, verbose,
                                         c_tol = c_tol, imp_func = imp_func, ...)
          if (length(recent_ems) != length(ems)) {
            leftover_imps <- imp_func(
              ems[duplicated(purrr::map_chr(ems, ~.$output_name))],
              cluster_gen$points, z, n = nth, ordered = TRUE)
            this_cutoff <- max(cluster_gen$cutoff,
                               sort(leftover_imps)[5*length(ranges)])
            points <- cluster_gen$points[leftover_imps <= this_cutoff,]
          }
          else {
            points <- cluster_gen$points
            this_cutoff <- cluster_gen$cutoff
          }
        }
        all_points <- rbind(all_points, points)
      }
      points <- all_points
    }
    if (nrow(points) >= n_points && this_cutoff == cutoff) {
      if (verbose & maybe_verbose) cat("Enough points generated from LHD - no need to apply other methods.\n") #nocov
      if (seek > 0) {
        if (verbose) cat("Searching for high-probability match points...\n") #nocov
        if (seek <= 1) seek <- floor(n_points*seek)
        extra_points <- seek_good(ems, seek, z, points, cutoff = cutoff, imp_func = imp_func, ...)
      }
      else extra_points <- NULL
      if (nrow(points) > n_points - seek) {
        if (verbose) cat("Selecting final points using maximin criterion...\n") #nocov
        points <- maximin_sample(points, n_points-seek, nms = names)
      }
      return(rbind(extra_points, points))
    }
  }
  else {
    point_imps <- imp_func(ems, plausible_set, z, n = nth, max_imp = Inf, ordered = TRUE)
    optimal_cut <- sort(point_imps)[min(length(point_imps)-1, floor(0.8*length(point_imps)), 5*length(ranges))]
    if (optimal_cut > cutoff && (optimal_cut - sort(point_imps)[1] < c_tol)) {
      if (verbose) cat("Point proposal seems to be asymptoting around implausibility", #nocov
                          round(optimal_cut, 3), "- terminating.\n") #nocov
      if (!is.null(ems$expectation)) {
        recent_ems <- ems$expectation[!duplicated(purrr::map_chr(ems$expectation, ~.$output_name))]
      }
      else
        recent_ems <- ems[!duplicated(purrr::map_chr(ems, ~.$output_name))]
      recent_imps <- do.call(
        'cbind.data.frame',
        purrr::map(
          recent_ems,
          ~.$implausibility(plausible_set, z[[.$output_name]])))
      recent_exps <- do.call(
        'cbind.data.frame',
        purrr::map(
          recent_ems, ~.$get_exp(plausible_set)))
      preflight(cbind(plausible_set, recent_exps), z)
      name <- value <- NULL
      plot_imps <- tidyr::pivot_longer(recent_imps, cols = everything())
      plot_imps$name <- factor(plot_imps$name, levels = names(recent_ems))
      if (verbose) { #nocov start
        print(ggplot(data = plot_imps, aes(x = name, y = value)) +
          geom_boxplot() +
            labs(title = "Implausibility Boxplot",
              x = "Output", y = "Implausibility"))
        cat("Inspect implausibility boxplot for problematic outputs,",
                  "and consider transforming them or",
                  "removing them from this wave.\n")
      } #nocov end
      if (!is.null(to_file)) #nocov start
        write.csv(plausible_set[point_imps <= optimal_cut,],
                  file = to_file, row.names = FALSE) #nocov end
      return(list(points = plausible_set[point_imps <= optimal_cut,],
                  cutoff = optimal_cut))
    }
    if (optimal_cut < cutoff) this_cutoff <- cutoff
    else this_cutoff <- round(optimal_cut, 3)
    points <- plausible_set[point_imps <= this_cutoff,]
  }
  if (length(ranges) == 1)
    points <- setNames(data.frame(temp = points), names(ranges))
  n_current <- nrow(points)
  if (is.null(nrow(points)) || nrow(points) == 0) {
    warning("No non-implausible points found from initial step.")
    return(points)
  }
  if (verbose) cat(n_current, " initial valid points generated for I=", #nocov
                            round(this_cutoff, 3), "\n", sep = "") #nocov
  if (!is.null(to_file)) write.csv(points, file = to_file, row.names = FALSE)
  if ("optical" %in% which_methods && nrow(points) < n_points) {
    if (verbose) cat("Performing optical depth sampling...\n") #nocov
    points <- op_depth_gen(ems, ranges, n_points, z, cutoff = this_cutoff,
                           nth = nth, plausible_set = points,
                           verbose = verbose, imp_func = imp_func, ...)
    if (verbose) cat("Optical depth sampling generated", #nocov
                             nrow(points)-n_current, "more points.\n") #nocov
    n_current <- nrow(points)
  }
  if (!is.null(to_file)) write.csv(points, file = to_file, row.names = FALSE)
  if ("line" %in% which_methods && nrow(points) < n_points) {
    if (verbose) cat("Performing line sampling...\n") #nocov
    points <- line_sample(ems, ranges, z, points,
                          cutoff = this_cutoff, nth = nth, imp_func = imp_func, ...)
    if (verbose) cat("Line sampling generated", #nocov
                        nrow(points)-n_current, "more points.\n") #nocov
    n_current <- nrow(points)
  }
  if ("slice" %in% which_methods) {
    if (verbose) cat("Performing slice sampling...\n") #nocov
    spoints <- slice_gen(ems, ranges, n_points, z, points,
                         this_cutoff, nth, imp_func = imp_func, ...)
    if (verbose) cat("Slice sampling generated", #nocov
                        nrow(spoints)-nrow(points), "more points.\n") #nocov
    points <- spoints
    n_current <- nrow(points)
  }
  if (!is.null(to_file)) write.csv(points, file = to_file, row.names = FALSE) #nocov
  if ("importance" %in% which_methods && nrow(points) < n_points) {
    if (verbose) cat("Performing importance sampling...\n") #nocov
    points <- importance_sample(ems, n_points, z, points,
                                this_cutoff, nth, to_file = to_file, imp_func = imp_func, ...)
    if (verbose) cat("Importance sampling generated", #nocov
                          nrow(points)-n_current, "more points.\n") #nocov
    n_current <- nrow(points)
  }
  if (this_cutoff - cutoff > i_tol) {
    if (!is.null(to_file)) write.csv(points, file = to_file, row.names = FALSE)
    if (length(which_methods) == 1 && which_methods == c('lhs')) {
      warning("Could not generate points to desired target implausibility.")
      return(points)
    }
    points <- generate_new_runs(ems, n_points, z,
                                which_methods[!which_methods %in% c('lhs')],
                                cutoff = cutoff, nth = nth,
                                plausible_set = points, verbose = verbose,
                                resample = 0, c_tol = c_tol, i_tol = i_tol,
                                to_file = to_file, chain.call = TRUE, ...)
  }
  else if (this_cutoff != cutoff) {
    if (verbose) #nocov start
      cat("Point implausibilities within tolerance.",
                  "Proposed points have maximum implausibility",
                  round(this_cutoff, 3), "\n") #nocov end
  }
  if (!is.null(points$cutoff)) {
    cutoff <- points$cutoff
    points <- points$points
  }
  if(!is.null(to_file)) write.csv(points, file = to_file, row.names = FALSE)
  if (("importance" %in% which_methods ||
       "line" %in% which_methods ||
       "slice" %in% which_methods) && resample > 0) {
    for (nsamp in 1:resample) {
      if (verbose) cat(paste("Resample", nsamp, "\n")) #nocov
      points <- maximin_sample(points, min(nrow(points), ceiling(n_points/2)), nms = names(ranges))
      n_current <- nrow(points)
      if ("line" %in% which_methods) {
        if (verbose) cat("Performing line sampling...\n") #nocov
        points <- line_sample(ems, ranges, z, points,
                              cutoff = cutoff, nth = nth, imp_func = imp_func, ...)
        if (verbose) cat("Line sampling generated", #nocov
                                 nrow(points)-n_current, "more points.\n") #nocov
        if (!is.null(to_file)) write.csv(points, file = to_file, #nocov
                                         row.names = FALSE) #nocov
        n_current <- nrow(points)
      }
      if ("slice" %in% which_methods) {
        if (verbose) cat("Performing slice sampling...\n") #nocov
        spoints <- slice_gen(ems, ranges, n_points, z, points,
                             cutoff, nth, imp_func = imp_func, ...)
        if (verbose) cat("Slice sampling generated", #nocov
                                 nrow(spoints)-nrow(points), "more points.\n") #nocov
        points <- spoints
        if (!is.null(to_file)) write.csv(points, file = to_file, #nocov
                                         row.names = FALSE) #nocov
        n_current <- nrow(points)
      }
      if ("importance" %in% which_methods && nrow(points) < n_points) {
        if (verbose) cat("Performing importance sampling...\n") #nocov
        points <- importance_sample(ems, n_points, z, points, cutoff,
                                    nth, to_file = to_file, imp_func = imp_func, ...)
        if (verbose) cat("Importance sampling generated", #nocov
                                 nrow(points)-n_current, "more points.\n") #nocov
        n_current <- nrow(points)
      }
    }
  }
  if (seek > 0) {
    if (verbose) cat("Searching for high-probability match points...\n") #nocov
    if (seek <= 1) seek <- floor(n_points*seek)
    extra_points <- seek_good(ems, seek, z, points, cutoff = cutoff, imp_func, ...)
  }
  else extra_points <- NULL
  if (nrow(points) > n_points - seek) {
    if (verbose) cat("Selecting final points using maximin criterion...\n") #nocov
    points <- maximin_sample(points, n_points-seek, nms = names(ranges))
  }
  chained <- list(...)[['chain.call']]
  if (!is.null(chained)) return(list(points = rbind(points, extra_points),
                                     cutoff = cutoff))
  if(!is.null(to_file)) write.csv(rbind(points, extra_points), #nocov
                                  file = to_file, row.names = FALSE) #nocov
  return(rbind(points, extra_points))
}

## LHS Generation
lhs_gen <- function(ems, ranges, n_points, z, cutoff = 3,
                    nth = 1, points.factor = 40, use_pca = FALSE, imp_func, ...) {
  if (use_pca) {
    in_range <- function(data, ranges) {
      apply(data, 1,
            function(x) all(
              purrr::map_lgl(
                seq_along(ranges),
                ~x[.] >= ranges[[.]][1] && x[.] <= ranges[[.]][2])))
    }
    if (!is.null(ems$expectation)) {
      train_pts <- ems$expectation[[1]]$in_data
      init_ranges <- ems$expectation[[1]]$ranges
    }
    else if (!is.null(ems$mode1)) {
      train_pts <- ems$mode1$expectation[[1]]$in_data
      init_ranges <- ems$mode1$expectation[[1]]$ranges
    }
    else {
      train_pts <- ems[[1]]$in_data
      init_ranges <- ems[[1]]$ranges
    }
    actual_points <- eval_funcs(scale_input, train_pts, init_ranges, FALSE)
    pca_points <- pca_transform(actual_points, actual_points)
    pca_ranges <- purrr::map(seq_len(ncol(pca_points)), ~range(pca_points[,.])) |> setNames(paste0("X", seq_len(ncol(pca_points))))
    pca_ranges <- purrr::map(pca_ranges, ~.*c(0.9, 1.1))
    temp_pts <- eval_funcs(scale_input,
                         setNames(
                           data.frame(2 * lhs::randomLHS(n_points * points.factor, length(pca_ranges)) - 0.5),
                         paste0("X", seq_along(pca_ranges))),
                         pca_ranges, FALSE)
    points <- data.frame(pca_transform(temp_pts, actual_points, FALSE)) |> setNames(names(ranges))
    points <- points[in_range(points, init_ranges),]
  }
  else {
    points <- eval_funcs(
      scale_input,
      setNames(
        data.frame(
          2 * (lhs::randomLHS(n_points * points.factor, length(ranges)) - 0.5)),
        names(ranges)), ranges, FALSE)
  }
  point_imps <- imp_func(ems, points, z, n = nth, max_imp = Inf, ordered = TRUE)
  required_points <- min(length(point_imps)-1, floor(0.8*length(point_imps)), 5*length(ranges))
  if (sum(point_imps <= cutoff) < required_points) {
    cutoff_current <- sort(point_imps)[required_points]
    if (sort(point_imps)[1] >= 20) {
      warning(paste("Parameter space has no points below implausibility 20;",
                    "terminating early. This may not indicate model inadequacy:",
                    "inspect results and re-run if applicable."))
      return(list(points = points[point_imps <= 0,], cutoff = 0))
    }
  }
  else
    cutoff_current <- cutoff
  return(list(points = points[point_imps <= cutoff_current,],
              cutoff = cutoff_current))
}

## LHS Generation with emulator clustering
lhs_gen_cluster <- function(ems, ranges, n_points, z, cutoff = 3, nth = 1,
                            verbose = FALSE,
                            points.factor = 10, c_tol = 0.1, imp_func, ...) {
  which_active <- setNames(
    data.frame(
      do.call(
        'rbind',
        purrr::map(ems, ~.$active_vars))), names(ranges))
  cluster_id <- tryCatch(
    Mclust(which_active, G = 1:2, verbose = FALSE)$classification,
    error = function(e) {
      warning("Cannot distinguish two clusters in emulator active variables.")
      return(NULL)
    }
  )
  if (is.null(cluster_id) || length(unique(cluster_id)) == 1)
    return(lhs_gen(ems, ranges, n_points, z,
                   cutoff, nth, points.factor, imp_func = imp_func, ...))
  c1 <- ems[cluster_id == 1]
  c2 <- ems[cluster_id == 2]
  p1 <- unique(do.call(c, purrr::map(c1, ~names(ranges)[.$active_vars])))
  p2 <- unique(do.call(c, purrr::map(c2, ~names(ranges)[.$active_vars])))
  pn <- intersect(p1, p2)
  if (length(union(p1, p2)) == length(pn) && all(union(p1, p2) %in% pn))
    return(lhs_gen(ems, ranges, n_points, z, cutoff, nth, imp_func = imp_func, ...))
  if (length(p1) > length(p2)) {
    c1 <- ems[cluster_id == 2]
    c2 <- ems[cluster_id == 1]
    p1 <- unique(do.call(c, purrr::map(c1, ~names(ranges)[.$active_vars])))
    p2 <- unique(do.call(c, purrr::map(c2, ~names(ranges)[.$active_vars])))
  }
  if (verbose) cat("Clusters determined. Cluster 1 has length", #nocov start
                      length(c1), "with", length(p1),
                      "active variables; cluster 2 has length",
                      length(c2), "with", length(p2),
                      "active variables -", length(pn), "shared.\n")
  if (verbose) cat("Proposing from clusters.\n") #nocov end
  lhs1 <- setNames(
    data.frame(2 * (lhs::randomLHS(n_points * 10, length(p1))-0.5)), p1)
  spare1 <- ranges[!names(ranges) %in% p1]
  for (i in names(spare1)) lhs1[[i]] <- runif(nrow(lhs1), -1, 1)
  lhs1 <- eval_funcs(scale_input, lhs1[,names(ranges)], ranges, FALSE)
  imps1 <- imp_func(c1, lhs1, z, n = nth, max_imp = Inf, ordered = TRUE)
  required_points <- 5*length(ranges)
  if (sum(imps1 <= cutoff) < required_points) {
    cutoff_current <- sort(imps1)[required_points]
    if (sort(imps1)[1] >= 20) {
      warning(paste("Parameter space has no points below implausibility 20;",
      "terminating early. This may not indicate model inadequacy:",
      "inspect results and re-run if applicable."))
      return(list(points = lhs1[imps1 <= 0,], cutoff = 0))
    }
  }
  else cutoff_current <- cutoff
  valid1 <- lhs1[imps1 <= cutoff_current,]
  LHS_augment <- function(df, p1, p2, ranges, n_points) {
    n_lhs <- ceiling(n_points * 10/nrow(df))
    new_lhs <- setNames(
      data.frame(2 * (lhs::randomLHS(n_lhs, length(p2)) - 0.5)), p2)
    new_lhs <- eval_funcs(scale_input, new_lhs, ranges[p2], FALSE)
    new_lhs <- new_lhs[rep(seq_len(nrow(new_lhs)), each = nrow(df)),]
    new_lhs[,setdiff(p1, intersect(p1, p2))] <- df[
      rep(seq_len(nrow(df)), n_lhs), setdiff(p1, intersect(p1, p2))]
    if (length(intersect(p1, p2)) > 0)
      for (i in intersect(p1, p2))
        new_lhs[,i] <- df[rep(seq_len(nrow(df)), n_lhs), i] +
      runif(nrow(new_lhs)) *
      (new_lhs[,i] - df[rep(seq_len(nrow(df)), n_lhs), i])
    if (length(setdiff(names(ranges), union(p1, p2))) > 0)
      for (i in setdiff(names(ranges), union(p1, p2)))
        new_lhs[,i] <- runif(nrow(new_lhs), ranges[[i]][1], ranges[[i]][2])
    return(new_lhs[,names(ranges)])
  }
  second_stage <- function(df, cutoff, p1, p2, ranges, n_points, ems, z, nth) {
    lhs2 <- LHS_augment(df, p1, p2, ranges, n_points)
    imps2 <- imp_func(ems, lhs2, z, n = nth, max_imp = Inf, ordered = TRUE)
    if (sum(imps2 <= cutoff) < required_points) {
      coff <- sort(imps2)[required_points]
      if (sort(imps2)[1] >= 20) {
        warning(paste("Parameter space has no points below implausibility 20;",
                      "terminating early. This may not indicate model inadequacy:",
                      "inspect results and re-run if applicable."))
        return(list(points = lhs2[imps2 <= 0, ], cutoff = 0))
      }
    }
    else coff <- cutoff
    return(list(points = lhs2[imps2 <= coff,], cutoff = coff))
  }
  valid2 <- second_stage(valid1, cutoff_current, p1, p2,
                         ranges, n_points, c2, z, nth)
  while (valid2$cutoff - cutoff_current > c_tol) {
    if (verbose) cat("Cutoff increase to", valid2$cutoff, #nocov
                          "- resampling from cluster 1.\n") #nocov
    cutoff_current <- valid2$cutoff
    valid1 <- lhs1[imps1 <= cutoff_current,]
    valid2 <- second_stage(valid1, cutoff_current, p1, p2,
                           ranges, n_points, c2, z, nth)
  }
  return(valid2)
}

## Line Sampling Function
line_sample <- function(ems, ranges, z, s_points, n_lines = 20,
                        ppl = 50, cutoff = 3, nth = 1, imp_func, ...) {
  in_range <- function(data, ranges) {
    apply(data, 1,
          function(x) all(
            purrr::map_lgl(
              seq_along(ranges),
              ~x[.] >= ranges[[.]][1] && x[.] <= ranges[[.]][2])))
  }
  if (nrow(s_points) < 2) return(s_points)
  n_lines <- min(nrow(s_points)*(nrow(s_points)-1)/2, n_lines)
  if (ppl %% 4 == 1) ppl <- ppl + 1
  s_lines <- lapply(1:(10*n_lines), function(x) {
    pts <- s_points[sample(nrow(s_points), 2),]
    pt_dist <- dist(pts)
    return(list(p1 = pts[1,], p2 = pts[2,], d = pt_dist))
  })
  s_lines <- s_lines[!duplicated(purrr::map_dbl(s_lines, ~.$d))]
  best_pts <- s_lines[order(purrr::map_dbl(s_lines, ~.$d),
                            decreasing = TRUE)][1:n_lines]
  # get_limits <- function(points) {
  #   point_dist <- sqrt(sum((points[[1]]-points[[2]])^2))
  #   range_dist <- sqrt(sum(purrr::map_dbl(ranges, ~(.[[2]]-.[[1]])^2)))
  #   dist_ratio <- range_dist/point_dist
  #   l_seq <- seq(1-dist_ratio, dist_ratio, length.out = 500)
  #   line_points <- setNames(
  #     do.call(
  #       'rbind.data.frame',
  #       purrr::map(
  #         l_seq,
  #         ~points[[1]] + .*(points[[2]]-points[[1]]))), names(ranges))
  #   is_in_range <- in_range(line_points, ranges)
  #   valid_ls <- l_seq[is_in_range]
  #   return(c(valid_ls[1], valid_ls[length(valid_ls)]))
  # }
  samp_pts <- lapply(best_pts, function(x) {
    #these_lims <- get_limits(x)
    tryCatch(
      {
        do.call(
          'rbind',
          lapply(seq(-x[[3]], x[[3]], length.out = ppl),
                 function(y) (x[[1]]+x[[2]])/2 + y*(x[[2]]-x[[1]])))
       },
      error = function(e) {
        return(NULL)
      }
    )
  })
  samp_pts <- samp_pts[!purrr::map_lgl(samp_pts, is.null)]
  samp_pts <- purrr::map(samp_pts, ~.[in_range(., ranges),])
  samp_pts <- samp_pts[!purrr::map_lgl(
    samp_pts, ~is.null(nrow(.)) || is.null(.) || nrow(.) == 0)]
  imps <- purrr::map(samp_pts, ~imp_func(ems, .,
                                                z, n = nth, cutoff = cutoff,
                                                ordered = TRUE))
  include_pts <- purrr::map(seq_along(samp_pts), function(x) {
    pts <- samp_pts[[x]]
    imp <- imps[[x]]
    included <- purrr::map_lgl(seq_along(imp), function(y) {
      if (!imp[y]) return(FALSE)
      if (y == 1 || y == length(imp)) return(TRUE)
      if (!imp[y+1] || !imp[y-1]) return(TRUE)
      return(FALSE)
    })
    return(pts[included,])
  })
  out_df <- rbind(s_points, do.call('rbind', include_pts))
  uniqueness <- row.names(unique(signif(out_df, 6)))
  return(out_df[uniqueness, ])
}

# Importance Sampling function
importance_sample <- function(ems, n_points, z, s_points, cutoff = 3,
                              nth = 1, distro = 'sphere', sd = NULL,
                              to_file = NULL, imp_func, ...) {
  if (nrow(s_points) >= n_points)
    return(s_points)
  m_points <- n_points - nrow(s_points)
  ranges <- getRanges(ems, FALSE)
  new_points <- s_points
  in_range <- function(data, ranges) {
    apply(
      data, 1,
      function(x) all(purrr::map_lgl(
        seq_along(ranges), ~x[.] >= ranges[[.]][1] && x[.] <= ranges[[.]][2])))
  }
  s_trafo <- sweep(
    sweep(
      s_points, 2,
      apply(s_points, 2, mean), "-"), 2, apply(s_points, 2, sd), "/")
  s_estruct <- eigen(cov(s_trafo))
  s_estruct$values <- purrr::map_dbl(
    s_estruct$values, function(x) {if(x < 1e-10) 1e-10 else x})
  pca_transform <- function(x, forward = TRUE) {
    if (forward) x <- sweep(
      sweep(
        x, 2,
        apply(s_points, 2, mean), "-"), 2, apply(s_points, 2, sd), "/")
    if ("data.frame" %in% class(x)) x <- data.matrix(x)
    if (forward) return(x %*% s_estruct$vectors %*%
                          diag(1/sqrt(s_estruct$values)))
    pre_traf <- x %*% diag(sqrt(s_estruct$values)) %*% t(s_estruct$vectors)
    return(sweep(
      sweep(
        pre_traf, 2,
        apply(s_points, 2, sd), "*"), 2, apply(s_points, 2, mean), "+"))
  }
  propose_points <- function(sp, sd, how_many = n_points) {
    sp_trafo <- pca_transform(sp)
    wp <- sp_trafo[sample(nrow(sp_trafo), how_many, replace = TRUE), ]
    if (distro == "normal")
      pp <- t(apply(wp, 1,
                    function(x) mvtnorm::rmvnorm(
                      1, mean = unlist(x, use.names = F),
                      sigma = sd)))
    else pp <- t(apply(
      wp, 1,
      function(x) runifs(1, length(s_points), unlist(x, use.names = F),
                         r = sd)))
    prop_points <- data.frame(pp)
    back_traf <- setNames(data.frame(pca_transform(pp, FALSE)), names(ranges))
    valid <- in_range(back_traf, ranges) &
      imp_func(ems, back_traf, z, n = nth, cutoff = cutoff, ordered = TRUE)
    if (distro == "normal") {
      tweights <- apply(
        prop_points, 1,
        function(x) 1/nrow(sp_trafo) * sum(
          apply(
            sp_trafo, 1,
            function(y) mvtnorm::dmvnorm(x, mean = y, sigma = sd))))
      min_w <- min(
        apply(
          sp_trafo, 1,
          function(x) 1/nrow(sp_trafo) * sum(
            apply(
              sp_trafo, 1,
              function(y) mvtnorm::dmvnorm(x, mean = y, sigma = sd)))))
      weights <- min_w/tweights
    }
    else
      weights <- apply(
        prop_points, 1,
        function(x) sum(apply(sp_trafo, 1, function(y) punifs(x, y, r = sd))))
    allow <- runif(length(weights)) < 1/weights
    accepted <- valid & allow
    return(back_traf[accepted,])
  }
  if (is.null(sd)) {
    if (distro == "normal")
      sd <- diag(2, length(ranges))
    else
      sd <- rep(2, length(ranges))
  }
  accept_rate <- NULL
  upper_accept <- 0.225
  lower_accept <- 0.075
  while ((is.null(accept_rate) ||
          accept_rate > upper_accept ||
          accept_rate < lower_accept) &&
         nrow(new_points) < n_points) {
    if (!is.null(accept_rate)) {
      if (accept_rate > upper_accept) sd <- sd * 1.1
      else sd <- sd * 0.9
    }
    how_many <- max(floor(n_points/4), 500)
    prop_points <- propose_points(s_points, sd, how_many)
    new_points <- rbind(new_points, prop_points)
    uniqueness <- row.names(unique(signif(new_points, 7)))
    new_points <- new_points[uniqueness,]
    if (!is.null(to_file)) #nocov
      write.csv(new_points, file = to_file, row.names = FALSE) #nocov
    accept_rate <- nrow(prop_points)/how_many
  }
  while (nrow(new_points) < n_points) {
    prop_points <- propose_points(s_points, sd,
                                  ceiling(1.5*m_points/accept_rate))
    new_points <- rbind(new_points, prop_points)
    uniqueness <- row.names(unique(signif(new_points, 7)))
    new_points <- new_points[uniqueness,]
    if (!is.null(to_file)) write.csv(new_points, #nocov
                                     file = to_file, row.names = FALSE) #nocov
  }
  return(new_points)
}

# Slice sampling point generation
slice_gen <- function(ems, ranges, n_points, z, points, cutoff = 3, nth = 1, pca = FALSE, imp_func, ...) {
  in_range <- function(data, ranges) {
    apply(
      data, 1,
      function(x) all(
        purrr::map_lgl(
          seq_along(ranges),
          ~x[.] >= ranges[[.]][1] && x[.] <= ranges[[.]][2])))
  }
  pca_transform <- function(x, points, forward = TRUE) {
    s_trafo <- sweep(
      sweep(
        points, 2,
        apply(points, 2, mean), "-"), 2, apply(points, 2, sd), "/")
    s_estruct <- eigen(cov(s_trafo))
    s_estruct$values <- purrr::map_dbl(
      s_estruct$values, function(x) {if(x < 1e-10) 1e-10 else x})
    if (forward) x <- sweep(
      sweep(
        x, 2,
        apply(points, 2, mean), "-"), 2, apply(points, 2, sd), "/")
    if ("data.frame" %in% class(x)) x <- data.matrix(x)
    if (forward)
      return(x %*% s_estruct$vectors %*% diag(1/sqrt(s_estruct$values)))
    pre_traf <- x %*% diag(sqrt(s_estruct$values)) %*% t(s_estruct$vectors)
    return(sweep(
      sweep(
        pre_traf, 2,
        apply(points, 2, sd), "*"), 2, apply(points, 2, mean), "+"))
  }
  make_slice <- function(points, ranges, indices) {
    new_coords <- vapply(ranges, function(x) runif(1, x[1], x[2]), numeric(1))
    old_values <- rep(0, nrow(points))
    for (i in seq_along(indices)) {
      old_values[i] <- points[i, indices[i]]
      points[i, indices[i]] <- new_coords[i]
    }
    return(list(p = points, o = old_values))
  }
  complete_points <- pca_base <- points
  if (pca) {
    points <- pca_transform(points, pca_base)
    pca_ranges <- purrr::map(seq_along(pca_base), ~c(-5, 5))
  }
  else
    pca_ranges <- ranges
  index_list <- rep(1, nrow(points))
  while(nrow(complete_points) < n_points) {
    range_list <- purrr::map(seq_along(index_list),
                             ~pca_ranges[[index_list[.]]])
    new_slice <- make_slice(points, range_list, index_list)
    points <- new_slice$p
    old_vals <- new_slice$o
    if (pca) {
      imps <- imp_func(ems,
                              setNames(
                                data.frame(
                                  matrix(
                                    pca_transform(points, pca_base, FALSE),
                                    nrow = nrow(points))), names(ranges)),
                              z, n = nth, cutoff = cutoff, ordered = TRUE)
      in_ranges <- in_range(pca_transform(points, pca_base, FALSE), ranges)
    }
    else {
      imps <- imp_func(ems, points, z, n = nth, cutoff = cutoff, ordered = TRUE)
      in_ranges <- in_range(points, ranges)
    }
    for (i in seq_along(imps)) {
      if (imps[i] && in_ranges[i]) {
        range_list[[index_list[i]]] <- pca_ranges[[index_list[i]]]
        if (index_list[i] == length(pca_ranges)) {
          if (pca)
            complete_points <- rbind(complete_points,
                                     setNames(
                                       data.frame(
                                         matrix(
                                           pca_transform(points[i,],
                                                         pca_base, FALSE),
                                           nrow = 1)), names(ranges)))
          else
            complete_points <- rbind(complete_points, points[i,])
          index_list[i] <- 1
        }
        else
          index_list[i] <- index_list[i]+1
      }
      else {
        if (points[i, index_list[i]] < old_vals[i])
          range_list[[i]][1] <- points[i, index_list[i]]
        else range_list[[i]][2] <- points[i, index_list[i]]
      }
    }
  }
  return(complete_points)
}

# Optical depth point generation
op_depth_gen <- function(ems, ranges, n_points, z, n.runs = 100, cutoff = 3,
                         nth = 1, plausible_set, verbose = interactive(), imp_func, ...) {
  get_depth <- function(p_set, v_name, nbins = 100) {
    output <- c()
    varseq <- seq(ranges[[v_name]][1],
                  ranges[[v_name]][2], length.out = (nbins+1))
    for (i in 1:(length(varseq)-1)) {
      odepth <- nrow(p_set[p_set[[v_name]]>=varseq[i] &
                             p_set[[v_name]]<varseq[i+1],])/nrow(p_set)
      output <- c(output, (varseq[i]+varseq[i+1])/2, odepth)
    }
    return(setNames(data.frame(t(matrix(output, nrow = 2))), c('bin', 'prob')))
  }
  out_stuff <- c()
  for (i in seq_along(ranges)) {
    probs <- get_depth(plausible_set, names(ranges)[i], ...)
    new_pts <- sample(
      probs$bin,
      n_points*10, prob = probs$prob, replace = TRUE) +
      runif(n_points*10, 0, probs$bin[2]-probs$bin[1])
    out_stuff <- c(out_stuff, new_pts)
  }
  df <- setNames(
    data.frame(
      matrix(out_stuff, nrow = n_points*10)), names(ranges))
  df <- df[imp_func(ems, df, z, n = nth, cutoff = cutoff, ordered = TRUE),]
  if (nrow(df) > n_points) {
    if(verbose) cat("Selecting final points using maximin criterion...\n") #nocov
    c_measure <- op <- NULL
    for (i in 1:100) {
      tp <- df[sample(nrow(df), n_points),]
      if (!"data.frame" %in% class(tp))
        tp <- setNames(data.frame(tp), names(ranges))
      measure <- min(dist(tp))
      if (is.null(c_measure) || measure > c_measure) {
        op <- tp
        c_measure <- measure
      }
    }
    df <- op
  }
  return(df)
}

## Good point generation
#'
#' @importFrom stats pnorm
seek_good <- function(ems, n_points, z, plausible_set, cutoff = 3,
                      distro = "norm", imp_func, ...) {
  dist_func <- get(paste0("p", distro))
  get_prob <- function(ems, points, targets) {
    for (i in seq_along(targets)) {
      if (!is.atomic(targets[[i]]))
        targets[[i]] <- c(targets[[i]]$val - 3*targets[[i]]$sigma,
                          targets[[i]]$val + 3*targets[[i]]$sigma)
    }
    em_exps <- do.call('cbind.data.frame',
                       purrr::map(ems, ~.$get_exp(points)))
    em_sds <- sqrt(do.call('cbind.data.frame',
                           purrr::map(ems, ~.$get_cov(points))))
    em_probs <- do.call('rbind.data.frame',
                        purrr::map(seq_len(nrow(em_exps)), function(x) {
      purrr::map_dbl(seq_along(em_exps[x,]), function(y) {
        pnorm(targets[[ems[[y]]$output_name]][2],
              em_exps[x,y], em_sds[x,y]) -
          pnorm(targets[[ems[[y]]$output_name]][1],
                em_exps[x,y], em_sds[x,y])
      })
    }))
    result <- apply(
      em_probs, 1,
      function(x) prod(
        purrr::map_dbl(
          unique(names(targets)),
          ~min(x[purrr::map_chr(ems, function(a) a$output_name) == .]))))
    return(result)
  }
  select_minimal <- function(data, first_index, how_many) {
    data_dist <- data.matrix(
      dist(data, upper = TRUE, diag = TRUE))[-first_index,]
    picked_rows <- c(first_index)
    while(length(picked_rows) < how_many) {
      new_point_dists <- purrr::map_dbl(row.names(data_dist), function(x) {
        temp_dd <- data_dist[!row.names(data_dist) %in% x,]
        dists <- apply(temp_dd[,c(picked_rows, as.numeric(x))], 1, min)
        return(max(dists))
      })
      next_row <- row.names(data_dist)[which.min(new_point_dists)]
      data_dist <- data_dist[!row.names(data_dist) %in% next_row,]
      picked_rows <- c(picked_rows, as.numeric(next_row))
    }
    return(data[picked_rows,])
  }
  point_set <- importance_sample(ems, max(20*n_points, 4*nrow(plausible_set)),
                                 z, plausible_set, cutoff = cutoff, imp_func = imp_func, ...)
  probs <- get_prob(ems, point_set, z)
  o_points <- point_set[order(probs, decreasing = TRUE),]
  keep_points <- o_points[1:(10*n_points),]
  row.names(keep_points) <- seq_len(nrow(keep_points))
  final_set <- select_minimal(keep_points, 1, n_points)
  return(final_set)
}
