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
# Helper function to obtain the most 'recent' emulators for each output.
obtain_recent <- function(ems, dup = FALSE, showReduce = TRUE) {
  if (!is.null(ems$expectation)) {
    em_indices <- duplicated(map_chr(ems$expectation, "output_name"))
    if (!dup) em_indices <- !em_indices
    ems <- list(
      expectation = ems$expectation[em_indices],
      variance = ems$variance[em_indices]
    )
  }
  else if (!is.null(ems$mode1)) {
    em_indices <- duplicated(map_chr(ems$mode1$expectation, "output_name"))
    if (!dup) em_indices <- !em_indices
    ems <- list(
      mode1 = list(
        expectation = ems$mode1$expectation[em_indices],
        variance = ems$mode1$variance[em_indices]
      ),
      mode2 = list(
        expectation = ems$mode2$expectation[em_indices],
        variance = ems$mode2$variance[em_indices]
      ),
      prop = ems$prop[[1]]
    )
  }
  else {
    em_indices <- duplicated(map_chr(ems, "output_name"))
    if (!dup) em_indices <- !em_indices
    ems <- ems[em_indices]
  }
  if (showReduce) return(list(ems = ems, reduced = any(em_indices == FALSE)))
  return(ems)
}
## Function for pca transformation based on a sample of points
pca_transform <- function(x, s_points, forward = TRUE) {
  if (is.data.frame(s_points)) s_points <- data.matrix(s_points)
  if (is.data.frame(x)) x <- data.matrix(x)
  s_trafo <- sweep(
    sweep(
      s_points, 2,
      apply(s_points, 2, mean), "-"
    ),
    2, apply(s_points, 2, sd), "/"
  )
  s_estruct <- eigen(cov(s_trafo))
  s_estruct$values <- map_dbl(s_estruct$values, ~ifelse(. < 1e-10, 1e-10, .))
  if (forward) {
    x <- sweep(
      sweep(
        x, 2, apply(s_points, 2, mean), "-"
      ),
      2, apply(s_points, 2, sd), "/"
    )
    return(x %*% s_estruct$vectors %*% diag(1/sqrt(s_estruct$values)))
  }
  pre_traf <- x %*% diag(sqrt(s_estruct$values)) %*% t(s_estruct$vectors)
  return(sweep(sweep(pre_traf, 2, apply(s_points, 2, sd), "*"), 2, apply(s_points, 2, mean), "+"))
}
## Function to check points are in range
in_range <- function(data, ranges) {
  apply(data, 1, function(x)
    all(map_lgl(seq_along(ranges),
                       ~x[.] >= ranges[[.]][1] && x[.] <= ranges[[.]][2])))
}
## Function to obtain a maximin sample from a set of proposed points
maximin_sample <- function(points, n, reps = 1000, nms) {
  c_measure <- op <- NULL
  opts <- map(1:reps, function(rep) {
    tp <- points[sample(nrow(points), n), , drop = FALSE]
    if (!is.data.frame(tp)) tp <- setNames(tp, nms)
    measure <- min(dist(tp))
    return(list(points = tp, value = measure))
  })
  return(opts[[which.max(map_dbl(opts, "value"))]]$points)
}

#' Generate Proposal Points
#'
#' Given a set of trained emulators, this finds the next set of points that will be
#' informative for a subsequent wave of emulation or, in the event that the
#' current wave is the last desired, a set of points that optimally span the
#' parameter region of interest. There are a number of different methods that can
#' be utilised, alone or in combination with one another, to generate the points.
#'
#' If the \code{method} argument contains \code{'lhs'}, a Latin hypercube is generated and
#' non-implausible points from this design are retained. If more points are accepted
#' than the next design requires, then points are subselected using a maximin argument.
#'
#' If \code{method} contains \code{'line'}, then line sampling is performed. Given an
#' already established collection of non-implausible points, rays are drawn between
#' pairs of points (selected so as to maximise the distance between them) and more
#' points are sampled along the rays. Points thus sampled are retained if they lie
#' near a boundary of the non-implausible space, or on the boundary of the parameter
#' region of interest.
#'
#' If \code{method} contains \code{'importance'}, importance sampling is performed.
#' Given a collection of non-implausible points, a mixture distribution of either
#' multivariate normal or uniform ellipsoid proposals around the current non-implausible
#' set are constructed. The optimal standard deviation (in the normal case) or radius
#' (in the ellipsoid case) is determined using a burn-in phase, and points are
#' proposed until the desired number of points have been found.
#'
#' If \code{method} contains \code{'slice'}, then slice sampling is performed. Given
#' a single known non-implausible point, a minimum enclosing hyperrectangle (perhaps
#' after transforming the space) is determined and points are sampled for each dimension
#' of the parameter space uniformly, shrinking the minimum enclosing hyperrectangle as
#' appropriate. This method is akin to to a Gibbs sampler.
#'
#' If \code{method} contains \code{'optical'}, then optical depth sampling is used.
#' Given a set of non-implausible points, an approximation of the one-dimensional
#' marginal distributions for each parameter can be determined. From these derived
#' marginals, points are sampled and subject to rejection as in the LHD sampling.
#'
#' For any sampling strategy, the parameters \code{ems}, \code{n_points}, and \code{z}
#' must be provided. All methods rely on a means of assessing point suitability, which
#' we refer to as an implausibility measure. By default, this uses nth-maximum implausibility
#' as provided by \code{\link{nth_implausible}}; a user-defined method can be provided
#' instead by supplying the function call to \code{opts[["accept_measure"]]}. Any
#' such function must take at least five arguments: the emulators, the points, the
#' targets, and a cutoff, as well as a \code{...} argument to ensure compatibility with
#' the default behaviour of the point proposal method. Note that, in accordance with
#' the default functionality of \code{\link{nth_implausible}}, if emulating more than
#' 10 outputs and an explicit \code{opts$nth} argument is not provided, then second-max
#' implausibility is used as the measure.
#'
#' The option \code{opts[["seek"]]} determines how many points should be chosen that
#' have a higher probability of matching targets, as opposed to not missing targets. Due
#' to the danger of such an approach if a representative space-filling design over the
#' space, this value should not be too high and should be used sparingly at early waves;
#' even at later waves, it is inadvisable to seek more than 10\% of the output points
#' using this metric. The default is \code{seek = 0}, and can be provided as either
#' a percentage of points desired (in the range [0,1]) or the fixed number of points.
#'
#' The default behaviour is as follows. A set of initial points are generated from a
#' large LHD; line sampling is performed to find the boundaries of the space; then importance
#' sampling is used to fill out the space. The proposed set of points are thinned and
#' both line and importance sampling are applied again; this resampling behaviour is
#' controlled by \code{opts[["resample"]]}, where \code{resample = n} indicates that
#' the proposal will be thinned and resampled from \code{n} times (resulting in \code{n+1}
#' proposal stages).
#'
#' In regions where the non-implausible space at a given cutoff value is very hard to find,
#' the point proposal will start at a higher cutoff where it can find a space-filling design.
#' Given such a design at a higher cutoff, it can subselect to a lower cutoff by demanding
#' some percentage of the proposed points are retained and repeat. This approach terminates
#' if the 'ladder' of cutoffs reaches the desired cutoff, or if the process asymptotes at
#' a particular higher cutoff. The opts \code{ladder_tolerance} and \code{cutoff_tolerance}
#' determine the minimum improvement required in consecutive cutoffs for the process to not
#' be considered to be asymptoting and the level of closeness to the desired cutoff at whihc
#' we are prepared to stop, respectively. For instance, setting \code{ladder_tolerance} to
#' 0.1 and \code{cutoff_tolerance} to 0.01, with a cutoff of 3, will terminate the process
#' if two consecutive cutoffs proposed are within 0.1 of each other, or when the points proposed
#' all have implausibility less than the 3.01.
#'
#' These methods may work slowly, or not at all, if the target space is extremely small in
#' comparison with the initial non-yet-ruled-out (NROY) space; it may also fail to give a
#' representative sample if the target space is formed of disconnected regions of different
#' volumes.
#'
#' @section Arguments within \code{opts}:
#'  \describe{
#'  \item{accept_measure}{A custom implausibility measure to be used.}
#'  \item{cluster}{Whether to try to apply emulator clustering.}
#'  \item{cutoff_tolerance}{Tolerance for an obtained cutoff to be similar enough to that desired.}
#'  \item{ladder_tolerance}{Tolerance with which to determine if the process is asymptoting.}
#'  \item{nth}{The level of nth implausibility to apply, if using the default implausibility.}
#'  \item{resample}{How many times to perform the resampling step once points are found.}
#'  \item{seek}{How many 'good' points should be sought: either as an integer or a ratio.}
#'  \item{to_file}{If output is to be written to file periodically, the file location.}
#'  \item{points.factor (LHS, Cluster LHS)}{How many more points than desired to sample.}
#'  \item{pca_lhs (LHS)}{Whether to apply PCA to the space before proposing.}
#'  \item{n_lines (Line)}{How many lines to draw.}
#'  \item{ppl (Line)}{The number of points to sample per line.}
#'  \item{imp_distro (Importance)}{The distribution to propose around points.}
#'  \item{imp_scale (Importance)}{The radius, or standard deviation, of proposed distributions.}
#'  \item{pca_slice (Slice)}{Whether to apply PCA to the space before slice sampling.}
#'  \item{seek_distro (Seek)}{The distribution to apply when looking for 'good' points.}
#'  }
#'
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom stats setNames runif dist cov
#' @importFrom utils write.csv stack
#' @importFrom purrr imap
#' @importFrom lhs randomLHS
#'
#' @param ems A list of \code{\link{Emulator}} objects, trained
#' on previous design points.
#' @param n_points The desired number of points to propose.
#' @param z The targets to match to.
#' @param method Which methods to use.
#' @param cutoff The value of the cutoff to use to assess suitability.
#' @param plausible_set An optional set of known non-implausible points, to avoid LHD sampling.
#' @param verbose Should progress statements be printed to the console?
#' @param opts A named list of opts as described.
#' @param ... Any parameters to pass via chaining to individual sampling functions (eg \code{distro}
#' for importance sampling or \code{ordering} for collecting emulators).
#'
#' @return A data.frame containing the set of new points upon which to run the model.
#'
#' @export
#'
#' @examples
#'  \donttest{ # Excessive runtime
#'   # A simple example that uses  number of the native and ... parameter opts.
#'   pts <- generate_new_design(SIREmulators$ems, 100, SIREmulators$targets,
#'   distro = 'sphere', opts = list(resample = 0))
#'   # Non-default methods
#'   pts_slice <- generate_new_design(SIREmulators$ems, 100, SIREmulators$targets,
#'   method = 'slice')
#'   ## Example using custom measure functionality
#'   custom_measure <- function(ems, x, z, cutoff, ...) {
#'   imps_df <- nth_implausible(ems, x, z, get_raw = TRUE)
#'   sorted_imps <- t(apply(imps_df, 1, sort, decreasing = TRUE))
#'   imps1 <- sorted_imps[,1] <= cutoff
#'   imps2 <- sorted_imps[,2] <= cutoff - 0.5
#'   constraint <- apply(x, 1, function(y) y[[1]] <= 0.4)
#'   return(imps1 & imps2 & constraint)
#'   }
#'   pts_custom <- generate_new_design(SIREmulators$ems, 100, SIREmulators$targets,
#'   opts = list(accept_measure = custom_measure))
#'  }
generate_new_design <- function(ems, n_points, z, method = "default", cutoff = 3,
                              plausible_set, verbose = interactive(),
                              opts = NULL, ...) {
  if (is.null(opts)) opts <- list(...)
  else {
    collected_opts <- unlist(c(opts, list(...)), recursive = FALSE)
    opts <- as.list(collected_opts[!duplicated(names(collected_opts))])
  }
  if (is.null(opts$accept_measure)) opts$accept_measure <- "default"
  if (!is.null(opts$cluster) && opts$cluster == 1) opts$cluster <- TRUE
  if (is.null(opts$cluster) || !is.logical(opts$cluster)) opts$cluster <- FALSE
  if (is.null(opts$use_collect)) opts$use_collect <- TRUE
  else tryCatch(is.logical(opts$use_collect), warning = function(e) {warning("Emulator collection setting not logical; setting to TRUE"); opts$use_collect <- TRUE})
  if (is.null(opts$cutoff_tolerance)) opts$cutoff_tolerance <- 0.01
  else tryCatch(opts$cutoff_tolerance <- as.numeric(opts$cutoff_tolerance), warning = function(e) {warning("Cutoff tolerance is not numeric; setting to 0.01"); opts$cutoff_tolerance <- 0.01})
  if (is.null(opts$ladder_tolerance)) opts$ladder_tolerance <- 0.1
  else tryCatch(opts$ladder_tolerance <- as.numeric(opts$ladder_tolerance), warning = function(e) {warning("Ladder tolerance is not numeric; setting to 0.1"); opts$ladder_tolerance <- 0.1})
  if (is.null(opts$nth)) opts$nth <- NA
  else tryCatch(opts$nth <- as.numeric(opts$nth), warning = function(e) {warning("Nth-implausibility nth is not numeric; setting to 1"); opts$nth <- 1})
  if (is.null(opts$resample)) opts$resample <- 1
  else tryCatch(opts$resample <- as.numeric(opts$resample), warning = function(e) {warning("Resample number is not numeric; setting to 1"); opts$resample <- 1})
  if (is.null(opts$seek)) opts$seek <- 0
  else tryCatch(opts$seek <- as.numeric(opts$seek), warning = function(e) {warning("Seek value is not numeric; setting to 0"); opts$seek <- 0})
  if (opts$seek > n_points) {
    warning("Value of seek is larger than the number of desired points; setting equal to the number of points.")
    opts$seek <- n_points
  }
  if (is.null(opts$to_file)) opts$to_file <- NA
  if (!is.na(opts$to_file)) { #nocov start
    tryCatch(
      write.csv(data.frame(), file = opts$to_file, row.names = FALSE),
      error = function(e) {warning("Cannot open directory provided in to_file; output will not be saved to an external file"); opts$to_file <- NA}
    )
  } #nocov end
  if (is.null(list(...)[["cutoff_info"]])) min_cutoff <- cutoff
  else min_cutoff <- list(...)[["cutoff_info"]][1]
  if (opts$use_collect)
    ems <- collect_emulators(ems, z, cutoff, ...)
  ranges <- getRanges(ems)
  if (is.na(opts$nth)) {
    if (!is.null(ems$expectation))
      nems <- length(unique(map_chr(ems$expectation, "output_name")))
    else if (!is.null(ems$mode1))
      nems <- length(unique(map_chr(ems$mode1$expectation, "output_name")))
    else
      nems <- length(unique(map_chr(ems, "output_name")))
    opts$nth <- ifelse(nems > 10, 2, 1)
  }
  if (is.character(opts$accept_measure) && opts$accept_measure == "default") imp_func <- function(ems, x, z, ...) nth_implausible(ems, x, z, n = opts$nth, ...)
  if (length(method) == 1 && method == "default") {
    method <- c('lhs', 'line', 'importance')
    user_select <- FALSE
  }
  else user_select <- TRUE
  possible_methods <- c('lhs', 'line', 'importance', 'slice', 'optical')
  which_methods <- possible_methods[possible_methods %in% method]
  n_current <- 0
  if (any(!method %in% possible_methods)) {
    possible_user_methods <- which(!method %in% possible_methods)
    for (i in possible_user_methods) {
      tryCatch({get(i); which_methods <- c(which_methods, i)},
               error = function(e) warning(paste("Method", i, "not found.")))
    }
  }
  if (missing(plausible_set) || 'lhs' %in% which_methods) {
    if (verbose) cat("Proposing from LHS...\n") #nocov
    if (!opts$cluster) {
      lh_gen <- lhs_gen(ems, ranges, max(n_points, 10*length(ranges)),
                        z, cutoff, verbose, opts)
      points <- lh_gen$points
      this_cutoff <- lh_gen$cutoff
    }
    else {
      cluster_ems <- obtain_recent(ems)
      cluster_gen <- lhs_gen_cluster(cluster_ems$ems, ranges,
                                     max(n_points, 10*length(ranges)),
                                     z, cutoff, verbose, opts)
      if (cluster_ems$reduced) {
        if (is.character(opts$accept_measure) && opts$accept_measure == "default") {
          leftover_imps <- imp_func(
            obtain_recent(ems, TRUE, FALSE),
            cluster_gen$points, z, ordered = TRUE
          )
          this_cutoff <- max(cluster_gen$cutoff,
                             sort(leftover_imps)[min(length(leftover_imps)-1, floor(0.8*leftover_imps), 5*length(ranges))])
          points <- cluster_gen$points[leftover_imps <= this_cutoff,]
        }
        else {
          i_bools <- rep(FALSE, nrow(cluster_gen$points))
          required_points <- min(nrow(cluster_gen$points)-1, floor(0.8*nrow(cluster_gen$points)), 5*length(ranges))
          c_current <- cluster_gen$cutoff
          while (sum(i_bools) < required_points) {
            false_indices <- which(!i_bools)
            imp_bools <- opts$accept_measure(obtain_recent(ems, TRUE, FALSE), cluster_gen$points[false_indices,], z, n = opts$nth, cutoff = c_current)
            i_bools[false_indices] <- i_bools[false_indices] | imp_bools
            c_current <- c_current + 0.5
          }
          this_cutoff <- c_current - 0.5
          points <- cluster_gen$points[i_bools,]
        }
      }
      else {
        points <- cluster_gen$points
        this_cutoff <- cluster_gen$cutoff
      }
    }
    meta_verbose <- TRUE
    if (this_cutoff == cutoff && nrow(points) >= 0.25 * n_points && !user_select) {
      if (verbose) cat("LHS has high yield; no other methods required.\n") #nocov
      meta_verbose <- FALSE
      all_points <- points
      while(nrow(all_points) < n_points) {
        if (verbose) cat("Proposing from LHS...\n") #nocov
        if (!opts$cluster) {
          lh_gen <- lhs_gen(ems, ranges, max(n_points, 10*length(ranges)),
                            z, cutoff, verbose, opts)
          points <- lh_gen$points
          this_cutoff <- lh_gen$cutoff
        }
        else {
          cluster_ems <- obtain_recent(ems)
          cluster_gen <- lhs_gen_cluster(cluster_ems$ems, ranges,
                                         max(n_points, 10*length(ranges)),
                                         z, cutoff, verbose, opts)
          if (cluster_ems$reduced) {
            if (is.character(opts$accept_measure) && opts$accept_measure == "default") {
              leftover_imps <- imp_func(
                obtain_recent(ems, TRUE, FALSE),
                cluster_gen$points, z, ordered = TRUE
              )
              this_cutoff <- max(cluster_gen$cutoff,
                                 sort(leftover_imps)[min(length(leftover_imps)-1, floor(0.8*leftover_imps), 5*length(ranges))])
              points <- cluster_gen$points[leftover_imps <= this_cutoff,]
            }
            else {
              i_bools <- rep(FALSE, nrow(cluster_gen$points))
              required_points <- min(nrow(cluster_gen$points)-1, floor(0.8*nrow(cluster_gen$points)), 5*length(ranges))
              c_current <- cluster_gen$cutoff
              while (sum(i_bools) < required_points) {
                false_indices <- which(!i_bools)
                imp_bools <- opts$accept_measure(obtain_recent(ems, TRUE, FALSE), cluster_gen$points[false_indices,], z, n = opts$nth, cutoff = c_current)
                i_bools[false_indices] <- i_bools[false_indices] | imp_bools
                c_current <- c_current + 0.5
              }
              this_cutoff <- c_current - 0.5
              points <- cluster_gen$points[i_bools,]
            }
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
    if (is.null(nrow(points))) points <- data.frame(matrix(points, ncol = 1)) |> setNames(names(ems[[1]]$ranges))
    if (nrow(points) >= n_points && this_cutoff == cutoff) {
      if (verbose && meta_verbose) cat("Enough points generated from LHD - no need to apply other methods.\n") #nocov
      if (opts$seek > 0) {
        if (verbose) cat("Searching for high probability match points...\n") #nocov
        if (opts$seek <= 1) opts$seek <- floor(n_points * opts$seek)
        extra_points <- seek_good(ems, z, points, cutoff, verbose, opts)
      }
      else extra_points <- NULL
      if (nrow(points) > n_points - opts$seek) {
        if (verbose) cat("Selecting final points using maximin criterion...\n") #nocov
        points <- maximin_sample(points, n_points - opts$seek, nms = names(ranges))
      }
      return(rbind(extra_points, points))
    }
  }
  else {
    plausible_set <- plausible_set[,names(ranges)]
    if (is.character(opts$accept_measure) && opts$accept_measure == "default") {
      point_imps <- imp_func(ems, plausible_set, z, max_imp = Inf, ordered = TRUE)
      optimal_cut <- sort(point_imps)[min(length(point_imps)-1, floor(0.8*length(point_imps)), 5*length(ranges))]
      is_asymptoting <- (optimal_cut > cutoff && (optimal_cut - sort(point_imps)[1] < opts$cutoff_tolerance))
    }
    else {
      if (!is.null(list(...)[["cutoff_info"]])) {
        c_details <- list(...)[["cutoff_info"]]
        cutoff_sequence <- c(cutoff, seq(mean(c(cutoff, c_details[1])), c_details[2], by = diff(c_details)/20+(c_details[1]-cutoff)/40))
      }
      else
        cutoff_sequence <- seq(cutoff, 20, by = 0.05)
      points_accept <- c(0)
      optimal_cut <- cutoff
      if (!verbose || !requireNamespace("progressr", quietly = TRUE)) {
        required_points <- min(max(1, nrow(plausible_set)-1, floor(0.8*nrow(plausible_set))), 5*length(ranges))
        which_match <- rep(FALSE, nrow(plausible_set))
        for (i in cutoff_sequence) {
          c_bools <- opts$accept_measure(ems, plausible_set[!which_match,], z, cutoff = i, n = opts$nth)
          which_match[!which_match] <- c_bools
          points_accept <- c(points_accept, sum(which_match))
          if (sum(which_match) >= required_points) {
            optimal_cut <- i
            break
          }
        }
      }
      else {
        progressr::with_progress({
          prog <- progressr::progressor(steps = length(cutoff_sequence))
          required_points <- min(max(1, nrow(plausible_set)-1, floor(0.8*nrow(plausible_set))), 5*length(ranges))
          which_match <- rep(FALSE, nrow(plausible_set))
          for (i in seq_along(cutoff_sequence)) {
            c_bools <- opts$accept_measure(ems, plausible_set[!which_match,], z, cutoff = cutoff_sequence[i], n = opts$nth)
            which_match[!which_match] <- c_bools
            points_accept <- c(points_accept, sum(which_match))
            if (sum(which_match) >= required_points) {
              optimal_cut <- cutoff_sequence[i]
              break
            }
            prog(message = sprintf("Implausibility ladder check step %g", i))
          }
        })
      }
      min_cutoff <- cutoff_sequence[which(points_accept != 0)[1]-1]
      is_asymptoting <- (optimal_cut > cutoff &&
                           (points_accept[length(points_accept)-1] == 0) || points_accept[length(points_accept)] >= n_points)
    }
    if (is_asymptoting && optimal_cut-cutoff > opts$cutoff_tolerance) {
      if (is.character(opts$accept_measure) && opts$accept_measure == "default") {
        if (verbose) cat("Point proposal seems to be asymptoting around cutoff", #nocov
                         round(optimal_cut, 3), "- terminating.\n") #nocov
        recent_ems <- obtain_recent(ems, FALSE, FALSE)
        if (!is.null(recent_ems$mode1))
          relev_ems <- list(mode1 = recent_ems$mode1$expectation,
                            mode2 = recent_ems$mode2$expectation)
        else if (!is.null(recent_ems$expectation))
          relev_ems <- recent_ems$expectation
        else
          relev_ems <- recent_ems
        if (is.null(relev_ems$mode1)) {
          recent_imps <- do.call("cbind.data.frame",
                                 map(
                                   relev_ems,
                                   ~.$implausibility(plausible_set, z[[.$output_name]])
                                 ))
          recent_exps <- do.call("cbind.data.frame",
                                 map(
                                   relev_ems, ~.$get_exp(plausible_set)
                                 ))
          preflight(cbind(plausible_set, recent_exps), z)
        }
        else {
          recent_imps1 <- do.call('cbind.data.frame',
                                  map(
                                    relev_ems$mode1,
                                    ~.$implausibility(plausible_set, z[[.$output_name]])
                                  ))
          recent_imps2 <- do.call('cbind.data.frame',
                                  map(
                                    relev_ems$mode2,
                                    ~.$implausibility(plausible_set, z[[.$output_name]])
                                  ))
          recent_imps <- pmin(recent_imps1, recent_imps2)
        }
        ind <- values <- NULL
        plot_imps <- stack(recent_imps)
        plot_imps$ind <- factor(plot_imps$ind, levels = names(relev_ems))
        if (verbose) { #nocov start
          print(ggplot(data = plot_imps, aes(x = ind, y = values)) +
                  geom_boxplot() +
                  labs(title = "Implausibility Boxplot"))
          cat("Inspect implausibility boxplot for problematic outputs,",
              "and consider transforming them or removing them from",
              "this wave.\n")
        } #nocov end
      }
      else {
        if (verbose) cat("Point proposal seems to be asymptoting around cutoff", #nocov
                         round(optimal_cut, 3), "- terminating.\n") #nocov
      }
      if (!is.na(opts$to_file)) { #nocov start
          if (is.character(opts$accept_measure) && opts$accept_measure == "default") {
            write.csv(plausible_set[point_imps <= optimal_cut,],
                      file = opts$to_file, row.names = FALSE)
            return(list(points = plausible_set[point_imps <= optimal_cut,],
                        cutoff = optimal_cut))
          }
          else {
            good_points <- plausible_set[opts$accept_measure(ems, plausible_set, z, n = opts$nth, cutoff = optimal_cut),]
            write.csv(good_points, file = opts$to_file, row.names = FALSE)
            return(list(points = good_points, cutoff = optimal_cut))
          }
        } #nocov end
      if (is.character(opts$accept_measure) && opts$accept_measure == "default")
        return(list(points = plausible_set[point_imps <= optimal_cut,],
                    cutoff = optimal_cut))
      good_points <- plausible_set[opts$accept_measure(ems, plausible_set, z, cutoff = optimal_cut),]
      return(list(points = good_points, cutoff = optimal_cut))
    }
    if (optimal_cut < cutoff) this_cutoff <- cutoff
    else this_cutoff <- round(optimal_cut, 3)
    if (is.character(opts$accept_measure) && opts$accept_measure == "default")
      points <- plausible_set[point_imps <= this_cutoff,]
    else
      points <- plausible_set[opts$accept_measure(ems, plausible_set, z, n = opts$nth, cutoff = this_cutoff),]
  }
  if (length(ranges) == 1)
    points <- setNames(data.frame(temp = points), names(ranges))
  n_current <- nrow(points)
  if (is.null(n_current) || n_current == 0) {
    warning("No non-implausible points found from initial step.")
    return(points)
  }
  if (verbose) cat(n_current, " initial valid points generated for I=", #nocov
                   round(this_cutoff, 3), "\n", sep = "") #nocov
  if (!is.na(opts$to_file)) write.csv(points, file = opts$to_file, row.names = FALSE)
  # if ("optical" %in% which_methods && nrow(points) < n_points) {
  #   if (verbose) cat("Performing optical depth sampling...\n") #nocov
  #   points <- op_depth_gen(ems, ranges, n_points, z, cutoff = this_cutoff,
  #                          plausible_set = points, verbose = verbose, opts = opts)
  #   if (verbose) cat("Optical depth sampling generated", #nocov
  #                    nrow(points)-n_current, "more points.\n") #nocov
  #   n_current <- nrow(points)
  # }
  # if (!is.na(opts$to_file)) write.csv(points, file = opts$to_file, row.names = FALSE)
  if ("line" %in% which_methods && nrow(points) < n_points) {
    if (verbose) cat("Performing line sampling...\n") #nocov
    points <- line_sample(ems, ranges, z, points, cutoff = this_cutoff, opts = opts)
    if (verbose) cat("Line sampling generated", #nocov
                     nrow(points) - n_current, "more points.\n") #nocov
    n_current <- nrow(points)
  }
  if ("slice" %in% which_methods) {
    if (verbose) cat("Performing slice sampling...\n") #nocov
    points <- slice_gen(ems, ranges, n_points, z, points, this_cutoff, opts)
    if (verbose) cat("Slice sampling generated", #nocov
                     nrow(points) - n_current, "more points.\n") #nocov
    n_current <- nrow(points)
  }
  if (!is.na(opts$to_file)) write.csv(points, file = opts$to_file, row.names = FALSE) #nocov
  if ("importance" %in% which_methods && nrow(points) < n_points) {
    if (verbose) cat("Performing importance sampling...\n") #nocov
    points <- importance_sample(ems, n_points, z, points, this_cutoff, opts)
    if (verbose) cat("Importance sampling generated", #nocov
                     nrow(points)-n_current, "more points.\n") #nocov
    n_current <- nrow(points)
  }
  if (this_cutoff - cutoff > opts$cutoff_tolerance) {
    if(!is.na(opts$to_file)) write.csv(points, file = opts$to_file, row.names = FALSE) #nocov
    if (length(which_methods) == 1 && which_methods == "lhs") {
      warning("Could not generate points to desired target cutoff")
      return(points)
    }
    new_opts <- opts
    new_opts$chain_call <- TRUE
    new_opts$resample <- 0
    points <- generate_new_design(ems, n_points, z, which_methods[!which_methods %in% c('lhs')],
                         cutoff = cutoff, plausible_set = points, verbose = verbose,
                         opts = new_opts, cutoff_info = c(min_cutoff, this_cutoff))
  }
  else if (this_cutoff != cutoff) {
    if (verbose) #nocov start
      cat("Point implausibilities within tolerance.",
          "Proposed points have maximum implausibility",
          round(this_cutoff, 4), "\n") #nocov end
  }
  if (!is.null(points$cutoff)) {
    cutoff <- points$cutoff
    points <- points$points
  }
  if (!is.na(opts$to_file)) write.csv(points, file = opts$to_file, row.names = FALSE)
  if (length(intersect(which_methods, c("importance", "line", "slice"))) > 0 && opts$resample > 0) {
    for (nsamp in 1:opts$resample) {
      if (verbose) cat("Resample", nsamp, "\n") #nocov
      points <- maximin_sample(points, min(nrow(points), ceiling(n_points/2)), nms = names(ranges))
      n_current <- nrow(points)
      if ("line" %in% which_methods) {
        if (verbose) cat("Performing line sampling...\n") #nocov
        points <- line_sample(ems, ranges, z, points, cutoff = cutoff, opts = opts)
        if (verbose) cat("Line sampling generated", #nocov
                         nrow(points)-n_current, "more points.\n") #nocov
        if (!is.na(opts$to_file)) write.csv(points, file = opts$to_file, row.names = FALSE) #nocov
        n_current <- nrow(points)
      }
      if ("slice" %in% which_methods) {
        if (verbose) cat("Performing slice sampling...\n") #nocov
        points <- slice_gen(ems, ranges, n_points, z, points, cutoff, opts)
        if (verbose) cat("Slice sampling generated", #nocov
                         nrow(points)-n_current, "more points.\n") #nocov
        if (!is.na(opts$to_file)) write.csv(points, file = opts$to_file, row.names = FALSE) #nocov
        n_current <- nrow(points)
      }
      if ("importance" %in% which_methods && nrow(points) < n_points) {
        if (verbose) cat("Performing importance sampling...\n") #nocov
        points <- importance_sample(ems, n_points, z, points, cutoff, opts)
        if (verbose) cat("Importance sampling generated", #nocov
                         nrow(points) - n_current, "more points.\n") #nocov
        n_current <- nrow(points)
      }
    }
  }
  if (opts$seek > 0) {
    if (verbose) cat("Searching for points with high matching probability...\n") #nocov
    if (opts$seek <= 1) opts$seek <- floor(n_points*opts$seek)
    extra_points <- seek_good(ems, z, points, cutoff, verbose, opts)
  }
  else extra_points <- NULL
  if (nrow(points) > n_points - opts$seek) {
    if (verbose) cat("Selecting final points using maximin criterion...\n") #nocov
    points <- maximin_sample(points, n_points - opts$seek, nms = names(ranges))
  }
  if (!is.null(opts$chain_call)) return(list(points = rbind(points, extra_points), cutoff = cutoff))
  if (!is.na(opts$to_file)) write.csv(rbind(points, extra_points), #nocov
                                         file = opts$to_file, row.names = FALSE)
  return(rbind(points, extra_points))
}
## LHD sampling function
lhs_gen <- function(ems, ranges, n_points, z, cutoff = 3, verbose, opts = NULL) {
  if (is.null(opts$points.factor)) opts$points.factor <- 40
  tryCatch(opts$points.factor <- as.numeric(opts$points.factor),
           warning = function(e) {warning("opts$points.factor is not numeric; setting to 40"); opts$points.factor <- 40})

  if (is.null(opts$pca_lhs) || !is.logical(opts$pca_lhs)) opts$pca_lhs <- FALSE
  if (opts$pca_lhs) {
    if (!is.null(ems$expectation)) {
      train_pts <- ems$expectation[[1]]$in_data
      init_ranges <- ems$expectation[[1]]$ranges
    }
    else if (!is.null(ems$mode1)) {
      init_ranges <- ems$mode1$expectation[[1]]$ranges
      train_pts <- unique(rbind(ems$mode1$expectation[[1]]$in_data, ems$mode2$expectation[[1]]$in_data)[,init_ranges])
    }
    else {
      train_pts <- ems[[1]]$in_data
      init_ranges <- ems[[1]]$ranges
    }
    actual_points <- eval_funcs(scale_input, train_pts, init_ranges, FALSE)
    pca_points <- pca_transform(actual_points, actual_points)
    pca_ranges <- map(seq_len(ncol(pca_points)), ~range(pca_points[,.])) |> setNames(paste0("X", seq_len(ncol(pca_points))))
    pca_ranges <- map(pca_ranges, ~.*c(0.9, 1.1))
    temp_pts <- eval_funcs(scale_input,
                           setNames(
                             data.frame(2 * randomLHS(n_points * opts$points.factor, length(pca_ranges)) - 0.5),
                             paste0("X", seq_along(pca_ranges))
                           ), pca_ranges, FALSE)
    points <- data.frame(pca_transform(temp_pts, actual_points, FALSE)) |> setNames(names(init_ranges))
    points <- points[in_range(points, init_ranges),]
  }
  else {
    points <- eval_funcs(scale_input, setNames(
      data.frame(2* (randomLHS(n_points * opts$points.factor, length(ranges)) - 0.5)),
    names(ranges)), ranges, FALSE)
  }
  if (is.character(opts$accept_measure) && opts$accept_measure == "default") {
    imp_func <- function(ems, x, z, ...) nth_implausible(ems, x, z, n = opts$nth, ...)
    point_imps <- imp_func(ems, points, z, max_imp = Inf, ordered = TRUE)
    required_points <- min(max(1, length(point_imps)-1, floor(0.8*length(point_imps))), 5*length(ranges))
    if (sum(point_imps <= cutoff) < required_points) {
      cutoff_current <- sort(point_imps)[required_points]
      if (sort(point_imps)[1] >= 20) {
        warning(paste("Parameter space has no points below implausibility 20;",
                      "terminating early. This may not indicate model inadequacy:",
                      "inspect results and re-run if applicable."))
        return(list(points = points[point_imps <= 0,], cutoff = 0))
      }
    }
    else cutoff_current <- cutoff
    return(list(points = points[point_imps <= cutoff_current,],
                cutoff = cutoff_current))
  }
  else {
    i_bools <- rep(FALSE, nrow(points))
    required_points <- min(max(1, nrow(points)-1, floor(0.8*nrow(points))), 5*length(ranges))
    c_current <- cutoff
    while (sum(i_bools) < required_points) {
      if (c_current > 20) {
        warning(paste("Parameter space has no points below implausibility 20;",
                      "terminating early. This may not indicate model inadequacy:",
                      "inspect results and re-run if applicable."))
        return(list(points = points[rep(FALSE, nrow(points)),], cutoff = 0))
      }
      false_indices <- which(!i_bools)
      imp_bools <- opts$accept_measure(ems, points[false_indices,], z, cutoff = c_current, n = opts$nth)
      i_bools[false_indices] <- i_bools[false_indices] | imp_bools
      c_current <- c_current + 0.5
    }
    return(list(points = points[i_bools,], cutoff = c_current - 0.5))
  }
}
## Cluster LHD sampling function
#'
#' @importFrom cluster daisy fanny
lhs_gen_cluster <- function(ems, ranges, n_points, z, cutoff = 3, verbose = FALSE, opts) {
  if (is.null(opts$points.factor)) opts$points.factor <- 40
  tryCatch(opts$points.factor <- as.numeric(opts$points.factor),
           warning = function(e) {warning("opts$points.factor is not numeric; setting to 40"); opts$points.factor <- 40})
  which_active <- setNames(
    data.frame(do.call('rbind',
                       map(ems, ~.$active_vars))),
    names(ranges)
  )
  cluster_id <- NULL
  tryCatch({
      dist_mat <- suppressWarnings(daisy(which_active, metric = 'gower'))
      single_clust <- fanny(dist_mat, k = 1)
      double_clust <- fanny(dist_mat, k = 2)
      if (single_clust$objective[["objective"]] < double_clust$objective[["objective"]])
        cluster_id <- c(single_clust$clustering, use.names = FALSE)
      else cluster_id <- c(double_clust$clustering, use.names = FALSE)
    },
    error = function(e) {
      warning("Cannot distinguish two clusters in emulator active variables.")
    }
  )
  # cluster_id <- tryCatch(
  #   Mclust(which_active, G = 1:2, verbose = FALSE)$classification,
  #   error = function(e) {
  #     warning("Cannot distinguish two clusters in emulator active variables.")
  #     return(NULL)
  #   }
  # )
  if (is.null(cluster_id) || length(unique(cluster_id)) == 1)
    return(lhs_gen(ems, ranges, n_points, z, cutoff, verbose, opts))
  c1 <- ems[cluster_id == 1]
  c2 <- ems[cluster_id == 2]
  p1 <- unique(do.call(c, map(c1, ~names(ranges)[.$active_vars])))
  p2 <- unique(do.call(c, map(c2, ~names(ranges)[.$active_vars])))
  pn <- intersect(p1, p2)
  if (length(union(p1, p2)) == length(pn) && all(union(p1, p2) %in% pn))
    return(lhs_gen(ems, ranges, n_points, z, cutoff, verbose, opts))
  if (length(p1) > length(p2)) {
    c1 <- ems[cluster_id == 2]
    c2 <- ems[cluster_id == 1]
    p1 <- unique(do.call(c, map(c2, ~names(ranges)[.$active_vars])))
    p2 <- unique(do.call(c, map(c1, ~names(ranges)[.$active_vars])))
  }
  if (verbose) cat("Clusters determined. Cluster 1 has length", #nocov start
                   length(c1), "with", length(p1),
                   "active variables; cluster 2 has length",
                   length(c2), "with", length(p2),
                   "active variables;", length(pn), "shared.\n")
  if (verbose) cat("Proposing from clusters...\n") #nocov end
  cluster_proposal <- function(ems, params, ranges, np, what_cutoff = cutoff) {
    req_points <- min(floor(np * opts$points.factor/2 - 1),
                      floor(0.4 * np * opts$points.factor),
                      5*length(ranges))
    prop_lhs <- setNames(
      data.frame(2 * (randomLHS(np * opts$points.factor, length(params)) - 0.5)),
      params
    )
    spare_p <- ranges[!names(ranges) %in% params]
    for (i in names(spare_p)) prop_lhs[[i]] <- runif(nrow(prop_lhs), -1, 1)
    prop_lhs <- eval_funcs(scale_input, prop_lhs[,names(ranges)], ranges, FALSE)
    if (is.character(opts$accept_measure) && opts$accept_measure == "default") {
      imps <- nth_implausible(ems, prop_lhs, z, n = max(1, opts$nth-1), max_imp = Inf, ordered = TRUE)
      if (sum(imps <= what_cutoff) < req_points) {
        cutoff_current <- sort(imps)[req_points]
        if (sort(imps)[1] >= 20) {
          warning(paste("Parameter space has no points below implausibility 20;",
                        "terminating early. This may not indicate model inadequacy:",
                        "inspect results and re-run if applicable."))
          return(list(points = prop_lhs[imps <= 0,], cutoff = 0))
        }
      }
      else cutoff_current <- what_cutoff
      valid <- prop_lhs[imps <= cutoff_current,]
    }
    else {
      i_bools <- rep(FALSE, nrow(prop_lhs))
      c_current <- what_cutoff
      while (sum(i_bools) < req_points) {
        if (c_current > 20) {
          warning(paste("Parameter space has no points below implausibility 20;",
                        "terminating early. This may not indicate model inadequacy:",
                        "inspect results and re-run if applicable."))
          return(list(points = prop_lhs[rep(FALSE, nrow(prop_lhs)),], cutoff = 0))
        }
        false_indices <- which(!i_bools)
        imp_bools <- opts$accept_measure(ems, prop_lhs[false_indices,], z, cutoff = c_current, n = opts$nth)
        i_bools[false_indices] <- i_bools[false_indices] | imp_bools
        c_current <- c_current + 0.5
      }
      cutoff_current <- c_current - 0.5
      valid <- prop_lhs[i_bools,]
    }
    return(list(points = valid, cutoff = cutoff_current))
  }
  clust_1_prop <- cluster_proposal(c1, p1, ranges, n_points)
  clust_2_prop <- cluster_proposal(c2, p2, ranges, n_points)
  if (nrow(clust_1_prop$points) == 0) return(clust_1_prop)
  if (nrow(clust_2_prop$points) == 0) return(clust_2_prop)
  if (clust_1_prop$cutoff != clust_2_prop$cutoff) {
    if (clust_1_prop$cutoff < clust_2_prop$cutoff)
      clust_1_prop <- cluster_proposal(c1, p1, ranges, n_points, clust_2_prop$cutoff)
    else
      clust_2_prop <- cluster_proposal(c2, p2, ranges, n_points, clust_1_prop$cutoff)
  }
  cutoff_val <- clust_1_prop$cutoff
  if (length(intersect(p1, p2)) == 0) {
    if (verbose) cat("No shared variables in clusters.\n")
    relev_1 <- clust_1_prop$points[,p1, drop = FALSE] |> setNames(p1)
    relev_2 <- clust_2_prop$points[,p2, drop = FALSE] |> setNames(p2)
    complete_df <- data.frame(cbind(relev_1[rep(seq_len(nrow(relev_1)),
                                     each = nrow(relev_2)),],
                         relev_2[rep(seq_len(nrow(relev_2)),
                                     nrow(relev_1)),])) |> setNames(c(p1, p2))
    if (length(setdiff(names(ranges), union(p1, p2))) != 0) {
      for (i in setdiff(names(ranges), union(p1, p2)))
        complete_df[,i] <- runif(nrow(complete_df), ranges[[i]][1], ranges[[i]][2])
    }
    complete_df <- complete_df[,names(ranges)]
    return(list(points = complete_df, cutoff = cutoff_val))
  }
  relev_1 <- data.frame(clust_1_prop$points[,setdiff(p1, pn), drop = FALSE]) |>
    setNames(setdiff(p1, pn))
  relev_2 <- data.frame(clust_2_prop$points[,setdiff(p2, pn), drop = FALSE]) |>
    setNames(setdiff(p2, pn))
  combine_relev <- rbind(clust_1_prop$points[,pn, drop = FALSE],
                         clust_2_prop$points[,pn, drop = FALSE]) |> setNames(pn)
  combine_ranges <- apply(combine_relev, 2, range)
  complete_df <- data.frame(cbind(relev_1[rep(seq_len(nrow(relev_1)),
                                   each = nrow(relev_2)),],
                       relev_2[rep(seq_len(nrow(relev_2)),
                                   nrow(relev_1)),])) |>
    setNames(c(names(relev_1), names(relev_2)))
  for (i in pn) {
    complete_df[,i] <- runif(nrow(complete_df), combine_ranges[1,i], combine_ranges[2,i])
  }
  complete_df <- complete_df[,names(ranges)]
  req_points <- min(floor(n_points * opts$points.factor/2 - 1),
                    floor(0.4 * n_points * opts$points.factor),
                    5*length(ranges))
  if (is.character(opts$accept_measure) && opts$accept_measure == "default") {
    imps <- nth_implausible(ems, complete_df, z, n = opts$nth, max_imp = Inf, ordered = TRUE)
    if (sum(imps <= cutoff_val) < req_points) {
      cutoff_current <- sort(imps)[req_points]
      if (sort(imps)[1] >= 20) {
        warning(paste("Parameter space has no points below implausibility 20;",
                      "terminating early. This may not indicate model inadequacy:",
                      "inspect results and re-run if applicable."))
        return(list(points = complete_df[imps <= 0,], cutoff = 0))
      }
    }
    else cutoff_current <- cutoff_val
    valid <- complete_df[imps <= cutoff_current,]
  }
  else {
    i_bools <- rep(FALSE, nrow(complete_df))
    c_current <- cutoff_val
    while (sum(i_bools) < req_points) {
      if (c_current > 20) {
        warning(paste("Parameter space has no points below implausibility 20;",
                      "terminating early. This may not indicate model inadequacy:",
                      "inspect results and re-run if applicable."))
        return(list(points = complete_df[rep(FALSE, nrow(complete_df)),], cutoff = 0))
      }
      false_indices <- which(!i_bools)
      imp_bools <- opts$accept_measure(ems, complete_df[false_indices,], z, cutoff = c_current, n = opts$nth)
      i_bools[false_indices] <- i_bools[false_indices] | imp_bools
      c_current <- c_current + 0.5
    }
    cutoff_current <- c_current - 0.5
    valid <- complete_df[i_bools,]
  }
  return(list(points = valid, cutoff = cutoff_current))
}

## Line sampling function
line_sample <- function(ems, ranges, z, s_points, cutoff = 3, opts) {
  if (is.null(opts$n_lines)) opts$n_lines <- 20
  tryCatch(opts$n_lines <- as.numeric(opts$n_lines), warning = function(e) {warning("opts$n_lines not numeric; setting to 20"); opts$n_lines <- 20})
  if (is.null(opts$ppl)) opts$ppl <- 50
  tryCatch(opts$ppl <- as.numeric(opts$ppl), warning = function(e) {warning("opts$ppl not numeric; setting to 50"); opts$ppl <- 50})
  if (nrow(s_points) < 2) return(s_points)
  n_lines <- min(nrow(s_points)*(nrow(s_points)-1)/2, opts$n_lines)
  if (opts$ppl %% 4 == 1) opts$ppl <- opts$ppl + 1
  s_lines <- lapply(1:(10*opts$n_lines), function(x) {
    pts <- s_points[sample(nrow(s_points), 2),]
    pt_dist <- dist(pts)
    return(list(p1 = pts[1,], p2 = pts[2,], d = pt_dist))
  })
  s_lines <- s_lines[!duplicated(map_dbl(s_lines, ~.$d))]
  best_pts <- s_lines[order(map_dbl(s_lines, ~.$d),
                            decreasing = TRUE)][1:opts$n_lines]
  samp_pts <- lapply(best_pts, function(x) {
    tryCatch(
      {
        do.call('rbind', lapply(seq(-x[[3]], x[[3]], length.out = opts$ppl),
                                function(y) (x[[1]]+x[[2]])/2 + y * (x[[2]]-x[[1]])))
      },
      error = function(e) {
        warning("Problem with line sampling; investigate.")
        return(NULL)
      }
    )
  })
  samp_pts <- samp_pts[!map_lgl(samp_pts, is.null)]
  samp_pts <- map(samp_pts, ~.[in_range(., ranges),])
  samp_pts <- samp_pts[!map_lgl(
    samp_pts, ~(is.null(.) || is.null(nrow(.)) || nrow(.) == 0)
  )]
  if (is.character(opts$accept_measure) && opts$accept_measure == "default")
    imps <- map(samp_pts,
                       ~nth_implausible(
                         ems, ., z,
                         n = opts$nth, cutoff = cutoff,
                         ordered = TRUE
                         ))
  else imps <- map(samp_pts, ~opts$accept_measure(ems, ., z, cutoff = cutoff, n = opts$nth))
  include_pts <- imap(samp_pts, function(x, i) {
    pts <- x
    imp <- imps[[i]]
    included <- map_lgl(seq_along(imp), function(y) {
      if (!imp[y]) return(FALSE)
      if (y == 1 || y == length(imp)) return(TRUE)
      if (!imp[y+1] || !imp[y-1]) return(TRUE)
      return(FALSE)
    })
    return(pts[included,])
  })
  out_df <- rbind(s_points, do.call('rbind', include_pts))
  uniqueness <- row.names(unique(signif(out_df, 7)))
  return(out_df[uniqueness,])
}

importance_sample <- function(ems, n_points, z, s_points, cutoff = 3, opts) {
  if (is.null(opts$imp_distro) || !opts$imp_distro %in% c('sphere', 'normal')) opts$imp_distro <- 'sphere'
  if (is.null(opts$imp_scale)) opts$imp_scale <- 2
  tryCatch(opts$imp_scale <- as.numeric(opts$imp_scale),
           warning = function(e) {warning("opts$imp_scale is not numeric; setting to 2"); opts$imp_scale <- 2})
  if (nrow(s_points) >= n_points)
    return(s_points)
  m_points <- n_points - nrow(s_points)
  ranges <- getRanges(ems, FALSE)
  new_points <- s_points
  propose_points <- function(sp, sd, how_many = n_points) {
    if (is.character(opts$accept_measure) && opts$accept_measure == "default")
      imp_func <- function(ems, x, z, ...) nth_implausible(ems, x, z, ...)
    else
      imp_func <- opts$accept_measure
    sp_trafo <- pca_transform(sp, s_points)
    wp <- sp_trafo[sample(nrow(sp_trafo), how_many, replace = TRUE), ]
    if (opts$imp_distro == "normal")
      pp <- t(apply(wp, 1, function(x)
        mvtnorm::rmvnorm(1, mean = unlist(x, use.names = FALSE), sigma = sd)))
    else
      pp <- t(apply(wp, 1, function(x)
        runifs(1, length(s_points), unlist(x, use.names = FALSE), r = sd)))
    prop_points <- data.frame(pp)
    back_traf <- setNames(data.frame(pca_transform(pp, s_points, FALSE)), names(ranges))
    valid <- in_range(back_traf, ranges) &
      imp_func(ems, back_traf, z, cutoff = cutoff, n = opts$nth, ordered = TRUE)
    if (opts$imp_distro == "normal") {
      tweights <- apply(prop_points, 1, function(x)
        1/nrow(sp_trafo) * sum(
          apply(sp_trafo, 1, function(y)
            mvtnorm::dmvnorm(x, mean = y, sigma = sd))
        ))
      min_w <- min(apply(sp_trafo, 1, function(x)
        1/nrow(sp_trafo) * sum(
          apply(sp_trafo, 1, function(y) mvtnorm::dmvnorm(x, mean = y, sigma = sd))
        )))
      weights <- min_w/tweights
    }
    else {
      weights <- apply(prop_points, 1, function(x)
        sum(apply(sp_trafo, 1, function(y)
          punifs(x, y, r = sd))))
    }
    allow <- runif(length(weights)) < 1/weights
    accepted <- valid & allow
    return(back_traf[accepted,])
  }
  if (opts$imp_distro == "normal" && !is.matrix(opts$imp_scale)) sd <- diag(2, length(ranges))
  else if (length(opts$imp_scale) == 1) sd <- rep(2, length(ranges))
  accept_rate <- NULL
  upper_accept <- 0.225
  lower_accept <- 0.1
  while((is.null(accept_rate) || accept_rate > upper_accept ||
         accept_rate < lower_accept) && nrow(new_points) < n_points) {
    if (!is.null(accept_rate)) {
      if (accept_rate > upper_accept) sd <- sd * 1.1
      else sd <- sd * 0.9
    }
    how_many <- max(floor(n_points/4), 1000)
    prop_points <- propose_points(s_points, sd, how_many)
    new_points <- rbind(new_points, prop_points)
    uniqueness <- row.names(unique(signif(new_points, 7)))
    new_points <- new_points[uniqueness,]
    if (!is.na(opts$to_file)) write.csv(new_points, file = opts$to_file, row.names = FALSE) #nocov
    accept_rate <- nrow(prop_points)/how_many
  }
  while(nrow(new_points) < n_points) {
    prop_points <- propose_points(s_points, sd, ceiling(1.5*m_points/accept_rate))
    new_points <- rbind(new_points, prop_points)
    uniqueness <- row.names(unique(signif(new_points, 7)))
    new_points <- new_points[uniqueness,]
    if (!is.na(opts$to_file)) write.csv(new_points, file = opts$to_file, row.names = FALSE) #nocov
  }
  return(new_points)
}

slice_gen <- function(ems, ranges, n_points, z, points, cutoff, opts) {
  if (is.null(opts$pca_slice) || !is.logical(opts$pca_slice)) opts$pca_slice <- FALSE
  if (is.character(opts$accept_measure) && opts$accept_measure == "default")
    imp_func <- function(ems, x, z, ...) nth_implausible(ems, x, z, ...)
  else
    imp_func <- opts$accept_measure
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
  if (opts$pca_slice) {
    points <- pca_transform(points, pca_base)
    pca_ranges <- map(seq_along(pca_base), ~c(-5,5))
  }
  else
    pca_ranges <- ranges
  index_list <- rep(1, nrow(points))
  while (nrow(complete_points) < n_points) {
    range_list <- map(seq_along(index_list),
                             ~pca_ranges[[index_list[.]]])
    new_slice <- make_slice(points, range_list, index_list)
    points <- new_slice$p
    old_vals <- new_slice$o
    if (opts$pca_slice) {
      imps <- imp_func(ems,
                       setNames(data.frame(matrix(
                         pca_transform(points, pca_base, FALSE),
                         nrow = nrow(points)
                         )), names(ranges)),
                       z, cutoff = cutoff, n = opts$nth, ordered = TRUE)
      in_ranges <- in_range(pca_transform(points, pca_base, FALSE), ranges)
    }
    else {
      imps <- imp_func(ems, points, z, cutoff = cutoff, n = opts$nth, ordered = TRUE)
      in_ranges <- in_range(points, ranges)
    }
    for (i in seq_along(imps)) {
      if (imps[i] && in_ranges[i]) {
        range_list[[index_list[i]]] <- pca_ranges[[index_list[i]]]
        if (index_list[i] == length(pca_ranges)) {
          if (opts$pca_slice)
            complete_points <- rbind(complete_points,
                                     setNames(data.frame(matrix(
                                       pca_transform(points[i,], pca_base, FALSE),
                                       nrow = 1
                                     )), names(ranges)))
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

## Good point generation
#'
#' @importFrom stats pnorm
seek_good <- function(ems, z, plausible_set, cutoff = 3, verbose,
                      opts) {
  if (!is.null(ems$expectation)) ems <- ems$expectation
  n_points <- opts$seek
  if (is.null(opts$seek_distro)) opts$seek_distro <- "norm"
  dist_func <- get(paste0("p", opts$seek_distro))
  get_prob <- function(ems, points, targets) {
    for (i in seq_along(targets)) {
      if (!is.atomic(targets[[i]]))
        targets[[i]] <- c(targets[[i]]$val - 3*targets[[i]]$sigma,
                          targets[[i]]$val + 3*targets[[i]]$sigma)
    }
    em_exps <- do.call('cbind.data.frame',
                       map(ems, ~.$get_exp(points)))
    em_sds <- sqrt(do.call('cbind.data.frame',
                           map(ems, ~.$get_cov(points))))
    em_probs <- do.call('rbind.data.frame',
                        map(seq_len(nrow(em_exps)), function(x) {
                          map_dbl(seq_along(em_exps[x,]), function(y) {
                            pnorm(targets[[ems[[y]]$output_name]][2],
                                  em_exps[x,y], em_sds[x,y]) -
                              pnorm(targets[[ems[[y]]$output_name]][1],
                                    em_exps[x,y], em_sds[x,y])
                          })
                        }))
    result <- apply(
      em_probs, 1,
      function(x) prod(
        map_dbl(
          unique(names(targets)),
          ~min(x[map_chr(ems, function(a) a$output_name) == .]))))
    return(result)
  }
  select_minimal <- function(data, first_index, how_many) {
    data_dist <- data.matrix(
      dist(data, upper = TRUE, diag = TRUE))[-first_index,]
    picked_rows <- c(first_index)
    while(length(picked_rows) < how_many) {
      new_point_dists <- map_dbl(row.names(data_dist), function(x) {
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
  if (nrow(plausible_set) > 5*n_points) {
    plausible_set <- maximin_sample(plausible_set, 5*n_points, nms = names(plausible_set))
  }
  point_set <- importance_sample(ems, 20*n_points,
                                 z,
                                 plausible_set,
                                 cutoff = cutoff, opts)
  if (verbose) cat("Generated candidate set...\n")
  probs <- get_prob(ems, point_set, z)
  o_points <- point_set[order(probs, decreasing = TRUE),]
  keep_points <- o_points[1:(10*n_points),]
  row.names(keep_points) <- seq_len(nrow(keep_points))
  if (verbose) cat("Subselected highest probability points. Thinning...\n")
  final_set <- select_minimal(keep_points, 1, n_points)
  return(final_set)
}

#' Generate Proposal Points
#'
#' This function is deprecated in favour of \code{\link{generate_new_design}}.
#'
#' Given a set of trained emulators, this finds the next set of points that will be
#' informative for a subsequent wave of emulation or, in the event that the
#' current wave is the last desired, a set of points that optimally span the
#' parameter region of interest. There are a number of different methods that can
#' be utilised, alone or in combination with one another, to generate the points.
#'
#' If the \code{method} argument contains \code{'lhs'}, a Latin hypercube is generated and
#' non-implausible points from this design are retained. If more points are accepted
#' than the next design requires, then points are subselected using a maximin argument.
#'
#' If \code{method} contains \code{'line'}, then line sampling is performed. Given an
#' already established collection of non-implausible points, rays are drawn between
#' pairs of points (selected so as to maximise the distance between them) and more
#' points are sampled along the rays. Points thus sampled are retained if they lie
#' near a boundary of the non-implausible space, or on the boundary of the parameter
#' region of interest.
#'
#' If \code{method} contains \code{'importance'}, importance sampling is performed.
#' Given a collection of non-implausible points, a mixture distribution of either
#' multivariate normal or uniform ellipsoid proposals around the current non-implausible
#' set are constructed. The optimal standard deviation (in the normal case) or radius
#' (in the ellipsoid case) is determined using a burn-in phase, and points are
#' proposed until the desired number of points have been found.
#'
#' If \code{method} contains \code{'slice'}, then slice sampling is performed. Given
#' a single known non-implausible point, a minimum enclosing hyperrectangle (perhaps
#' after transforming the space) is determined and points are sampled for each dimension
#' of the parameter space uniformly, shrinking the minimum enclosing hyperrectangle as
#' appropriate. This method is akin to to a Gibbs sampler.
#'
#' If \code{method} contains \code{'optical'}, then optical depth sampling is used.
#' Given a set of non-implausible points, an approximation of the one-dimensional
#' marginal distributions for each parameter can be determined. From these derived
#' marginals, points are sampled and subject to rejection as in the LHD sampling.
#'
#' For any sampling strategy, the parameters \code{ems}, \code{n_points}, and \code{z}
#' must be provided. All methods rely on a means of assessing point suitability, which
#' we refer to as an implausibility measure. By default, this uses nth-maximum implausibility
#' as provided by \code{\link{nth_implausible}}; a user-defined method can be provided
#' instead by supplying the function call to \code{opts[["accept_measure"]]}. Any
#' such function must take at least five arguments: the emulators, the points, the
#' targets, and a cutoff, as well as a \code{...} argument to ensure compatibility with
#' the default behaviour of the point proposal method.
#'
#' The option \code{opts[["seek"]]} determines how many points should be chosen that
#' have a higher probability of matching targets, as opposed to not missing targets. Due
#' to the danger of such an approach if a representative space-filling design over the
#' space, this value should not be too high and should be used sparingly at early waves;
#' even at later waves, it is inadvisable to seek more than 10\% of the output points
#' using this metric. The default is \code{seek = 0}, and can be provided as either
#' a percentage of points desired (in the range [0,1]) or the fixed number of points.
#'
#' The default behaviour is as follows. A set of initial points are generated from a
#' large LHD; line sampling is performed to find the boundaries of the space; then importance
#' sampling is used to fill out the space. The proposed set of points are thinned and
#' both line and importance sampling are applied again; this resampling behaviour is
#' controlled by \code{opts[["resample"]]}, where \code{resample = n} indicates that
#' the proposal will be thinned and resampled from \code{n} times (resulting in \code{n+1}
#' proposal stages).
#'
#' In regions where the non-implausible space at a given cutoff value is very hard to find,
#' the point proposal will start at a higher cutoff where it can find a space-filling design.
#' Given such a design at a higher cutoff, it can subselect to a lower cutoff by demanding
#' some percentage of the proposed points are retained and repeat. This approach terminates
#' if the 'ladder' of cutoffs reaches the desired cutoff, or if the process asymptotes at
#' a particular higher cutoff. The opts \code{ladder_tolerance} and \code{cutoff_tolerance}
#' determine the minimum improvement required in consecutive cutoffs for the process to not
#' be considered to be asymptoting and the level of closeness to the desired cutoff at whihc
#' we are prepared to stop, respectively. For instance, setting \code{ladder_tolerance} to
#' 0.1 and \code{cutoff_tolerance} to 0.01, with a cutoff of 3, will terminate the process
#' if two consecutive cutoffs proposed are within 0.1 of each other, or when the points proposed
#' all have implausibility less than the 3.01.
#'
#' These methods may work slowly, or not at all, if the target space is extremely small in
#' comparison with the initial non-yet-ruled-out (NROY) space; it may also fail to give a
#' representative sample if the target space is formed of disconnected regions of different
#' volumes.
#'
#' @section Arguments within \code{opts}:
#'  \describe{
#'  \item{accept_measure}{A custom implausibility measure to be used.}
#'  \item{cluster}{Whether to try to apply emulator clustering.}
#'  \item{cutoff_tolerance}{Tolerance for an obtained cutoff to be similar enough to that desired.}
#'  \item{ladder_tolerance}{Tolerance with which to determine if the process is asymptoting.}
#'  \item{nth}{The level of nth implausibility to apply, if using the default implausibility.}
#'  \item{resample}{How many times to perform the resampling step once points are found.}
#'  \item{seek}{How many 'good' points should be sought: either as an integer or a ratio.}
#'  \item{to_file}{If output is to be written to file periodically, the file location.}
#'  \item{points.factor (LHS, Cluster LHS)}{How many more points than desired to sample.}
#'  \item{pca_lhs (LHS)}{Whether to apply PCA to the space before proposing.}
#'  \item{n_lines (Line)}{How many lines to draw.}
#'  \item{ppl (Line)}{The number of points to sample per line.}
#'  \item{imp_distro (Importance)}{The distribution to propose around points.}
#'  \item{imp_scale (Importance)}{The radius, or standard deviation, of proposed distributions.}
#'  \item{pca_slice (Slice)}{Whether to apply PCA to the space before slice sampling.}
#'  \item{seek_distro (Seek)}{The distribution to apply when looking for 'good' points.}
#'  }
#'
#' @param ems A list of \code{\link{Emulator}} objects, trained
#' on previous design points.
#' @param n_points The desired number of points to propose.
#' @param z The targets to match to.
#' @param method Which methods to use.
#' @param cutoff The value of the cutoff to use to assess suitability.
#' @param plausible_set An optional set of known non-implausible points, to avoid LHD sampling.
#' @param verbose Should progress statements be printed to the console?
#' @param opts A named list of opts as described.
#' @param ... Any parameters to pass via chaining to individual sampling functions (eg \code{distro}
#' for importance sampling or \code{ordering} for collecting emulators).
#'
#' @return A data.frame containing the set of new points upon which to run the model.
#'
#' @export
#'
generate_new_runs <- function(ems, n_points, z, method = "default", cutoff = 3,
                                plausible_set, verbose = interactive(),
                                opts = NULL, ...) {
  .Deprecated("generate_new_design")
  generate_new_design(ems, n_points, z, method, cutoff, plausible_set,
                      verbose, opts, ...)
}
