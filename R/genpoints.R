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
#'
#' The option \code{seek} determines how many points should be chosen that have a higher
#' probability of matching targets, as opposed to not missing targets. Due to the danger of
#' such an approach in terms of obtaining a representative space-filling design over the
#' space, this value should not be too high: a rough guide is that it should be no larger
#' than 10% of the desired number of points. The default is \code{seek = 0}.
#'
#' The default behaviour is as follows. A set of initial points are generated from an LHD;
#' line sampling is performed to find the boundaries; and finally this collection of points
#' is augmented to the desired number of points by importance sampling using uniform
#' spherical proposals.
#'
#' These methods may not work if the target space is very small comparative to the current
#' not-yet-ruled-out space, or it may miss small disconnected regions of parameter space.
#'
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom stats setNames runif dist cov
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
#' @param seek How many 'good' points to search for
#' @param ... Any parameters to pass to individual sampling functions, eg \code{distro} for importance sampling.
#'
#' @return A data.frame containing the set of new points to run the model at.
#'
#' @export
#'
#' @examples
#' # A simple example that runs quickly with some passed parameters to subroutines
#' pts <- generate_new_runs(sample_emulators$ems, 20, sample_emulators$targets,
#'  measure.method = 'maximin', distro = 'sphere', resample = 0)
#' \donttest{
#'  pts_optical <- generate_new_runs(sample_emulators$ems, 100, sample_emulators$targets,
#'   method = c('optical'), plausible_set = pts)
#'  pts_slice <- generate_new_runs(sample_emulators$ems, 100, sample_emulators$targets,
#'   method = c('slice'))
#'  pts_no_importance <- generate_new_runs(sample_emulators$ems, 100, sample_emulators$targets,
#'   method = c('line'))
#' }
generate_new_runs <- function(ems, n_points, z, method = c('lhs', 'line', 'importance'), cutoff = 3, nth = 1, plausible_set, verbose = FALSE, cluster = FALSE, resample = 1, seek = 0, ...) {
  ranges <- if ("Emulator" %in% class(ems)) ems$ranges else ems[[1]]$ranges
  if (n_points < 10*length(ranges)) n_points <- 10*length(ranges)
  possible_methods <- c('lhs', 'line', 'importance', 'slice', 'optical')
  which_methods <- possible_methods[possible_methods %in% method]
  n_current <- 0
  if (any(!method %in% possible_methods)) warning(paste("Unrecognised method name(s)", method[!method %in% possible_methods], "ignored."))
  if (missing(plausible_set) || 'lhs' %in% which_methods) {
    if (verbose) print("Performing Latin Hypercube sampling...")
    if (cluster) {
      points <- lhs_gen_cluster(ems, n_points, z, cutoff, nth, verbose = verbose, ...)
      cutoff_current <- cutoff
      while(nrow(points) < 5*length(ranges)) {
        cutoff_current <- cutoff_current + 1
        points <- lhs_gen_cluster(ems, n_points, z, cutoff, nth, verbose = verbose, ...)
      }
      while(cutoff_current > cutoff) {
        plaus_points <- points[nth_implausible(ems, points, z, cutoff = cutoff_current),]
        if (nrow(plaus_points) == 0) break
        points <- generate_new_runs(ems, n_points, z, method = which_methods[!which_methods %in% c('lhs')], cutoff = cutoff_current, nth = nth, plausible_set = plaus_points, verbose = verbose, resample = resample - 1, ...)
        cutoff_current <- cutoff_current - 0.5
      }
      if (cutoff_current != cutoff) {
        if (verbose) print(paste("Cannot reach implausibility cutoff", cutoff, "- points returned will have implausibility no bigger than", cutoff_current))
        points <- points[sample(1:nrow(points), floor(nrow(points)/2)),]
        cutoff <- cutoff_current
      }
      else {
        points <- points[nth_implausible(ems, points, z, cutoff = cutoff),]
      }
    }
    else {
      points <- eval_funcs(scale_input, setNames(data.frame(2 * (lhs::randomLHS(n_points * 10, length(ranges))-0.5)), names(ranges)), ranges, FALSE)
      point_imps <- nth_implausible(ems, points, z)
      required_points <- 5*length(ranges)
      if (sum(point_imps <= cutoff) < required_points) {
        cutoff_current <- sort(point_imps)[required_points]
        if (cutoff_current == 20) {
          warning("Parameter space has no points below implausibility 20: terminating early.")
          return(points[point_imps < 0, ])
        }
        if (length(which_methods) == 1)
          return(points[nth_implausible(ems, points, z, cutoff = cutoff_current),])
        while(cutoff_current > cutoff) {
          plaus_points <- points[point_imps <= cutoff_current,]
          if (nrow(plaus_points) == 0) {
            if (verbose) print("No plausible points in the set; cannot proceed.")
            break
          }
          if (verbose) print(paste0("Proposing at implausibility I=", round(cutoff_current, 4)))
          points <- generate_new_runs(ems, n_points, z, method = which_methods[!which_methods %in% c('lhs')], cutoff = cutoff_current, nth = nth, plausible_set = plaus_points, verbose = verbose, resample = resample - 1, ...)
          point_imps <- nth_implausible(ems, points, z)
          cutoff_temp <- max(cutoff, sort(point_imps)[required_points])
          if (cutoff_temp == cutoff_current) break
          cutoff_current <- cutoff_temp
        }
        if (cutoff_current != cutoff) {
          if (verbose) print(paste("Cannot reach implausibility cutoff", cutoff, "- points returned will have implausibility no bigger than", cutoff_current))
          points <- points[sample(1:nrow(points), floor(nrow(points)/2)),]
          cutoff <- cutoff_current
        }
        else {
          points <- points[nth_implausible(ems, points, z, cutoff = cutoff),]
        }
      }
      else {
        points <- points[point_imps <= cutoff,]
        if (!"data.frame" %in% class(points)) points <- setNames(data.frame(points), names(ranges))
      }
    }
  }
  else {
    points <- plausible_set[nth_implausible(ems, plausible_set, z, n = nth, cutoff = cutoff),]
  }
  n_current <- nrow(points)
  if (verbose) print(paste0(n_current, " points generated from LHS at I=", cutoff))
  if (nrow(points) == 0) {
    warning("No non-implausible points found from initial step.")
    return(points)
  }
  if ("slice" %in% which_methods || nrow(points) == 1) {
    if (verbose) print("Performing slice sampling...")
    starter_point <- points[sample(nrow(points), 1),]
    if (nrow(points) == 1)
      spoints <- slice_gen(ems, floor(n_points)/10, z, cutoff, nth, x1 = starter_point, ...)
    else
      spoints <- slice_gen(ems, n_points-nrow(points), z, cutoff, nth, x1 = starter_point, ...)
    if (verbose) print(paste("Slice sampling generated", nrow(spoints)-1, "more points."))
    points <- rbind(points, spoints[-1,])
    n_current <- nrow(points)
  }
  if ("optical" %in% which_methods) {
    if (verbose) print("Performing optical depth sampling...")
    points <- op_depth_gen(ems, n_points, z, cutoff = cutoff, nth = nth, plausible_set = points, verbose = verbose, ...)
    if (verbose) print(paste("Optical depth sampling generated", nrow(points)-n_current, "more points."))
    n_current <- nrow(points)
  }
  if ("line" %in% which_methods) {
    if (verbose) print("Performing line sampling...")
    points <- line_sample(ems, z, points, cutoff = cutoff, nth = nth, ...)
    if (verbose) print(paste("Line sampling generated", nrow(points)-n_current, "more points."))
    n_current <- nrow(points)
  }
  if ("importance" %in% which_methods) {
    if (verbose) print("Performing importance sampling...")
    points <- importance_sample(ems, n_points, z, points, cutoff, nth, ...)
    if (verbose) print(paste("Importance sampling generated", nrow(points)-n_current, "more points."))
    n_current <- nrow(points)
  }
  if (nrow(points) > n_points) {
    if (verbose) print("Selecting final points using maximin criterion...")
    c_measure <- op <- NULL
    for (i in 1:1000) {
      tp <- points[sample(nrow(points), n_points),]
      if (!"data.frame" %in% class(tp)) tp <- setNames(data.frame(tp), names(ranges))
      measure <- min(dist(tp))
      if (is.null(c_measure) || measure > c_measure) {
        op <- tp
        c_measure <- measure
      }
    }
    points <- op
  }
  if (("importance" %in% which_methods || "line" %in% which_methods) && resample > 0) {
    for (nsamp in 1:resample) {
      points <- points[sample(1:nrow(points), floor(n_points)/2),]
      n_current <- nrow(points)
      if ("line" %in% which_methods) {
        if (verbose) print("Resampling: Performing line sampling...")
        points <- line_sample(ems, z, points, cutoff = cutoff, nth = nth, ...)
        if (verbose) print(paste("Line sampling generated", nrow(points)-n_current, "more points."))
        n_current <- nrow(points)
      }
      if ("importance" %in% which_methods) {
        if (verbose) print("Resampling: Performing importance sampling...")
        points <- importance_sample(ems, n_points, z, points, cutoff, nth, ...)
        if (verbose) print(paste("Importance sampling generated", nrow(points)-n_current, "more points."))
        n_current <- nrow(points)
      }
      if (seek > 0) {
        if (verbose) print("Looking for high-probability match points...")
        extra_points <- seek_good(ems, seek, z, points, cutoff = cutoff, ...)
      }
      if (nrow(points) > n_points-seek) {
        if (verbose) print("Selecting final points using maximin criterion...")
        c_measure <- op <- NULL
        for (i in 1:1000) {
          tp <- points[sample(nrow(points), n_points-seek),]
          if (!"data.frame" %in% class(tp)) tp <- setNames(data.frame(tp), names(ranges))
          measure <- min(dist(tp))
          if (is.null(c_measure) || measure > c_measure) {
            op <- tp
            c_measure <- measure
          }
        }
        points <- op
      }
    }
  }
  if (seek > 0) points <- rbind(points, extra_points)
  return(points)
}

# LHS Sampling function
lhs_gen <- function(ems, n_points, z, cutoff = 3, nth = 1, measure.method = 'V_optimal', n.runs = 100, verbose = TRUE, ...) {
  ranges <- if("Emulator" %in% class(ems)) ems$ranges else ems[[1]]$ranges
  current_trace <- out_points <- NULL
  if (!measure.method %in% c('maximin', 'V_optimal')) {
    warning("Optimality method not recognised. Changing to V_optimal.")
    measure.method <- "V_optimal"
  }
  new_points <- eval_funcs(scale_input, setNames(data.frame(2 * (lhs::randomLHS(n_points * 10, length(ranges))-0.5)), names(ranges)), ranges, FALSE)
  if (!missing(z)) new_points <- new_points[nth_implausible(ems, new_points, z, n = nth, cutoff = cutoff),]
  if (!"data.frame" %in% class(new_points)) new_points <- setNames(data.frame(new_points), names(ranges))
  if (nrow(new_points) < n_points) {
    if (nrow(new_points) == 0) return(new_points)
    return(setNames(data.frame(new_points), names(ranges)))
  }
  for (i in 1:n.runs) {
    test_points <- new_points[sample(nrow(new_points), n_points),]
    if (!"data.frame" %in% class(test_points)) test_points <- setNames(data.frame(test_points), names(ranges))
    if (measure.method == "maximin") {
      measure <- min(dist(test_points))
      if (is.null(current_trace) || measure > current_trace) {
        out_points <- test_points
        current_trace <- measure
      }
    }
    else if (measure.method == "V_optimal") {
      if ("Emulator" %in% class(ems)) measure <- sum(ems$get_cov(test_points))
      else measure <- mean(purrr::map_dbl(seq_along(ems), ~sum(ems[[.]]$get_cov(test_points))))
      if (is.null(current_trace) || measure < current_trace) {
        out_points <- test_points
        current_trace <- measure
      }
    }
  }
  return(setNames(data.frame(out_points), names(ranges)))
}

lhs_gen_cluster <- function(ems, n_points, z, cutoff = 3, nth = 1, previous_ems = NULL, verbose = FALSE, ...) {
  if ("Emulator" %in% class(ems)) return(lhs_gen(ems, n_points, z, cutoff, nth))
  ranges <- ems[[1]]$ranges
  which_active <- setNames(data.frame(do.call('rbind', purrr::map(ems, ~.$active_vars))), names(ranges))
  cluster_id <- Mclust(which_active, G = 1:2, verbose = FALSE)$classification
  c1 <- ems[cluster_id == 1]
  c2 <- ems[cluster_id == 2]
  p1 <- unique(do.call(c, purrr::map(c1, ~names(ranges)[.$active_vars])))
  p2 <- unique(do.call(c, purrr::map(c2, ~names(ranges)[.$active_vars])))
  pn <- intersect(p1, p2)
  if (verbose) print(paste("Clusters determined. Cluster 1 is length", length(c1), "with", length(p1), "active variables; cluster 2 is length", length(c2), "with", length(p2), "parameters -", length(pn), "shared."))
  if (length(union(p1, p2)) == length(intersect(p1, p2)) && all(union(p1,p2) == intersect(p1,p2))) return(lhs_gen(ems, n_points, z, cutoff, nth))
  if (length(p1) == length(ranges) && length(p2) == length(ranges)) return(lhs_gen(ems, n_points, z, cutoff, nth))
  if (length(p1) > length(p2)) {
    c1 <- ems[cluster_id == 2]
    c2 <- ems[cluster_id == 1]
    p1 <- unique(do.call(c, purrr::map(c1, ~names(ranges)[.$active_vars])))
    p2 <- unique(do.call(c, purrr::map(c2, ~names(ranges)[.$active_vars])))
  }
  ## Proposing from c1
  lhs1 <- setNames(data.frame(2 * (lhs::randomLHS(n_points * 10, length(p1))-0.5)), p1)
  spare1 <- ranges[!names(ranges) %in% p1]
  for (i in 1:length(spare1)) {
    lhs1[[names(spare1)[i]]] <- runif(nrow(lhs1), -1, 1)
  }
  lhs1 <- lhs1[,names(ranges)]
  lhs1 <- eval_funcs(scale_input, lhs1, ranges, FALSE)
  valid1 <- lhs1[nth_implausible(c1, lhs1, z, cutoff = cutoff, n = nth),]
  if (nrow(valid1) == 0) return(valid1)
  if(verbose) print(paste("Proposing from cluster 1:", nrow(valid1), "points deemed acceptable."))
  ## Proposing from c2, with valid1 as a 'prior'
  n_lhs <- ceiling(n_points * 10/nrow(valid1))
  lhs2 <- setNames(data.frame(2 * (lhs::randomLHS(n_lhs, length(p2))-0.5)), p2)
  lhs2 <- eval_funcs(scale_input, lhs2, ranges[p2], FALSE)
  lhs2 <- lhs2[rep(seq_len(nrow(lhs2)), each = nrow(valid1)),]
  lhs2[,setdiff(p1, pn)] <- valid1[rep(seq_len(nrow(valid1)), n_lhs), setdiff(p1, pn)]
  if (length(pn) > 0) {
    for (i in pn)
      lhs2[,i] <- valid1[rep(seq_len(nrow(valid1)), n_lhs), i] + runif(nrow(lhs2)) * (lhs2[,i] - valid1[rep(seq_len(nrow(valid1)), n_lhs), i])
  }
  if (length(setdiff(names(ranges), union(p1, p2))) > 0) {
    for (i in setdiff(names(ranges), union(p1, p2)))
      lhs2[,i] <- runif(nrow(lhs2), ranges[[i]][1], ranges[[i]][2])
  }
  lhs2 <- lhs2[,names(ranges)]
  valid2 <- lhs2[nth_implausible(c2, lhs2, z, cutoff = cutoff, n = nth),]
  if (nrow(valid2) == 0) return(valid2)
  ## Put these back into c1, and then into previous waves (if any)
  valid2 <- valid2[nth_implausible(c1, valid2, z, cutoff = cutoff, n = nth),]
  if(verbose) print(paste("Proposing from cluster 2:", nrow(valid2), "points deemed acceptable."))
  if (!is.null(previous_ems) && nrow(valid2) != 0)
    final_out <- valid2[nth_implausible(previous_ems, valid2, z, cutoff = cutoff, n = nth),]
  else
    final_out <- valid2
  if (nrow(final_out) > n_points) final_out <- final_out[sample(1:nrow(final_out), n_points),]
  return(final_out)
}

# Line Sampling function
line_sample <- function(ems, z, s_points, nlines = 20, ppl = 50, cutoff = 3, nth = 1, ...) {
  ranges <- if("Emulator" %in% class(ems)) ems$ranges else ems[[length(ems)]]$ranges
  in_range <- function(data, ranges) {
    apply(data, 1, function(x) all(purrr::map_lgl(seq_along(ranges), ~x[.] >= ranges[[.]][1] && x[.] <= ranges[[.]][2])))
  }
  if (nrow(s_points) < 2) return(s_points)
  nlines <- min(nrow(s_points)*(nrow(s_points)-1)/2, nlines)
  if (ppl %% 4 == 1) ppl <- ppl + 1
  s_lines <- lapply(1:(10*nlines), function(x) {
    pts <- s_points[sample(nrow(s_points), 2),]
    pt_dist <- dist(pts)
    return(list(p1 = pts[1,], p2 = pts[2,], d = pt_dist))
  })
  best_pts <- s_lines[order(purrr::map_dbl(s_lines, ~.$d), decreasing = TRUE)][1:nlines]
  get_limits <- function(points) {
    point_dist <- sqrt(sum((points[[1]]-points[[2]])^2))
    range_dist <- sqrt(sum(purrr::map_dbl(ranges, ~(.[[2]]-.[[1]])^2)))
    dist_ratio <- range_dist/point_dist
    l_seq <- seq(1-dist_ratio, dist_ratio, length.out = 1000)
    line_points <- setNames(data.frame(do.call('rbind', purrr::map(l_seq, ~points[[1]] + .*(points[[2]]-points[[1]])))), names(ranges))
    is_in_range <- in_range(line_points, ranges)
    valid_ls <- l_seq[is_in_range]
    return(c(valid_ls[1], valid_ls[length(valid_ls)]))
  }
  samp_pts <- do.call('rbind', lapply(best_pts, function(x) {
    these_lims <- get_limits(x)
    do.call('rbind', lapply(seq(these_lims[1], these_lims[2], length.out = ppl), function(y) (x[[1]]+x[[2]])/2 + y*(x[[2]]-x[[1]])))
  }))
  imps <- nth_implausible(ems, samp_pts, z, n = nth, cutoff = cutoff)
  in_rg <- in_range(samp_pts, ranges)
  include_points <- purrr::map_lgl(seq_along(imps), function(x) {
    if ((x %% ppl == 0 || x %% ppl == 1) && imps[x] && in_rg[x]) return(TRUE)
    if (imps[x] && in_rg[x] && (!imps[x-1] || !imps[x+1])) return(TRUE)
    return(FALSE)
  })
  out_df <- rbind(s_points, samp_pts[include_points,])
  uniqueness <- row.names(unique(signif(out_df, 6)))
  return(out_df[uniqueness,])
}

# Importance Sampling function
importance_sample <- function(ems, n_points, z, s_points, cutoff = 3, nth = 1, distro = 'sphere', sd = NULL, ...) {
  if (nrow(s_points) >= n_points)
    return(s_points)
  m_points <- n_points - nrow(s_points)
  ranges <- if("Emulator" %in% class(ems)) ems$ranges else ems[[length(ems)]]$ranges
  new_points <- s_points
  in_range <- function(data, ranges) {
    apply(data, 1, function(x) all(purrr::map_lgl(seq_along(ranges), ~x[.] >= ranges[[.]][1] && x[.] <= ranges[[.]][2])))
  }
  s_trafo <- sweep(sweep(s_points, 2, apply(s_points, 2, mean), "-"), 2, apply(s_points, 2, sd), "/")
  s_estruct <- eigen(cov(s_trafo))
  s_estruct$values <- purrr::map_dbl(s_estruct$values, function(x) {if(x < 1e-10) 1e-10 else x})
  pca_transform <- function(x, forward = TRUE) {
    if (forward) x <- sweep(sweep(x, 2, apply(s_points, 2, mean), "-"), 2, apply(s_points, 2, sd), "/")
    if ("data.frame" %in% class(x)) x <- data.matrix(x)
    if (forward) return(x %*% s_estruct$vectors %*% diag(1/sqrt(s_estruct$values)))
    pre_traf <- x %*% diag(sqrt(s_estruct$values)) %*% t(s_estruct$vectors)
    return(sweep(sweep(pre_traf, 2, apply(s_points, 2, sd), "*"), 2, apply(s_points, 2, mean), "+"))
  }
  propose_points <- function(sp, sd, how_many = n_points) {
    sp_trafo <- pca_transform(sp)
    wp <- sp_trafo[sample(nrow(sp_trafo), how_many, replace = TRUE), ]
    if (distro == "normal") pp <- t(apply(wp, 1, function(x) mvtnorm::rmvnorm(1, mean = unlist(x, use.names = F), sigma = sd)))
    else pp <- t(apply(wp, 1, function(x) runifs(1, length(s_points), unlist(x, use.names = F), r = sd)))
    prop_points <- setNames(data.frame(pp), names(ranges))
    back_traf <- setNames(data.frame(pca_transform(pp, FALSE)), names(ranges))
    valid <- in_range(back_traf, ranges) & nth_implausible(ems, back_traf, z, n = nth, cutoff = cutoff)
    if (distro == "normal") {
      tweights <- apply(prop_points, 1, function(x) 1/nrow(sp_trafo) * sum(apply(sp_trafo, 1, function(y) mvtnorm::dmvnorm(x, mean = y, sigma = sd))))
      min_w <- min(apply(sp_trafo, 1, function(x) 1/nrow(sp_trafo) * sum(apply(sp_trafo, 1, function(y) mvtnorm::dmvnorm(x, mean = y, sigma = sd)))))
      weights <- min_w/tweights
    }
    else
      weights <- apply(prop_points, 1, function(x) sum(apply(sp_trafo, 1, function(y) punifs(x, y, r = sd))))
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
  upper_accept <- 0.125
  lower_accept <- 0.075
  while ((is.null(accept_rate) || accept_rate > upper_accept || accept_rate < lower_accept) && nrow(new_points) < n_points) {
    if (!is.null(accept_rate)) {
      if (accept_rate > upper_accept) sd <- sd * 1.1
      else sd <- sd * 0.9
    }
    how_many <- max(n_points, 250)
    prop_points <- propose_points(s_points, sd, how_many)
    new_points <- rbind(new_points, prop_points)
    uniqueness <- row.names(unique(signif(new_points, 7)))
    new_points <- new_points[uniqueness,]
    accept_rate <- nrow(prop_points)/how_many
  }
  while (nrow(new_points) < n_points) {
    prop_points <- propose_points(s_points, sd, ceiling(1.5*m_points/accept_rate))
    new_points <- rbind(new_points, prop_points)
    uniqueness <- row.names(unique(signif(new_points, 7)))
    new_points <- new_points[uniqueness,]
  }
  return(new_points)
}

# Slice sampling point generation
slice_gen <- function(ems, n_points, z, cutoff = 3, nth = 1, x1, ...) {
  ranges <- if ("Emulator" %in% class(ems)) ems$ranges else ems[[1]]$ranges
  out_points <- x1
  x0 <- x_new <- as.numeric(c(out_points, use.names = FALSE))
  for (i in 1:n_points) {
    for (j in 1:length(ranges)) {
      xl <- ranges[[j]][1]
      xr <- ranges[[j]][2]
      repeat {
        x_new[j] <- runif(1, xl, xr)
        if (nth_implausible(ems, setNames(data.frame(matrix(x_new, nrow = 1)), names(ranges)), z, n = nth, cutoff = cutoff)) break
        if(x_new[j] < x0[j]) xl <- x_new[j] else xr <- x_new[j]
      }
    }
    out_points <- rbind(out_points, x_new)
    x0 <- x_new
  }
  return(setNames(data.frame(out_points), names(ranges)))
}

# Optical depth point generation
op_depth_gen <- function(ems, n_points, z, n.runs = 100, cutoff = 3, nth = 1, plausible_set, verbose = TRUE, ...) {
  ranges <- if ("Emulator" %in% class(ems)) ems$ranges else ems[[1]]$ranges
  get_depth <- function(p_set, v_name, nbins = 100) {
    output <- c()
    varseq <- seq(ranges[[v_name]][1], ranges[[v_name]][2], length.out = (nbins+1))
    for (i in 1:(length(varseq)-1)) {
      odepth <- nrow(p_set[p_set[[v_name]]>=varseq[i] & p_set[[v_name]]<varseq[i+1],])/nrow(p_set)
      output <- c(output, (varseq[i]+varseq[i+1])/2, odepth)
    }
    return(setNames(data.frame(t(matrix(output, nrow = 2))), c('bin', 'prob')))
  }
  out_stuff <- c()
  for (i in 1:length(ranges)) {
    probs <- get_depth(plausible_set, names(ranges)[i], ...)
    new_pts <- sample(probs$bin, n_points*10, prob = probs$prob, replace = TRUE) + runif(n_points*10, 0, probs$bin[2]-probs$bin[1])
    out_stuff <- c(out_stuff, new_pts)
  }
  df <- setNames(data.frame(matrix(out_stuff, nrow = n_points*10)), names(ranges))
  df <- df[nth_implausible(ems, df, z, n = nth, cutoff = cutoff),]
  if (nrow(df) > n_points) {
    if(verbose) print("Selecting final points using maximin criterion...")
    c_measure <- op <- NULL
    for (i in 1:100) {
      tp <- df[sample(nrow(df), n_points),]
      if (!"data.frame" %in% class(tp)) tp <- setNames(data.frame(tp), names(ranges))
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
seek_good <- function(ems, n_points, z, plausible_set, cutoff = 3, distro = "norm", ...) {
  dist_func <- get(paste0("p", distro))
  get_prob <- function(ems, points, targets) {
    for (i in 1:length(targets)) {
      if (!is.atomic(targets[[i]])) targets[[i]] <- c(targets[[i]]$val - 3*targets[[i]]$sigma, targets[[i]]$val + 3*targets[[i]]$sigma)
    }
    em_exps <- data.frame(do.call('cbind', purrr::map(ems, ~.$get_exp(points))))
    em_sds <- sqrt(data.frame(do.call('cbind', purrr::map(ems, ~.$get_cov(points)))))
    em_probs <- data.frame(do.call('rbind', purrr::map(1:nrow(em_exps), function(x) {
      purrr::map_dbl(1:length(em_exps[x,]), function(y) {
        pnorm(targets[[ems[[y]]$output_name]][2], em_exps[x,y], em_sds[x,y]) - pnorm(targets[[ems[[y]]$output_name]][1], em_exps[x,y], em_sds[x,y])
      })
    })))
    result <- apply(em_probs, 1, function(x) prod(purrr::map_dbl(unique(names(targets)), ~min(x[purrr::map_chr(ems, function(a) a$output_name) == .]))))
    return(result)
  }
  select_minimal <- function(data, first_index, how_many) {
    data_dist <- data.matrix(dist(data, upper = TRUE, diag = TRUE))[-first_index,]
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
  point_set <- importance_sample(ems, max(20*n_points, 4*nrow(plausible_set)), z, plausible_set, cutoff = cutoff, ...)
  probs <- get_prob(ems, point_set, z)
  o_points <- point_set[order(probs, decreasing = TRUE),]
  keep_points <- o_points[1:(10*n_points),]
  row.names(keep_points) <- 1:nrow(keep_points)
  final_set <- select_minimal(keep_points, 1, n_points)
  return(final_set)
}
