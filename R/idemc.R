obtain_clusters <- function(x) { #nocov start
  clusters <- suppressWarnings(purrr::map(2:6, ~fanny(x, .)))
  which_clust <- which.max(purrr::map_dbl(clusters, ~.$silinfo$avg.width))
  points_clustered <- purrr::map(1:(which_clust+1), function(i) {
    x[clusters[[which_clust]]$clustering == i,]
  })
  cluster_data <- purrr::map(points_clustered, function(a) {
    return(list(mu = apply(a, 2, mean),
                sigma = var(a),
                siginv = tryCatch(
                  chol2inv(chol(var(a))),
                  error = function(e) {
                    return(MASS::ginv(var(a)))
                  }
                  )
                ))
  })
  return(list(
    clusters = cluster_data,
    total = list(
      mu = apply(x, 2, mean), sigma = var(x),
      siginv = tryCatch(
        chol2inv(chol(var(x))),
        error = function(e) {
          return(MASS::ginv(var(x)))
        })
    ))
  )
}
predict_cluster <- function(x, clusters) {
  vals <- do.call('rbind', purrr::map(clusters, function(c) {
    point_diffs <- data.matrix(sweep(x, 2, c$mu))
    return(diag(point_diffs %*% c$siginv %*% t(point_diffs)))
  }))
  if (is.null(dim(vals))) return(which.min(vals))
  return(apply(vals, 2, which.min))
}

idemc_step <- function(ems, targets, points, point_imps, ladder, clusters,
                       order_active, ranges, M = 10, pm = 0.9, w = 0.8) {
  input_names <- names(ems[[1]]$ranges)
  in_range <- function(x, ranges) {
    is_in <- purrr::map_lgl(seq_along(ranges), ~x[[.]] >= ranges[[.]][1] &&
                              x[[.]] <= ranges[[.]][2])
    return(all(is_in))
  }
  proposal <- points
  prop_imps <- point_imps
  n <- length(points)
  nmutate <- M
  ncross <- ceiling((n+1)/2)
  nex <- n+1
  do_mutate <- runif(1)
  if (do_mutate <= pm) {
    ## Do the mutation
    for (j in seq_along(points)) {
      if (j == 1) {
        proposal[[j]] <- setNames(
          data.frame(
            matrix(
              purrr::map_dbl(
                clusters[[j]]$ranges,
                ~runif(1, .[[1]], .[[2]])
              ), nrow = 1
            )
          ), input_names
        )
        prop_imps[[j]] <- nth_implausible(ems, proposal[[j]], targets, max_imp = Inf)
      }
      else {
        for (i in 1:M) {
          which_clust <- predict_cluster(points[[j]], clusters[[j]]$clusters)
          this_clust <- clusters[[j]]$clusters[[which_clust]]
          clust_mean <- this_clust$mu
          clust_var <- this_clust$sigma
          all_mean <- clusters[[j]]$total$mu
          all_var <- clusters[[j]]$total$sigma
          if (runif(1) < w)
            new_point <- mvtnorm::rmvnorm(1, mean = clust_mean, sigma = clust_var)
          else
            new_point <- mvtnorm::rmvnorm(1, mean = all_mean, sigma = all_var)
          new_clust <- predict_cluster(new_point, clusters[[j]]$clusters)
          new_var <- clusters[[j]]$clusters[[new_clust]]$sigma
          new_imp <- nth_implausible(ems, setNames(data.frame(matrix(new_point, nrow = 1)),
                                                   input_names),
                                     targets, max_imp = Inf)
          if (new_imp <= ladder[[j]] && in_range(new_point, ranges)) {
            dp <- unlist(proposal[[j]])
            np <- new_point[1,]
            accept_num <- (w * mvtnorm::dmvnorm(dp, np, clust_var) +
                              (1-w) * mvtnorm::dmvnorm(dp, np, all_var))
            accept_denom <- (w * mvtnorm::dmvnorm(np, dp, new_var) +
                                (1-w) * mvtnorm::dmvnorm(np, dp, all_var))
            prob_accept <- min(accept_num/accept_denom, 1)
            if (!is.numeric(prob_accept)) prob_accept <- 0
            accepted <- if (runif(1) <= prob_accept) new_imp else NULL
          }
          else accepted <- NULL
          if (!is.null(accepted)) {
            proposal[[j]] <- setNames(data.frame(matrix(new_point, nrow = 1)),
                                      input_names)
            prop_imps[[j]] <- new_imp
          }
        }
      }
    }
  }
  ## Do the crossover
  else {
    for (i in seq_len(ncross)) {
      index1 <- sample(1:n, 1, prob = purrr::map_dbl(1:n, ~2*./(n*(n+1))))
      if (n == 2)
        index2 <- ifelse(index1 == 1, 2, 1)
      else
        index2 <- sample((1:n)[-index1], 1,
                         prob = purrr::map_dbl((1:n)[-index1],
                                               ~2*(n+1-.)/((n-2)*(n+1)+2*index1)))
      x1 <- proposal[[index1]][,order_active]
      x2 <- proposal[[index2]][,order_active]
      c_point <- sample(1:(length(proposal[[1]])-1), 1)
      y1 <- unlist(c(x1[1:c_point], x2[(c_point+1):length(x2)]),
                   use.names = FALSE)
      y2 <- unlist(c(x2[1:c_point], x1[(c_point+1):length(x1)]),
                   use.names = FALSE)
      point_df <- rbind.data.frame(x1, x2) |> setNames(input_names[order_active])
      point_df <- point_df[,input_names]
      imps <- nth_implausible(ems, point_df, targets, max_imp = Inf)
      if (imps[[1]] <= ladder[[index1]] &&
          imps[[2]] <= ladder[[index2]] &&
          in_range(point_df[1,], ranges) &&
          in_range(point_df[2,], ranges)) {
        proposal[[index1]] <- point_df[1,]
        proposal[[index2]] <- point_df[2,]
        prop_imps[[index1]] <- imps[[1]]
        prop_imps[[index2]] <- imps[[2]]
      }
    }
  }
  ## Do the exchange
  for (i in seq_len(nex)) {
    index1 <- sample(n, 1)
    if (index1 == 1) index2 <- 2
    else if (index1 == n) index2 <- n-1
    else index2 <- index1 + sample(c(-1, 1), 1)
    if (index1 > index2) {
      index1 <- index1-1
      index2 <- index2+1
    }
    if (prop_imps[[index1]] <= ladder[[index2]]) {
      temp_point <- proposal[[index1]]
      proposal[[index1]] <- proposal[[index2]]
      proposal[[index2]] <- temp_point
      temp_imp <- prop_imps[[index1]]
      prop_imps[[index1]] <- prop_imps[[index2]]
      prop_imps[[index2]] <- temp_imp
    }
  }
  return(list(points = proposal, imps = prop_imps))
}

#' IDEMC Point Generation
#'
#' Performs Implausibility-driven Evolutionary Monte Carlo
#'
#' This method for generating points is focused on finding non-implausible
#' regions that are either extremely small relative to the initial parameter
#' domain, or have interesting structure (particularly disconnected structure)
#' that would potentially be overlooked by more standard point generation methods.
#' The method is robust but computationally intensive, compared to normal methods,
#' and should not be used as a default - see \code{\link{generate_new_runs}} for
#' less computationally expensive methods.
#'
#' The IDEMC method operates on an 'implausibility ladder', in the vein of common
#' annealing methods. Each rung of the ladder is characterised by the implausibility
#' threshold, and determinations are made about the structure of the points in each
#' rung using clustering. One step of the evolutionary algorithm can consist of the
#' following steps:
#'
#' Mutation. A point is modified using a random-walk proposal, which can be a global
#' move or a within-cluster move. Within-cluster moves are chosen with probability
#' \code{w}. The move is retained if the new point satisfies the implausibility
#' constraints of the rung.
#'
#' Crossover. Points are re-organised in descending order of how active each variable
#' is for the emulated outputs, and two different rungs are selected randomly. The
#' points are 'mixed' using one-point real crossover at a random crossover point,
#' producing two new points. The move is retained if both new points satisfy the
#' relevant implausibility constraints of their rung.
#'
#' Exchange. Two adjacent rungs are chosen and their points are swapped. The move
#' is retained if the higher-implausibility rung is appropriate for being in the
#' lower implausibility rung.
#'
#' At a given step, one of mutation or crossover is performed, with probability
#' of mutation being chosen determined by \code{pm}. If mutation is chosen, then
#' \code{M} mutation moves are performed; else \code{(n+1)/2} crossover moves are
#' performed on the \code{n} rungs. Exchange is always perfomed and \code{n+1} such
#' moves are performed.
#'
#' The choice of 'implausibility ladder' and clusters can be determined a priori,
#' or else this function performs a burn-in phase to determine them. Points are
#' generated using the idemc steps at the current rungs, and the next ladder rung
#' implausibility is chosen by requiring that a proportion \code{p} of points from
#' the previous rung are accepted in the new one. At each stage, \code{s} idemc
#' steps are performed. Once the final rung has implausibility no larger than the
#' desired \code{cutoff}, a final set of \code{sn} idemc steps are performed across
#' all rungs to determine final clusters.
#'
#' @importFrom cluster fanny
#'
#' @param ems The emulators to evaluate implausibility on
#' @param N The desired number of final points to generate
#' @param targets The target values for the emulated outputs
#' @param cutoff The desired implausibility cutoff of the final proposal
#' @param s The number of points to generate at intermediate burn-in steps
#' @param sn The number of points to generate at the final burn-in stage
#' @param p The proportion of space that should remain between ladder rungs
#' @param thin The thinning factor: a factor T means N*T points are generated to obtain N
#' @param pm The probability that an idemc step will use mutation moves
#' @param w The probability of local random walks in the mutation step
#' @param M The number of mutations to perform in an IDEMC step
#' @param detailed If TRUE, points proposed at every rung will be returned
#' @param verbose Should information about burn-in be displayed during the process?
#' @param get_burnt If TRUE, the procedure stops after burn-in, returning seeding for a full IDEMC proposal
#' @param burnt If provided, this is assumed to be the result of a burn-in (or a priori analysis)
#'
#' @return Either a list of data.frames, one per rung, or a single data.frame of points.
#'
#' @seealso \code{\link{generate_new_runs}} for more standard point generation methods
#'
#' @references
#'     Vernon & Williamson (2013) <arXiv:1309.3520>
#'
#' @export
#'
#' @examples
#'  \donttest{
#'   idemc_res <- idemc(SIREmulators$ems, 200, SIREmulators$targets, s = 100, p = 0.3)
#'  }
#'
idemc <- function(ems, N, targets, cutoff = 3, s = max(500, ceiling(N/5)),
                  sn = s, p = 0.4, thin = 1, pm = 0.9, w = 0.8, M = 10,
                  detailed = FALSE, verbose = interactive(),
                  get_burnt = FALSE, burnt = NULL) {
  ems <- collect_emulators(ems)
  ranges <- getRanges(ems, FALSE)
  order_active <- order(
    apply(do.call('rbind', purrr::map(ems, "active_vars")), 2, sum),
          decreasing = TRUE)
  if (is.null(burnt) ||
      is.null(burnt$points) ||
      is.null(burnt$clusters) ||
      is.null(burnt$imps) ||
      is.null(burnt$ladder)) {
    # Burn-in stage one: identify all rungs
    points <- setNames(do.call('cbind.data.frame',
                               purrr::map(ranges, ~runif(s, .[[1]], .[[2]]))),
                       names(ranges))
    clusters_list <- list(list(ranges = ranges))
    ladder = c(Inf)
    imps <- nth_implausible(ems, points, targets)
    imp_cutoff <- sort(imps)[floor(p*s)+1]
    orig_imp <- max(imps)
    point_imps <- c(Inf)
    points_list <- list(points[nrow(points),])
    while (imp_cutoff > cutoff) {
      if (verbose) cat(imp_cutoff, "\n")
      next_points <- points[imps <= imp_cutoff,]
      next_cluster <- obtain_clusters(next_points)
      points_list[[length(points_list)+1]] <- next_points[nrow(next_points),]
      all_points <- list()
      for (i in seq_along(points_list)) {
        df <- data.frame(matrix(0, nrow = s+1, ncol = length(ranges))) |>
          setNames(names(ranges))
        df[1,] <- points_list[[i]]
        all_points[[length(all_points)+1]] <- df
      }
      clusters_list[[length(clusters_list)+1]] <- next_cluster
      ladder <- c(ladder, imp_cutoff)
      point_imps <- c(point_imps, imps[imps <= imp_cutoff][sum(imps <= imp_cutoff)])
      all_imps <- matrix(0, nrow = s+1, ncol = length(points_list))
      all_imps[1,] <- point_imps
      if (!verbose || !requireNamespace("progressr", quietly = TRUE)) {
        for (i in seq_len(s)) {
          these_points <- purrr::map(all_points, ~.[i,])
          idemc_res <- idemc_step(ems, targets, these_points, all_imps[i,],
                                  ladder, clusters_list, order_active, ranges,
                                  M = M, pm = pm, w = w)
          all_imps[i+1,] <- idemc_res$imps
          for (j in seq_along(idemc_res$points)) {
            all_points[[j]][i+1,] <- idemc_res$points[[j]]
          }
        }
      }
      else {
        progressr::with_progress({
          prog <- progressr::progressor(steps = s)
          for (i in seq_len(s)) {
            these_points <- purrr::map(all_points, ~.[i,])
            idemc_res <- idemc_step(ems, targets, these_points, all_imps[i,],
                                    ladder, clusters_list, order_active, ranges,
                                    M = M, pm = pm, w = w)
            all_imps[i+1,] <- idemc_res$imps
            for (j in seq_along(idemc_res$points)) {
              all_points[[j]][i+1,] <- idemc_res$points[[j]]
            }
            prog(message = sprintf("IDEMC generation (%g rungs) step %g", length(all_points), i))
          }
        })
      }
      points <- all_points[[length(all_points)]]
      clusters_list <- purrr::map(seq_along(all_points), function(i) {
        if (i == 1) return(list(ranges = ranges))
        return(obtain_clusters(all_points[[i]]))
      })
      imps <- all_imps[,ncol(all_imps)]
      imp_cutoff <- sort(imps)[floor(p*s)+1]
      points_list <- purrr::map(all_points, ~.[nrow(.),])
    }
    if (imp_cutoff < cutoff) imp_cutoff <- cutoff
    if (verbose) cat(imp_cutoff, "\n")
    next_points <- points[imps <= imp_cutoff,]
    next_cluster <- obtain_clusters(next_points)
    points_list[[length(points_list)+1]] <- next_points[nrow(next_points),]
    all_points <- list()
    for (i in seq_along(points_list)) {
      df <- data.frame(matrix(0, nrow = sn+1, ncol = length(ranges))) |>
        setNames(names(ranges))
      df[1,] <- points_list[[i]]
      all_points[[length(all_points)+1]] <- df
    }
    clusters_list[[length(clusters_list)+1]] <- next_cluster
    ladder <- c(ladder, imp_cutoff)
    point_imps <- c(point_imps, imps[imps <= imp_cutoff][sum(imps <= imp_cutoff)])
    all_imps <- matrix(0, nrow = sn+1, ncol = length(points_list))
    all_imps[1,] <- point_imps
    # Burn-in stage two: generate points to determine clusters
    if (sn != s) {
      if (!verbose || !requireNamespace("progressr", quietly = TRUE)) {
        for (i in 1:sn) {
          these_points <- purrr::map(all_points, ~.[i,])
          idemc_res <- idemc_step(ems, targets, these_points, all_imps[i,],
                                  ladder, clusters_list, order_active, ranges,
                                  M = M, pm = pm, w = w)
          all_imps[i+1,] <- idemc_res$imps
          for (j in seq_along(idemc_res$points))
            all_points[[j]][i+1,] <- idemc_res$points[[j]]
        }
      }
      else {
        progressr::with_progress({
          prog <- progressr::progressor(steps = sn)
          for (i in 1:sn) {
            these_points <- purrr::map(all_points, ~.[i,])
            idemc_res <- idemc_step(ems, targets, these_points, all_imps[i,],
                                    ladder, clusters_list, order_active, ranges,
                                    M = M, pm = pm, w = w)
            all_imps[i+1,] <- idemc_res$imps
            for (j in seq_along(idemc_res$points))
              all_points[[j]][i+1,] <- idemc_res$points[[j]]
            prog(message = sprintf("IDEMC burn-in Final %g", i))
          }
        })
      }
      points <- all_points[[length(all_points)]]
      points_list <- purrr::map(all_points, ~.[nrow(.),])
      clusters_list <- purrr::map(seq_along(all_points), function(i) {
        if (i == 1) return(list(ranges = ranges))
        return(obtain_clusters(all_points[[i]]))
      })
      point_imps <- all_imps[nrow(all_imps),]
    }
    # Burn-in complete. Generate actual points
    if (verbose) cat("Burn-in completed: implausibility ladder is (",
                     round(orig_imp, 3), "; ", paste0(round(ladder[-1], 3),
                                                      collapse = "; "), ").\n",
                     sep = "")
    if (get_burnt)
      return(
        list(points = points_list,
             clusters = clusters_list,
             imps = point_imps,
             ladder = ladder)
      )
  }
  else {
    points_list <- burnt$points
    clusters_list <- burnt$clusters
    point_imps <- burnt$imps
    ladder <- burnt$ladder
  }
  if (verbose) cat("Performing final point generation...\n")
  all_points <- list()
  for (i in seq_along(points_list)) {
    df <- data.frame(matrix(0, nrow = N*thin+1, ncol = length(ranges))) |>
      setNames(names(ranges))
    df[1,] <- points_list[[i]]
    all_points[[length(all_points)+1]] <- df
  }
  all_imps <- matrix(0, nrow = N*thin+1, ncol = length(points_list))
  all_imps[1,] <- point_imps
  if (!verbose || !requireNamespace("progressr", quietly = TRUE)) {
    for (i in seq_len(N*thin)) {
      these_points <- purrr::map(all_points, ~.[i,])
      idemc_res <- idemc_step(ems, targets, these_points, all_imps[i,],
                              ladder, clusters_list, order_active, ranges,
                              M = M, pm = pm, w = w)
      all_imps[i+1,] <- idemc_res$imps
      for (j in seq_along(idemc_res$points))
        all_points[[j]][i+1,] <- idemc_res$points[[j]]
    }
  }
  else {
    progressr::with_progress({
      prog <- progressr::progressor(steps = N*thin)
      for (i in seq_len(N*thin)) {
        these_points <- purrr::map(all_points, ~.[i,])
        idemc_res <- idemc_step(ems, targets, these_points, all_imps[i,],
                                   ladder, clusters_list, order_active, ranges,
                                   M = M, pm = pm, w = w)
        all_imps[i+1,] <- idemc_res$imps
        for (j in seq_along(idemc_res$points))
          all_points[[j]][i+1,] <- idemc_res$points[[j]]
        prog(message = sprintf("IDEMC Generation Step %g", i))
      }
    })
  }
  all_points <- purrr::map(all_points, ~.[-1,])
  all_imps <- all_imps[-1]
  if (thin > 1) {
    all_points <- purrr::map(all_points, ~.[seq(thin, N*thin, by = thin),])
    all_imps <- all_imps[seq(thin, N*thin, by = thin),]
  }
  if (detailed) return(all_points)
  return(all_points[[length(all_points)]])
} #nocov end

