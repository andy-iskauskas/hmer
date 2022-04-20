idemc_step <- function(ems, targets, points, point_imps, ladder, clusters, #nocov start
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
    for (j in seq_along(points)) {
      if (j == 1) {
        proposal[[j]] <- setNames(
          data.frame(
            matrix(
              purrr::map_dbl(
                clusters[[j]]$ranges,
                ~runif(1, .[[1]], .[[2]])), nrow = 1)), input_names)
        prop_imps[[j]] <- nth_implausible(ems, proposal[[j]],
                                          targets, max_imp = Inf)
      }
      else {
        for (i in 1:M) {
          which_clust <- predict(clusters[[j]]$cluster,
                                 points[[j]])$classification
          clust_mean <- clusters[[j]]$cluster$parameters$mean[,which_clust]
          clust_var <- clusters[[j]]$cluster$parameters$variance$sigma[,,which_clust]
          all_mean <- clusters[[j]]$total$mean
          all_var <- clusters[[j]]$total$sigma
          if (runif(1) < w)
            new_point <- mvtnorm::rmvnorm(1, mean = clust_mean,
                                          sigma = clust_var)
          else
            new_point <- mvtnorm::rmvnorm(1, mean = all_mean, sigma = all_var)
          new_clust <- predict(clusters[[j]]$cluster, new_point)$classification
          new_var <- clusters[[j]]$cluster$parameters$variance$sigma[,,new_clust]
          new_imp <- nth_implausible(ems,
                                     setNames(
                                       data.frame(matrix(new_point, nrow = 1)),
                                       input_names),
                                     targets, max_imp = Inf)
          if (new_imp <= ladder[[j]] && in_range(new_point, ranges)) {
            accept_num <- (w * mclust::dmvnorm(points[[j]],
                                               mean = new_point,
                                               sigma = clust_var) +
                             (1-w)*mclust::dmvnorm(points[[j]],
                                                   mean = new_point,
                                                   sigma = all_var))
            accept_denom <- (w * mclust::dmvnorm(new_point,
                                                 mean = points[[j]],
                                                 sigma = new_var) +
                               (1-w)*mclust::dmvnorm(new_point,
                                                     mean = points[[j]],
                                                     sigma = all_var))
            prob_accept <- min(accept_num/accept_denom, 1)
            if(!is.numeric(prob_accept)) prob_accept <- -1
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
    for (i in 1:ncross) {
      index1 <- sample(1:n, 1, prob = purrr::map_dbl(1:n, ~2*./(n*(n+1))))
      if (n == 2)
        index2 <- ifelse(index1 == 1, 2, 1)
      else
        index2 <- sample((1:n)[-index1], 1,
                         prob = purrr::map_dbl((1:n)[-index1],
                                               ~(n+1-.)/(n*(n+1)/2+index1-n-1)))
      x1 <- proposal[[index1]][,order_active]
      x2 <- proposal[[index2]][,order_active]
      c_point <- sample(1:(length(proposal[[1]])-1), 1)
      y1 <- unlist(c(x1[1:c_point], x2[(c_point+1):length(x2)]),
                   use.names = FALSE)
      y2 <- unlist(c(x2[1:c_point], x1[(c_point+1):length(x1)]),
                   use.names = FALSE)
      imps <- nth_implausible(ems,
                              setNames(data.frame(rbind(y1, y2)), input_names),
                              targets, max_imp = Inf)
      if (imps[[1]] <= ladder[[index1]] &&
          imps[[2]] <= ladder[[index2]] &&
          in_range(y1, ranges[order_active]) &&
          in_range(y2, ranges[order_active])) {
        proposal[[index1]] <- setNames(data.frame(matrix(y1, nrow = 1)),
                                       order_active)[,input_names]
        proposal[[index2]] <- setNames(data.frame(matrix(y2, nrow = 1)),
                                       order_active)[,input_names]
        prop_imps[[index1]] <- imps[[1]]
        prop_imps[[index2]] <- imps[[2]]
      }
    }
  }
  ## Do the exchanges
  for (i in 1:nex) {
    index1 <- sample(n, 1)
    if (index1 == 1) index2 <- 2
    else if (index1 == n) index2 <- n-1
    else index2 <- index1 + sample(c(-1, 1), 1)
    if (index1 > index2) {
      temp <- index1
      index1 <- index2
      index2 <- temp
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
#' This method for generating points is focused on finding non-implausible regions
#' that are either extremely small relative to the initial space, are have interesting
#' structure (particularly disconnected structure) that would potentially be overlooked
#' by more standard point generation methods. The method is robust but computationally
#' intensive, and should not be used as a default - for more standard methods of finding
#' points, see \code{\link{generate_new_runs}} and the methods therein.
#'
#' The IDEMC method operates on an `implausibility ladder', in the spirit of annealing
#' methods. Each `rung' of the ladder is characterised by within-cluster and overall
#' variance. The stages performed in one step of the evolutionary algorithm are as follows:
#'
#' Mutation: A point is modified using a process akin to a random walk step. The parameters
#' that determine the walk can be a global step (determined by the second-order quantities of
#' the entire rung) or a within-cluster step (where the point's cluster is determined and the
#' second-order quantities are drawn from that particular cluster). Local, within-cluster,
#' moves are chosen with probability \code{w}. The move is retained if the new point satisfies
#' the constraints of its rung.
#'
#' Crossover: Points are reorganised in descending order of how active their variables are in
#' the emulated outputs, and two different rungs are selected. The points are 'mixed' using a
#' one-point crossover: given a randomly selected index k and two points x1, x2, the new points
#' are y1 = (x11, x12, ..., x1k, x2(k+1), ... x2n) and similarly for y2. The move is retained
#' if both new points satisfy the constraints on their respective rungs. Choices of rung where
#' the first is a later (more restrictive) rung are favoured.
#'
#' Exchange: Two adjacent rungs are chosen and their points swapped wholesale. The move is
#' retained if the point coming from the less restrictive rung satisfies the constraints of the
#' rung it moves to.
#'
#' At a given step, only one of mutation or crossover is performed: the probability of performing
#' mutation is given by \code{pm}. If mutation is chosen, \code{M} such moves are performed on
#' each rung; if crossover is chosen, then \code{(n+1)/2} such moves are performed across the
#' \code{n} rungs. Exchange is always performed and \code{n+1} such moves are performed.
#'
#' The choice of `implausibility ladder' and clusters has a large bearing on the results. This
#' function performs a `burn-in' to determine a reasonable ladder by starting with a uniform
#' sample across the whole space and defining the next rung by demanding that a percentage
#' (determined by \code{p}) of the original points satisfy the constraint of this new rung.
#' The IDEMC process is performed on these two rungs to generate \code{s} points, from which
#' the process is repeated. Once the desired implausibility has been reached, \code{sn} steps
#' of the algorithm are performed on all rungs to determine final clusters.
#'
#' @importFrom stats predict
#'
#' @param ems The emulators to evaluate implausibility on
#' @param N The desired number of points to generate
#' @param targets The target values for the outputs
#' @param cutoff The desired final implausibility cutoff
#' @param s The number of points to generate at each intermediate burn-in stage
#' @param sn The number of points to generate in the final burn-in stage
#' @param p The proportion of space that should remain in a move along the ladder
#' @param thin The thinning factor: a factor T means that N*T points are generated to get N
#' @param pm The probability that mutation is chosen in an IDEMC step
#' @param w The probability of local random walk moves in the mutation step
#' @param M The number of mutations to perform in an IDEMC step, if chosen
#' @param detailed If TRUE, points from every ladder rung are returned
#' @param verbose Should information about burn-in be displayed during the process?
#'
#' @return Either a list of points (for each rung), or a single data.frame from the last rung.
#'
#' @seealso \code{\link{generate_new_runs}} for more standard point generation methods
#'
#' @references
#'  Vernon & Williamson (2013) <arXiv:1309.3520>
#'
#' @examples
#'  \donttest{
#'   idemc_res <- idemc(SIREmulators$ems, 500, SIREmulators$targets, s = 250, p = 0.3)
#'  }
#'
#' @export
idemc <- function(ems, N, targets, cutoff = 3, s = max(500, ceiling(N/5)),
                  sn = s, p = 0.4, thin = 1, pm = 0.9, w = 0.8,
                  M = 10, detailed = FALSE, verbose = interactive()) {
  ems <- collect_emulators(ems)
  ranges <- getRanges(ems, FALSE)
  order_active <- names(ranges)[order(
    apply(
      do.call(
        'rbind', purrr::map(ems, ~.$active_vars)), 2, sum), decreasing = TRUE)]
  ## Burn-in
  points <- setNames(
    do.call(
      'cbind.data.frame',
      purrr::map(ranges, ~runif(s, .[[1]], .[[2]]))), names(ranges))
  clusters_list <- list(list(ranges = ranges))
  ladder <- c(Inf)
  imps <- nth_implausible(ems, points, targets, max_imp = Inf)
  imp_cutoff <- sort(imps)[floor(p*s)+1]
  orig_imp <- max(imps)
  point_imps <- c(Inf)
  points_list <- list(points[nrow(points),])
  while(imp_cutoff > cutoff) {
    if (verbose) cat(imp_cutoff, "\n")
    next_points <- points[imps <= imp_cutoff,]
    next_cluster <- Mclust(next_points, G = 1:4,
                           control = emControl(tol = 1e-3), verbose = FALSE)
    next_mean <- apply(next_points, 2, mean)
    next_var <- var(next_points)
    points_list[[length(points_list)+1]] <- next_points[nrow(next_points),]
    all_points <- list()
    for (i in seq_along(points_list)) {
      df <- setNames(
        data.frame(
          matrix(0, nrow = s+1, ncol = length(ranges))), names(ranges))
      df[1,] <- points_list[[i]]
      all_points[[length(all_points)+1]] <- df
    }
    clusters_list[[length(clusters_list)+1]] <- list(cluster = next_cluster,
                                                     total = list(
                                                       mean = next_mean,
                                                       sigma = next_var))
    ladder <- c(ladder, imp_cutoff)
    point_imps <- c(point_imps,
                    imps[imps <= imp_cutoff][sum(imps <= imp_cutoff)])
    all_imps <- matrix(0, nrow = s+1, ncol = length(points_list))
    all_imps[1,] <- point_imps
    for (i in 1:s) {
      these_points <- purrr::map(all_points, ~.[i,])
      idemc_result <- idemc_step(ems, targets, these_points, all_imps[i,],
                                 ladder, clusters_list, order_active, ranges,
                                 M = M, pm = pm, w = w)
      all_imps[i+1,] <- idemc_result$imps
      for (j in seq_along(idemc_result$points)) {
        all_points[[j]][i+1,] <- idemc_result$points[[j]]
      }
    }
    points <- all_points[[length(all_points)]]
    clusters_list <- purrr::map(seq_along(all_points), function(i) {
      if (i == 1) return(list(ranges = ranges))
      return(list(cluster = Mclust(all_points[[i]], G = 1:4,
                                   control = emControl(tol = 1e-3),
                                   verbose = FALSE),
                  total = list(mean = apply(all_points[[i]], 2, mean),
                               sigma = var(all_points[[i]]))))
    })
    imps <- all_imps[,ncol(all_imps)]
    imp_cutoff <- sort(imps)[floor(p*s)+1]
    points_list <- purrr::map(all_points, ~.[nrow(.),])
  }
  if (imp_cutoff < cutoff) imp_cutoff <- cutoff
  if (verbose) cat(imp_cutoff, "\n")
  next_points <- points[imps <= imp_cutoff,]
  next_cluster <- Mclust(next_points, G = 1:4, verbose = FALSE)
  next_mean <- apply(next_points, 2, mean)
  next_var <- var(next_points)
  points_list[[length(points_list)+1]] <- next_points[nrow(next_points),]
  all_points <- list()
  for (i in seq_along(points_list)) {
    df <- setNames(data.frame(matrix(0, nrow = sn+1,
                                     ncol = length(ranges))), names(ranges))
    df[1,] <- points_list[[i]]
    all_points[[length(all_points)+1]] <- df
  }
  clusters_list[[length(clusters_list)+1]] <- list(cluster = next_cluster,
                                                   total = list(
                                                     mean = next_mean,
                                                     sigma = next_var))
  ladder <- c(ladder, imp_cutoff)
  point_imps <- c(point_imps, imps[imps <= imp_cutoff][sum(imps <= imp_cutoff)])
  all_imps <- matrix(0, nrow = sn+1, ncol = length(points_list))
  all_imps[1,] <- point_imps
  for (i in 1:sn) {
    these_points <- purrr::map(all_points, ~.[i,])
    idemc_result <- idemc_step(ems, targets, these_points, all_imps[i,],
                               ladder, clusters_list, order_active, ranges,
                               M = M, pm = pm, w = w)
    all_imps[i+1,] <- idemc_result$imps
    for (j in seq_along(idemc_result$points)) {
      all_points[[j]][i+1,] <- idemc_result$points[[j]]
    }
  }
  points <- all_points[[length(all_points)]]
  imps <- all_imps[,ncol(all_imps)]
  points_list <- purrr::map(all_points, ~.[nrow(.),])
  ## Burn in done: generate points
  if (verbose) cat("Completed burn-in: implausibility ladder is (",
               round(orig_imp, 3), "; ", paste0(round(ladder[-1], 3),
                                                collapse = "; "), ").\n",
               sep = "")
  clusters_list <- purrr::map(seq_along(all_points), function(i) {
    if (i == 1) return(list(ranges = ranges))
    return(list(cluster = Mclust(all_points[[i]], G = 1:4,
                                 control = emControl(tol = 1e-3),
                                 verbose = FALSE),
                total = list(mean = apply(all_points[[i]], 2, mean),
                             sigma = var(all_points[[i]]))))
  })
  if (verbose) cat("Performing final point generation...\n")
  all_points <- list()
  for (i in seq_along(points_list)) {
    df <- setNames(
      data.frame(
        matrix(0, nrow = N*thin+1, ncol = length(ranges))), names(ranges))
    df[1,] <- points_list[[i]]
    all_points[[length(all_points)+1]] <- df
  }
  point_imps <- all_imps[nrow(all_imps),]
  all_imps <- matrix(0, nrow = N*thin+1, ncol = length(points_list))
  all_imps[1,] <- point_imps
  for (i in 1:(N*thin)) {
    these_points <- purrr::map(all_points, ~.[i,])
    idemc_result <- idemc_step(ems, targets, these_points, all_imps[i,],
                               ladder, clusters_list, order_active, ranges,
                               M = M, pm = pm, w = w)
    all_imps[i+1,] <- idemc_result$imps
    for (j in seq_along(idemc_result$points))
      all_points[[j]][i+1,] <- idemc_result$points[[j]]
  }
  all_points <- purrr::map(all_points, ~.[-1,])
  all_imps <- all_imps[-1,]
  if (thin > 1) {
    all_points <- purrr::map(all_points, ~.[seq(thin, N*thin, by = thin),])
    all_imps <- all_imps[seq(thin, N*thin, by = thin),]
  }
  if (detailed) return(all_points)
  return(all_points[[length(all_points)]])
} #nocov end
