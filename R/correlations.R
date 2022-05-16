#' Distance function
#'
#' Calculates Euclidean distances between points in two data.frames
#'
#' If the number of points is moderate, use the \code{dist} function (to take advantage of speed);
#' when the number of points is large this can try to allocate a large object so we instead use
#' direct calculation using \code{colSums}.
#'
#' @param df1 The first data.frame
#' @param df2 The second data.frame
#'
#' @return The distance matrix
#'
#' @keywords internal
#' @noRd
#'
get_dist <- function(df1, df2) {
  if (ncol(df1) != ncol(df2)) stop("Data frames do not have the same dimension.")
  d <- ncol(df1)
  p <- max(nrow(df1), nrow(df2))
  if (d < (833*p^2-198400*p+144350000)/(55550*p+6910000) || p > 2500) {
    return(apply(df1, 1, function(a) sqrt(colSums((a-t(df2))^2))))
  }
  all_dists <- as.matrix(dist(rbind(df1, df2)))
  dists <- all_dists[((nrow(df1)+1):(nrow(df1)+nrow(df2))), seq_len(nrow(df1))]
  row.names(dists) <- colnames(dists) <- NULL
  return(dists)
}


#' Exponential squared correlation function
#'
#' For points \code{x}, \code{xp} and a correlation length \code{theta}, gives the exponent
#' of the squared distance between \code{x} and \code{xp}, weighted by \code{theta} squared.
#'
#' @param x A data.frame of rows corresponding to position vectors
#' @param xp A data.frame of rows corresponding to position vectors
#' @param hp The hyperparameter theta (correlation length)
#'
#' @return The exponential-squared correlation between x and xp.
#' @export
#'
#' @references Rasmussen & Williams (2005) <ISBN: 9780262182539>
#' @examples
#' exp_sq(data.frame(a=1), data.frame(a=2), list(theta = 0.1))
#' #> 3.720076e-44
#' exp_sq(data.frame(a=1,b=2,c=-1),data.frame(a=1.5,b=2.9,c=-0.7), list(theta = 0.2))
#' #> 3.266131e-13
exp_sq <- function(x, xp, hp) {
  if (is.null(hp$theta)) stop("Correlation length theta must be specified.")
  dists <- get_dist(x/hp$theta, xp/hp$theta)^2
  return(exp(-dists))
}

exp_sq_d <- function(x, xp, hp, xi, xpi = NULL) {
  if (is.null(hp$theta)) stop("Correlation length theta must be specified.")
  diff_1 <- outer(x[,xi, drop = FALSE], xp[,xi, drop = FALSE], "-")[,,,1]
  if (is.null(nrow(diff_1))) diff_1 <- t(diff_1)
  if (is.null(xpi)) return(-2 * t(diff_1) * exp_sq(x, xp, hp)/hp$theta^2)
  diff_2 <- if (is.null(xpi)) diff_1 else outer(x[,xpi, drop = FALSE],
                                                xp[,xpi, drop = FALSE],
                                                "-")[,,,1]
  if (is.null(nrow(diff_2))) diff_2 <- t(diff_2)
  return(2/hp$theta^2 * (if (xi == xpi) 1 else 0) * exp_sq(x, xp, hp) -
           4/hp$theta^4 * t(diff_1) * t(diff_2) * exp_sq(x, xp, hp))
}

#' Matern correlation function
#'
#' For points \code{x}, \code{xp}, and a pair of hyperparameters \code{nu} and \code{theta}, gives
#' the Matern correlation between the two points.
#'
#' At present, only half-integer arguments for nu are supported.
#'
#' @param x A data.frame of rows corresponding to position vectors
#' @param xp A data.frame of rows corresponding to position vectors
#' @param hp The hyperparameters nu (smoothness) and theta (correlation length), as a named list
#'
#' @return The Matern correlation between x and xp.
#' @export
#'
#' @references Rasmussen & Williams (2005) <ISBN: 9780262182539>
#' @examples
#' matern(data.frame(a=1), data.frame(a=2), list(nu = 1.5, theta = 0.1))
#' #> 5.504735e-07
#' matern(data.frame(a=1,b=2,c=-1),data.frame(a=1.5,b=2.9,c=-0.7), list(nu = 1.5, theta = 0.2))
#' #> 0.0009527116
matern <- function(x, xp, hp) {
  if (floor(hp$nu*2) != hp$nu*2 || floor(hp$nu) == hp$nu)
    stop("Matern hyperparameter nu must be half-integer.")
  if (is.null(hp$theta))
    stop("No correlation length theta specified.")
  p <- hp$nu-0.5
  d <- get_dist(x, xp)
  exp(-sqrt(2*p+1)*d/hp$theta) * factorial(p)/factorial(2*p) *
    Reduce('+', purrr::map(0:p, ~factorial(p+.)/(factorial(.)*factorial(p-.)) *
                             (2*sqrt(2*p+1)*d/hp$theta)^(p-.)))
}

matern_d <- function(x, xp, hp, xi, xpi = NULL) {
  if (floor(hp$nu*2) != hp$nu*2 || floor(hp$nu) == hp$nu)
    stop("Matern hyperparameter nu must be half-integer.")
  if (floor(hp$nu) < 1)
    stop("This correlation function is not differentiable.")
  if (floor(hp$nu) < 2 && !is.null(xpi))
    stop("This correlation function is not twice differentiable.")
  if (is.null(hp$theta))
    stop("No correlation length theta specified")
  p <- hp$nu-0.5
  inner_arg <- sqrt(2*p+1) * get_dist(x, xp)/hp$theta
  diff_1 <- outer(xp[,xi, drop = FALSE], x[,xi, drop = FALSE], "-")[,,,1]
  if (is.null(nrow(diff_1))) diff_1 <- t(diff_1)
  if (is.null(xpi)) {
    non_sum <- -4*hp$nu/hp$theta^2 * diff_1 * factorial(p)/factorial(2*p) *
      exp(-inner_arg)
    sum <- Reduce('+', purrr::map(0:(p-1),
                                  ~factorial(p-1+.)/(factorial(.) *
                                                       factorial(p-1-.)) *
                                    (2*inner_arg)^(p-1-.)))
    return(non_sum*sum)
  }
  extra_bit <- if(xi == xpi) 4*hp$nu/hp$theta^2 * factorial(p)/factorial(2*p) *
    exp(-inner_arg) * Reduce('+',
                             purrr::map(0:(p-1),
                                        ~factorial(p-1+.)/(factorial(.) *
                                                             factorial(p-1-.)) *
                                          (2*inner_arg)^(p-1-.))) else 0
  diff_2 <- outer(xp[,xpi, drop = FALSE], x[,xpi, drop = FALSE], "-")[,,,1]
  if (is.null(nrow(diff_2))) diff_2 <- t(diff_2)
  non_sum <- -16*hp$nu^2/hp$theta^4 * diff_1 * diff_2 *
    factorial(p)/factorial(2*p) * exp(-inner_arg)
  sum <- Reduce("+", purrr::map(0:(p-2),
                                ~factorial(p-2+.)/
                                  (factorial(.)*factorial(p-2-.)) *
                                  (2*inner_arg)^(p-2-.)))
  return(extra_bit+non_sum*sum)
}

#' Ornstein-Uhlenbeck correlation function
#'
#' For points \code{x}, \code{xp}, and a hyperparameter \code{theta}, gives
#' the Ornstein-Uhlenbeck correlation between the two points.
#'
#' This correlation function can be seen as a specific case of the Matern correlation function
#' when nu = 1/2.
#'
#' @param x A data.frame of rows corresponding to position vectors
#' @param xp A data.frame of rows corresponding to position vectors
#' @param hp The hyperparameter theta (correlation length) in a named list
#'
#' @return The Ornstein-Uhlenbeck correlation between x and xp.
#' @export
#'
#' @references Rasmussen & Williams (2005) <ISBN: 9780262182539>
#' @examples
#' orn_uhl(data.frame(a=1), data.frame(a=2), list(theta = 0.1))
#' #> 4.539993e-05
#' orn_uhl(data.frame(a=1,b=2,c=-1),data.frame(a=1.5,b=2.9,c=-0.7), list(theta = 0.2))
#' #> 0.00469197
#' orn_uhl(data.frame(a=1,b=1,c=1), data.frame(a=1.2,b=0.9,c=0.6), list(theta = 0.2)) ==
#'  matern(data.frame(a=1,b=1,c=1), data.frame(a=1.2,b=0.9,c=0.6), list(theta = 0.2, nu = 0.5)) #> TRUE
orn_uhl <- function(x, xp, hp) {
  if (is.null(hp$theta))
    stop("No correlation length theta specified")
  dists <- get_dist(x, xp)
  exp(-dists/hp$theta)
}

#' Gamma-exponential correlation function
#'
#' For points \code{x}, \code{xp}, and a pair of hyperparameters \code{gamma} and \code{theta},
#' gives the gamma-exponential correlation between the two points.
#'
#' The gamma-exponential correlation function, for d = |x-x'|, is given by
#' \eqn{\exp(-(d/\theta)^\gamma)}. Gamma must be between 0 (exclusive) and 2 (inclusive).
#'
#' @param x A data.frame of rows corresponding to position vectors
#' @param xp A data.frame of rows corresponding to position vectors
#' @param hp The hyperparameters theta (correlation length) and gamma (exponent), as a named list
#'
#' @return The gamma-exponential correlation between x and xp.
#'
#' @references Rasmussen & Williams (2005) <ISBN: 9780262182539>
#' @export
#'
#' @examples
#' gamma_exp(data.frame(a=1), data.frame(a=2), list(gamma = 1.5, theta = 0.1))
#' #> 1.846727e-14
#' gamma_exp(data.frame(a=1,b=2,c=-1),data.frame(a=1.5,b=2.9,c=-0.7), list(gamma = 1.3, theta = 0.2))
#' #> 0.0001399953
gamma_exp <- function(x, xp, hp) {
  if (hp$gamma > 2 || hp$gamma <= 0)
    stop("Gamma hyperparameter must be between 0 (exclusive) and 2 (inclusive)")
  if (is.null(hp$theta))
    stop("No correlation length theta specified")
  dists <- get_dist(x, xp)
  exp(-(dists/hp$theta)^hp$gamma)
}

#' Rational Quadratic correlation function
#'
#' For points \code{x}, \code{xp}, and a pair of hyperparameters \code{alpha} and \code{theta},
#' gives the rational quadratic correlation between the two points.
#'
#' This correlation function, for d = |x-x'|, has the form
#' \eqn{(1+d^2/(2\alpha\theta^2))^{-\alpha}}, and can be seen as a superposition of exponential-squared
#' correlation functions.
#'
#' @param x A data.frame of rows corresponding to position vectors
#' @param xp A data.frame of rows corresponding to position vectors
#' @param hp The hyperparameters alpha (exponent and scale) and theta (correlation length)
#'
#' @return The rational quadratic correlation between x and xp.
#' @export
#'
#' @references Rasmussen & Williams (2005) <ISBN: 9780262182539>
#' @examples
#' rat_quad(data.frame(a=1), data.frame(a=2), list(alpha = 1.5, theta = 0.1))
#' #> 0.004970797
#' rat_quad(data.frame(a=1,b=2,c=-1),data.frame(a=1.5,b=2.9,c=-0.7), list(alpha = 1.5, theta = 0.2))
#' #> 0.02904466
rat_quad <- function(x, xp, hp) {
  if (is.null(hp$theta))
    stop("No correlation length theta specified")
  if (is.null(hp$alpha))
    stop("No reciprocal power alpha specified")
  dists <- get_dist(x, xp)^2
  (1+dists/(2*hp$alpha*hp$theta^2))^(-hp$alpha)
}

rat_quad_d <- function(x, xp, hp, xi, xpi = NULL) {
  if (is.null(hp$theta))
    stop("No correlation length theta specified")
  if (is.null(hp$alpha))
    stop("No reciprocal power alpha specified")
  dists <- get_dist(x, xp)^2
  diff_1 <- outer(xp[,xi, drop = FALSE], x[,xi, drop = FALSE], "-")[,,,1]
  if (is.null(nrow(diff_1))) diff_1 <- t(diff_1)
  if (is.null(xpi))
    return(-diff_1/hp$theta^2 * (1+dists/(2*hp$alpha*hp$theta^2))^(-hp$alpha-1))
  diff_2 <- outer(xp[,xpi, drop = FALSE], x[,xpi, drop = FALSE], "-")[,,,1]
  if (is.null(nrow(diff_2))) diff_2 <- t(diff_2)
  if(xi == xpi)
    extra_bit <- (1+dists/(2*hp$alpha*hp$theta^2))^(-hp$alpha-1)
  else
    extra_bit <- 0
  return(-(hp$alpha+1)/hp$alpha * diff_1*diff_2/hp$theta^4 *
           (1+dists/(2*hp$alpha*hp$theta^2))^(-hp$alpha-2) + extra_bit)
}

Correlator <- R6::R6Class(
  "Correlator",
  public = list(
    hyper_p = NULL,
    corr_type = NULL,
    corr_name = NULL,
    nugget = NULL,
    initialize = function(corr = 'exp_sq', hp = list(theta = 0.1), nug = 0) {
      self$corr_name <- corr
      self$corr_type <- get(corr)
      self$hyper_p <- hp
      self$nugget <- nug
      invisible(self)
    },
    get_corr = function(x, xp = NULL, actives = TRUE, use.nugget = TRUE) {
      if (is.null(xp) && nrow(x) == 1) return(1)
      if (is.null(xp)) return(self$get_corr(x, x, actives, use.nugget))
      x <- data.matrix(x)
      xp <- data.matrix(xp)
      active <- self$corr_type(x[,actives, drop = FALSE],
                               xp[,actives, drop = FALSE], self$hyper_p)
      if (!use.nugget) return(active)
      extra <- as.matrix(dist(rbind(x, xp)))[(nrow(x)+1):(nrow(x)+nrow(xp)),
                                             seq_len(nrow(x))]
      extra[extra < 1e-10] <- 1
      extra[extra != 1] <- 0
      return((1 - self$nugget) * active + self$nugget * extra)
    },
    get_corr_d = function(x, xp = NULL, p1, p2 = NULL,
                          actives = rep(TRUE, ncol(x))) {
      if (is.null(xp)) return(self$get_corr_d(x, x, p1, p2, actives))
      x <- data.matrix(x)
      xp <- data.matrix(xp)
      if (!actives[p1] || (!is.null(p2) && !actives[p2]))
        return(matrix(0, nrow = nrow(xp), ncol = nrow(x)))
      active_p1 <- sum(actives[1:p1])
      active_p2 <- if (!is.null(p2)) sum(actives[1:p2]) else NULL
      d_func <- get(paste0(self$corr_name, "_d"))
      return(d_func(x[,actives, drop = FALSE], xp[,actives, drop = FALSE],
                    self$hyper_p, active_p1, active_p2))
    },
    set_hyper_p = function(new_hp, nug = self$nugget) {
      new_corr <- self$clone()
      if (is.null(names(new_hp)))
        new_hp <- setNames(new_hp, names(self$hyper_p))
      new_hp <- purrr::map(new_hp, ~.)
      new_corr$hyper_p <- new_hp
      new_corr$nugget <- nug
      return(new_corr)
    },
    get_hyper_p = function() {
      return(self$hyper_p)
    },
    print = function(prepend = NULL, ...) {
      cat(prepend, "Correlation type:", self$corr_name, "\n")
      cat(prepend, "Hyperparameters: ", paste0(names(self$hyper_p), ": ",
                                               round(unlist(self$hyper_p), 4),
                                               collapse = "; "), "\n")
      cat(prepend, "Nugget term:", round(self$nugget, 4), "\n")
    }
  )
)
