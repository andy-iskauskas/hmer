#' Exponential squared correlation function
#'
#' For points \code{x}, \code{xp} and a correlation length \code{theta}, gives the exponent
#' of the squared distance between \code{x} and \code{xp}, weighted by \code{theta} squared.
#'
#' @param x A numeric position vector
#' @param xp A numeric position vector
#' @param hp The hyperparameter theta (correlation length)
#'
#' @return The exponential-squared correlation between x and xp.
#' @export
#' @examples
#' exp_sq(1,2, list(theta = 0.1))
#' #> 3.720076e-44
#' exp_sq(c(1,2,-1),c(1.5,2.9,-0.7), list(theta = 0.2))
#' #> 3.266131e-13
exp_sq <- function(x, xp, hp) {
  return(exp(-sum((x-xp)^2/hp$theta^2)))
}

exp_sq_d <- function(x, xp, hp, xi, xpi = NULL) {
  if (is.null(xpi)) return(-2 * (x[xi] - xp[xi]) * exp_sq(x, xp, hp)/hp$theta^2)
  return(2/hp$theta^2 * (if (xi == xpi) 1 else 0) * exp_sq(x, xp, hp) - 4/hp$theta^4 * (x[xi] - xp[xi]) * (x[xpi] - xp[xpi]) * exp_sq(x, xp, hp))
}

#' Matern correlation function
#'
#' For points \code{x}, \code{xp}, and a pair of hyperparameters \code{nu} and \code{rho}, gives
#' the Matern correlation between the two points.
#'
#' At present, only half-integer arguments for nu are supported.
#'
#' @param x A numeric position vector
#' @param xp A numeric position vector
#' @param hp The hyperparameters nu (smoothness) and theta (correlation length), as a named list
#'
#' @return The Matern correlation between x and xp.
#' @export
#' @examples
#' matern(1, 2, list(nu = 1.5, theta = 0.1))
#' #> 5.504735e-07
#' matern(c(1,2,-1), c(1.5, 2.9,-0.7), list(nu = 1.5, theta = 0.2))
#' #> 0.0009527116
matern <- function(x, xp, hp) {
  if (floor(hp$nu*2) != hp$nu*2 || floor(hp$nu) == hp$nu) stop("Matern hyperparameter nu must be half-integer.")
  p <- hp$nu-0.5
  d <- sqrt(sum((x-xp)^2))
  exp(-sqrt(2*p+1)*d/hp$theta) * factorial(p)/factorial(2*p) * sum(purrr::map_dbl(0:p, ~factorial(p+.)/(factorial(.)*factorial(p-.)) * (2*sqrt(2*p+1)*d/hp$theta)^(p-.)))
}

matern_d <- function(x, xp, hp, xi, xpi = NULL) {
  if (floor(hp$nu*2) != hp$nu*2 || floor(hp$nu) == hp$nu) stop("Matern hyperparameter nu must be half-integer.")
  if (floor(hp$nu) < 1) stop("This correlation function is not differentiable.")
  if (floor(hp$nu) < 2 && !is.null(xpi)) stop("This correlation function is not twice differentiable.")
  p <- hp$nu-0.5
  inner_arg <- sqrt(2*p+1) * sqrt(sum((x-xp)^2))/hp$theta
  if (is.null(xpi)) {
    non_sum <- -4*hp$nu/hp$theta^2 * (x[xi] - xp[xi]) * factorial(p)/factorial(2*p) * exp(-inner_arg)
    sum <- sum(purrr::map_dbl(0:(p-1), ~factorial(p-1+.)/(factorial(.) * factorial(p-1-.)) * (2*inner_arg)^(p-1-.)))
    return(non_sum*sum)
  }
  extra_bit <- if(xi == xpi) 4*hp$nu/hp$theta^2 * factorial(p)/factorial(2*p) * exp(-inner_arg) * sum(purrr::map_dbl(0:(p-1), ~factorial(p-1+.)/(factorial(.) * factorial(p-1-.)) * (2*inner_arg)^(p-1-.))) else 0
  non_sum <- -16*hp$nu^2/hp$theta^4 * (x[xi]-xp[xi]) * (x[xpi]-xp[xpi]) * factorial(p)/factorial(2*p) * exp(-inner_arg)
  sum <- sum(purrr::map_dbl(0:(p-2), ~factorial(p-2+.)/(factorial(.)*factorial(p-2-.)) * (2*inner_arg)^(p-2-.)))
  return(extra_bit+non_sum*sum)
}

#' Ornstein-Uhlenbeck correlation function
#'
#' For points \code{x}, \code{xp}, and a pair of hyperparameters \code{nu} and \code{rho}, gives
#' the Ornstein-Uhlenbeck correlation between the two points.
#'
#' This correlation function can be seen as a specific case of the Matern correalation function
#' when nu = 1/2
#'
#' @param x A numeric position vector
#' @param xp A numeric position vector
#' @param hp The hyperparameter theta (correlation length) in a named list
#'
#' @return The Ornstein-Uhlenbeck correlation between x and xp.
#' @export
#' @examples
#' orn_uhl(1, 2, list(theta = 0.1))
#' #> 4.539993e-05
#' orn_uhl(c(1,2,-1), c(1.5, 2.9,-0.7), list(theta = 0.2))
#' #> 0.00469197
orn_uhl <- function(x, xp, hp) {
  exp(-sqrt(sum((x-xp)^2))/hp$theta)
}

#' Gamma-exponential correlation function
#'
#' For points \code{x}, \code{xp}, and a pair of hyperparameters \code{nu} and \code{rho}, gives
#' the gamma-exponential correlation between the two points.
#'
#' The gamma-exponential correlation function, for d = |x-x'|, is given by
#' exp(-(d/l)^gamma). Gamma must be between 0 (exclusive) and 2 (inclusive).
#'
#' @param x A numeric position vector
#' @param xp A numeric position vector
#' @param hp The hyperparameters theta (correlation length) and gamma (exponent), as a named list
#'
#' @return The gamma-exponential correlation between x and xp.
#' @export
#' @examples
#' gamma_exp(1, 2, list(gamma = 1.5, theta = 0.1))
#' #> 1.846727e-14
#' gamma_exp(c(1,2,-1), c(1.5, 2.9,-0.7), list(gamma = 1.3, theta = 0.2))
#' #> 6.020197e-05
gamma_exp <- function(x, xp, hp) {
  if (hp$gamma > 2 || hp$gamma <= 0) stop("Gamma hyperparameter must be between 0 (exclusive) and 2 (inclusive)")
  exp(-(sum((x-xp)^2)/hp$theta)^hp$gamma)
}

#' Rational Quadratic correlation function
#'
#' For points \code{x}, \code{xp}, and a pair of hyperparameters \code{nu} and \code{rho}, gives
#' the rational quadratic correlation between the two points.
#'
#' This correlation function, for d = |x-x'|, has the form
#' (1+d^2/(2 alpha theta^2))^(-alpha), and can be seen as a superposition of exponential-squared
#' correlation functions.
#'
#' @param x A numeric position vector
#' @param xp A numeric position vector
#' @param hp The hyperparameters alpha (exponent and scale) and theta (correlation length)
#'
#' @return The rational quadratic correlation between x and xp.
#' @export
#' @examples
#' rat_quad(1, 2, list(alpha = 1.5, rho = 0.1))
#' #> 0.110858
#' rat_quad(c(1,2,-1), c(1.5, 2.9,-0.7), list(alpha = 1.5, rho = 0.2))
#' #> 0.2194183
rat_quad <- function(x, xp, hp) {
  (1+sum((x-xp)^2)/(2*hp$alpha*hp$theta^2))^(-hp$alpha)
}

rat_quad_d <- function(x, xp, hp, xi, xpi = NULL) {
  if (is.null(xpi)) return(-(x[xi]-xp[xi])/hp$theta^2 * (1+sum((x-xp)^2)/(2*hp$alpha*hp$theta^2))^(-hp$alpha-1))
  extra_bit <- if(xi == xpi) (1+sum((x-xp)^2)/(2*hp$alpha*hp$theta^2))^(-hp$alpha-1) else 0
  return(-(hp$alpha+1)/hp$alpha * (x[xi]-xp[xi])*(x[xpi]-xp[xpi])/hp$theta^4 * (1+sum((x-xp)^2)/(2*hp$alpha*hp$theta^2))^(-hp$alpha-2) + extra_bit)
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
      if (is.null(xp)) return(1)
      active <- self$corr_type(x[actives], xp[actives], self$hyper_p)
      extra <- if (sum((x-xp)^2) < 1e-10 && use.nugget) 1 else 0
      return((1 - self$nugget) * active + self$nugget * extra)
    },
    get_corr_d = function(x, xp = NULL, p1, p2 = NULL, actives = rep(TRUE, length(x))) {
      if (!actives[p1] || (!is.null(p2) && !actives[p2]) || is.null(xp)) return(0)
      active_p1 <- sum(actives[1:p1])
      active_p2 <- if (!is.null(p2)) sum(actives[1:p2]) else NULL
      d_func <- get(paste0(self$corr_name, "_d"))
      return(d_func(x[actives], xp[actives], self$hyper_p, active_p1, active_p2))
    },
    set_hyper_p = function(new_hp, nug = self$nugget) {
      new_corr <- self$clone()
      if (is.null(names(new_hp))) new_hp <- setNames(new_hp, names(self$hyper_p))
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
      cat(prepend, "Hyperparameters: ", paste0(names(self$hyper_p), ": ", round(unlist(self$hyper_p), 4), collapse = "; "), "\n")
      cat(prepend, "Nugget term:", round(self$nugget, 4), "\n")
    }
  )
)
