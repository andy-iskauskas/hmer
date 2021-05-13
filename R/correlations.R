#' Exponential squared correlation function
#'
#' For points \code{x}, \code{xp} and a correlation length \code{theta}, gives the exponent
#' of the squared distance between \code{x} and \code{xp}, weighted by \code{theta} squared.
#'
#' @param x A numeric position vector
#' @param xp A numeric position vector
#' @param hp The hyperparameter theta (correlation length)
#'
#' @return The exponental-squared correlation between x and xp.
#' @export
#' @examples
#' exp_sq(1,2, list(theta = 0.1))
#' #> 3.720076e-44
#' exp_sq(c(1,2,-1),c(1.5,2.9,-0.7), list(theta = 0.2))
#' #> 3.266131e-13
exp_sq <- function(x, xp, hp) {
  return(exp(-sum((x-xp)^2/hp$theta^2)))
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
    get_corr = function(x, xp = NULL, actives = TRUE) {
      if (is.null(xp)) return(1)
      active <- self$corr_type(x[actives], xp[actives], self$hyper_p)
      extra <- if (sum((x-xp)^2) < 1e-10) 1 else 0
      return((1 - self$nugget) * active + self$nugget * extra)
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
