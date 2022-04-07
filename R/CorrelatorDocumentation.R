# ------------------------
# Correlator Documentation
# ------------------------
#' @title Correlation Structure
#'
#' @description Creates a correlation structure, with the necessary specifications.
#'
#'     The correlator has three main elements: the type of correlator, the associated
#'     hyperparameters, and the nugget term. The nugget term is broadly separate from
#'     the other two parameters, being type-independent.
#'
#' @name Correlator
#'
#' @section Constructor: \code{Correlator$new(corr, hp, nug)}
#'
#' @section Arguments:
#'
#'     \code{corr} The type of correlation function. This is provided as a string which
#'     corresponds exactly with a function - the function should take three arguments
#'     \code{x, xp, hp}. This gives a correlation function u(x, xp) defined by
#'     hyperparameters \code{hp}. For a simple example, see \code{\link{exp_sq}}.
#'
#'     \code{hp} The associated hyperparameters needed to define the correlation
#'     structure, as a named list. In the case of \code{exp_sq}, this is a list of
#'     one element, \code{list(theta)}.
#'
#'     \code{nug} The size of the nugget term. In situations where not all variables
#'     are active, the main part of u(x) operates only on the active parts, xA. The
#'     presence of the nugget term accounts for the fact that points at the same
#'     position in the active space need not be at the same position in the full space.
#'
#'     By default, \code{Correlator$new()} initialises with \code{corr = exp_sq},
#'     \code{hp = list(theta = 0.1)}, and \code{nug = 0}.
#'
#' @section Accessor Methods:
#'
#'     \code{get_corr(x, xp = NULL, actives = TRUE)} Returns the correlation
#'     between two points. If \code{xp} is \code{NULL}, then this is correlation
#'     between a set of points and themselves (i.e. 1 on the diagonal). All variables
#'     are assumed to be active unless otherwise stated in \code{actives}.
#'
#'     \code{get_hyper_p()} Returns the list of hyperparameters.
#'
#'     \code{print()} Produces a summary of the correlation structure specification.
#'
#' @section Object Methods:
#'
#'     \code{set_hyper_p(hp, nugget)} Modifies the hyperparameter and/or nugget
#'     terms. Returns a new \code{Correlator} object.
#'
#' @section Options for Correlations:
#'
#'     The default choice (and that supported by other functions in this package, particularly
#'     emulator_from_data) for the correlation structure is exponential-squared, due to the
#'     useful properties it possesses. However, one can manually instantiate a Correlator with
#'     a different underlying structure. Built-in alternatives are as follows, as well as whether
#'     a form exists for its derivative:
#'
#'     \describe{
#'      \item{\code{\link{matern}}}{the Matérn function (derivative exists)}
#'      \item{\code{\link{orn_uhl}}}{the Ornstein-Uhlenbeck function (no derivative)}
#'      \item{\code{\link{rat_quad}}}{the rational quadratic function (derivative exists)}
#'      }
#'
#'     One more function, \code{\link{gamma_exp}}, is available but not directly supported
#'     by \code{emulator_from_data}, for example, due to its very limited suitability to
#'     emulating model outputs. However, this can be used as a test case for writing one's
#'     own correlation functions and using them with \code{emulator_from_data}.
#'
#'     A user-defined correlation function can be provided to the Correlator: the requirements
#'     are that the function accept data.matrix objects as its first and second arguments,
#'     and accept a named list of hyperparameters as its third argument, and return a matrix
#'     of correlations between rows of the data.matrices. If a derivative also
#'     exists, it should take the same name as the correlation function with "_d" appended to
#'     it, and the directions to differentiate with respect to should come after the
#'     hyperparameter argument. For example, the rational quadratic functions have the form
#'
#'      \code{rat_quad(x1, x2, hp = list(alpha, theta))}
#'
#'      \code{rat_quad_d(x1, x2, hp = list(alpha, theta), dx1, dx2)}
#'
#'     If defining a custom correlation function, care should be taken with hyperparameter
#'     estimation - see \code{\link{emulator_from_data}} examples for details.
#'
#' @export
#'
#' @examples
#' test_corr <- Correlator$new(nug = 0.1)
#' test_corr
#' point1 <- data.frame(a = 0.1, b = 0.2, c = 0.3)
#' point2 <- data.frame(a = 0.15, b = 0.18, c = 0.295)
#' test_corr$get_corr(point1) #> 1
#' test_corr$get_corr(point1, point2) #> 0.6717557
#' test_corr$get_corr(point1, point2, actives = c(TRUE, TRUE, FALSE)) #> 0.6734372
#'
#' new_corr <- test_corr$set_hyper_p(list(theta = 0.5), nug = 0.01)
#' new_corr$get_corr(point1, point2) #> 0.9784845
#' new_corr$get_corr(point1, point2, actives = c(TRUE, TRUE, FALSE)) #> 0.9785824
#'
#' mat_corr <- Correlator$new('matern', list(nu = 1.5, theta = 0.5))
#' mat_corr$get_corr(data.frame(a = c(1, 0.9), b = c(4, 4.2)))
NULL
