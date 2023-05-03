# -----------------------------------
# Hierarchical Emulator Documentation
# -----------------------------------
#' @title Hierarchical Bayes Linear Emulator
#'
#' @description Creates a univariate emulator with hierarchical structure.
#'
#'     This object does not differ extensively from the standard \code{\link{Emulator}} object, so
#'     most of the functionality will not be listed here: the main difference is that
#'     it allows for the variance structure of the emulator to be modified by a higher
#'     order object. The typical usage is to create a variance emulator, whose predictions
#'     inform the behaviour of a mean emulator with regard to a stochastic process.
#'
#' @name HierarchicalEmulator
#'
#' @section Constructor: \code{HierarchicalEmulator$new(basis_f, beta, u, ranges, ...)}
#'
#' @section Arguments:
#'
#'    For details of shared arguments, see \code{\link{Emulator}}.
#'
#'    \code{s_diag} The function that modifies the structure of the Bayes Linear adjustment.
#'
#'    \code{samples} A numeric vector that indicates how many replicates each of the training
#'    points has.
#'
#'    \code{em_type} Whether the emulator is emulating a mean surface or a variance surface.
#'
#' @section Constructor Details:
#'
#'    See \code{\link{Emulator}}: the constructor structure is the same save for the
#'    new arguments discussed above.
#'
#' @section Accessor Methods:
#'
#'    \code{get_exp(x, samps = NULL)} Similar in form to the normal Emulator method; the
#'    \code{samps} argument allows the estimation of summary statistics derived from
#'    multiple realisations.
#'
#'    \code{get_cov(x, xp = NULL, full = FALSE, samps = NULL)} Differences here are in
#'    line with those described in \code{get_exp}.
#'
#' @section Object Methods:
#'
#'    Identical to those of \code{\link{Emulator}}: the one internal difference is that
#'    \code{adjust} returns a HierarchicalEmulator rather than a standard one.
#'
#' @references Goldstein & Vernon (2016), in preparation
#' @export
#'
#' @examples
#'  h_em <- emulator_from_data(BirthDeath$training, c('Y'),
#'   list(lambda = c(0, 0.08), mu = c(0.04, 0.13)), emulator_type = "variance")
#'  names(h_em) # c("expectation', 'variance')

NULL
