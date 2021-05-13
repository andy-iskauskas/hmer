#' Emulator Standardised Errors
#'
#' Finds and plots emulator standardised errors.
#'
#' For an emulator of a simulator function \code{f(x)} and a validation set \code{X}, this
#' finds the standardised errors in the form \code{(f(x)-E[f(x)])/sqrt(Var[f(x)])}, where
#' \code{E[f(x)]} is the emulator expectation, and \code{Var[f(x)]} the emulator uncertainty,
#' at each point \code{x} in \code{X}.
#'
#' If \code{plt = TRUE}, then the histogram of standardised errors is provided. Regardless,
#' the set of points that fail this diagnostic are returned as a data.frame. If the output
#' target is provided as a (val, sigma) pair, then the 'failing' points are checked for
#' relevance: if they are a long way from the target region, then they are discounted.
#'
#' @importFrom graphics hist
#'
#' @param emulator An \code{\link{Emulator}} object.
#' @param validation The validation data.frame, with all inputs and outputs.
#' @param targets The output targets (to check the relevance of failed points).
#' @param plt Should the histogram be produced?
#' @param ... Dummy parameters for compatibility with the diagnostic wrapper
#'
#' @return A data.frame of points that fail the diagnostic.
#' @export
#'
#' @examples
#' standard_errors(sample_emulators$ems$nS, GillespieValidation, sample_emulators$targets)
#' ## An empty data.frame.
#'
standard_errors <- function(emulator, validation, targets = NULL, plt = TRUE, ...) {
  input_points <- validation[,names(emulator$ranges)]
  output_points <- validation[,emulator$output_name]
  errors <- (emulator$get_exp(input_points) - output_points)/sqrt(emulator$get_cov(input_points))
  if (plt) {
    hist(errors, xlab = "Standardised Error", main = emulator$output_name)
  }
  emulator_invalid <- abs(errors) > 3
  if (!is.null(targets)) {
    this_target <- targets[[emulator$output_name]]
    point_invalid <- ((output_points < this_target$val - 6*this_target$sigma) | (output_points > this_target$val + 6*this_target$sigma))
    emulator_invalid <- (!point_invalid & emulator_invalid)
  }
  return(input_points[emulator_invalid,])
}

#' Emulator Comparison Diagnostic
#'
#' Produces a comparison diagnostic of emulator and simulator prediction.
#'
#' The emulator output \code{E[f(x)]} is plotted against the simulator output \code{f(x)},
#' for points in a validation set. Error bars are determined by the emulator standard
#' deviation \code{sqrt(Var[f(x)])}. Points whose emulator expectation bound does not
#' contain the actual simulator output are deemed to fail this diagnostic.
#'
#' If \code{plt = TRUE}, then the equivalent plot is provided. As with
#' \code{\link{standard_errors}}, if a target is provided to the diagnostic, then points
#' whose emulator bound does not encompass the simulator result but lie far from the
#' desired output are discounted.
#'
#' @importFrom graphics abline arrows plot
#'
#' @param emulator An \code{\link{Emulator}} object.
#' @param validation The validation data.frame containing all inputs and outputs.
#' @param targets The output targets (to check the relevance of failed points).
#' @param sd The number of allowed standard deviations (default 3).
#' @param plt Should the plot be produced?
#' @param ... Dummy parameters for compatibility with the diagnostic wrapper.
#'
#' @return A data.frame of points that fail the diagnostic.
#' @export
#'
#' @examples
#' comparison_diag(sample_emulators$ems$nS, GillespieValidation, sample_emulators$targets)
#' ## An empty data.frame.
#' comparison_diag(sample_emulators$ems$nS, GillespieValidation, sample_emulators$targets,
#'  sd = 0.5)
#' ## A data.frame containing one point.
comparison_diag <- function(emulator, validation, targets = NULL, sd = 3, plt = T, ...) {
  input_points <- validation[,names(emulator$ranges)]
  output_points <- validation[,emulator$output_name]
  emulator_exp <- emulator$get_exp(input_points)
  emulator_unc <- sd * sqrt(emulator$get_cov(input_points))
  em_ranges <- range(c(emulator_exp + emulator_unc, emulator_exp - emulator_unc))
  emulator_invalid <- (output_points > emulator_exp + emulator_unc) | (output_points < emulator_exp - emulator_unc)
  if (!is.null(targets)) {
    this_target <- targets[[emulator$output_name]]
    point_invalid <- ((output_points < this_target$val - 6*this_target$sigma) | (output_points > this_target$val + 6*this_target$sigma))
    emulator_invalid <- (!point_invalid & emulator_invalid)
  }
  if (plt) {
    plot(output_points, emulator_exp, pch = 16, col = ifelse(emulator_invalid, 'red', 'black'),
         xlim = range(output_points), ylim = range(em_ranges), xlab = 'f(x)', ylab = 'E[f(x)',
         panel.first = c(abline(a=0, b=1, col = 'green')),
         main = emulator$output_name)
    for (i in seq_along(input_points[,1])) {
      if (emulator_unc[[i]] < 1e-8) next
      arrows(x0 = output_points[[i]], y0 = emulator_exp[[i]] - emulator_unc[[i]],
             x1 = output_points[[i]], y1 = emulator_exp[[i]] + emulator_unc[[i]],
             col = ifelse(emulator_invalid[[i]], 'red', 'blue'),
             code = 3, angle = 90, length = .1)
    }
  }
  which_invalid <- input_points[emulator_invalid,]
  return(which_invalid)
}

#' Classification Diagnostic
#'
#' Checks for emulator misclassifications.
#'
#' Both the emulator implausibility and the simulator implausibility can be computed
#' and compared. Points for which the emulator implausibility is outside the desired cutoff
#' but for which the simulator implausibility is not are misclassified points.
#'
#' If \code{plt = TRUE}, then the emulator and simulator implausibilities are plotted
#' against one another. Misclassified points (i.e. those in the bottom right 'quadrant')
#' are highlighted in red.
#'
#' @importFrom graphics abline plot
#'
#' @param emulator An \code{\link{Emulator}} object.
#' @param validation The validation data.frame containing all inputs and outputs.
#' @param targets The targets, as (val, sigma) pairs.
#' @param cutoff The cutoff for the implausibility measures.
#' @param plt Should the plot be produced?
#' @param ... Dummy parameters for compatibility with the diagnostic wrapper.
#'
#' @return A data.frame of points that fail the diagnostic check.
#' @export
#'
#' @examples
#' classification_diag(sample_emulators$ems$nS, GillespieValidation, sample_emulators$targets)
#' ## An empty data.frame.
#' classification_diag(sample_emulators$ems$nS, GillespieValidation, sample_emulators$targets,
#'  cutoff = 1)
#' ## An empty data.frame.
classification_diag <- function(emulator, validation, targets, cutoff = 3, plt = T, ...) {
  input_points <- validation[,names(emulator$ranges)]
  output_points <- validation[,emulator$output_name]
  this_target <- targets[[emulator$output_name]]
  if (is.numeric(this_target)) z <- list(val = this_target, sigma = 0.001)
  else z <- this_target
  em_imp <- emulator$implausibility(input_points, z)
  sim_imp <- purrr::map_dbl(output_points, ~sqrt((z$val-.)^2/z$sigma^2))
  misclass <- (em_imp > cutoff) & (sim_imp <= cutoff)
  if (plt) {
    plot(em_imp, sim_imp, pch = 16,
         col = ifelse(misclass, 'red', 'black'), xlab = "Emulator Implausibility", ylab = "Simulator Implausibility",
         main = emulator$output_name, panel.first = c(abline(h = cutoff), abline(v = cutoff)))
  }
  which_invalid <- input_points[misclass,]
  return(which_invalid)
}

#' Emulator Diagnostics
#'
#' Performs the standard set of validation diagnostics on emulators.
#'
#' All the diagnostics here assume the existence of a validation (or 'hold-out') dataset
#' from the simulator. The presence of a set of targets is optional for some checks but
#' mandatory for others: the appropriate warnings will be given in the event that some
#' checks cannot be applied.
#'
#' The current options for diagnostics (with the codes for \code{which_diag}) are:
#'
#'   Standardised Errors (se)
#'
#'   Comparison Diagnostics (cd)
#'
#'   Classification Errors (ce)
#'
#'   All of the above (all)
#'
#' For details of each of the tests, see the help files for \code{\link{standard_errors}},
#' \code{\link{comparison_diag}} and \code{\link{classification_diag}} respectively.
#'
#' @importFrom graphics par
#'
#' @param emulators A list of \code{\link{Emulator}} objects.
#' @param validation The validation set, containing all inputs and outputs.
#' @param targets The list of observations for the outputs
#' @param which_diag Which diagnostics should be performed (see description)
#' @param ... Any additional parameters to pass to the diagnostics (eg sd, cutoff, ...)
#'
#' @return A data.frame containing points that failed one or more diagnostic tests.
#'
#' @export
#'
#' @examples
#' validation_diagnostics(sample_emulators$ems, GillespieValidation, sample_emulators$targets)
#' # Empty data.frame and a 3x3 set of plots
#' validation_diagnostics(sample_emulators$ems, GillespieValidation, sample_emulators$targets,
#'  c('ce','cd'))
#' # Empty data.frame and a 3x2 set of plots
#' validation_diagnostics(sample_emulators$ems, GillespieValidation, sample_emulators$targets,
#'  cutoff = 2, sd = 2)
#' # Data.frame with one point in it.
validation_diagnostics <- function(emulators, validation, targets = NULL, which_diag = 'all', ...) {
  if ("Emulator" %in% class(emulators)) {
    emulators <- setNames(list(emulators), emulators$output_name)
    if (!is.null(targets)) {
      if (!is.null(targets$val)) targets <- setNames(list(targets), emulators[[1]]$output_name)
    }
  }
  fail_point_list <- list()
  if (length(which_diag) == 1 && which_diag == 'all') actual_diag <- c('se', 'cd', 'ce')
  else {
    actual_diag <- which_diag[which_diag %in% c('se', 'cd', 'ce')]
    if (length(which_diag) != length(actual_diag)) warning(paste("Unrecognised diagnostics:", paste0(which_diag[!which_diag %in% c('se', 'cd', 'ce')], collapse = ", "), "\n\tValid diagnostic labels are cd, se, ce or all."))
  }
  mf <- length(actual_diag)
  if (('ce' %in% actual_diag) && is.null(targets)) {
    warning("No targets provided. Cannot perform classification diagnostics.")
    actual_diag <- actual_diag[-which(actual_diag == 'ce')]
  }
  op <- par(mfrow = c(3, mf))
  for (i in 1:length(emulators)) {
    if ('cd' %in% actual_diag) fail_point_list[[length(fail_point_list)+1]] <- comparison_diag(emulators[[i]], validation, targets, ...)
    if ('ce' %in% actual_diag) fail_point_list[[length(fail_point_list)+1]] <- classification_diag(emulators[[i]], validation, targets, ...)
    if ('se' %in% actual_diag) fail_point_list[[length(fail_point_list)+1]] <- standard_errors(emulators[[i]], validation, targets, ...)
  }
  par(op)
  failed_points <- unique(do.call('rbind', fail_point_list))
  return(failed_points)
}
