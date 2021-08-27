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
  errors <- (emulator$get_exp(input_points) - output_points)/sqrt(emulator$get_cov(input_points) + emulator$disc$internal^2 + emulator$disc$external^2)
  if (plt) {
    hist(errors, xlab = "Standardised Error", main = emulator$output_name)
  }
  emulator_invalid <- abs(errors) > 3
  if (!is.null(targets)) {
    this_target <- targets[[emulator$output_name]]
    if (is.atomic(this_target)) point_invalid <- ((output_points < this_target[1]-diff(this_target)/2) | (output_points > this_target[2]+diff(this_target)/2))
    else point_invalid <- ((output_points < this_target$val - 6*this_target$sigma) | (output_points > this_target$val + 6*this_target$sigma))
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
  emulator_unc <- sd * sqrt(emulator$get_cov(input_points) + emulator$disc$internal^2 + emulator$disc$external^2)
  em_ranges <- range(c(emulator_exp + emulator_unc, emulator_exp - emulator_unc))
  emulator_invalid <- (output_points > emulator_exp + emulator_unc) | (output_points < emulator_exp - emulator_unc)
  if (!is.null(targets)) {
    this_target <- targets[[emulator$output_name]]
    if (is.atomic(this_target)) point_invalid <- ((output_points < this_target[1]-diff(this_target)/2) | (output_points > this_target[2]+diff(this_target)/2))
    else point_invalid <- ((output_points < this_target$val - 6*this_target$sigma) | (output_points > this_target$val + 6*this_target$sigma))
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
#' @param targets The targets to compare suitability against.
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
  if (is.atomic(this_target)) {
    t_cutoff <- 1 * cutoff/3
    sim_imp <- abs(output_points - rep(mean(this_target), length(output_points)))/rep(diff(this_target)/2, length(output_points))
  }
  else {
    t_cutoff <- cutoff
    sim_imp <- purrr::map_dbl(output_points, ~sqrt((this_target$val-.)^2/this_target$sigma^2))
  }
  em_imp <- emulator$implausibility(input_points, this_target)
  misclass <- (em_imp > cutoff) & (sim_imp <= t_cutoff)
  if (plt) {
    plot(em_imp, sim_imp, pch = 16,
         col = ifelse(misclass, 'red', 'black'), xlab = "Emulator Implausibility", ylab = "Simulator Implausibility",
         main = emulator$output_name, panel.first = c(abline(h = t_cutoff), abline(v = cutoff)))
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
#' # data.frame of failed points and a 3x3 set of plots
#' validation_diagnostics(sample_emulators$ems, GillespieValidation, sample_emulators$targets,
#'  c('ce','cd'))
#' # data.frame and a 3x2 set of plots
#' validation_diagnostics(sample_emulators$ems, GillespieValidation, sample_emulators$targets,
#'  cutoff = 2, sd = 2)
validation_diagnostics <- function(emulators, validation, targets = NULL, which_diag = 'all', ...) {
  if ("Emulator" %in% class(emulators))
    emulators <- setNames(list(emulators), emulators$output_name)
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

#' Predictive Error Plots
#'
#' Plots the predictive error with respect to a variety of quantities.
#'
#' The choice of errors to plot is controlled by \code{errtype}, and can be one
#' of four things: normal, corresponding to the regular standardised errors; eigen,
#' corresponding to the errors after reordering given by the eigendecomposition
#' of the emulator covariance matrix; chol, similarly deriving errors after Cholesky
#' decomposition; and cholpivot, deriving the errors after pivoted Cholesky decomposition.
#'
#' What the errors are plotted with respect to is controlled by \code{xtype}. The options
#' are index, which plots them in their order in the validation set; em, which plots errors
#' with respect to the emulator prediction at that point; and any named parameter of the
#' model, which plots with respect to the values of that parameter.
#'
#' Finally, the plot type is controlled by \code{plottype}: this can be one of normal,
#' which plots the errors; or qq, which produces a Q-Q plot of the errors.
#'
#' The default output is to plot the standardised errors (with no decomposition)
#' against the ordering in the validation set; i.e. \code{errtype = "normal"},
#' \code{xtype = "index"}, \code{plottype = "normal"}.
#'
#' Some combinations are not permitted, as the output would not be meaningful. Errors
#' arising from an eigendecomposition cannot be plotted against either emulator prediction
#' or a particular parameter (due to the transformation induced by the eigendecomposition);
#' Q-Q plots cannoot be plotted for a non-decomposed set of errors, as the correlation
#' between errors makes it much harder to interpret.
#'
#' @importFrom stats qqnorm qqline
#'
#' @param em The emulator to perform diagnostics on
#' @param validation The validation set of points with output(s)
#' @param errtype The type of individual error to be plotted.
#' @param xtype The value to plot against
#' @param plottype Whether to plot a standard or Q-Q plot.
#'
#' @return The relevant plot.
#'
#' @references L. Bastos & A. O'Hagan: Diagnostics for Gaussian Process Emulators.
#'
#' @export
#'
#' @examples
#' individual_errors(sample_emulators$ems$nS, GillespieValidation)
#' individual_errors(sample_emulators$ems$nS, GillespieValidation, "chol", "em")
#' individual_errors(sample_emulators$ems$nS, GillespieValidation, "eigen", plottype = "qq")
#' individual_errors(sample_emulators$ems$nS, GillespieValidation, "cholpivot", xtype = "aSI")
#'
individual_errors <- function(em, validation, errtype = "normal", xtype = "index", plottype = "normal") {
  if (!errtype %in% c("normal", "eigen", "chol", "cholpivot"))
    stop("Error type not recognised (options are normal, eigen, chol, or cholpivot).")
  if (!xtype %in% c("index", "em", names(em$ranges)))
    stop("x-axis measure not recognised (options are index, em, or a parameter name.")
  if (!plottype %in% c("normal", "qq"))
    stop("Plot type not recognised (options are normal or qq).")
  if (plottype == "qq" && errtype == "normal") {
    warning("Not meaningful to create Q-Q plot with untransformed errors. Changing to pivoted Cholesky.")
    errtype <- "cholpivot"
  }
  if (xtype %in% c(names(em$ranges), 'em') && errtype == "eigen") {
    warning("Not meaningful to plot parameter or emulator prediction against eigendecomposed errors. Changing to pivoted Cholesky.")
    errtype <- "cholpivot"
  }
  points <- validation[,names(em$ranges)]
  outputs <- validation[,em$output_name]
  em_pred <- em$get_exp(points)
  em_cov <- em$get_cov(points, full = TRUE)
  if (errtype == "normal")
    indiv_errors <- (outputs - em_pred)/sqrt(diag(em_cov))
  else {
      if (errtype == "eigen")
        G <- eigen(em_cov)$vectors %*% diag(sqrt(eigen(em_cov)$values))
      else if (errtype == "chol")
        G <- t(chol(em_cov))
      else {
        Q <- chol(em_cov, pivot = TRUE)
        P <- diag(1, nrow = nrow(Q))[attr(Q, 'pivot'),]
        G <- P %*% t(Q)
      }
      G_inv <- tryCatch(
        chol2inv(chol(G)),
        error = function(e) {MASS::ginv(G)}
      )
      indiv_errors <- G_inv %*% (outputs - em_pred)
  }
  if (xtype == "index") x_vals <- 1:length(outputs)
  else if (xtype == "em") x_vals <- em_pred
  else x_vals <- validation[,xtype]
  if (plottype == "normal") {
    x_lab <- switch(xtype, "index" = "Index", "em" = "Emulator Prediction", xtype)
    appended <- switch(errtype, "normal" = "", "eigen" = " (eigendecomposition)", "chol" = " (Choleksy decomposition)", "cholpivot" = " (pivoted Cholesky decomposition)")
    plot(x_vals, indiv_errors, pch = 16, xlab = x_lab, ylab = "Error", main = paste0("Errors against ", x_lab, appended), panel.first = abline(h = c(-2,2), lty = 2))
  }
  if (plottype == "qq") {
    appended <- switch(errtype, "eigen" = " (eigendecomposition)", "chol" = " (Choleksy decomposition)", "cholpivot" = " (pivoted Cholesky decomposition)")
    qqnorm(indiv_errors, pch = 16, main = paste0("Q-Q plot for Errors", appended), panel.first = qqline(indiv_errors, col = "steelblue"))
  }
  return(NULL)
}

#' Summary Statistics for Emulators
#'
#' Generates measures for emulator quality
#'
#' A couple of summary statistics can be generated for emulators, based on their
#' prediction errors on a validation set. This function produces the test statistic
#' for a comparison to a relevant chi-squared distribution, and the similar test
#' statistic for an F-distribution. In both cases, the expectation and standard
#' deviation of the underlying distribution are also provided.
#'
#' The output of this function is a logical vector stating whether the derived
#' value lies within 3-sigma of the expected value. In systems where errors are
#' expected to be correlated, higher weight should be given to the Mahalanobis
#' measure than the chi-squared measure. Any anomalous results can be investigated
#' in more depth using the \code{\link{individual_errors}} function.
#'
#' @param em The emulator to test
#' @param validation The validation set, consisting of points and output(s)
#'
#' @return Whether the observed value lies within 3-sigma of the expected value..
#'
#' @export
#'
#' @examples
#'  summary_diagnostic(sample_emulators$ems$nR, GillespieValidation)
summary_diagnostic <- function(em, validation) {
  points <- validation[,names(em$ranges)]
  outputs <- validation[,em$output_name]
  m <- nrow(validation)
  n <- nrow(em$in_data)
  q <- length(em$basis_f)
  indiv_errs <- outputs - em$get_exp(points)
  chi_sq_measure <- sum(indiv_errs^2/em$get_cov(points))
  print(paste("Chi-squared:", round(chi_sq_measure,4), "against mean", m, "with standard deviation", round(sqrt(2*m),4)))
  cov_mat <- em$get_cov(points, full = TRUE)
  cov_inv <- tryCatch(
    chol2inv(chol(cov_mat)),
    error = function(e) {MASS::ginv(cov_mat)}
  )
  mahal_measure <- t(indiv_errs) %*% cov_inv %*% indiv_errs
  print(paste("Mahalanobis:", round(mahal_measure,4), "against mean", m, "with standard deviation", round(sqrt(2*m*(m+n-q-2)/(n-q-4)),4)))
  #return(list(chi = chi_sq_measure, mahalanobis = as.numeric(mahal_measure)))
  chi_valid <- (abs(chi_sq_measure - m)/sqrt(2*m) <= 3)
  mahal_valid <- (abs(mahal_measure - m)/sqrt(2*m*(m+n-q-2)/(n-q-4)) <= 3)
  return(c(chi_valid, mahal_valid))
}
