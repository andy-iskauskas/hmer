#' Data Cleaning for Stochastic Emulators
#'
#' Reforms the data to allow for meaningful diagnostics on variance emulators.
#'
#' Takes a data.frame of points (with replicates), and generates a set of
#' meaningful statistics based on them and dependent on the type of emulator.
#'
#' @param data The data to reform
#' @param in_names The names of the input parameters
#' @param out_name The output value
#' @param is.variance Is the emulator a mean or variance emulator?
#' @param boots The number of samples for bootstrapping (if required)
#'
#' @importFrom stats sd
#'
#' @keywords internal
#' @noRd
#'
#' @return A data.frame consisting of the summary stats (mean, sd, n)
clean_data <- function(data, in_names, out_name, is.variance = FALSE, boots = 1000) {
  unique_points <- unique(data[,in_names])
  summary_stats <- purrr::map(1:nrow(unique_points), function(i) {
    relev <- data[purrr::map_lgl(1:nrow(data), ~all(data[.,in_names] == unique_points[i,])),out_name]
    if (!is.variance)
      return(c(mean(relev), sd(relev)^2/length(relev), length(relev)))
    var_mean <- sd(relev)^2
    bootstrap <- purrr::map_dbl(1:boots, function(a) {
      repsamp <- sample(relev, length(relev), replace = TRUE)
      sd(repsamp)^2
    })
    return(c(var_mean, sd(bootstrap)^2, length(relev)))
  })
  return(setNames(data.frame(do.call('rbind', summary_stats)), c('mean', 'var', 'n')))
}

k_fold_measure <- function(em, target = NULL, k = 1, ...) {
  train_data <- setNames(data.frame(cbind(eval_funcs(scale_input, data.frame(em$in_data), em$ranges, FALSE), em$out_data)), c(names(em$ranges), em$output_name))
  ordering <- 1:nrow(train_data)
  if (k > 1) {
    if (nrow(train_data) %% k != 0) k <- 1
    else {
      ordering <- sample(1:nrow(train_data), nrow(train_data))
      loo_data <- purrr::map(1:(nrow(train_data)/k), function(x) {
        relev_em <- em$o_em$clone()
        relev_pts <- ordering[((x-1)*k+1):(k*x)]
        adjust_data <- train_data[-relev_pts,]
        adj_em <- relev_em$adjust(adjust_data, em$output_name)
        outp <- cbind(cbind(train_data[relev_pts,], adj_em$get_exp(train_data[relev_pts, names(em$ranges)])), adj_em$get_cov(train_data[relev_pts, names(em$ranges)]))
        if (!is.null(target)) outp <- c(outp, adj_em$implausibility(train_data[relev_pts, names(em$ranges)], target))
        return(outp)
      })
    }
  }
  if (k == 1) {
    loo_data <- purrr::map(1:nrow(train_data), function(x) {
      relev_em <- em$o_em$clone()
      adjust_data <- train_data[-x,]
      adj_em <- relev_em$adjust(adjust_data, em$output_name)
      outp <- c(train_data[x,], adj_em$get_exp(train_data[x,names(em$ranges)]), adj_em$get_cov(train_data[x,names(em$ranges)]))
      if (!is.null(target)) outp <- c(outp, adj_em$implausibility(train_data[x, names(em$ranges)], target))
      return(outp)
    })
  }
  if (!is.null(target)) nms <- c(names(em$ranges), em$output_name, "E", "V", "I")
  else nms <- nms <- c(names(em$ranges), em$output_name, "E", "V")
  out_df <- setNames(data.frame(apply(data.frame(do.call('rbind', loo_data)), 2, unlist)), nms)
  return(out_df[order(as.numeric(row.names(out_df))),])
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
#' @param emulator The emulator to test
#' @param validation The validation set, consisting of points and output(s)
#' @param verbose Should statistics be printed out?
#'
#' @return Whether the observed value lies within 3-sigma of the expected value.
#'
#' @export
#'
#' @examples
#'  summary_diag(sample_emulators$ems$nR, GillespieValidation)
summary_diag <- function(emulator, validation, verbose = FALSE) {
  points <- validation[,names(emulator$ranges)]
  outputs <- validation[,emulator$output_name]
  m <- nrow(validation)
  n <- nrow(emulator$in_data)
  q <- length(emulator$basis_f)
  indiv_errs <- outputs - emulator$get_exp(points)
  chi_sq_measure <- sum(indiv_errs^2/emulator$get_cov(points))
  print(paste("Chi-squared:", round(chi_sq_measure,4), "against mean", m, "with standard deviation", round(sqrt(2*m),4)))
  cov_mat <- emulator$get_cov(points, full = TRUE)
  cov_inv <- tryCatch(
    chol2inv(chol(cov_mat)),
    error = function(e) {MASS::ginv(cov_mat)}
  )
  mahal_measure <- t(indiv_errs) %*% cov_inv %*% indiv_errs
  if (verbose) print(paste("Mahalanobis:", round(mahal_measure,4), "against mean", m, "with standard deviation", round(sqrt(2*m*(m+n-q-2)/(n-q-4)),4)))
  chi_valid <- (abs(chi_sq_measure - m)/sqrt(2*m) <= 3)
  mahal_valid <- (abs(mahal_measure - m)/sqrt(2*m*(m+n-q-2)/(n-q-4)) <= 3)
  return(c(chi_valid, mahal_valid))
}

#' Emulator Regression Residuals
#'
#' Plots the emulator residuals.
#'
#' An emulator is composed of two parts: a global regression surface, and a local
#' correlation structure. It can sometimes be informative to examine the residuals
#' of the regression surface on the training set, to determine the extent to which
#' the regression surface is being 'corrected for' by the correlation structure.
#'
#' @param emulator The emulator to consider.
#' @param histogram Should a histogram or a scatter plot be shown? Default: FALSE
#' @param ... Any additional arguments (used internally)
#'
#' @return A set of residuals, standardised by the regression surface residual standard error.
#' @export
#'
#' @examples
#' residual_diag(sample_emulators$ems$nS)
#' residual_diag(sample_emulators$ems$nI, TRUE)
#'
residual_diag <- function(emulator, histogram = FALSE, ...) {
  in_points <- eval_funcs(scale_input, data.frame(emulator$in_data), emulator$ranges, FALSE)
  if (is.numeric(emulator$u_sigma))
    standardised_residuals <- emulator$model$residuals/emulator$u_sigma
  else
    standardised_residuals <- emulator$model$residuals/apply(in_points, 1, emulator$u_sigma)
  if (histogram) hist(standardised_residuals, xlab = "Standadized Residual", ylab = "Frequency", main = paste(emulator$output_name, "Residual Histogram"))
  else plot(standardised_residuals, xlab = "Point", ylab = "Standardized Residual", pch = 16, main = emulator$output_name, panel.first = abline(h = c(-3, 0, 3), lty = c(2, 1, 2)))
  return(in_points[abs(standardised_residuals) > 3,])
}

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
#' If a validation dataset is not provided, then cross-validation is performed using the
#' training set. By default, leave-one-out cross-validation is used; the parameter \code{k}
#' can be provided as an optional argument to this function to perform k-fold cross validation
#' (where k must be a multiple of the dataset size).
#'
#' @importFrom graphics hist
#'
#' @param emulator An \code{\link{Emulator}} object.
#' @param targets The output targets (to check the relevance of failed points).
#' @param validation The validation data.frame, with all inputs and outputs.
#' @param plt Should the histogram be produced?
#' @param ... Dummy parameters for compatibility with the diagnostic wrapper
#'
#' @return A data.frame of points that fail the diagnostic.
#' @export
#'
#' @examples
#' standard_errors(sample_emulators$ems$nS, sample_emulators$targets, GillespieValidation)
#' ## An empty data.frame.
#' standard_errors(sample_emulators$ems$nS, sample_emulators$targets)
#' # leave-one-out cross-validation
#'
standard_errors <- function(emulator, targets = NULL, validation = NULL, plt = TRUE, ...) {
  if ("Hierarchical" %in% class(emulator)) {
    if (is.null(validation)) stop("LOO diagnostics not (yet) supported for hierarchical emulators.")
    cleaned <- list(...)[['cleaned']]
    if (!is.null(cleaned)) cleaned_data <- cleaned
    else cleaned_data <- clean_data(validation, names(emulator$ranges), emulator$output_name, emulator$em_type == "variance")
    input_points <- unique(validation[,names(emulator$ranges)])
    em_exps <- emulator$get_exp(input_points)
    em_vars <- emulator$get_cov(input_points)
    output_points <- cleaned_data$mean
    num <- em_exps - cleaned_data$mean
    denom <- sqrt(em_vars + cleaned_data$var + emulator$disc$internal^2 + emulator$disc$external^2)
    errors <- num/denom
  }
  else {
    if (is.null(validation)) {
      dat <- k_fold_measure(emulator, ...)
      input_points <- dat[,names(emulator$ranges)]
      output_points <- dat[,emulator$output_name]
      errors <- (dat[,"E"] - output_points)/sqrt(dat[,"V"] + emulator$disc$internal^2 + emulator$disc$external^2)
    }
    else {
      input_points <- validation[,names(emulator$ranges)]
      output_points <- validation[,emulator$output_name]
      errors <- (emulator$get_exp(input_points) - output_points)/sqrt(emulator$get_cov(input_points) + emulator$disc$internal^2 + emulator$disc$external^2)
    }
  }
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
#' If a validation dataset is not provided, then cross-validation is performed using the
#' training set. By default, leave-one-out cross-validation is used; the parameter \code{k}
#' can be provided as an optional argument to this function to perform k-fold cross validation
#' (where k must be a multiple of the dataset size).
#'
#' @importFrom graphics abline arrows plot
#'
#' @param emulator An \code{\link{Emulator}} object.
#' @param targets The output targets (to check the relevance of failed points).
#' @param validation The validation data.frame containing all inputs and outputs.
#' @param sd The number of allowed standard deviations (default 3).
#' @param plt Should the plot be produced?
#' @param ... Dummy parameters for compatibility with the diagnostic wrapper.
#'
#' @return A data.frame of points that fail the diagnostic.
#' @export
#'
#' @examples
#' comparison_diag(sample_emulators$ems$nS, sample_emulators$targets, GillespieValidation)
#' ## An empty data.frame.
#' comparison_diag(sample_emulators$ems$nS, sample_emulators$targets, GillespieValidation,
#'  sd = 0.5)
#' ## A data.frame containing one point.
comparison_diag <- function(emulator, targets = NULL, validation = NULL, sd = 3, plt = T, ...) {
  if ("Hierarchical" %in% class(emulator)) {
    if (is.null(validation)) stop("LOO diagnostics not (yet) supported for hierarchical emulators.")
    cleaned <- list(...)[['cleaned']]
    if (!is.null(cleaned)) cleaned_data <- cleaned
    else cleaned_data <- clean_data(validation, names(emulator$ranges), emulator$output_name, emulator$em_type == "variance")
    input_points <- unique(validation[,names(emulator$ranges)])
    emulator_exp <- emulator$get_exp(input_points)
    em_vars <- emulator$get_cov(input_points)
    emulator_unc <- sd * sqrt(em_vars + cleaned_data$var + emulator$disc$internal^2 + emulator$disc$external^2)
    output_points <- cleaned_data$mean
  }
  else {
    if (is.null(validation)) {
      dat <- k_fold_measure(emulator, ...)
      input_points <- dat[,names(emulator$ranges)]
      output_points <- dat[,emulator$output_name]
      emulator_exp <- dat[,"E"]
      emulator_unc <- sd * sqrt(dat[,"V"] + emulator$disc$internal^2 + emulator$disc$external^2)
    }
    else {
      input_points <- validation[,names(emulator$ranges)]
      output_points <- validation[,emulator$output_name]
      emulator_exp <- emulator$get_exp(input_points)
      emulator_unc <- sd * sqrt(emulator$get_cov(input_points) + emulator$disc$internal^2 + emulator$disc$external^2)
    }
  }
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
#' @param targets The targets to compare suitability against.
#' @param validation The validation data.frame containing all inputs and outputs.
#' @param cutoff The cutoff for the implausibility measures.
#' @param plt Should the plot be produced?
#' @param ... Dummy parameters for compatibility with the diagnostic wrapper.
#'
#' @return A data.frame of points that fail the diagnostic check.
#' @export
#'
#' @examples
#' classification_diag(sample_emulators$ems$nS, sample_emulators$targets, GillespieValidation)
#' ## An empty data.frame.
#' classification_diag(sample_emulators$ems$nS, sample_emulators$targets, GillespieValidation,
#'  cutoff = 1)
#' ## An empty data.frame.
classification_diag <- function(emulator, targets, validation = NULL, cutoff = 3, plt = T, ...) {
  this_target <- targets[[emulator$output_name]]
  if ("Hierarchical" %in% class(emulator)) {
    if (is.null(validation)) stop("LOO diagnostics not (yet) supported for hierarchical emulators.")
    input_points <- unique(validation[,names(emulator$ranges)])
    cleaned <- list(...)[['cleaned']]
    if (!is.null(cleaned)) output_points <- cleaned$mean
    else output_points <- clean_data(validation, names(emulator$ranges), emulator$output_name, emulator$em_type == "validation")$mean
    em_imp <- emulator$implausibility(input_points, this_target)
  }
  else {
    if (is.null(validation)) {
      dat <- k_fold_measure(emulator, this_target, ...)
      input_points <- dat[,names(emulator$ranges)]
      output_points <- dat[,emulator$output_name]
      em_imp <- dat[,"I"]
    }
    else {
      input_points <- validation[,names(emulator$ranges)]
      output_points <- validation[,emulator$output_name]
      em_imp <- emulator$implausibility(input_points, this_target)
    }
  }
  if (is.atomic(this_target)) {
    t_cutoff <- 1 * cutoff/3
    sim_imp <- abs(output_points - rep(mean(this_target), length(output_points)))/rep(diff(this_target)/2, length(output_points))
  }
  else {
    t_cutoff <- cutoff
    sim_imp <- purrr::map_dbl(output_points, ~sqrt((this_target$val-.)^2/this_target$sigma^2))
  }
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
#' All the diagnostics here can be performed with or without a validation (or 'holdout')
#' set of data. The presence of a set of targets is optional for some checks but
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
#'   Regression Residuals (re)
#'
#'   All of the above (all)
#'
#' For details of each of the tests, see the help files for \code{\link{standard_errors}},
#' \code{\link{comparison_diag}}, \code{\link{classification_diag}} and
#' \code{\link{residual_diag}} respectively.
#'
#' @importFrom graphics par
#'
#' @param emulators A list of \code{\link{Emulator}} objects.
#' @param targets The list of observations for the outputs
#' @param validation The validation set, containing all inputs and outputs.
#' @param which_diag Which diagnostics should be performed (see description)
#' @param ... Any additional parameters to pass to the diagnostics (eg sd, cutoff, ...)
#'
#' @return A data.frame containing points that failed one or more diagnostic tests.
#'
#' @export
#'
#' @examples
#' validation_diagnostics(sample_emulators$ems, sample_emulators$targets, GillespieValidation)
#' # data.frame of failed points and a 3x3 set of plots
#' validation_diagnostics(sample_emulators$ems, sample_emulators$targets, GillespieValidation,
#'  c('ce','cd'))
#' # data.frame and a 3x2 set of plots
#' validation_diagnostics(sample_emulators$ems, sample_emulators$targets, GillespieValidation,
#'  cutoff = 2, sd = 2)
#' # k-fold (with k = 3)
#' validation_diagnostics(sample_emulators$ems, sample_emulators$targets, k = 3)
validation_diagnostics <- function(emulators, targets = NULL, validation = NULL, which_diag = c('se', 'cd', 'ce'), ...) {
  if ("Emulator" %in% class(emulators))
    emulators <- setNames(list(emulators), emulators$output_name)
  if (is.null(validation) && any(purrr::map_lgl(emulators, ~"Hierarchical" %in% class(.)))) {
    stop("LOO diagnostics not (yet) supported for hierarchical emulators.")
  }
  cleaning_dat <- purrr::map(emulators, function(x) {
    if ("Hierarchical" %in% class(x)) valid_dat <- clean_data(validation, names(x$ranges), x$output_name, x$em_type == "validation")
    else valid_dat <- NULL
    valid_dat
  })
  fail_point_list <- list()
  if (length(which_diag) == 1 && which_diag == 'all') actual_diag <- c('se', 'cd', 'ce', 're')
  else {
    actual_diag <- which_diag[which_diag %in% c('se', 'cd', 'ce', 're')]
    if (length(which_diag) != length(actual_diag)) warning(paste("Unrecognised diagnostics:", paste0(which_diag[!which_diag %in% c('se', 'cd', 'ce')], collapse = ", "), "\n\tValid diagnostic labels are cd, se, ce or all."))
  }
  mf <- length(actual_diag)
  if (('ce' %in% actual_diag) && is.null(targets)) {
    warning("No targets provided. Cannot perform classification diagnostics.")
    actual_diag <- actual_diag[-which(actual_diag == 'ce')]
  }
  op <- par(mfrow = c(3, mf))
  for (i in 1:length(emulators)) {
    if ('cd' %in% actual_diag) fail_point_list[[length(fail_point_list)+1]] <- comparison_diag(emulators[[i]], targets, validation, cleaned = cleaning_dat[[i]], ...)
    if ('ce' %in% actual_diag) fail_point_list[[length(fail_point_list)+1]] <- classification_diag(emulators[[i]], targets, validation, cleaned = cleaning_dat[[i]], ...)
    if ('se' %in% actual_diag) fail_point_list[[length(fail_point_list)+1]] <- standard_errors(emulators[[i]], targets, validation, cleaned = cleaning_dat[[i]], ...)
    if ('re' %in% actual_diag) fail_point_list[[length(fail_point_list)+1]] <- residual_diag(emulators[[i]], ...)
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
#' Q-Q plots are not plotted for a non-decomposed set of errors, as the correlation
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
        P <- diag(1, nrow = nrow(Q))[,attr(Q, 'pivot')]
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
  output <- data.frame(variable = x_vals, error = indiv_errors)
  return(output)
}
