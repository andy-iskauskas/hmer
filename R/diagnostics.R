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
#' @importFrom rlang hash
#'
#' @keywords internal
#' @noRd
#'
#' @return A data.frame consisting of the summary stats (mean, sd, n)
clean_data <- function(data, in_names, out_name,
                       is.variance = FALSE, boots = 1000) {
  unique_points <- unique(data[,in_names])
  uids <- apply(unique_points, 1, hash)
  summary_stats <- purrr::map(uids, function(i) {
    relev <- data[apply(data[,in_names], 1, hash) == i, out_name]
    if (!is.variance)
      return(c(mean(relev), sd(relev)^2/length(relev), length(relev)))
    var_mean <- sd(relev)^2
    bootstrap <- purrr::map_dbl(1:boots, function(a) {
      repsamp <- sample(relev, length(relev), replace = TRUE)
      sd(repsamp)^2
    })
    return(c(var_mean, sd(bootstrap)^2, length(relev)))
  })
  out_data <- setNames(
    cbind.data.frame(
      unique_points, do.call('rbind', summary_stats)),
    c(in_names, 'mean', 'var', 'n'))
  return(out_data[out_data$n > 1, ])
}

# A function for calculating k-fold cross-validation statistics.
k_fold_measure <- function(em, target = NULL, k = 1, ...) {
  train_data <- setNames(
    cbind.data.frame(
      eval_funcs(
        scale_input, data.frame(em$in_data), em$ranges, FALSE),
      em$out_data), c(names(em$ranges), em$output_name))
  ordering <- seq_len(nrow(train_data))
  if (k > 1) {
    if (nrow(train_data) %% k != 0) k <- 1
    else {
      ordering <- sample(seq_len(nrow(train_data)), nrow(train_data))
      loo_data <- purrr::map(1:(nrow(train_data)/k), function(x) {
        relev_em <- em$o_em$clone()
        relev_pts <- ordering[((x-1)*k+1):(k*x)]
        adjust_data <- train_data[-relev_pts,]
        adj_em <- relev_em$adjust(adjust_data, em$output_name)
        outp <- cbind(
          cbind(
            train_data[relev_pts,],
            adj_em$get_exp(train_data[relev_pts, names(em$ranges)])),
          adj_em$get_cov(train_data[relev_pts, names(em$ranges)]))
        if ("Hierarchical" %in% class(em))
          outp <- cbind(
            outp, purrr::map_dbl(relev_pts,
                                 ~em$s_diag(train_data[., names(em$ranges)],
                                            1)))
        if (!is.null(target))
          outp <- cbind(
            outp, adj_em$implausibility(train_data[relev_pts, names(em$ranges)],
                                        target))
        return(c(outp, use.names = FALSE))
      })
    }
  }
  if (k == 1) {
    loo_data <- purrr::map(seq_len(nrow(train_data)), function(x) {
      relev_em <- em$o_em$clone()
      adjust_data <- train_data[-x,]
      adj_em <- relev_em$adjust(adjust_data, em$output_name)
      outp <- c(train_data[x,], adj_em$get_exp(train_data[x,names(em$ranges)]),
                adj_em$get_cov(train_data[x,names(em$ranges)]))
      if ("Hierarchical" %in% class(em))
        outp <- c(outp, em$s_diag(train_data[x, names(em$ranges)], 1))
      if (!is.null(target))
        outp <- c(outp, adj_em$implausibility(train_data[x, names(em$ranges)],
                                              target))
      return(c(outp, use.names = FALSE))
    })
  }
  res_names <- c(names(em$ranges), em$output_name, "E", "V")
  if ("Hierarchical" %in% class(em)) res_names <- c(res_names, "Vest")
  if (!is.null(target)) res_names <- c(res_names, "I")
  out_df <- setNames(
    data.frame(apply(do.call('rbind.data.frame', loo_data), 2, unlist)),
    res_names)
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
#' @references Bastos & O'Hagan (2009) <doi:10.1198/TECH.2009.08019>
#' @family diagnostic functions
#' @export
#'
#' @examples
#'  summary_diag(SIREmulators$ems$nR, SIRSample$validation)
summary_diag <- function(emulator, validation, verbose = interactive()) {
  if ("EmProto" %in% class(emulator))
    stop("summary_diag not applicable for Proto_emulator objects.")
  points <- validation[,names(emulator$ranges)]
  outputs <- validation[,emulator$output_name]
  m <- nrow(validation)
  n <- nrow(emulator$in_data)
  q <- length(emulator$basis_f)
  indiv_errs <- outputs - emulator$get_exp(points)
  chi_sq_measure <- sum(indiv_errs^2/emulator$get_cov(points))
  if(verbose) cat("Chi-squared:", round(chi_sq_measure,4), #nocov start
                          "against mean", m,
                          "with standard deviation", round(sqrt(2*m),4), "\n") #nocov end
  cov_mat <- emulator$get_cov(points, full = TRUE)
  cov_inv <- tryCatch(
    chol2inv(chol(cov_mat)),
    error = function(e) {MASS::ginv(cov_mat)}
  )
  mahal_measure <- t(indiv_errs) %*% cov_inv %*% indiv_errs
  if (verbose) cat("Mahalanobis:", round(mahal_measure,4), #nocov start
                           "against mean", m, "with standard deviation",
                           round(sqrt(2*m*(m+n-q-2)/(n-q-4)),4), "\n") #nocov end
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
#' the regression surface is being `corrected for' by the correlation structure.
#'
#' @param emulator The emulator to consider.
#' @param histogram Should a histogram or a scatter plot be shown? Default: FALSE
#' @param ... Any additional arguments (used internally)
#'
#' @return A set of residuals, standardised by the regression surface residual standard error.
#'
#' @family diagnostic functions
#' @export
#'
#' @examples
#' residual_diag(SIREmulators$ems$nS)
#' residual_diag(SIREmulators$ems$nI, TRUE)
#'
residual_diag <- function(emulator, histogram = FALSE, ...) {
  if ("EmProto" %in% class(emulator))
    stop("residual_diag not applicable for Proto_emulator objects.")
  in_points <- eval_funcs(scale_input, data.frame(emulator$in_data),
                          emulator$ranges, FALSE)
  if (is.numeric(emulator$u_sigma))
    standardised_residuals <- emulator$model$residuals/emulator$u_sigma
  else
    standardised_residuals <- emulator$model$residuals/
      apply(in_points, 1, emulator$u_sigma)
  if (histogram) hist(standardised_residuals, xlab = "Standadized Residual", #nocov start
                      ylab = "Frequency",
                      main = paste(emulator$output_name, "Residual Histogram")) #nocov end
  else plot(standardised_residuals, xlab = "Point",
            ylab = "Standardized Residual", pch = 16,
            main = emulator$output_name,
            panel.first = abline(h = c(-3, 0, 3), lty = c(2, 1, 2)))
  return(in_points[abs(standardised_residuals) > 3,])
}

#' Percentage of Space Removed
#'
#' For a wave of emulators, estimates the proportion of space removed at this wave.
#'
#' Given a collection of emulators corresponding to a wave, we can look at an estimate of
#' the proportion of points from previous waves that will be accepted at this wave, either
#' on an emulator-by-emulator basis (to see which outputs are most restrictive) or as an all-wave
#' determination.
#'
#' Naturally, such a statement will be an estimate of the restriction on the full space (which will
#' become more unreliable for higher dimensions), but it can give an order-of-magnitude statement,
#' or useful comparators between different emulators in a wave.
#'
#' If no points are provided, the training points for the emulators are used. For best results, a
#' good number of points should be given: typically one should consider using as many points as one
#' knows to be in the NROY space (including any validation points, if accessible).
#'
#' @param ems The emulators to compute over, as a list
#' @param targets The output target values
#' @param points The points to test against
#' @param ppd If no points are provided and uniform grid is wanted, the number of
#'            points per parameter dimension.
#' @param cutoff The cutoff value for implausibility
#' @param individual If true, gives emulator-by-emulator results; otherwise works with maximum implausibility
#'
#' @return A numeric corresponding to the proportions of points removed.
#'
#' @seealso \code{\link{space_removed}} for a visualisation of the space removal.
#' @export
#'
#' @examples
#'  space_removal(SIREmulators$ems, SIREmulators$targets,
#'   rbind(SIRSample$training, SIRSample$validation))
#'  space_removal(SIREmulators$ems, SIREmulators$targets,
#'   rbind(SIRSample$training, SIRSample$validation), individual = FALSE)
space_removal <- function(ems, targets, points = NULL, ppd = NULL,
                          cutoff = 3, individual = TRUE) {
  if (is.null(points)) {
    if (is.null(ppd)) {
      if ("EmProto" %in% class(ems[[1]]))
        stop("Cannot access emulator training points as basis for space removal calculation for proto_ems.")
      input_points <- eval_funcs(scale_input, data.frame(ems[[1]]$in_data),
                                 ems[[1]]$ranges, FALSE)
      output_points <- setNames(do.call('cbind.data.frame',
                                        purrr::map(ems, ~.$out_data)),
                                purrr::map_chr(ems, ~.$output_name))
      points <- data.frame(cbind(input_points, output_points))
    }
    else {
      maxpoints <- 50000
      ranges <- ems[[1]]$ranges
      if (ppd^length(ranges) > maxpoints) {
        points <- setNames(do.call('cbind.data.frame',
                                   purrr::map(ranges, ~runif(maxpoints, .[[1]], .[[2]]))),
                           names(ranges))
      }
      else {
        points <- setNames(
          expand.grid(purrr::map(ranges, ~seq(.[[1]], .[[2]], length.out = ppd))),
          names(ranges))
      }
    }
  }
  if (individual) {
    imps <- setNames(
      do.call(
        'cbind.data.frame',
        purrr::map(ems, ~.$implausibility(points, targets[[.$output_name]],
                                          cutoff = cutoff))),
      purrr::map_chr(ems, ~.$output_name))
    return(1-apply(imps, 2, sum)/nrow(points))
  }
  nth_imps <- nth_implausible(ems, points, targets, cutoff = cutoff)
  return(1-sum(nth_imps)/nrow(points))
}

#' Diagnostic Tests for Emulators
#'
#' Given an emulator, return a diagnostic measure.
#'
#' An emulator's suitability can be checked in a number of ways. This function combines all
#' current diagnostics available in the package, returning a context-dependent data.frame
#' containing the results.
#'
#' Comparison Diagnostics (cd): Given a set of points, the emulator expectation and variance
#' are calculated. This gives a predictive range for the input point, according to the
#' emulator. We compare this against the actual value given by the simulator: points whose
#' emulator prediction is further away from the simulator prediction are to be investigated.
#' This 'distance' is given by \code{stdev}, and an emulator prediction correspondingly
#' should not be further away from the simulator value than stdev*uncertainty.
#'
#' Classification Error (ce): Given a set of targets, the emulator can determine implausibility
#' of a point with respect to the relevant target, accepting or rejecting it as appropriate.
#' We can define a similar `implausibility' function for the simulator: the combination of
#' the two rejection schemes gives four classifications of points. Any point where the
#' emulator would reject the point but the simulator would not should be investigated.
#'
#' Standardized Error (se): The known value at a point, combined with the emulator expectation
#' and uncertainty, can be combined to provide a standardized error for a point. This error
#' should not be too large, in general. but the diagnostic is more useful when looking at
#' a collection of such measures, where systematic bias or over/underconfidence can be seen.
#'
#' Which of the diagnostics is performed can be controlled by the \code{which_diag} argument.
#' If performing classification error diagnostics, a set of targets must be provided; for all
#' diagnostics, a validation (or holdout) set can be provided. If no such set is given, then
#' the emulator diagnostics are performed with respect to its training points, using k-fold
#' cross-validation.
#'
#' @param emulator An object of class Emulator
#' @param targets If desired, the target values for the output(s) of the system
#' @param validation If provided, the emulator is tested against the outputs of these points
#' @param which_diag Which diagnostic measure to use (choosing from cd, ce, se above)
#' @param stdev For `cd', a measure of the allowed distance from prediction and reality
#' @param cleaned Internal for stochastic emulators
#' @param warn Should a warning be shown if ce is chosen and no targets provided?
#' @param kfold Mainly internal: pre-computed k-fold diagnostic results for output
#' @param ... Any other parameters to be passed through to subfunctions.
#'
#' @return A data.frame consisting of the input points, output values, and diagnostic measures.
#' @export
#'
#' @family diagnostic functions
#' @seealso validation_diagnostics
#'
#' @examples
#'  # Use the simple SIR model via SIREmulators
#'  get_diagnostic(SIREmulators$ems$nS, validation = SIRSample$validation)
#'  # Classification error fails without the set of targets
#'  get_diagnostic(SIREmulators$ems$nI, SIREmulators$targets, SIRSample$validation, 'ce')
#'  # No validation set: k-fold cross-validation will be used.
#'  get_diagnostic(SIREmulators$ems$nR, which_diag = 'se')
#'
get_diagnostic <- function(emulator, targets = NULL, validation = NULL,
                           which_diag = 'cd', stdev = 3, cleaned = NULL,
                           warn = TRUE, kfold = NULL, ...) {
  if (is.null(validation) && "EmProto" %in% class(emulator))
    stop("Proto_emulator object requires validation set.")
  if (is.null(targets) && which_diag == 'ce')
    stop("Targets must be provided for classification error diagnostics.")
  if (is.null(validation)) {
    if (is.null(kfold)) {
      if (which_diag == "ce"){
        dat <- k_fold_measure(emulator, targets[[emulator$output_name]], ...)
      }
      else dat <- k_fold_measure(emulator, ...)
    }
    else dat <- kfold
    input_points <- dat[,names(emulator$ranges)]
    output_points <- dat[,emulator$output_name]
    if (which_diag == 'ce')
      em_imp <- dat[,"I"]
    else {
      em_exps <- dat[,"E"]
      em_vars <- dat[,"V"]
      if ("Hierarchical" %in% class(emulator))
        point_vars <- dat[,"Vest"]
    }
  }
  else {
    if ("Hierarchical" %in% class(emulator)) {
      if (!is.null(cleaned)) cleaned_data <- cleaned
      else cleaned_data <- clean_data(validation, names(emulator$ranges),
                                      emulator$output_name,
                                      emulator$em_type == "variance")
      input_points <- cleaned_data[,names(emulator$ranges)]
      output_points <- cleaned_data$mean
      if (which_diag == "ce")
        em_imp <- emulator$implausibility(input_points,
                                          targets[[emulator$output_name]])
      else {
        em_exps <- emulator$get_exp(input_points)
        em_vars <- emulator$get_cov(input_points)
        point_vars <- cleaned_data$var/(cleaned_data$n-1)
      }
    }
    else {
      input_points <- validation[,names(emulator$ranges)]
      output_points <- validation[,emulator$output_name]
      if (which_diag == "ce")
        em_imp <- emulator$implausibility(input_points,
                                          targets[[emulator$output_name]])
      else {
        em_exps <- emulator$get_exp(input_points)
        em_vars <- emulator$get_cov(input_points)
      }
    }
  }
  if (which_diag == 'se') {
    num <- em_exps - output_points
    if ("Hierarchical" %in% class(emulator))
      denom <- sqrt(em_vars + point_vars +
                      emulator$disc$internal^2 + emulator$disc$external^2)
    else
      denom <- sqrt(em_vars +
                      emulator$disc$internal^2 + emulator$disc$external^2)
    errors <- num/denom
    out_data <- setNames(
      do.call(
        'cbind.data.frame',
        list(input_points, output_points, errors)),
      c(names(input_points), emulator$output_name, 'error'))
  }
  if (which_diag == 'cd') {
    if (is.null(stdev)) stdev <- 3
    if ("Hierarchical" %in% class(emulator))
      emulator_unc <- stdev * sqrt(em_vars + point_vars +
                                     emulator$disc$internal^2 +
                                     emulator$disc$external^2)
    else
      emulator_unc <- stdev * sqrt(em_vars +
                                     emulator$disc$internal^2 +
                                     emulator$disc$external^2)
    out_data <- setNames(
      do.call(
        'cbind.data.frame',
        list(input_points, output_points, em_exps, emulator_unc)),
      c(names(input_points), emulator$output_name, 'exp', 'unc'))
  }
  if (which_diag == 'ce') {
    this_target <- targets[[emulator$output_name]]
    if (is.atomic(this_target))
      sim_imp <- abs(output_points -
                       rep(mean(this_target), length(output_points)))/
        rep(diff(this_target)/2, length(output_points))
    else
      sim_imp <- purrr::map_dbl(output_points,
                                ~sqrt((this_target$val-.)^2/
                                        this_target$sigma^2))
    out_data <- setNames(
      do.call(
        'cbind.data.frame',
        list(input_points, output_points, em_imp, sim_imp)),
      c(names(input_points), emulator$output_name, 'em', 'sim'))
  }
  return(out_data)
}

#' Diagnostic Analysis for Emulators
#'
#' Produces summary and plots for diagnostics
#'
#' Given diagnostic information (almost certainly provided from \code{\link{get_diagnostic}}),
#' we can plot the results and highlight the points that are worthy of concern or further
#' consideration. Each diagnostic available has a plot associated with it which can be produced
#' here:
#'
#' Standardized Error: A histogram of standardized errors. Outliers should be considered, as well
#' as whether very many points have either large or small errors.
#'
#' Comparison Diagnostics: Error bars around points, corresponding to emulator prediction plus or
#' minus emulator uncertainty. A green line indicates where the emulator and simulator prediction
#' would be in complete agreement: error bars that do not overlap with this line (coloured red) are
#' to be considered. Where targets are provided, the colouration is limited only to points where
#' the simulator prediction would be close to the targets.
#'
#' Classification Error: A point plot comparing emulator implausibility to simulator
#' implausibility, sectioned into regions horizontally and vertically by \code{cutoff}. Points
#' that lie in the lower right quadrant (i.e. emulator would reject; simulator would not) should
#' be considered.
#'
#' This function takes a data.frame that contains the input points, simulator values and, depending
#' on the diagnostic, a set of summary measures. It returns a data.frame of any points that failed
#' the diagnostic.
#'
#' We may also superimpose the target bounds on the comparison diagnostics, to get a sense of
#' where it is most important that the emulator and simulator agree. The \code{target_viz}
#' argument controls this, and has three options: 'interval' (a horizontal interval); 'solid'
#' (a solid grey box whose dimensions match the target region in both vertical and horizontal
#' extent); and 'hatched' (similar to solid, but a semi-transparent box with hatching inside).
#' Any such vizualisation has extent equal to the target plus/mius 4.5 times the target
#' uncertainty. By default, \code{target_viz = NULL}, indicating that no superposition is shown.
#'
#' @importFrom graphics abline arrows hist rect
#' @importFrom grDevices rgb
#'
#' @param in_data The data to perform the analysis on
#' @param output_name The name of the output emulated
#' @param targets If required or desired, the targets for the system outputs
#' @param plt Whether or not to plot the analysis
#' @param cutoff The implausibility cutoff for diagnostic `ce'
#' @param target_viz How to show the targets on the diagnostic plots
#' @param ... Any other parameters to pass to subfunctions
#'
#' @return A data.frame of failed points
#'
#' @family diagnostic functions
#'
#' @references Jackson (2018) <http://etheses.dur.ac.uk/12826>
#' @export
#'
#' @seealso \code{\link{get_diagnostic}}
analyze_diagnostic <- function(in_data, output_name, targets = NULL,
                               plt = interactive(), cutoff = 3, target_viz = NULL, ...) {
  if (!is.null(target_viz))
    if (!target_viz %in% c("interval", "solid", "hatched")) target_viz <- NULL
  output_points <- in_data[,output_name]
  input_points <- in_data[, !names(in_data) %in% c(output_name,
                                                   'error', 'em',
                                                   'sim', 'exp', 'unc')]
  if (!is.null(in_data$error)) {
    if (plt) {
      h1 <- hist(in_data$error, plot = FALSE)
      plot(h1, xlab = "Standardised Error", main = output_name)
    }
    emulator_invalid <- abs(in_data$error) > 3
    if (!is.null(targets)) {
      this_target <- targets[[output_name]]
      if (is.atomic(this_target))
        point_invalid <- ((output_points < this_target[1]-diff(this_target)/2) |
                            (output_points > this_target[2]+diff(this_target)/2))
      else
        point_invalid <- ((output_points < this_target$val - 6*this_target$sigma) |
                            (output_points > this_target$val + 6*this_target$sigma))
      if (plt) {
        errors_restricted <- in_data[!point_invalid, 'error']
        h2 <- hist(errors_restricted, breaks = h1$breaks, plot = FALSE)
        plot(h2, add = TRUE, col = 'blue')
      }
      emulator_invalid <- (!point_invalid & emulator_invalid)
    }
  }
  if (!is.null(in_data$exp)) {
    em_ranges <- range(c(in_data$exp + in_data$unc, in_data$exp - in_data$unc))
    emulator_invalid <- (output_points > in_data$exp + in_data$unc) |
      (output_points < in_data$exp - in_data$unc)
    if (!is.null(targets)) {
      this_target <- targets[[output_name]]
      if (is.atomic(this_target)) {
        point_invalid <- ((output_points < this_target[1]-diff(this_target)/2) |
                            (output_points > this_target[2]+diff(this_target)/2))
        panlims <- c(this_target[1]-diff(this_target)/4,
                   this_target[2]+diff(this_target)/4)
      }
      else {
        point_invalid <- ((output_points < this_target$val - 6*this_target$sigma) |
                            (output_points > this_target$val + 6*this_target$sigma))
        panlims <- c(this_target$val - 4.5*this_target$sigma,
                     this_target$val + 4.5*this_target$sigma)
      }
      emulator_invalid <- (!point_invalid & emulator_invalid)
    }
    else panlims <- NULL
    if (plt) {
      plot(output_points, in_data$exp, pch = 16,
           col = ifelse(emulator_invalid, 'red', 'black'),
           xlim = range(output_points), ylim = range(em_ranges),
           xlab = 'f(x)', ylab = 'E[f(x)]',
           panel.first = c(
             if (!is.null(target_viz)) {
               if (target_viz == "interval")
                 arrows(x0 = panlims[1], x1 = panlims[2],
                        y0 = min(in_data$exp)+0.05*diff(range(in_data$exp)),
                        y1 = min(in_data$exp)+0.05*diff(range(in_data$exp)),
                        length = 0.05, code = 3, angle = 90)
               else if (target_viz == "solid")
                 rect(xleft = panlims[1], xright = panlims[2],
                      ybottom = panlims[1], ytop = panlims[2],
                      col = rgb(40, 40, 40, 51, maxColorValue = 255))
               else
                 rect(xleft = panlims[1], xright = panlims[2],
                      ybottom = panlims[1], ytop = panlims[2],
                      col = rgb(40, 40, 40, 102, maxColorValue = 255),
                      density = 15, angle = 135)
             }
             else NULL,
               abline(a = 0, b = 1, col = 'green')),
           main = output_name)
      for (i in seq_along(input_points[,1])) {
        if (in_data$unc[[i]] < 1e-6) next
        suppressWarnings(arrows(x0 = output_points[[i]],
               y0 = in_data$exp[[i]] - in_data$unc[[i]],
               x1 = output_points[[i]],
               y1 = in_data$exp[[i]] + in_data$unc[[i]],
               col = ifelse(emulator_invalid[[i]], 'red', 'blue'),
               code = 3, angle = 90, length = 0.05))
      }
    }
  }
  if (!is.null(in_data$em)) {
    if (is.null(targets)) stop("Require target to check classification error.")
    if (is.null(cutoff)) cutoff <- 3
    this_target <- targets[[output_name]]
    t_cutoff <- if (is.atomic(this_target)) 1 * cutoff/3 else cutoff
    emulator_invalid <- (in_data$em > cutoff) & (in_data$sim <= t_cutoff)
    if (plt) {
      plot(in_data$em, in_data$sim, pch = 16,
           col = ifelse(emulator_invalid, 'red', 'black'),
           xlab = "Emulator Implausibility", ylab = "Simulator Implausibility",
           main = output_name,
           panel.first = c(abline(h = t_cutoff), abline(v = cutoff)))
    }
  }
  return(input_points[emulator_invalid,])
}

#' Emulator Diagnostics
#'
#' Performs the standard set of validation diagnostics on emulators.
#'
#' All the diagnostics here can be performed with or without a validation (or `holdout')
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
#'   All of the above (all)
#'
#' For details of each of the tests, see the help file for
#' \code{\link{get_diagnostic}}.
#'
#' @importFrom graphics par
#'
#' @param emulators A list of \code{\link{Emulator}} objects.
#' @param targets The list of observations for the outputs
#' @param validation The validation set, containing all inputs and outputs.
#' @param which_diag Which diagnostics should be performed (see description)
#' @param analyze Should plotting and/or failing points be returned?
#' @param diagnose For bimodal systems, should the expectation or variance be considered?
#' @param ... Any additional parameters to pass to the diagnostics (eg sd, cutoff, ...)
#'
#' @return A data.frame containing points that failed one or more diagnostic tests.
#'
#' @family diagnostic functions
#' @export
#'
#' @examples
#' validation_diagnostics(SIREmulators$ems, SIREmulators$targets, SIRSample$validation)
#' # data.frame of failed points (empty) and a 3x3 set of plots
#' validation_diagnostics(SIREmulators$ems, SIREmulators$targets, SIRSample$validation,
#'  c('ce','cd'))
#' # empty data.frame and a 3x2 set of plots
#' validation_diagnostics(SIREmulators$ems, SIREmulators$targets, SIRSample$validation,
#'  cutoff = 2, sd = 2)
#' # k-fold (with k = 3)
#' validation_diagnostics(SIREmulators$ems, SIREmulators$targets, k = 3)
validation_diagnostics <- function(emulators, targets = NULL,
                                   validation = NULL,
                                   which_diag = c('cd', 'ce', 'se'),
                                   analyze = TRUE,
                                   diagnose = "expectation", ...) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  if ("Emulator" %in% class(emulators))
    emulators <- setNames(list(emulators), emulators$output_name)
  if ("EmProto" %in% unlist(purrr::map(emulators, class), use.names = FALSE) && is.null(validation))
    stop("Proto_emulator objects require a validation set.")
  if (length(which_diag) == 1 && which_diag == 'all') actual_diag <- c('cd', 'ce', 'se')
  else {
    actual_diag <- which_diag[which_diag %in% c('cd', 'ce', 'se')]
    if (length(actual_diag) != length(which_diag))
      warning(paste("Unrecognised diagnostics:",
                    paste0(which_diag[!which_diag %in% actual_diag],
                           collapse = ","),
                    "\n\tValid diagnostics labels are 'cd', 'se', 'ce', or 'all'."))
  }
  actual_diag <- unique(actual_diag)
  if (('ce' %in% actual_diag) && is.null(targets)) {
    warning("No targets provided; cannot perform classification diagnostics.")
    actual_diag <- actual_diag[-which(actual_diag == 'ce')]
  }
  if (!is.null(validation))
    validation <- validation[apply(validation, 1, function(x) all(!is.na(x))),]
  if (!is.null(emulators$mode1) && !is.null(emulators$mode2)) {
    if (diagnose == "variance") {
      m1_ems <- emulators$mode1$variance
      m2_ems <- emulators$mode2$variance
    }
    else {
      m1_ems <- emulators$mode1$expectation
      m2_ems <- emulators$mode2$expectation
    }
    if (!is.null(validation)) {
      v_outputs <- validation[,
                              unique(c(purrr::map_chr(m1_ems, ~.$output_name),
                                       purrr::map_chr(m2_ems, ~.$output_name))), drop = FALSE]
      v_class <- fanny(suppressWarnings(daisy(v_outputs)), k = 2)$clustering
      #v_class <- Mclust(v_outputs, G = 1:2, verbose = FALSE)$classification
      valid_one <- validation[v_class == 1,]
      valid_two <- validation[v_class == 2,]
      cleaning_dat1 <- purrr::map(m1_ems, function(x) {
        if ("Hierarchical" %in% class(x))
          valid_dat <- clean_data(valid_one,
                                  names(x$ranges),
                                  x$output_name, x$em_type == "variance")
        else valid_dat <- NULL
        valid_dat
      })
      cleaning_dat2 <- purrr::map(m1_ems, function(x) {
        if ("Hierarchical" %in% class(x))
          valid_dat <- clean_data(valid_two,
                                  names(x$ranges),
                                  x$output_name, x$em_type == "variance")
        else valid_dat <- NULL
        valid_dat
      })
      match_modes <- purrr::map_dbl(m1_ems, function(x) {
        outputs1 <- c(x$out_data, cleaning_dat1[[x$output_name]]$mean)
        outputs2 <- c(x$out_data, cleaning_dat2[[x$output_name]]$mean)
        if (sd(outputs1) < sd(outputs2)) return(1)
        return(2)
      })
      cluster1 <- purrr::map(
        seq_along(m1_ems),
        ~list(emulator = m1_ems[[.]],
              validation = list(valid_one, valid_two)[[match_modes[.]]]))
      cluster2 <- purrr::map(
        seq_along(m2_ems),
        ~list(emulator = m2_ems[[.]],
              validation = list(valid_two, valid_one)[[match_modes[.]]]))
      res_one <- setNames(unlist(purrr::map(cluster1, function(x) {
        suppressWarnings(
          validation_diagnostics(x$emulator, targets, x$validation,
                                 which_diag, analyze = FALSE, ...))
      }), recursive = FALSE), purrr::map_chr(cluster1, ~.$emulator$output_name))
      res_two <- setNames(unlist(purrr::map(cluster2, function(x) {
        suppressWarnings(
          validation_diagnostics(x$emulator, targets, x$validation,
                                 which_diag, analyze = FALSE, ...))
      }), recursive = FALSE), purrr::map_chr(cluster2, ~.$emulator$output_name))
    }
    else {
      res_one <- setNames(
        unlist(
          purrr::map(
            m1_ems,
            ~suppressWarnings(
              validation_diagnostics(., targets, validation,
                                     which_diag, analyze = FALSE, ...))),
          recursive = FALSE), purrr::map_chr(m1_ems, ~.$output_name))
      res_two <- setNames(
        unlist(
          purrr::map(
            m2_ems,
            ~suppressWarnings(
              validation_diagnostics(., targets, validation,
                                     which_diag, analyze = FALSE, ...))),
          recursive = FALSE), purrr::map_chr(m2_ems, ~.$output_name))
    }
    diag_results <- setNames(
      purrr::map(unique(c(names(res_one), names(res_two))), function(x) {
      if (!is.null(res_one[[x]]) && !is.null(res_two[[x]]))
        return(purrr::map(seq_along(res_one[[x]]), ~rbind(res_one[[x]][[.]],
                                                          res_two[[x]][[.]])))
      if (is.null(res_one[[x]])) return(res_two[[x]])
      return(res_one[[x]])
    }), unique(c(names(res_one), names(res_two))))
  }
  else {
    if (!is.null(emulators$expectation)) {
      emulators <- emulators[[ifelse(diagnose == "variance", "variance", "expectation")]]
    }
    if (!is.null(validation)) {
      cleaning_dat <- purrr::map(emulators, function(x) {
        if ("Hierarchical" %in% class(x))
          valid_dat <- clean_data(validation, names(x$ranges),
                                  x$output_name, x$em_type == "variance")
        else valid_dat <- NULL
        valid_dat
      })
      kfolded <- NULL
    }
    else {
      cleaning_dat <- NULL
      kfolded <- setNames(
        purrr::map(
          emulators,
          function(x) k_fold_measure(x, targets[[x$output_name]], ...)),
        purrr::map_chr(emulators, ~.$output_name))
    }
    diag_results <- purrr::map(
      emulators,
      function(x) purrr::map(actual_diag, function(y)
        get_diagnostic(x, targets, validation, y,
                       cleaned = cleaning_dat[[x$output_name]],
                       kfold = kfolded[[x$output_name]], ...)))
  }
  if (!analyze) return(diag_results)
  fail_point_list <- list()
  mf <- length(actual_diag)
  n_row <- list(...)[['row']]
  if (!is.null(n_row)) op <- par(mfrow = c(n_row, mf))
  else op <- par(mfrow = c(3, mf))
  for (i in seq_along(diag_results)) {
    for (j in seq_along(diag_results[[i]])) {
      fail_point_list[[length(fail_point_list)+1]] <- analyze_diagnostic(
        diag_results[[i]][[j]], names(diag_results)[[i]], targets, ...)
    }
  }
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
#' @references Bastos & O'Hagan (2009) <doi:10.1198/TECH.2009.08019>
#'
#' @family diagnostic functions
#' @export
#'
#' @examples
#' i1 <- individual_errors(SIREmulators$ems$nS, SIRSample$validation)
#' i2 <- individual_errors(SIREmulators$ems$nS, SIRSample$validation, "chol", "em")
#' i3 <- individual_errors(SIREmulators$ems$nS, SIRSample$validation, "eigen", plottype = "qq")
#' i4 <- individual_errors(SIREmulators$ems$nS, SIRSample$validation, "cholpivot", xtype = "aSI")
#'
individual_errors <- function(em, validation, errtype = "normal",
                              xtype = "index", plottype = "normal") {
  if ("EmProto" %in% class(em))
    stop("individual_errors not applicable for Proto_emulator objects.")
  if (!errtype %in% c("normal", "eigen", "chol", "cholpivot"))
    stop(paste("Error type not recognised",
    "(options are normal, eigen, chol, or cholpivot)."))
  if (!xtype %in% c("index", "em", names(em$ranges)))
    stop(paste("x-axis measure not recognised",
    "(options are index, em, or a parameter name."))
  if (!plottype %in% c("normal", "qq"))
    stop("Plot type not recognised (options are normal or qq).")
  if (plottype == "qq" && errtype == "normal") {
    warning(
      paste("Not meaningful to create Q-Q plot with untransformed errors.",
            "Changing to pivoted Cholesky."))
    errtype <- "cholpivot"
  }
  if (xtype %in% c(names(em$ranges), 'em') && errtype == "eigen") {
    warning(paste(
      "Not meaningful to plot parameter or emulator prediction",
      "against eigendecomposed errors. Changing to pivoted Cholesky."))
    errtype <- "cholpivot"
  }
  points <- validation[,names(em$ranges)]
  outputs <- validation[,em$output_name]
  em_pred <- em$get_exp(points)
  em_cov <- em$get_cov(points, full = TRUE)
  ## Ensure positive-definiteness
  em_cov_struct <- eigen(em_cov)
  em_cov_eval <- em_cov_struct$values
  em_cov_eval[em_cov_eval < 0] <- 1e-6
  em_cov <- em_cov_struct$vectors %*% diag(em_cov_eval, nrow = length(em_cov_eval)) %*% t(em_cov_struct$vectors)
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
  if (xtype == "index") x_vals <- seq_along(outputs)
  else if (xtype == "em") x_vals <- em_pred
  else x_vals <- validation[,xtype]
  if (plottype == "normal") {
    x_lab <- switch(xtype, "index" = "Index",
                    "em" = "Emulator Prediction", xtype)
    appended <- switch(errtype, "normal" = "",
                       "eigen" = " (eigendecomposition)",
                       "chol" = " (Choleksy decomposition)",
                       "cholpivot" = " (pivoted Cholesky decomposition)")
    plot(x_vals, indiv_errors, pch = 16, xlab = x_lab, ylab = "Error",
         main = paste0("Errors against ", x_lab, appended),
         panel.first = abline(h = c(-2,2), lty = 2))
  }
  if (plottype == "qq") {
    appended <- switch(errtype, "eigen" = " (eigendecomposition)",
                       "chol" = " (Choleksy decomposition)",
                       "cholpivot" = " (pivoted Cholesky decomposition)")
    qqnorm(indiv_errors, pch = 16,
           main = paste0("Q-Q plot for Errors", appended),
           panel.first = qqline(indiv_errors, col = "steelblue"))
  }
  output <- data.frame(variable = x_vals, error = indiv_errors)
  return(output)
}

#' Classification Diagnostics
#'
#' Shorthand function for diagnostic test `ce'.
#'
#' For details of the function, see \code{\link{get_diagnostic}} and for the plot
#' see \code{\link{analyze_diagnostic}}.
#'
#' @param emulator The emulator in question
#' @param targets The output targets
#' @param validation The validation set
#' @param cutoff The implausibility cutoff
#' @param plt Whether to plot or not
#'
#' @return A data.frame of failed points
#'
#' @family diagnostic functions
#' @export
#'
#' @references Jackson (2018) <http://etheses.dur.ac.uk/12826>
#'
#' @seealso \code{\link{get_diagnostic}}, \code{\link{analyze_diagnostic}},
#'  \code{\link{validation_diagnostics}}
#'
classification_diag <- function(emulator, targets, validation, cutoff = 3,
                                plt = interactive()) {
  analyze_diagnostic(get_diagnostic(emulator, targets, validation, 'ce'),
                     emulator$output_name, targets, plt, cutoff)
}

#' Comparison Diagnostics
#'
#' Shorthand function for diagnostic test `cd'.
#'
#' For details of the function, see \code{\link{get_diagnostic}} and for the plot
#' see \code{\link{analyze_diagnostic}}.
#'
#' @param emulator The emulator in question
#' @param targets The output targets
#' @param validation The validation set
#' @param sd The range of uncertainty allowed
#' @param plt Whether to plot or not
#'
#' @return A data.frame of failed points
#'
#' @family diagnostic functions
#' @export
#'
#' @references Jackson (2018) <http://etheses.dur.ac.uk/12826>
#'
#' @seealso \code{\link{get_diagnostic}}, \code{\link{analyze_diagnostic}},
#'  \code{\link{validation_diagnostics}}
#'
comparison_diag <- function(emulator, targets, validation, sd = 3,
                            plt = interactive()) {
  analyze_diagnostic(get_diagnostic(emulator, targets, validation, 'cd', sd),
                     emulator$output_name, targets, plt)
}

#' Standardized Error Diagnostics
#'
#' Shorthand function for diagnostic test `se'.
#'
#' For details of the function, see \code{\link{get_diagnostic}} and for the plot
#' see \code{\link{analyze_diagnostic}}.
#'
#' @param emulator The emulator in question
#' @param targets The output targets
#' @param validation The validation set
#' @param plt Whether to plot or not
#'
#' @return A data.frame of failed points
#'
#' @family diagnostic functions
#' @export
#'
#' @references Jackson (2018) <http://etheses.dur.ac.uk/12826>
#'
#' @seealso \code{\link{get_diagnostic}}, \code{\link{analyze_diagnostic}},
#'  \code{\link{validation_diagnostics}}
#'
standard_errors <- function(emulator, targets = NULL, validation = NULL,
                            plt = interactive()) {
  analyze_diagnostic(get_diagnostic(emulator, targets, validation, 'se'),
                     emulator$output_name, targets, plt)
}
