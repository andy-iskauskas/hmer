#' Output Plotting
#'
#' A simple diagnostic plot to determine relationships between emulated outputs and
#' each individual input dimension. If \code{points} is not provided, then a sample
#' of points across the full space is created by uniform sampling.
#'
#' @importFrom graphics plot par
#' @importFrom stats setNames
#'
#' @param ems A set of \code{\link{Emulator}} objects.
#' @param points A set of points at which to evaluate the emulator expectation
#'
#' @return The dependency plots.
#' @export
#'
#' @examples
#'  behaviour_plot(sample_emulators$ems, GillespieSIR)
#'  behaviour_plot(sample_emulators$ems$nS)
behaviour_plot <- function(ems, points) {
  if ("Emulator" %in% class(ems)) ems <- setNames(list(ems), ems$output_name)
  if (missing(points)) {
    ranges <- ems[[1]]$ranges
    points <- setNames(data.frame(do.call('cbind', purrr::map(ranges, ~runif(1000, .[1], .[2])))), names(ranges))
  }
  in_names <- names(ems[[1]]$ranges)
  out_names <- purrr::map_chr(ems, ~.$output_name)
  data <- setNames(data.frame(cbind(points, purrr::map(ems, ~.$get_exp(points)))), c(in_names, out_names))
  op <- par(mfrow = c(length(out_names), length(in_names)))
  for (i in 1:length(ems)) {
    for (j in 1:length(in_names)) {
      plot(data[,in_names[j]], data[,out_names[i]], pch = 16, xlab = in_names[[j]], ylab = out_names[[i]])
    }
  }
  par(op)
}

#' Space Removal Diagnostics
#'
#' Finds the proportion of space removed as a function of implausibility cut-off and of one of
#' structural discrepancy, emulator variance, or correlation hyperparameter(s).
#'
#' The reduction in space is found by evaluating a p^d regular grid, where p is chosen by
#' \code{ppd} and d is the dimension of the input psace. Larger values of p will give a more
#' accurate reflection of the space removed, at a corresponding computational cost. For the
#' purpose of quick-and-dirty diagnostics, \code{ppd = 5} is sufficient: the default is 10.
#'
#' The parameter \code{modified} can be one of three strings: \code{'disc'} corresponding
#' to model uncertainty (i.e. structural discrepancy); \code{'var'} corresponding to global
#' emulator variance (as given by \code{Emulator$u_sigma}), and \code{'hp'} corresponding to
#' the hyperparameters of the emulator correlation structure. In the first case, the
#' implausibilities are recalculated for each inflation value; in the other two cases the
#' emulators are retrained. For this reason, the \code{'var'} and \code{'hp'} options are
#' computationally more intensive. The default is \code{'disc'}.
#'
#' The inflationary/deflationary values are chosen by \code{u_mod}: the default is to take
#' 80%, 90%, 100%, 110%, and 120% of the original value as the variation. The proportion of
#' points deemed non-implausible is checked at a set of implausibility cutoffs defined by
#' \code{intervals}, and a plot is returned showing the relevant data.
#'
#' @import ggplot2
#' @importFrom stats setNames
#' @importFrom viridis scale_colour_viridis
#'
#' @param ems The \code{\link{Emulator}} objects.
#' @param targets The corresponding targets to match to.
#' @param ppd The number of points per input dimension to sample at.
#' @param u_mod The proportional values by which to inflate/deflate the relevant statistic.
#' @param intervals The interval values of the implausibility cutoff at which to evaluate.
#' @param modified The statistic to modify: disc, var or hp (see above)
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' space_removed(sample_emulators$ems, sample_emulators$targets, ppd = 5)
#' space_removed(sample_emulators$ems$nS, sample_emulators$targets$nS,
#'  ppd = 5, u_mod = seq(0.75, 1.25, by = 0.25), intervals = seq(2, 6, by = 0.1))
space_removed <- function(ems, targets, ppd = 10, u_mod = seq(0.8, 1.2, by = 0.1), intervals = seq(0, 10, length.out = 200), modified = 'disc') {
  value <- variable <- NULL
  if ("Emulator" %in% class(ems)) {
    ems <- setNames(list(ems), ems$output_name)
    if (!is.null(targets$val)) targets <- setNames(list(targets), ems[[1]]$output_name)
    else targets <- setNames(list(targets[[ems[[1]]$output_name]]), ems[[1]]$output_name)
  }
  ranges <- ems[[1]]$ranges
  z_vals <- purrr::map_dbl(targets, ~.$val)
  z_sigs <- purrr::map_dbl(targets, ~.$sigma)
  ptgrid <- setNames(expand.grid(purrr::map(ranges, ~seq(.[[1]], .[[2]], length.out = ppd))), names(ranges))
  imp_array <- array(0, dim = c(length(intervals), length(u_mod)))
  if (!modified %in% c('disc', 'var', 'hp')) {
    warning("Unrecognised vary parameter. Setting to structural discrepancy (disc).")
    modified = 'disc'
  }
  if (modified == 'disc') {
    exps <- do.call('cbind', purrr::map(ems, ~.$get_exp(ptgrid)))
    vars <- do.call('cbind', purrr::map(ems, ~.$get_cov(ptgrid)))
    for (i in u_mod) {
      imps <- abs(sweep(exps, 2, z_vals, "-"))/sqrt(sweep(vars, 2, (i*z_sigs)^2, "+"))
      m_imps <- apply(imps, 1, max)
      cutoff <- purrr::map_dbl(intervals, ~1-length(m_imps[m_imps <= .])/length(m_imps))
      imp_array[, match(i, u_mod)] <- cutoff
    }
  }
  else {
    for (i in u_mod) {
      if (modified == 'var')
        ems <- purrr::map(ems, ~.$set_sigma(i*.$u_sigma))
      else
        ems <- purrr::map(ems, ~.$set_hyperparams(purrr::map(.$corr$hyper_p, ~i*.)))
      imps <- nth_implausible(ems, ptgrid, targets)
      imp_array[, match(i, u_mod)] <- purrr::map_dbl(intervals, ~1-length(imps[imps <= .])/length(imps))
    }
  }
  df <- setNames(data.frame(imp_array), u_mod)
  df$cutoff <- intervals
  df_melt <- reshape2::melt(df, id.vars = 'cutoff')
  title <- switch(modified, 'disc' = 'structural discrepancy', 'var' = 'variance inflation', 'hp' = 'hyperparameter inflation')
  subtitle <- switch(modified, 'disc' = '% Structural\nDiscrepancy', 'var' = '% Variance\nInflation', 'hp' = '% Hyperparameter\nInflation')
  g <- ggplot(data = df_melt, aes(x = cutoff, y = value, group = variable, colour = variable)) +
    geom_line(lwd = 1.5) +
    viridis::scale_color_viridis(discrete = TRUE, option = 'cividis', labels = function(b) {paste0(round(as.numeric(b)*100, 0), "%")}) +
    scale_x_continuous("Implausibility cut-off", labels = function(b) round(b, 1)) +
    scale_y_continuous("Removed", labels = function(b) paste0(round(b*100, 0), "%")) +
    labs (title = paste("Space removed as a function of implausibilty cut-off and", title), colour = subtitle, x = "Cut-off", y = "% Removed") +
    theme_minimal()
  return(g)
}

#' Validation Set Diagnostics and Implausibility
#'
#' Creates pairs plots on the set of validation points of diagnostic suitability and
#' implausibility.
#'
#' The plots are organised as follows:
#'
#' a) Emulated versus simulated output (lower diagonal). This is similar in spirit to
#' \code{\link{comparison_diag}}: the plotted points are their location in the input
#' space and the points are coloured by the emulator prediction's deviation from the
#' simulator value.
#'
#' b) Implausibility (upper diagonal). The points are again plotted based on their
#' location in input space, but their colouration is now based on the implausibility
#' of the point.
#'
#' If \code{ranges} is provided, then the plotting region is created relative to these
#' ranges. This can be useful if on later waves of a history match and the plotting is
#' to be done relative to the original input space, rather than the (reduced) parameter
#' space upon which the emulators have been trained.
#'
#' @importFrom GGally ggpairs wrap
#' @importFrom rlang quo_get_expr
#' @import ggplot2
#'
#' @param ems The \code{\link{Emulator}} object(s).
#' @param points The set of validation points to plot.
#' @param targets The set of targets to match to.
#' @param ranges If provided, this gives the plotting region (see above).
#' @param nth The level of maximum implausibility to plot.
#' @param cb Whether or not the colour scheme should be colourblind friendly.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#'  validation_pairs(sample_emulators$ems, GillespieValidation, sample_emulators$targets)
#'  wider_ranges <- purrr::map(sample_emulators$ems[[1]]$ranges, ~.*c(-2, 2))
#'  validation_pairs(sample_emulators$ems, GillespieValidation,
#'   sample_emulators$targets, ranges = wider_ranges, cb = TRUE)
validation_pairs <- function(ems, points, targets, ranges, nth = 1, cb = FALSE) {
  if ("Emulator" %in% class(ems)) {
    ems <- setNames(list(ems), ems$output_name)
    if (!is.null(targets$val)) targets <- setNames(list(targets), ems[[1]]$output_name)
    else targets <- setNames(list(targets[[ems[[1]]$output_name]]), ems[[1]]$output_name)
  }
  if (missing(ranges)) ranges <- ems[[1]]$ranges
  em_exp <- data.frame(purrr::map(ems, ~.$get_exp(points)))
  em_var <- data.frame(purrr::map(ems, ~.$get_cov(points)))
  sim_vals <- points[,purrr::map_chr(ems, ~.$output_name)]
  diag_vals <- setNames(cbind(points[,names(ranges)], apply(abs(em_exp - sim_vals)/sqrt(em_var), 1, max)), c(names(ranges), 'd'))
  diag_vals$imp <- nth_implausible(ems, points, targets, n = nth)
  colour_breaks <- c(0, 0.3, 0.7, 1, 1.3, 1.7, 2, 2.3, 2.7, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 15, 100)
  colour_names <- c(0, '', '', 1, '', '', 2, '', '', 3, '', '', '', 5, '', '', '', 10, 15, '')
  plotfun <- function(data, mapping) {
    ggplot(data = data, mapping = mapping) +
      geom_point(cex = 2) +
      xlim(ranges[[rlang::quo_get_expr(mapping$x)]]) +
      ylim(ranges[[rlang::quo_get_expr(mapping$y)]])
  }
  cols <- if(cb) colourblindcont else redgreencont
  g <- ggpairs(diag_vals, columns = 1:length(ranges), aes(colour = diag_vals['d']), legend = c(1,2),
               title = "Emulator Diagnostics (lower) and Emulator Implausibility (upper)",
               lower = list(continuous = wrap(plotfun), mapping = aes(colour = diag_vals[,'d'])),
               upper = list(continuous = wrap(plotfun), mapping = aes(colour = diag_vals[,'imp'])),
               diag = 'blank', progress = FALSE) +
    scale_colour_gradient2(low = cols$low, mid = cols$mid, high = cols$high, midpoint = 3, breaks = colour_breaks, name = "Scale", labels = colour_names) +
    theme(legend.position = 'right') +
    theme_minimal()
  return(g)
}
