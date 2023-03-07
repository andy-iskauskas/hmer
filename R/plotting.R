# Plotting emulator expectation
exp_plot <- function(em, plotgrid = NULL, ppd = 30) {
  ranges <- em$ranges
  if (is.null(plotgrid)) {
    plotgrid <- setNames(
      expand.grid(
        seq(ranges[[1]][1], ranges[[1]][2], length.out = ppd),
        seq(ranges[[2]][1], ranges[[2]][2], length.out = ppd)),
      names(ranges)[1:2])
    for (i in 3:length(ranges)) {
      plotgrid[[names(ranges)[i]]] <- sum(ranges[[i]])/2
    }
  }
  em_exp <- em$get_exp(plotgrid[,names(ranges)])
  grid_data <- setNames(cbind(plotgrid[,1:2], em_exp),
                        c(names(plotgrid)[1:2],"E"))
  g <- tryCatch(
    {
      ## To account for rounding problems in ggplot's pretty_isoband_levels
      rng <- range(grid_data$E)
      if (diff(rng) == 0) warning("All output values are identical.")
      ac <- signif(diff(rng), 1)/10
      rng[1] <- floor(rng[1]/ac)*ac
      rng[2] <- ceiling(rng[2]/ac)*ac
      bks <- seq(rng[1], rng[2], length.out = 26)
      x_pos <- as.integer(factor(grid_data[,1],
                                 levels = sort(unique(grid_data[,1]))))
      y_pos <- as.integer(factor(grid_data[,2],
                                 levels = sort(unique(grid_data[,2]))))
      raster <- matrix(NA_real_, nrow = max(y_pos), ncol = max(x_pos))
      raster[cbind(y_pos, x_pos)] <- grid_data$E
      ibs <- isoband::isobands(x = sort(unique(grid_data[,1])),
                               y = sort(unique(grid_data[,2])),
                               z = raster, levels_low = bks[-length(bks)],
                               levels_high = bks[-1])
      int_lo <- gsub(":.*$", "", names(ibs))
      int_hi <- gsub("^[^:]*:", "", names(ibs))
      lab_lo <- format(as.numeric(int_lo), digits = 3, trim = TRUE)
      lab_hi <- format(as.numeric(int_hi), digits = 3, trim = TRUE)
      lab_check <- sprintf("(%s, %s]", lab_lo, lab_hi)
      if (length(unique(lab_check)) == 1)
        warning("Can't produce contours natively due to internal accuracy issues.")
      bns <- min(length(unique(lab_check)), 25)
      ggplot(data = grid_data, aes(x = grid_data[,1],
                                   y = grid_data[,2])) +
        geom_contour_filled(aes(z = grid_data[,'E']),
                            bins = bns, colour = 'black') +
        viridis::scale_fill_viridis(discrete = TRUE, option = "magma",
                                    name = "exp",
                                    guide = guide_legend(ncol = 1))
    },
    warning = function(w) { #nocov start
      exp_breaks <- seq(min(em_exp) - diff(range(em_exp))/(2 * 23),
                        max(em_exp) + diff(range(em_exp))/(2*23),
                        length.out = 25)
      exp_breaks <- unique(signif(exp_breaks, 10))
      if (length(exp_breaks) == 1) exp_breaks <- c(min(em_exp)-1e-6,
                                                   max(em_exp)+1e-6)
      intervals <- findInterval(em_exp, exp_breaks)
      fake_breaks <- seq_along(exp_breaks)
      ggplot(data = grid_data, aes(x = grid_data[,1],
                                   y = grid_data[,2])) +
        geom_contour_filled(aes(z = intervals), breaks = fake_breaks,
                            colour = 'black') +
        viridis::scale_fill_viridis(discrete = TRUE, option = "magma",
                                    name = "exp",
                                    guide = guide_legend(ncol = 1),
                                    labels = function(b)
                                      {signif(
                                        exp_breaks[as.numeric(
                                          stringr::str_extract(b, "\\d+"))],
                                        6)})
    } #nocov end
  )
  if (is.null(em$em_type)) extra_ident <- NULL
  else if (em$em_type == "mean") extra_ident <- "Mean"
  else if (em$em_type == "variance") extra_ident <- "Variance"
  else extra_ident <- em$em_type
  g <- g +
    labs(title = paste(em$output_name, extra_ident, "Emulator Expectation"),
         x = names(grid_data)[1], y = names(grid_data)[2]) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_minimal()
  return(g)
}

# Plotting emulator variance
var_plot <- function(em, plotgrid = NULL, ppd = 30, sd = FALSE) {
  ranges <- em$ranges
  if (is.null(plotgrid)) {
    plotgrid <- setNames(
      expand.grid(
        seq(
          ranges[[1]][1],
          ranges[[1]][2],
          length.out = ppd),
        seq(
          ranges[[2]][1],
          ranges[[2]][2],
          length.out = ppd)),
      names(ranges)[1:2])
    for (i in 3:length(ranges)) {
      plotgrid[[names(ranges)[i]]] <- sum(ranges[[i]])/2
    }
  }
  if(sd)
    em_cov <- sqrt(em$get_cov(plotgrid[,names(ranges)]))
  else
    em_cov <- em$get_cov(plotgrid[,names(ranges)])
  grid_data <- setNames(cbind(plotgrid[,1:2], em_cov),
                        c(names(plotgrid)[1:2], "V"))
  g <- tryCatch(
    {
      ## To account for rounding problems in ggplot's pretty_isoband_levels
      rng <- range(grid_data$V)
      if (diff(rng) == 0) warning("All output values are identical.")
      ac <- signif(diff(rng), 1)/10
      rng[1] <- floor(rng[1]/ac)*ac
      rng[2] <- ceiling(rng[2]/ac)*ac
      bks <- seq(rng[1], rng[2], length.out = 26)
      x_pos <- as.integer(factor(grid_data[,1],
                                 levels = sort(unique(grid_data[,1]))))
      y_pos <- as.integer(factor(grid_data[,2],
                                 levels = sort(unique(grid_data[,2]))))
      raster <- matrix(NA_real_, nrow = max(y_pos), ncol = max(x_pos))
      raster[cbind(y_pos, x_pos)] <- grid_data$V
      ibs <- isoband::isobands(x = sort(unique(grid_data[,1])),
                               y = sort(unique(grid_data[,2])),
                               z = raster, levels_low = bks[-length(bks)],
                               levels_high = bks[-1])
      int_lo <- gsub(":.*$", "", names(ibs))
      int_hi <- gsub("^[^:]*:", "", names(ibs))
      lab_lo <- format(as.numeric(int_lo), digits = 3, trim = TRUE)
      lab_hi <- format(as.numeric(int_hi), digits = 3, trim = TRUE)
      lab_check <- sprintf("(%s, %s]", lab_lo, lab_hi)
      if (length(unique(lab_check)) == 1)
        warning("Can't produce contours natively due to internal accuracy issues.")
      bns <- min(length(unique(lab_check)), 25)
      ggplot(data = grid_data, aes(x = grid_data[,1], y = grid_data[,2])) +
        geom_contour_filled(aes(z = grid_data[,'V']), bins = bns,
                            colour = 'black') +
        viridis::scale_fill_viridis(discrete = TRUE, option = "plasma",
                                    name = if(sd) "sd" else "var",
                                    guide = guide_legend(ncol = 1))
    },
    warning = function(w) {
      cov_breaks <- seq(
        min(em_cov) - diff(range(em_cov))/(2 * 23),
        max(em_cov) + diff(range(em_cov))/(2*23),
        length.out = 25)
      cov_breaks <- unique(signif(cov_breaks, 10))
      if (length(cov_breaks) == 1) cov_breaks <- c(min(em_cov)-1e-6,
                                                   max(em_cov)+1e-6)
      intervals <- findInterval(em_cov, cov_breaks)
      fake_breaks <- seq_along(cov_breaks)
      ggplot(data = grid_data, aes(x = grid_data[,1], y = grid_data[,2])) +
        geom_contour_filled(aes(z = intervals), breaks = fake_breaks,
                            colour = 'black') +
        viridis::scale_fill_viridis(discrete = TRUE, option = "plasma",
                                    name = if (sd) "sd" else "var",
                                    guide = guide_legend(ncol = 1),
                                    labels = function(b)
                                      {signif(
                                        cov_breaks[as.numeric(
                                          stringr::str_extract(b, "\\d+"))],
                                        6)})
    }
  )
  if (is.null(em$em_type)) extra_ident <- NULL
  else if (em$em_type == "mean") extra_ident <- "Mean"
  else if (em$em_type == "variance") extra_ident <- "Variance"
  else extra_ident <- em$em_type
  g <- g +
    labs(title = paste(em$output_name, extra_ident, "Emulator",
                       (if(sd) "Standard Deviation" else "Variance")),
         x = names(grid_data)[1], y = names(grid_data)[2]) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_minimal()
  return(g)
}

# Plotting emulator implausibility
imp_plot <- function(em, z, plotgrid = NULL, ppd = 30, cb = FALSE, nth = NULL,
                     imp_breaks = NULL) {
  if (!is.null(nth)) {
    if (nth == 1) ns <- ""
    else if (nth == 2) ns <- "Second"
    else if (nth == 3) ns <- "Third"
    else ns <- paste0(nth, "th")
    if (!is.null(em$expectation)) ranges <- em$expectation[[1]]$ranges
    else if (!is.null(em$mode1)) ranges <- em$mode1$expectation[[1]]$ranges
    else ranges <- em[[1]]$ranges
  }
  else ranges <- em$ranges
  if (is.null(plotgrid)) {
    plotgrid <- setNames(
      expand.grid(
        seq(
          ranges[[1]][1],
          ranges[[1]][2],
          length.out = ppd),
        seq(ranges[[2]][1],
            ranges[[2]][2],
            length.out = ppd)),
      names(ranges)[1:2])
    for (i in 3:length(ranges)) {
      plotgrid[[names(ranges)[i]]] <- sum(ranges[[i]])/2
    }
  }
  if (is.null(imp_breaks)) {
    imp_breaks <- c(0, 0.3, 0.7, 1, 1.3, 1.7, 2, 2.3, 2.7,
                    3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 15, Inf)
    imp_names <- c(0, '', '', 1, '', '', 2, '', '', 3, '',
                  '', '', 5, '', '', '', 10, 15, '')
  } else {
    if (!is.numeric(imp_breaks)) stop("imp_breaks not a numeric vector.")
    if (any(is.na(imp_breaks))) stop("imp_breaks cannot contain missing values.")
    if (any(diff(imp_breaks) <= 0)) stop("imp_breaks must be in ascending order.")
    if (length(imp_breaks) != 20) stop("imp_breaks must be of length 20.")
    imp_names <- as.character(imp_breaks)
  }
  if (!is.null(nth)) {
    if (!"Emulator" %in% class(em)) {
      em_imp <- nth_implausible(em, plotgrid[names(ranges)],
                                z, n = nth, max_imp = 99)
    }
    else
      stop(paste("Not all required parameters",
                 "(emulator list, target list, nth)",
                 "passed for nth maximum implausibility."))
  }
  else{
    em_imp <- em$implausibility(plotgrid[names(ranges)], z)
  }
  included <- c(purrr::map_lgl(imp_breaks[-1], ~any(em_imp < .)), TRUE)
  grid_data <- setNames(
    cbind(plotgrid[,1:2], em_imp), c(names(plotgrid)[1:2],"I"))
  col_scale <- if(cb) colourblind else redgreen
  g <- ggplot(data = grid_data, aes(x = grid_data[,1], y = grid_data[,2]))
  if (length(unique(em_imp)) == 1) {
    col <- col_scale[which.min(abs(imp_breaks - c(unique(em_imp))))]
    g <- g + geom_raster(aes(fill = grid_data[,"I"])) +
      scale_fill_gradient(low = col, high = col, name = "I")
  }
  else
    g <- g + geom_contour_filled(aes(z = grid_data[,"I"]),
                                 colour = 'black',
                                 breaks = imp_breaks[included]) +
    scale_fill_manual(values = col_scale[included], name = "I",
                      labels = imp_names[included],
                      guide = guide_legend(reverse = TRUE))
  if (is.null(em$em_type)) extra_ident <- NULL
  else if (em$em_type == "mean") extra_ident <- "Mean"
  else if (em$em_type == "variance") extra_ident <- "Variance"
  else extra_ident <- em$em_type
  g <- g + labs(title = paste(
    em$output_name,
    extra_ident,
    (if (is.null(nth)) "Emulator Implausibility" else paste(ns, "Maximum Implausibility"))),
    x = names(grid_data)[1], y = names(grid_data)[2]) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_minimal()
  return(g)
}

#' Plot Emulator Outputs
#'
#' A function for plotting emulator expectations, variances, and implausibilities
#'
#' Given a single emulator, or a set of emulators, the emulator statistics can be plotted
#' across a two-dimensional slice of the parameter space. Which statistic is plotted is
#' determined by \code{plot_type}: options are `exp', `var', `sd', `imp', and `nimp', which
#' correspond to expectation, variance, standard deviation, implausibility, and nth-max
#' implausibility.
#'
#' By default, the slice varies in the first two parameters of the emulators, and all other
#' parameters are taken to be fixed at their mid-range values. This behaviour can be changed
#' with the \code{params} and \code{fixed_vals} parameters (see examples).
#'
#' If the statistic is `exp', `var' or `sd', then the minimal set of parameters to pass to this
#' function are \code{ems} (which can be a list of emulators or a single one) and \code{plot_type}.
#' If the statistic is `imp' or `nimp', then the \code{targets} must be supplied - it is not
#' necessary to specify the individual target for a single emulator plot. If the statistic is
#' `nimp', then the level of maximum implausibility can be chosen with the parameter \code{nth}.
#'
#' Implausibility plots are typically coloured from green (low implausibility) to red (high
#' implausibility): a colourblind-friendly option is available and can be turned on by setting
#' \code{cb = TRUE}.
#'
#' The granularity of the plot is controlled by the \code{ppd} parameter, determining the number
#' of points per dimension in the grid. For higher detail, at the expense of longer computing
#' time, increase this value. The default is 30.
#'
#' @import ggplot2
#' @importFrom viridis scale_fill_viridis
#' @importFrom GGally ggally_text ggally_blank ggmatrix
#'
#' @param ems An \code{\link{Emulator}} object, or a list thereof.
#' @param plot_type The statistic to plot (see description or examples).
#' @param ppd The number of points per plotting dimension
#' @param targets If required, the targets from which to calculate implausibility
#' @param cb A boolean representing whether a colourblind-friendly plot is produced.
#' @param params Which two input parameters should be plotted?
#' @param fixed_vals For fixed input parameters, the values they are held at.
#' @param nth If plotting nth maximum implausibility, which level maximum to plot.
#' @param imp_breaks If plotting nth maximum implausibility, defines the levels at
#'                   which to draw contours.
#'
#' @return A ggplot object, or collection thereof.
#'
#' @family visualisation tools
#' @export
#'
#' @examples
#'  # Reducing ppd to 10 for speed.
#'  emulator_plot(SIREmulators$ems, ppd = 10)
#'  emulator_plot(SIREmulators$ems$nS, ppd = 10)
#'  emulator_plot(SIREmulators$ems, plot_type = 'var', ppd = 10, params = c('aIR', 'aSR'))
#'  \donttest{
#'     emulator_plot(SIREmulators$ems, plot_type = 'imp', ppd = 10,
#'      targets = SIREmulators$targets,
#'      fixed_vals = list(aSR = 0.02))
#'     emulator_plot(SIREmulators$ems, plot_type = 'nimp', cb = TRUE,
#'      targets = SIREmulators$targets, nth = 2, ppd = 10)
#'  }
#'
emulator_plot <- function(ems, plot_type = 'exp', ppd = 30, targets = NULL,
                          cb = FALSE, params = NULL, fixed_vals = NULL,
                          nth = 1, imp_breaks = NULL) {
  if ("Emulator" %in% class(ems)){
    ranges <- ems$ranges
    single_em <- TRUE
  }
  else {
    if (!is.null(ems$expectation)) ranges <- ems$expectation[[1]]$ranges
    else if (!is.null(ems$mode1)) ranges <- ems$mode1$expectation[[1]]$ranges
    else ranges <- ems[[1]]$ranges
    single_em <- FALSE
  }
  if (!is.null(ems$expectation)) ems <- ems$expectation
  if (is.null(params) ||
      length(params) != 2 ||
      any(!params %in% names(ranges))) p_vals <- c(1,2)
  else p_vals <- which(names(ranges) %in% params)
  plotgrid <- setNames(
    expand.grid(
      seq(
        ranges[[p_vals[1]]][1],
        ranges[[p_vals[1]]][2],
        length.out = ppd),
      seq(
        ranges[[p_vals[2]]][1],
        ranges[[p_vals[2]]][2],
        length.out = ppd)),
    names(ranges)[p_vals])
  if (!is.null(fixed_vals) && all(names(fixed_vals) %in% names(ranges))) {
    for (i in seq_along(fixed_vals))
      plotgrid[[names(fixed_vals)[i]]] <- fixed_vals[[names(fixed_vals)[i]]]
    used_names <- c(names(ranges)[p_vals], names(fixed_vals))
  }
  else used_names <- names(ranges)[p_vals]
  used_names <- unique(used_names)
  if (length(used_names) < length(ranges)) {
    unused_nms <- names(ranges)[which(!names(ranges) %in% used_names)]
    for (i in seq_along(unused_nms))
      plotgrid[[unused_nms[i]]] <- sum(ranges[[unused_nms[i]]])/2
  }
  get_plot <- function(em) {
    if (plot_type == 'exp') return(exp_plot(em, plotgrid, ppd))
    if (plot_type == 'var') return(var_plot(em, plotgrid, ppd))
    if (plot_type == 'sd') return(var_plot(em, plotgrid, ppd, sd = TRUE))
    if (plot_type == 'imp') {
      if (is.null(targets))
        stop("Cannot plot implausibility without target value.")
    }
    if (!is.null(targets$val)) return(imp_plot(em, targets, plotgrid, ppd, cb, NULL, imp_breaks))
    else return(imp_plot(em, targets[[em$output_name]], plotgrid, ppd, cb, NULL, imp_breaks))
  }
  if (single_em) return(get_plot(ems))
  if (plot_type == 'nimp') return(imp_plot(ems, targets, plotgrid, ppd, cb, nth, imp_breaks))
  else {
    plotlist <- purrr::map(ems, get_plot)
    replacement_function <- function(plots, title = NULL) {
      titles <- purrr::map_chr(
        plots,
        ~sub("(.*) Emulator (Expectation|Variance|Implausibility)",
             "\\1", .$labels$title))
      create_name_plot <- function(name) {
        return(ggally_text(label = name, colour = 'black') + theme_void())
      }
      plot_cols <- ceiling(sqrt(length(plots)))
      plot_rows <- ceiling(length(plots)/plot_cols)
      n_empty <- plot_cols^2 - length(plots)
      plot_list <- list()
      for (i in 1:(2*plot_rows)) {
        for (j in 1:plot_cols) {
          if (i %% 2 == 1) {
            in_bounds <- plot_cols*(i-1)/2+j
            if (in_bounds <= length(titles)) {
              plot_list[[length(plot_list)+1]] <- create_name_plot(titles[in_bounds])
            }
            else plot_list[[length(plot_list)+1]] <- ggally_blank()
          }
          else {
            in_bounds <- plot_cols*(i-2)/2+j
            if (in_bounds <= length(plots)) {
              plot_list[[length(plot_list)+1]] <- plots[[in_bounds]]
            }
            else plot_list[[length(plot_list)+1]] <- ggally_blank()
          }
        }
      }
      ggmatrix(plot_list, ncol = plot_cols, nrow = 2*plot_rows,
               xlab = plots[[1]]$labels$x, ylab = plots[[1]]$labels$y,
               title = title, yProportions = rep(c(0.05, 1), plot_rows),
               progress = FALSE)
    }
    if (plot_type == "exp") plot_title <- "Emulator Expectations"
    else if (plot_type == "var") plot_title <- "Emulator Variances"
    else if (plot_type == "imp") plot_title <- "Emulator Implausibilities"
    else plot_title <- NULL
    return(replacement_function(plotlist, plot_title))
  }
}

#' Emulator Expectation Against Target Outputs
#'
#' Plots emulator expectation across the parameter space, with comparison to the corresponding
#' target values (with appropriate uncertainty).
#'
#' If a \code{points} data.frame is not provided, then points are sampled uniformly from the
#' input region. Otherwise, the provided points are used: for example, if a representative
#' sample of the current NROY space is available.
#'
#' @importFrom ggplot2 ggplot aes labs geom_line geom_point geom_errorbar
#'
#' @param ems The \code{\link{Emulator}} objects.
#' @param targets A named list of observations, given in the usual form.
#' @param points A list of points at which the emulators should be evaluated.
#' @param npoints If no points are provided, the number of input points to evaluate at.
#'
#' @return A ggplot object
#'
#' @family visualisation tools
#' @export
#'
#' @examples
#'  output_plot(SIREmulators$ems, SIREmulators$targets)
#'  output_plot(SIREmulators$ems, SIREmulators$targets, points = SIRSample$training)
output_plot <- function(ems, targets, points = NULL, npoints = 1000) {
  ranges <- if ("Emulator" %in% class(ems)) ems$ranges else ems[[1]]$ranges
  if (is.null(points)) {
    points <- data.frame(purrr::map(ranges, ~runif(npoints, .[1], .[2])))
  }
  em_exp <- setNames(
    data.frame(purrr::map(ems, ~.$get_exp(points))), names(targets))
  em_exp$run <- seq_len(nrow(points))
  em_exp <- pivot_longer(em_exp, cols = !'run')
  for (i in seq_along(targets))
  {
    if (!is.atomic(targets[[i]]))
      targets[[i]] <- c(targets[[i]]$val - 3 * targets[[i]]$sigma,
                        targets[[i]]$val + 3*targets[[i]]$sigma)
  }
  target_data <- data.frame(label = names(targets),
                            mn = purrr::map_dbl(targets, ~.[1]),
                            md = purrr::map_dbl(targets, mean),
                            mx = purrr::map_dbl(targets, ~.[2]))
  name <- value <- run <- mn <- md <- mx <- label <- NULL
  em_exp$name <- factor(em_exp$name, levels = names(targets))
  ggplot(data = em_exp, aes(x = name, y = value)) +
    geom_line(colour = 'purple', aes(group = run), linewidth = 1) +
    geom_point(data = target_data, aes(x = label, y = md), size = 2) +
    geom_errorbar(data = target_data, aes(x = label, y = md, ymin = mn,
                                          ymax = mx), width = .1, linewidth = 1.25) +
    labs(title = "Emulator Runs versus Observations")
}

#' Plot Lattice of Emulator Implausibilities
#'
#' Plots a set of projections of the full-dimensional input space.
#'
#' The plots are:
#'
#' One dimensional optical depth plots (diagonal);
#'
#' Two dimensional optical depth plots (lower triangle);
#'
#' Two dimensional minimum implausibility plots (upper triangle).
#'
#' The optical depth is calculated as follows. A set of points is constructed across the
#' full d-dimensional parameter space, and implausibility is calculated at each point.
#' The points are collected into groups based on their placement in a projection to a
#' one- or two-dimensional slice of the parameter space. For each group, the proportion
#' of non-implausible points is calculated, and this value in [0,1] is plotted. The
#' minimum implausibility plots are similar, but with minimum implausibility calculated
#' rather than proportion of non-implausible points.
#'
#' The \code{maxpoints} argument is used as a cutoff for if a regular ppd grid would
#' result in a very large number of points. If this is the case, then \code{maxpoints} points
#' are sampled uniformly from the region instead of regularly spacing them.
#'
#' @importFrom stats xtabs
#' @importFrom GGally ggmatrix grab_legend
#' @import ggplot2
#'
#' @param ems The \code{\link{Emulator}} objects in question.
#' @param targets The corresponding target values.
#' @param ppd The number of points to sample per dimension.
#' @param cb Whether or not a colourblind-friendly plot should be produced.
#' @param cutoff The cutoff value for non-implausible points.
#' @param maxpoints The limit on the number of points to be evaluated.
#' @param imp_breaks If plotting nth maximum implausibility, defines the levels at
#'                   which to draw contours.
#' @param contour Logical determining whether to plot implausibility contours or not.
#'
#' @return A ggplot object
#'
#' @family visualisation tools
#' @export
#'
#' @references Bower, Goldstein & Vernon (2010) <doi:10.1214/10-BA524>
#'
#' @examples
#' \donttest{
#'  plot_lattice(SIREmulators$ems, SIREmulators$targets, ppd = 10)
#'  plot_lattice(SIREmulators$ems$nS, SIREmulators$targets)
#' }
plot_lattice <- function(ems, targets, ppd = 20, cb = FALSE,
                         cutoff = 3, maxpoints = 5e4, imp_breaks = NULL, 
                         contour = TRUE) {
  ems <- collect_emulators(ems)
  ranges <- if ("Emulator" %in% class(ems)) ems$ranges else ems[[1]]$ranges
  if (ppd^length(ranges) > maxpoints) {
    point_grid <- setNames(
      data.frame(
        do.call('cbind',
                purrr::map(ranges, ~runif(maxpoints, .[[1]], .[[2]])))),
      names(ranges))
    nbins <- 19
  }
  else {
    dim_bounds <- purrr::map(ranges, ~seq(.[[1]], .[[2]], length.out = ppd+1))
    dim_unif <- purrr::map(dim_bounds,
                           ~purrr::map_dbl(1:(length(.)-1),
                                           function(i) mean(.[i:(i+1)])))
    point_grid <- expand.grid(dim_unif)
  }
  point_grid$I <- nth_implausible(ems, point_grid, targets)
  one_dim <- function(data, parameter) {
    param_seq <- seq(
      ranges[[parameter]][1],
      ranges[[parameter]][2], length.out = ppd + 1)
    collection <- purrr::map(1:ppd, function(x) {
      valid_points <- data[data[,parameter] >= param_seq[x] &
                             data[,parameter] <= param_seq[x+1],]
      how_many_valid <- if (nrow(valid_points) == 0)
        NA
      else
        sum(valid_points$I <= cutoff)/nrow(valid_points)
      return(c(param_seq[x], how_many_valid))
    })
    setNames(do.call('rbind.data.frame', collection), c(parameter, 'op'))
  }
  two_dim <- function(data, parameters, op = FALSE) {
    param_seqs <- purrr::map(ranges[parameters],
                             ~seq(.[[1]], .[[2]], length.out = ppd + 1))
    param_list <- list()
    for (i in 1:ppd) {
      for (j in 1:ppd) {
        valid_points <- data[data[,parameters[1]] >= param_seqs[[1]][i] &
                               data[,parameters[1]] <= param_seqs[[1]][i+1] &
                               data[,parameters[2]] >= param_seqs[[2]][j] &
                               data[,parameters[2]] <= param_seqs[[2]][j+1],]
        if (op) {
          rel_stat <- if (nrow(valid_points) == 0)
            NA
          else
            sum(valid_points$I <= cutoff)/nrow(valid_points)
        }
        else {
          rel_stat <- if (nrow(valid_points) == 0) NA else min(valid_points$I)
        }
        param_list[[length(param_list)+1]] <- c(param_seqs[[1]][i],
                                                param_seqs[[2]][j], rel_stat)
      }
    }
    return(setNames(do.call('rbind.data.frame', param_list),
                    c(parameters, if(op) "op" else "imp")))
  }
  cols <- if(cb) colourblind else redgreen
  if (is.null(imp_breaks)) {
    imp_breaks <- c(0, 0.3, 0.7, 1, 1.3, 1.7, 2, 2.3, 2.7,
                    3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 15, Inf)
    imp_names <- c(0, '', '', 1, '', '', 2, '', '', 3, '',
                   '', '', 5, '', '', '', 10, 15, '')
  } else {
    if (!is.numeric(imp_breaks)) stop("imp_breaks not a numeric vector.")
    if (any(is.na(imp_breaks))) stop("imp_breaks cannot contain missing values.")
    if (any(diff(imp_breaks) <= 0)) stop("imp_breaks must be in ascending order.")
    if (length(imp_breaks) != 20) stop("imp_breaks must be of length 20.")
    imp_names <- as.character(imp_breaks)
  }
  parameter_combinations <- expand.grid(
    names(ranges),
    names(ranges),
    stringsAsFactors = FALSE)
  plot_list <- purrr::map(seq_len(nrow(parameter_combinations)), function(x) {
    parameters <- unlist(parameter_combinations[x,], use.names = FALSE)
    if (parameters[1] == parameters[2]) {
      pt <- one_dim(point_grid, parameters[1])
      g <- suppressMessages(ggplot(mapping = aes(x = pt[,1], y = pt[,2])) +
        geom_smooth(colour = 'black', se = FALSE) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 1)))
    }
    else if (which(
      names(ranges) == parameters[1]) >
      which(names(ranges) == parameters[2])) {
      pt <- two_dim(point_grid, parameters)
      g <- ggplot(mapping = aes(x = pt[,1], y = pt[,2], z = pt[,3])) +
        geom_contour_filled(breaks = imp_breaks, colour = ifelse(contour, "black", NA)) +
        scale_fill_manual(values = cols, labels = imp_names,
                          name = "Min. I",
                          guide = guide_legend(reverse = TRUE), drop = FALSE) +
        scale_y_continuous(expand = c(0,0))
    }
    else {
      pt <- two_dim(point_grid, parameters, TRUE)
      g <- ggplot(mapping = aes(x = pt[,1], y = pt[,2], fill = pt[,3])) +
        geom_raster(interpolate = TRUE) +
        scale_fill_gradient(low = 'black', high = 'white',
                            breaks = seq(0, 1, by = 0.1),
                            name = "Op. Depth") +
        scale_y_continuous(expand = c(0,0))
    }
    return(g + scale_x_continuous(expand = c(0,0)) + theme_minimal())
  })
  x <- y <- z <- NULL
  pointless_data <- expand.grid(x = 1:10, y = 1:10)
  pointless_data$fill <- seq(0, 1, length.out = 100)
  pointless_data$z <- seq(0, 20, length.out = 100)
  pointless_plot <- ggplot(data = pointless_data, aes(x = x, y = y)) +
    geom_line(aes(colour = fill)) +
    geom_contour_filled(aes(z = z), breaks = imp_breaks) +
    scale_fill_manual(values = cols, labels = imp_names,
                      name = "Min. Imp", drop = FALSE) +
    scale_colour_gradient(low = 'black', high = 'white',
                          breaks = seq(0, 1, by = 0.1), name = "Op. Depth") +
    guides(fill = guide_legend(order = 1, reverse = TRUE))
  return(ggmatrix(plot_list, length(ranges), length(ranges),
                  title = "Minimum Implausibility and Optical Depth",
                  xAxisLabels = names(ranges), yAxisLabels = names(ranges),
                  showYAxisPlotLabels = FALSE,
                  legend = grab_legend(pointless_plot)))
}

#' Active variable plotting
#'
#' For a set of emulators, demonstrate which variables are active.
#'
#' Each emulator has a list of `active' variables; those which contribute in an appreciable way
#' to its regression surface. It can be instructive to examine the differences in active variables
#' for a give collection of emulators. The plot here produces an nxp grid for n emulators in p
#' inputs; a square is blacked out if that variable does not contribute to that output.
#'
#' Both the outputs and inputs can be restricted to collections of interest, if desired, with the
#' optional \code{output_names} and \code{input_names} parameters.
#'
#' @param ems The list of emulators to consider
#' @param output_names The names of the outputs to include in the plot, if not all
#' @param input_names The names of the inputs to include in the plot, if not all
#' @return A ggplot object corresponding to the plot
#'
#' @family visualisation tools
#' @export
#'
#' @examples
#'  plot_actives(SIREmulators$ems)
#'  # Remove the nR output and aIR input from the plot
#'  plot_actives(SIREmulators$ems, c('nS', 'nI'), c('aSI', 'aSR'))
#'  # Note that we can equally restrict the emulator list...
#'  plot_actives(SIREmulators$ems[c('nS', 'nI')], input_names = c('aSI', 'aSR'))
plot_actives <- function(ems, output_names = NULL, input_names = NULL) {
  if ("Emulator" %in% class(ems)) {
    ems <- list(ems)
  }
  in_names <- names(ems[[1]]$ranges)
  active_list <- setNames(
    data.frame(
      do.call('rbind', purrr::map(ems, ~.$active_vars))), in_names)
  if (!is.null(input_names)) active_list <- active_list[, input_names, drop = FALSE]
  if (!is.null(output_names))
    active_list <- active_list[row.names(active_list) %in% output_names, , drop = FALSE]
  if (nrow(active_list) == 0 || length(active_list) == 0)
    stop("No inputs/outputs to plot.")
  pivoted <- pivot_longer(active_list, cols = everything(), names_to = "Var2")
  pivoted$Var1 <- rep(row.names(active_list), each = length(active_list))
  pivoted$value <- factor(pivoted$value, levels = c("FALSE", "TRUE"))
  pivoted$Var1 <- factor(pivoted$Var1, levels = purrr::map(ems, ~.$output_name))
  pivoted$Var2 <- factor(pivoted$Var2, levels = names(ems[[1]]$ranges))
  Var1 <- Var2 <- value <- NULL
  g <- ggplot(data = pivoted, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(colour = 'black') +
    scale_fill_manual(values = c('black', 'green'), drop = FALSE,
                      labels = c("FALSE", "TRUE"), name = "Active?") +
    labs(title = "Active variables", x = "Parameter", y = "Output") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(g)
}

#' Plot proposed points
#'
#' A wrapper around R's base plot to show proposed points
#'
#' Given a set of points proposed from emulators at a given wave, it's often useful to look at
#' how they are spread and where in parameter space they tend to lie relative to the original
#' ranges of the parameters. This function provides pairs plots of the parameters, with the
#' bounds of the plots calculated with respect to the parameter ranges provided.
#'
#' @param points The points to plot
#' @param ranges The parameter ranges
#' @param p_size The size of the plotted points (passed to \code{cex})
#'
#' @return The corresponding pairs plot
#'
#' @family visualisation tools
#' @export
#' @examples
#'  plot_wrap(SIRSample$training[,1:3], SIREmulators$ems[[1]]$ranges)
#'
plot_wrap <- function(points, ranges = NULL, p_size = 0.5) { #nocov start
  if (is.null(ranges))
    boundary_points <- setNames(
      do.call(
        'cbind.data.frame',
        purrr::map(names(points),
                   ~c(min(points[,.]), max(points[,.])))), names(points))
  else
    boundary_points <- setNames(
      do.call(
        'rbind.data.frame',
        purrr::map(
          1:2,
          ~purrr::map_dbl(ranges, function(x) x[[.]]))), names(ranges))
  plot(rbind(points, boundary_points),
       pch = 16, cex = p_size,
       col = c(rep('black', nrow(points)), 'white', 'white'))
} #nocov end
