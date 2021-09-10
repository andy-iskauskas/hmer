# Plotting emulator expectation
exp_plot <- function(em, plotgrid = NULL, ppd = 30) {
  ranges <- em$ranges
  if (is.null(plotgrid)) {
    plotgrid <- setNames(expand.grid(seq(ranges[[1]][1], ranges[[1]][2], length.out = ppd), seq(ranges[[2]][1], ranges[[2]][2], length.out = ppd)), names(ranges)[1:2])
    for (i in 3:length(ranges)) {
      plotgrid[[names(ranges)[i]]] <- sum(ranges[[i]])/2
    }
  }
  em_exp <- em$get_exp(plotgrid[,names(ranges)])
  grid_data <- setNames(cbind(plotgrid[,1:2], em_exp), c(names(plotgrid)[1:2],"E"))
  g <- tryCatch(
    {
      ## To account for rounding problems in ggplot's pretty_isoband_levels
      rng <- range(grid_data$E)
      if (diff(rng) == 0) warning("All output values are identical.")
      ac <- signif(diff(rng), 1)/10
      rng[1] <- floor(rng[1]/ac)*ac
      rng[2] <- ceiling(rng[2]/ac)*ac
      bks <- seq(rng[1], rng[2], length.out = 26)
      x_pos <- as.integer(factor(grid_data[,1], levels = sort(unique(grid_data[,1]))))
      y_pos <- as.integer(factor(grid_data[,2], levels = sort(unique(grid_data[,2]))))
      raster <- matrix(NA_real_, nrow = max(y_pos), ncol = max(x_pos))
      raster[cbind(y_pos, x_pos)] <- grid_data$E
      ibs <- isoband::isobands(x = sort(unique(grid_data[,1])), y = sort(unique(grid_data[,2])), z = raster, levels_low = bks[-length(bks)], levels_high = bks[-1])
      int_lo <- gsub(":.*$", "", names(ibs))
      int_hi <- gsub("^[^:]*:", "", names(ibs))
      lab_lo <- format(as.numeric(int_lo), digits = 3, trim = TRUE)
      lab_hi <- format(as.numeric(int_hi), digits = 3, trim = TRUE)
      lab_check <- sprintf("(%s, %s]", lab_lo, lab_hi)
      if (length(unique(lab_check)) == 1) warning("Can't produce contours natively due to internal accuracy issues.")
      bns <- min(length(unique(lab_check)), 25)
      ggplot(data = grid_data, aes(x = grid_data[,1], y = grid_data[,2])) + geom_contour_filled(aes(z = grid_data[,'E']), bins = bns, colour = 'black') + viridis::scale_fill_viridis(discrete = TRUE, option = "magma", name = "exp")
    },
    warning = function(w) {
      exp_breaks <- seq(min(em_exp) - diff(range(em_exp))/(2 * 23), max(em_exp) + diff(range(em_exp))/(2*23), length.out = 25)
      exp_breaks <- unique(signif(exp_breaks, 10))
      if (length(exp_breaks) == 1) exp_breaks <- c(min(em_exp)-1e-6, max(em_exp)+1e-6)
      intervals <- findInterval(em_exp, exp_breaks)
      fake_breaks <- 1:length(exp_breaks)
      ggplot(data = grid_data, aes(x = grid_data[,1], y = grid_data[,2])) + geom_contour_filled(aes(z = intervals), breaks = fake_breaks, colour = 'black') + viridis::scale_fill_viridis(discrete = TRUE, option = "magma", name = "exp", labels = function(b) {signif(exp_breaks[as.numeric(stringr::str_extract(b, "\\d+"))], 6)})
    }
  )
  g <- g +
    labs(title = paste(em$output_name, "Emulator Expectation"), x = names(grid_data)[1], y = names(grid_data)[2]) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_minimal()
  return(g)
}

# Plotting emulator variance
var_plot <- function(em, plotgrid = NULL, ppd = 30, sd = FALSE) {
  ranges <- em$ranges
  if (is.null(plotgrid)) {
    plotgrid <- setNames(expand.grid(seq(ranges[[1]][1], ranges[[1]][2], length.out = ppd), seq(ranges[[2]][1], ranges[[2]][2], length.out = ppd)), names(ranges)[1:2])
    for (i in 3:length(ranges)) {
      plotgrid[[names(ranges)[i]]] <- sum(ranges[[i]])/2
    }
  }
  if(sd)
    em_cov <- sqrt(em$get_cov(plotgrid[,names(ranges)]))
  else
    em_cov <- em$get_cov(plotgrid[,names(ranges)])
  grid_data <- setNames(cbind(plotgrid[,1:2], em_cov), c(names(plotgrid)[1:2], "V"))
  g <- tryCatch(
    {
      ## To account for rounding problems in ggplot's pretty_isoband_levels
      rng <- range(grid_data$V)
      if (diff(rng) == 0) warning("All output values are identical.")
      ac <- signif(diff(rng), 1)/10
      rng[1] <- floor(rng[1]/ac)*ac
      rng[2] <- ceiling(rng[2]/ac)*ac
      bks <- seq(rng[1], rng[2], length.out = 26)
      x_pos <- as.integer(factor(grid_data[,1], levels = sort(unique(grid_data[,1]))))
      y_pos <- as.integer(factor(grid_data[,2], levels = sort(unique(grid_data[,2]))))
      raster <- matrix(NA_real_, nrow = max(y_pos), ncol = max(x_pos))
      raster[cbind(y_pos, x_pos)] <- grid_data$V
      ibs <- isoband::isobands(x = sort(unique(grid_data[,1])), y = sort(unique(grid_data[,2])), z = raster, levels_low = bks[-length(bks)], levels_high = bks[-1])
      int_lo <- gsub(":.*$", "", names(ibs))
      int_hi <- gsub("^[^:]*:", "", names(ibs))
      lab_lo <- format(as.numeric(int_lo), digits = 3, trim = TRUE)
      lab_hi <- format(as.numeric(int_hi), digits = 3, trim = TRUE)
      lab_check <- sprintf("(%s, %s]", lab_lo, lab_hi)
      if (length(unique(lab_check)) == 1) warning("Can't produce contours natively due to internal accuracy issues.")
      bns <- min(length(unique(lab_check)), 25)
      ggplot(data = grid_data, aes(x = grid_data[,1], y = grid_data[,2])) + geom_contour_filled(aes(z = grid_data[,'V']), bins = bns, colour = 'black') + viridis::scale_fill_viridis(discrete = TRUE, option = "plasma", name = if(sd) "sd" else "var")
    },
    warning = function(w) {
      cov_breaks <- seq(min(em_cov) - diff(range(em_cov))/(2 * 23), max(em_cov) + diff(range(em_cov))/(2*23), length.out = 25)
      cov_breaks <- unique(signif(cov_breaks, 10))
      if (length(cov_breaks) == 1) cov_breaks <- c(min(em_cov)-1e-6, max(em_cov)+1e-6)
      intervals <- findInterval(em_cov, cov_breaks)
      fake_breaks <- 1:length(cov_breaks)
      ggplot(data = grid_data, aes(x = grid_data[,1], y = grid_data[,2])) + geom_contour_filled(aes(z = intervals), breaks = fake_breaks, colour = 'black') + viridis::scale_fill_viridis(discrete = TRUE, option = "plasma", name = if (sd) "sd" else "var", labels = function(b) {signif(cov_breaks[as.numeric(stringr::str_extract(b, "\\d+"))], 6)})
    }
  )
  g <- g +
    labs(title = paste(em$output_name, "Emulator", (if(sd) "Standard Deviation" else "Variance")), x = names(grid_data)[1], y = names(grid_data)[2]) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_minimal()
  return(g)
}

# Plotting emulator variance
imp_plot <- function(em, z, plotgrid = NULL, ppd = 30, cb = FALSE, nth = NULL) {
  if (!is.null(nth)) {
    if (nth == 1) ns <- ""
    else if (nth == 2) ns <- "Second"
    else if (nth == 3) ns <- "Third"
    else ns <- paste0(nth, "th")
    ranges <- em[[1]]$ranges
  }
  else ranges <- em$ranges
  if (is.null(plotgrid)) {
    ranges <- em$ranges
    plotgrid <- setNames(expand.grid(seq(ranges[[1]][1], ranges[[1]][2], length.out = ppd), seq(ranges[[2]][1], ranges[[2]][2], length.out = ppd)), names(ranges)[1:2])
    for (i in 3:length(ranges)) {
      plotgrid[[names(ranges)[i]]] <- sum(ranges[[i]])/2
    }
  }
  imp_breaks <- c(0, 0.3, 0.7, 1, 1.3, 1.7, 2, 2.3, 2.7, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 15, 100)
  imp_names <- c(0, '', '', 1, '', '', 2, '', '', 3, '', '', '', 5, '', '', '', 10, 15, '')
  if (!is.null(nth)) {
    if (!"Emulator" %in% class(em)) {
      em_imp <- nth_implausible(em, plotgrid[names(em[[1]]$ranges)], z, n = nth, max_imp = 99)
    }
    else stop("Not all required parameters (emulator list, target list, nth) passed for nth maximum implausibility.")
  }
  else{
    em_imp <- em$implausibility(plotgrid[names(ranges)], z)
  }
  included <- c(purrr::map_lgl(imp_breaks[-1], ~any(em_imp < .)), TRUE)
  grid_data <- setNames(cbind(plotgrid[,1:2], em_imp), c(names(plotgrid)[1:2],"I"))
  col_scale <- if(cb) colourblind else redgreen
  g <- ggplot(data = grid_data, aes(x = grid_data[,1], y = grid_data[,2]))
  if (length(unique(em_imp)) == 1) {
    col <- col_scale[which.min(abs(imp_breaks - c(unique(em_imp))))]
    g <- g + geom_raster(aes(fill = grid_data[,"I"])) + scale_fill_gradient(low = col, high = col, name = "I")
  }
  else
    g <- g + geom_contour_filled(aes(z = grid_data[,"I"]), colour = 'black', breaks = imp_breaks[included]) + scale_fill_manual(values = col_scale[included], name = "I", labels = imp_names[included], guide = guide_legend(reverse = TRUE))
  g <- g + labs(title = paste(em$output_name, (if (is.null(nth)) "Emulator Implausibility" else paste(ns, "Maximum Implausibility"))), x = names(grid_data)[1], y = names(grid_data)[2]) +
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
#' determined by \code{var_name}: options are 'exp', 'var', 'sd', 'imp', and 'nimp', which
#' correspond to expectation, variance, standard deviation, implausibility, and nth-max
#' implausibility.
#'
#' By default, the slice varies in the first two parameters of the emulators, and all other
#' parameters are taken to be fixed at their mid-range values. This behaviour can be changed
#' with the \code{params} and \code{fixed_vals} parameters (see examples).
#'
#' If the statistic is 'exp', 'var' or 'sd', then the minimal set of parameters to pass to this
#' function are \code{ems} (which can be a list of emulators or a single one) and \code{var_name}.
#' If the statistic is 'imp' or 'nimp', then the \code{targets} must be supplied - it is not
#' necessary to specify the individual target for a single emulator plot. If the statistic is
#' 'nimp', then the level of maximum implausiblility can be chosen with the parameter \code{nth}.
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
#' @importFrom cowplot plot_grid
#'
#' @param ems An \code{\link{Emulator}} object, or a list thereof.
#' @param var_name The statistic to plot (see description or examples).
#' @param ppd The number of points per plotting dimension
#' @param targets If required, the targets from which to calculate implausibility
#' @param cb A boolean representing whether a colourblind-friendly plot is produced.
#' @param params Which two input parameters should be plotted?
#' @param fixed_vals For fixed input parameters, the values they are held at.
#' @param nth If plotting nth maximum implausibility, which level maximum to plot.
#'
#' @return A ggplot object, or collection thereof.
#' @export
#'
#' @examples
#'  # Reducing ppd to 10 for speed.
#'  emulator_plot(sample_emulators$ems, ppd = 10)
#'  emulator_plot(sample_emulators$ems$nS, ppd = 10)
#'  emulator_plot(sample_emulators$ems, var_name = 'var', ppd = 10, params = c('aIR', 'aSR'))
#'  emulator_plot(sample_emulators$ems, var_name = 'imp', ppd = 10,
#'   targets = sample_emulators$targets,
#'   fixed_vals = list(aSR = 0.02))
#'  emulator_plot(sample_emulators$ems, var_name = 'nimp', cb = TRUE,
#'   targets = sample_emulators$targets, nth = 2, ppd = 10)
#'
emulator_plot <- function(ems, var_name = 'exp', ppd = 30, targets = NULL, cb = FALSE, params = NULL, fixed_vals = NULL, nth = 1) {
  if ("Emulator" %in% class(ems)){
    ranges <- ems$ranges
    single_em <- TRUE
  }
  else {
    ranges <- ems[[1]]$ranges
    single_em <- FALSE
  }
  if (is.null(params) || length(params) != 2 || any(!params %in% names(ranges))) p_vals <- c(1,2)
  else p_vals <- which(names(ranges) %in% params)
  plotgrid <- setNames(expand.grid(seq(ranges[[p_vals[1]]][1], ranges[[p_vals[1]]][2], length.out = ppd), seq(ranges[[p_vals[2]]][1], ranges[[p_vals[2]]][2], length.out = ppd)), names(ranges)[p_vals])
  if (!is.null(fixed_vals) && all(names(fixed_vals) %in% names(ranges))) {
    for (i in 1:length(fixed_vals)) plotgrid[[names(fixed_vals)[i]]] <- fixed_vals[[names(fixed_vals)[i]]]
    used_names <- c(names(ranges)[p_vals], names(fixed_vals))
  }
  else used_names <- names(ranges)[p_vals]
  used_names <- unique(used_names)
  if (length(used_names) < length(ranges)) {
    unused_nms <- names(ranges)[which(!names(ranges) %in% used_names)]
    for (i in 1:length(unused_nms)) plotgrid[[unused_nms[i]]] <- sum(ranges[[unused_nms[i]]])/2
  }
  get_plot <- function(em) {
    if (var_name == 'exp') return(exp_plot(em, plotgrid, ppd))
    if (var_name == 'var') return(var_plot(em, plotgrid, ppd))
    if (var_name == 'sd') return(var_plot(em, plotgrid, ppd, sd = TRUE))
    if (var_name == 'imp') {
      if (is.null(targets)) stop("Cannot plot implausibility without target value.")
    }
    if (!is.null(targets$val)) return(imp_plot(em, targets, plotgrid, ppd, cb))
    else return(imp_plot(em, targets[[em$output_name]], plotgrid, ppd, cb))
  }
  if (single_em) return(get_plot(ems))
  if (var_name == 'nimp') return(imp_plot(ems, targets, plotgrid, ppd, cb, nth))
  else {
    plotlist <- purrr::map(ems, get_plot)
    for (i in 1:length(plotlist)) {
      plotlist[[i]]$layers[[1]]$show.legend = FALSE
    }
    return(cowplot::plot_grid(plotlist = plotlist, ncol = ceiling(sqrt(length(ems)))))
  }
}

#' Emulator Expectation Against Target Outputs
#'
#' Plots emulator expectation across the parameter space, with comparison to the corresponding
#' target values (with appropriate uncertainty).
#'
#' If a \code{points} data.frame is not provided, then points are sampled uniformly from the
#' input region. Otherwise, the provided points are used: for example, if a representative
#' sample current NROY space is available.
#'
#' @importFrom ggplot2 ggplot aes labs geom_line geom_point geom_errorbar
#' @importFrom reshape2 melt
#'
#' @param ems The \code{\link{Emulator}} objects.
#' @param targets A named list of observations, given in the usual form.
#' @param points A list of points at which the emulators should be evaluated.
#' @param npoints If no points are provided, the number of input points to evaluate at.
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#'  output_plot(sample_emulators$ems, sample_emulators$targets)
#'  output_plot(sample_emulators$ems, sample_emulators$targets, points = GillespieSIR)
output_plot <- function(ems, targets, points = NULL, npoints = 1000) {
  ranges <- if ("Emulator" %in% class(ems)) ems$ranges else ems[[1]]$ranges
  if (is.null(points)) {
    points <- data.frame(purrr::map(ranges, ~runif(npoints, .[1], .[2])))
  }
  em_exp <- setNames(data.frame(purrr::map(ems, ~.$get_exp(points))), names(targets))
  em_exp$run <- 1:nrow(points)
  em_exp <- reshape2::melt(em_exp, id.vars = 'run')
  for (i in 1:length(targets))
  {
    if (!is.atomic(targets[[i]])) targets[[i]] <- c(targets[[i]]$val - 3 * targets[[i]]$sigma, targets[[i]]$val + 3*targets[[i]]$sigma)
  }
  target_data <- data.frame(label = names(targets), mn = purrr::map_dbl(targets, ~.[1]), md = purrr::map_dbl(targets, mean), mx = purrr::map_dbl(targets, ~.[2]))
  variable <- value <- run <- mn <- md <- mx <- label <- NULL
  ggplot(data = em_exp, aes(x = variable, y = value)) +
    geom_line(colour = 'purple', aes(group = run), lwd = 1) +
    geom_point(data = target_data, aes(x = label, y = md), size = 2) +
    geom_errorbar(data = target_data, aes(x = label, y = md, ymin = mn, ymax = mx), width = .1, size = 1.25) +
    labs(title = "Emulator Runs versus Observations")
}

#' Plot Lattice of Emulator Implausibilities
#'
#' Plots a set of projections of the full-dimensional input space. The plots are:
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
#' @importFrom cowplot plot_grid get_legend
#' @importFrom stats xtabs
#' @import ggplot2
#'
#' @param ems The \code{\link{Emulator}} objects in question.
#' @param targets The corresponding target values.
#' @param ppd The number of points to sample per dimension.
#' @param cb Whether or not a colourblind-friendly plot should be produced.
#' @param cutoff The cutoff value for non-implausible points.
#' @param maxpoints The limit on the number of points to be evaluated.
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \donttest{
#'  plot_lattice(sample_emulators$ems, sample_emulators$targets, ppd = 10)
#'  plot_lattice(sample_emulators$ems$nS, sample_emulators$targets)
#' }
plot_lattice <- function(ems, targets, ppd = 20, cb = FALSE, cutoff = 3, maxpoints = 5e4) {
  get_count <- function(df, nbins, nms) {
    out_df <- setNames(expand.grid(1:nbins, 1:nbins), nms)
    out_df$Freq <- 0
    for (i in 1:nbins) {
      for (j in 1:nbins)
        out_df[out_df[,nms[1]]==i & out_df[,nms[2]]==j, "Freq"] <- sum(apply(df, 1, function(x) x[nms[1]]==i && x[nms[2]]==j))
    }
    x_form <- as.formula(paste("Freq ~", nms[1], "+", nms[2], sep = ""))
    return(xtabs(x_form, data = out_df))
  }
  ranges <- if("Emulator" %in% class(ems)) ems$ranges else ems[[1]]$ranges
  if (ppd^length(ranges) > maxpoints) {
    plotgrid <- setNames(data.frame(do.call('cbind', purrr::map(ranges, ~runif(maxpoints, .[[1]], .[[2]])))), names(ranges))
  }
  else {
    dim_unif <- purrr::map(ranges, ~seq(.[[1]], .[[2]], length.out = ppd))
    plotgrid <- expand.grid(dim_unif)
  }
  imps <- nth_implausible(ems, plotgrid, targets)
  df <- setNames(data.frame(cbind(plotgrid, imps)), c(names(ranges), "I"))
  bins <- if (ppd^length(ranges) > maxpoints) 20 else ppd
  df_filtered <- subset(df, I <= cutoff)
  one_dim <- setNames(data.frame(do.call('cbind', purrr::map(seq_along(ranges), function(x) {
    binpoints <- seq(ranges[[x]][1]-1e-6, ranges[[x]][2]+1e-6, length.out = bins+1)
    intervals <- findInterval(df[,names(ranges)[[x]]], binpoints)
    tot_in_bin <- c(tabulate(intervals, bins))
    tot_in_bin[tot_in_bin == 0] <- 1
    tot_NROY <- c(tabulate(findInterval(df_filtered[,names(ranges)[[x]]], binpoints), bins))
    props <- tot_NROY/tot_in_bin
    return(props[intervals])
  }))), names(ranges))
  variable_combs <- unlist(purrr::map(outer(names(ranges), names(ranges), paste)[upper.tri(outer(names(ranges), names(ranges), paste))],
                                      ~strsplit(., " ")), recursive = FALSE)
  two_dim <- purrr::map(variable_combs, function(x) {
    binpoints <- purrr::map(seq_along(x), ~seq(ranges[[x[.]]][1]-1e-6, ranges[[x[.]]][2]+1e-6, length.out = bins+1))
    intervals <- setNames(data.frame(do.call('cbind', purrr::map(seq_along(x), ~findInterval(df[,x[.]], binpoints[[.]])))), x)
    tot_in_bin <- table(intervals)
    NROY_class <- setNames(data.frame(do.call('cbind', purrr::map(seq_along(x), ~findInterval(df_filtered[,x[.]], binpoints[[.]])))), x)
    tot_NROY <- get_count(NROY_class, bins, x)
    tot_in_bin[tot_in_bin == 0] <- 1
    props <- tot_NROY/tot_in_bin
    return(apply(intervals, 1, function(x) props[x[1], x[2]]))
  })
  mins <- purrr::map(variable_combs, function(x) {
    binpoints <- purrr::map(seq_along(x), ~seq(ranges[[x[.]]][1]-1e-6, ranges[[x[.]]][2]+1e-6, length.out = bins+1))
    intervals <- setNames(data.frame(do.call('cbind', purrr::map(seq_along(x), ~findInterval(df[,x[.]], binpoints[[.]])))), x)
    minimp <- matrix(rep(0, bins^2), nrow = bins)
    for (i in 1:bins) {
      for (j in 1:bins)
        minimp[i,j] <- min(df[intervals[,1] == i & intervals[,2] == j, "I"])
    }
    return(apply(intervals, 1, function(x) minimp[x[1], x[2]]))
  })
  full_df <- setNames(data.frame(cbind(one_dim, two_dim, mins)),
                      c(paste0(names(ranges), 'op'), paste0(purrr::map_chr(variable_combs, ~paste(., collapse = "")), "op"), paste0(purrr::map_chr(variable_combs, ~paste(., collapse = "")), "min")))
  if (ppd^length(ranges) > maxpoints) {
    grid_df <- expand.grid(purrr::map(ranges, ~seq(.[[1]], .[[2]], length.out = bins)))
    ordering <- apply(grid_df, 1, function(x) which.min(apply(plotgrid[,names(ranges)], 1, function(y) sum((x-y)^2))))
    full_df <- data.frame(cbind(grid_df, full_df[ordering,]))
  }
  else {
    full_df <- data.frame(cbind(expand.grid(purrr::map(ranges, ~seq(.[[1]], .[[2]], length.out = bins))), full_df))
  }
  u_mat <- l_mat <- matrix(rep("", length(ranges)^2), nrow = length(ranges))
  comb_names <- purrr::map(variable_combs, paste, collapse = "")
  u_mat[upper.tri(u_mat)] <- paste0(comb_names, 'min')
  l_mat[upper.tri(l_mat)] <- paste0(comb_names, 'op')
  l_mat <- t(l_mat)
  name_mat <- matrix(paste(u_mat, l_mat, sep = ""), nrow = length(ranges))
  diag(name_mat) <- paste0(names(ranges), 'op')
  cols <- if(cb) colourblind else redgreen
  imp_breaks <- c(0, 0.3, 0.7, 1, 1.3, 1.7, 2, 2.3, 2.7, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 15, 100)
  imp_names <- c(0, '', '', 1, '', '', 2, '', '', 3, '', '', '', 5, '', '', '', 10, 15, '')
  plot_list <- unlist(purrr::map(seq_along(1:nrow(name_mat)), function(x) {
    purrr::map(seq_along(1:ncol(name_mat)), function(y) {
      if (x == y)
        g <- ggplot(data = full_df, aes(x = full_df[,names(ranges)[x]])) +
          geom_smooth(aes(y = full_df[,name_mat[x,x]]), colour = 'black')
      else {
        g <- ggplot(data = full_df, aes(x = full_df[,names(ranges)[y]], y = full_df[,names(ranges)[x]]))
        if (x < y)
          g <- g + geom_contour_filled(aes(z = full_df[,name_mat[x,y]]), breaks = imp_breaks, colour = 'black') +
            scale_fill_manual(values = cols, labels = imp_names, name = "Min. I", guide = guide_legend(reverse = TRUE), drop = FALSE)
        else
          g <- g + geom_raster(aes(fill = full_df[,name_mat[x,y]]), interpolate = TRUE) +
            scale_fill_gradient(low = 'black', high = 'white', breaks = seq(0, 1, by = 0.1), name = "Op. Depth")
      }
      xlab <- ylab <- NULL
      if (x == nrow(name_mat)) xlab <- names(ranges)[y]
      if (y == 1) ylab <- names(ranges)[x]
      return(g + labs(x = xlab, y = ylab) + theme_minimal())
    })
  }), recursive = FALSE)
  legend_list <- list(cowplot::get_legend(plot_list[[2]]), cowplot::get_legend(plot_list[[length(ranges)+1]]))
  for (i in 1:length(plot_list)) plot_list[[i]] <- plot_list[[i]] + theme(legend.position = 'none')
  return(suppressMessages(cowplot::plot_grid(cowplot::plot_grid(plotlist = plot_list), cowplot::plot_grid(plotlist = legend_list, ncol = 1, rel_heights = c(1, 1/length(ranges))), rel_widths = c(1, 0.1))))
}
