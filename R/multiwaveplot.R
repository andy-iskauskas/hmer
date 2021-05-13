#' Multiple Wave Point Plotting
#'
#' Given multiple waves of points, produces pairs plots.
#'
#' Subsequent waves are overlaid on the same pairs plots, to determine the
#' evolution of the non-implausible region. One-dimensional density plots
#' are also created on the diagonal .
#'
#' @importFrom GGally ggpairs
#'
#' @param waves The list of data.frames, one for each set of points at that wave.
#' @param input_names The input names to be plotted.
#' @param surround If true, points are surrounded by black boundaries.
#' @param p_size The size of the points. Smaller values are better for high-dimensional spaces.
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#'  wave_points(GillespieMultiWaveData, c('aSI', 'aIR', 'aSR'))
#'  wave_points(GillespieMultiWaveData, c('aSI', 'aIR', 'aSR'), TRUE, 0.8)
wave_points <- function(waves, input_names, surround = FALSE, p_size = 1.5) {
  wave <- NULL
  out_list <- list()
  for (i in 0:(length(waves)-1)) {
    out_list[[i+1]] <- setNames(cbind(waves[[i+1]][,input_names], rep(i, nrow(waves[[i+1]]))), c(input_names, 'wave'))
  }
  total_data <- do.call('rbind', out_list)
  total_data$wave <- factor(total_data$wave)
  plotfun <- function(data, mapping) {
    g <- ggplot(data = data, mapping = mapping) +
      geom_point(cex = p_size)
    if(surround) g <- g + geom_point(cex = p_size, pch = 1, colour = 'black')
    return(g)
  }
  pal <- viridis::viridis(length(waves), option = "D", direction = -1)
  g <- ggpairs(total_data, columns = 1:length(input_names), aes(colour = wave),
          lower = list(continuous = plotfun),
          upper = 'blank', title = "Wave Points Location", progress = FALSE) +
    scale_colour_manual(values = pal) +
    scale_fill_manual(values = alpha(pal, 0.5)) +
    theme_bw()
  return(g)
}

#' Multiple Wave Output Plotting
#'
#' Given multiple waves of points, produces pairs plots of the outputs.
#'
#' This function operates in a similar fashion to \code{\link{wave_points}} - the main
#' difference is that the output values are plotted. Consequently, the set of targets is required
#' to overlay the region of interest onto the plot.
#'
#' @importFrom GGally ggpairs ggally_densityDiag
#'
#' @param waves The list of data.frame, one for each set of outputs at that wave.
#' @param targets The output targets.
#' @param output_names The outputs to plot.
#' @param surround As in \code{\link{wave_points}}.
#' @param p_size As in \code{\link{wave_points}}.
#' @param l_wid The width of the lines that create the target boxes.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#'  wave_values(GillespieMultiWaveData, sample_emulators$targets, surround = TRUE, p_size = 1)
#'  wave_values(GillespieMultiWaveData, sample_emulators$targets, c('nS', 'nI'), l_wid = 0.8)
wave_values <- function(waves, targets, output_names = names(targets), surround = FALSE, p_size = 1.5, l_wid = 1.5) {
  if (!is.null(targets[[1]]$val))
    targets <- purrr::map(targets, ~c(.$val - 3*.$sigma, .$val + 3*.$sigma))
  wave <- NULL
  out_list <- list()
  for (i in 0:(length(waves)-1)) {
    out_list[[i+1]] <- setNames(cbind(waves[[i+1]][,output_names], rep(i, nrow(waves[[i+1]]))), c(output_names, 'wave'))
  }
  total_data <- do.call('rbind', out_list)
  total_data$wave <- factor(total_data$wave)
  pal <- viridis::viridis(length(waves), option = "D", direction = -1)
  lfun <- function(data, mapping, targets, zoom = F) {
    xname <- rlang::quo_get_expr(mapping$x)
    yname <- rlang::quo_get_expr(mapping$y)
    g <- ggplot(data = data, mapping = mapping) +
      geom_point(cex = p_size) +
      scale_colour_manual(values = pal) +
      theme_bw() +
      geom_rect(xmin = targets[[xname]][1], xmax = targets[[xname]][2], ymin = targets[[yname]][1], ymax = targets[[yname]][2], colour = 'red', fill = NA, lwd = l_wid)
    if (zoom) {
      xrange <- (targets[[xname]][2]-targets[[xname]][1])/2
      yrange <- (targets[[yname]][2]-targets[[yname]][1])/2
      g <- g + coord_cartesian(xlim = c(targets[[xname]][1] - xrange, targets[[xname]][2] + xrange), ylim = c(targets[[yname]][1] - yrange, targets[[yname]][2] + yrange)) + theme_void()
    }
    if (surround)
      g <- g + geom_point(cex = p_size, pch = 1, colour = 'black')
    return(g)
  }
  dfun <- function(data, mapping, targets) {
    xname = rlang::as_string(rlang::quo_get_expr(mapping$x))
    g <- ggally_densityDiag(data = data, mapping = mapping, alpha = 0.4) +
      geom_vline(xintercept = targets[[xname]], colour = 'red', lwd = l_wid) +
      scale_fill_manual(values = pal) +
      theme_bw()
    return(g)
  }
  ggpairs(total_data, columns = 1:length(output_names), mapping = aes(colour = wave, group = wave),
          lower = list(continuous = wrap(lfun, targets = targets)),
          diag = list(continuous = wrap(dfun, targets = targets)),
          upper = list(continuous = wrap(lfun, targets = targets, zoom = TRUE)), progress = FALSE)
}

#' Plot simulator outputs for multiple waves
#'
#' Plots the simulator results for points at successive waves.
#'
#' The values plotted are the outputs from the simulator; the points passed to it are the
#' points suggested by that wave of emulators. By default, wave 0 is included. A colour
#' scheme is chosen outright for all invocations of this function: it is a 10-colour
#' palette. If more waves are required, then an alternative palette should be selected.
#'
#' @import ggplot2
#'
#' @param wave_points The set of wave points, as a list of data.frames
#' @param z The set of target values for each output
#' @param zero_in Is wave zero included? Default: TRUE
#' @param palette If a larger palette is required, it should be supplied here.
#' @param wave_numbers Which waves to plot. If not supplied, all waves are plotted.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#'  simulator_plot(GillespieMultiWaveData, sample_emulators$targets)
#'  simulator_plot(GillespieMultiWaveData[2:4], sample_emulators$targets,
#'   zero_in = FALSE, wave_numbers = c(1,3))
#'
simulator_plot <- function(wave_points, z, zero_in = TRUE, palette = NULL, wave_numbers = seq(ifelse(zero_in, 0, 1), length(wave_points)-ifelse(zero_in, 1, 0))) {
  variable <- value <- run <- wave <- val <- sigma <- NULL
  output_names <- names(z)
  sim_runs <- do.call('rbind', purrr::map(wave_numbers, ~data.frame(wave_points[[.+ifelse(zero_in, 1, 0)]][,output_names], wave = .)))
  sim_runs$run <- 1:length(sim_runs[,1])
  melted <- reshape2::melt(sim_runs, id.vars = c('run', 'wave'))
  melted$wave = as.factor(melted$wave)
  if (is.null(palette)) pal <- viridisLite::viridis(10, option = 'plasma', direction = -1)
  else pal <- palette
  pal <- pal[seq_along(pal) %in% (wave_numbers+ifelse(zero_in, 1, 0))]
  obs <- data.frame(variable = names(z), val = purrr::map_dbl(z, ~.$val), sigma = purrr::map_dbl(z, ~.$sigma))
  g <- ggplot(data = melted, aes(x = variable, y = value)) +
    geom_line(aes(group = run, colour = wave)) +
    scale_colour_manual(values = pal) +
    geom_point(data = obs, aes(x = variable, y = val)) +
    geom_errorbar(data = obs, aes(y = val, ymax = val + 3*sigma, ymin = val - 3*sigma), width = 0.1, size = 1.25) +
    labs(title = "Simulator evaluations at wave points") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(g)
}

#' Emulator Variance across waves
#'
#' Plots the emulator variance for each output across emulator waves.
#'
#' It is instructive to look at the change in emulator variance over successive
#' waves, rather than across successive outputs. This function provides a means
#' of doing so quickly for each emulator output.
#'
#' A 2d slice is taken across the input space, where mid-range values of any
#' non-plotted parameters are fixed. The emulator variance (or standard
#' deviation) is then calculated across the two-dimensional subspace for each
#' wave, and for each output.
#'
#' @import ggplot2
#' @importFrom patchwork wrap_plots plot_spacer
#' @importFrom purrr %>%
#' @importFrom viridis scale_fill_viridis
#'
#' @param waves A list of lists of \code{\link{Emulator}} objects, corresponding to the waves
#' @param output_names The list of desired outputs to be plotted
#' @param plot_dirs The (two) input parameters to be plotted.
#' @param wave_numbers A numeric vector of which waves to plot.
#' @param ppd The number of grid points per plotting dimension.
#' @param sd Should the standard deviation be plotted instead of the variance? Default: FALSE
#'
#' @return The ggplot object(s).
#'
#' @export
#'
#' @examples
#'  outputs <- c('nS', 'nI', 'nR')
#'  wave_variance(GillespieMultiWaveEmulators, outputs, ppd = 5)
#'  wave_variance(GillespieMultiWaveEmulators, c('nI', 'nR'),
#'   plot_dirs = c('aIR', 'aSR'), ppd = 5, sd = TRUE)
#'
wave_variance <- function(waves, output_names, plot_dirs = names(waves[[1]][[1]]$ranges)[1:2], wave_numbers = 1:length(waves), ppd = 20, sd = FALSE) {
  concurrent_plots <- max(3, length(waves[wave_numbers]))
  for (i in 1:length(waves)) {
    if (is.null(names(waves[[i]]))) {
      if (length(waves[[i]]) == length(output_names)) names(waves[[i]]) <- output_names
      else stop("One or more waves is missing an output emulation. Please specify emulator names explicitly.")
    }
  }
  if (length(plot_dirs) != 2) stop("Two input directions must be specified.")
  variable <- value <- NULL
  main_ranges <- waves[[1]][[output_names[1]]]$ranges
  grid_ranges <- setNames(data.frame(purrr::map(names(main_ranges), function(x) {
    if (x %in% plot_dirs) main_ranges[[x]]
    else rep((main_ranges[[x]][1]+main_ranges[[x]][2])/2, 2)
  })), names(main_ranges))
  on_grid <- setNames(expand.grid(seq(grid_ranges[[plot_dirs[1]]][[1]], grid_ranges[[plot_dirs[1]]][[2]], length.out = ppd),
                                  seq(grid_ranges[[plot_dirs[2]]][[1]], grid_ranges[[plot_dirs[2]]][[2]], length.out = ppd)),
                      plot_dirs)
  for (i in names(grid_ranges)) {
    if (!i %in% plot_dirs) on_grid[,i] <- grid_ranges[[i]][1]
  }
  on_grid <- on_grid[,names(main_ranges)]
  output <- purrr::map(output_names, ~setNames(cbind(on_grid, data.frame(purrr::map(waves[wave_numbers], function(x) {
    if (!sd) x[[.]]$get_cov(on_grid)
    else sqrt(x[[.]]$get_cov(on_grid))
  }))), c(names(main_ranges), wave_numbers)))
  melted_frames <- setNames(purrr::map(output, ~reshape2::melt(., id.vars = names(main_ranges))), output_names)
  for (i in 1:length(melted_frames)) melted_frames[[i]]$output <- output_names[i]
  for (i in 1:ceiling(length(melted_frames)/concurrent_plots)) {
    current_end <- min(concurrent_plots*i, length(melted_frames))
    dat <- melted_frames[(concurrent_plots*(i-1)+1):current_end]
    dat <- do.call('rbind', dat)
    plot_bins <- round(seq(0, max(dat$value), length.out = 25))
    g <- ggplot(data = dat, aes(x = dat[,plot_dirs[1]], y = dat[,plot_dirs[2]], group = variable)) +
      geom_contour_filled(aes(z = value), colour = 'black', breaks = plot_bins) +
      scale_fill_viridis(discrete = TRUE, option = 'plasma', labels = plot_bins, name = ifelse(sd, 'SD[f(x)]', 'Var[f(x)]')) +
      facet_grid(rows = vars(variable), cols = vars(output), labeller = labeller(variable = function(x) paste("Wave", x))) +
      labs(title = paste(ifelse(sd, "Standard Deviation", "Variance"), 'across waves'), x = plot_dirs[1], y = plot_dirs[2]) +
      theme_minimal()
    print(g)
  }
}
