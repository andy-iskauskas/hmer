#' Multiple Wave Point Plotting
#'
#' Given multiple waves of points, produces pairs plots
#'
#' Subsequent waves are overlaid on the same pairs plots, to determine the
#' evolution of the non-implausible region. One-dimensional density plots
#' are also created on the diagonal.
#'
#' @importFrom GGally ggpairs
#'
#' @param waves The list of data.frames, one for each set of points at that wave.
#' @param input_names The input names to be plotted.
#' @param surround If true, points are surrounded by black boundaries.
#' @param p_size The size of the points. Smaller values are better for high-dimensional spaces.
#' @param zero_in Is a wave 0 included in the waves list?
#' @param wave_numbers Which waves to plot
#' @param ... Optional parameters (not to be used directly)
#'
#' @return A ggplot object
#'
#' @family visualisation tools
#' @export
#'
#' @examples
#'  wave_points(SIRMultiWaveData, c('aSI', 'aIR', 'aSR'))
#'  \donttest{
#'      wave_points(SIRMultiWaveData, c('aSI', 'aIR', 'aSR'), TRUE, 0.8)
#'      # For many plots, it may be helpful to manually modify the font size
#'      wave_points(SIRMultiWaveData, c('aSI', 'aIR', 'aSR')) +
#'       ggplot2::theme(text = ggplot2::element_text(size = 5))
#'  }
wave_points <- function(waves, input_names, surround = FALSE, p_size = 1.5,
                        zero_in = TRUE,
                        wave_numbers = ifelse(
                          zero_in, 0, 1):
                          (length(waves)-ifelse(zero_in, 1, 0)), ...) {
  wave <- NULL
  out_list <- list()
  wave_names <- wave_numbers
  if(zero_in) wave_numbers <- wave_numbers + 1
  for (i in 1:length(wave_numbers)) {
      out_list[[i]] <- setNames(
        cbind(waves[[wave_numbers[i]]][,input_names],
              rep(wave_names[i], nrow(waves[[wave_numbers[i]]]))), c(input_names, 'wave'))
  }
  total_data <- do.call('rbind', out_list)
  total_data$wave <- factor(total_data$wave)
  plotfun <- function(data, mapping) {
    g <- ggplot(data = data, mapping = mapping) +
      geom_point(cex = p_size)
    if(surround) g <- g + geom_point(cex = p_size, pch = 1, colour = 'black')
    return(g)
  }
  pal <- viridis::viridis(max(length(waves), max(wave_numbers)),
                          option = "D", direction = -1)[wave_numbers]
  g <- ggpairs(total_data, columns = seq_along(input_names), aes(colour = wave),
          lower = list(continuous = plotfun),
          upper = 'blank',
          title = "Wave Points Location", progress = FALSE, legend = 1) +
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
#' To ensure that the wave numbers provided in the legend match, one should provide waves
#' as a list of data.frames with the earliest wave at the start of the list.
#'
#' The parameters \code{which_wave} and \code{upper_scale} control the level of `zoom' on
#' each of the lower-triangular and upper-triangular plots, respectively. For the lower
#' plots, \code{which_wave} determines which of the provided waves is to be used to determine
#' the output ranges to plot with respect to: generally, higher \code{which_wave} values
#' result in a more zoomed-in plot. For the upper plots, \code{upper_scale} determines the
#' plot window via a multiple of the target bounds: higher values result in a more zoomed-out
#' plot. If not provided, these default to \code{which_wave=0} (or 1 if no wave 0 is given)
#' and \code{upper_scale = 1}. If the value provided to \code{which_wave} does not correspond
#' to a provided wave (or one explicitly not included in \code{wave_numbers}), it defaults to
#' the closest available wave to the value of \code{which_wave}.
#'
#' If \code{ems} is provided, it should follow the same structure as \code{waves}: at the very
#' least, it should contain all emulators trained over the course of the waves. The emulator
#' predictions for a target are made by the emulator for that target whose ranges are the
#' smallest such that contain the point.
#'
#' @importFrom GGally ggpairs ggally_densityDiag
#'
#' @param waves The list of data.frames, one for each set of outputs at that wave.
#' @param targets The output targets.
#' @param output_names The outputs to plot.
#' @param ems If provided, plots the emulator expectations and 3-standard deviations.
#' @param surround As in \code{\link{wave_points}}.
#' @param restrict Should the plotting automatically restrict to failing target windows?
#' @param p_size As in \code{\link{wave_points}}.
#' @param l_wid The width of the lines that create the target boxes.
#' @param zero_in Is a wave 0 included in the waves list?
#' @param wave_numbers Which waves to plot.
#' @param which_wave Scaling for lower plots (see description)
#' @param upper_scale Scaling for upper plots (ibid)
#' @param ... Optional parameters (not to be used directly)
#'
#' @return A ggplot object.
#'
#' @family visualisation tools
#' @export
#'
#' @examples
#'  wave_values(SIRMultiWaveData, SIREmulators$targets, surround = TRUE, p_size = 1)
#'  \donttest{
#'    wave_values(SIRMultiWaveData, SIREmulators$targets, c('nS', 'nI'), l_wid = 0.8)
#'      wave_values(SIRMultiWaveData, SIREmulators$targets, l_wid = 0.8,
#'       wave_numbers = c(0, 1, 3), which_wave = 2, upper_scale =  1.5)
#'      # For many plots, it may be helpful to manually modify the font size
#'      wave_values(SIRMultiWaveData, SIREmulators$targets) +
#'       ggplot2::theme(text = ggplot2::element_text(size = 5))
#'  }
wave_values <- function(waves, targets, output_names = names(targets),
                        ems = NULL, surround = FALSE,
                        restrict = FALSE, p_size = 1.5, l_wid = 1.5,
                        zero_in = TRUE,
                        wave_numbers = ifelse(zero_in, 0, 1):
                          (length(waves)-ifelse(zero_in, 0, 1)),
                        which_wave = ifelse(zero_in, 0, 1),
                        upper_scale = 1, ...) {
  waves <- tryCatch(
    purrr::map(waves, data.frame),
    error = function(e) {
      stop(paste("Cannot convert wave objects to data.frame:", e))
    }
  )
  if (!is.null(ems)) ems <- collect_emulators(ems)
  if (!which_wave %in% wave_numbers) {
    which_wave <- wave_numbers[which.min(abs(wave_numbers - which_wave))]
  }
  for (i in seq_along(targets)) {
    if (!is.atomic(targets[[i]]))
      targets[[i]] <- c(targets[[i]]$val - 3*targets[[i]]$sigma,
                        targets[[i]]$val + 3*targets[[i]]$sigma)
  }
  in_range <- function(data, ranges) {
    apply(
      data, 1,
      function(x) all(
        purrr::map_lgl(
          seq_along(ranges),
          ~x[.] >= ranges[[.]][1] && x[.] <= ranges[[.]][2])))
  }
  wave <- NULL
  out_list <- list()
  var_list <- list()
  for (i in ifelse(zero_in, 0, 1):(length(waves)-ifelse(zero_in, 1, 0))) {
    if (i %in% wave_numbers) {
      if (!is.null(ems)) {
        collected_ems <- collect_emulators(ems)
        which_ems_fit <- do.call(
          'cbind.data.frame',
          purrr::map(
            collected_ems,
            ~in_range(waves[[i+1]][,names(.$ranges)], .$ranges)))
        which_em_use <- setNames(
          do.call(
            'cbind.data.frame',
            purrr::map(
              names(targets),
              ~apply(
                t(
                  apply(
                    which_ems_fit, 1,
                    function(x) names(x) == . & x)), 1,
                function(x) which(x)[1]))), names(targets))
        em_exps <- do.call(
          'cbind.data.frame',
          purrr::map(ems, ~.$get_exp(waves[[i+1]])))
        em_vars <- do.call(
          'cbind.data.frame',
          purrr::map(ems, ~sqrt(.$get_cov(waves[[i+1]]))))
        out_exps <- setNames(
          do.call(
            'rbind.data.frame',
            purrr::map(
              seq_len(nrow(which_em_use)),
              function(x)
                purrr::map_dbl(
                  seq_along(which_em_use),
                  function(y) em_exps[x,y]))), names(targets))
        out_vars <- setNames(
          do.call(
            'rbind.data.frame',
            purrr::map(
              seq_len(nrow(which_em_use)),
              function(x) purrr::map_dbl(
                seq_along(which_em_use),
                function(y) em_vars[x,y]))), names(targets))
        out_list[[i+1]] <- setNames(
          cbind(
            out_exps, rep(i, nrow(waves[[i+1]]))), c(output_names, 'wave'))
        var_list[[i+1]] <- setNames(
          cbind(
            out_vars,
            rep(i,
                nrow(waves[[i+1]]))),
          c(paste0(output_names, "V", sep = ""), 'wave'))
      }
      else{
        out_list[[i+1]] <- setNames(
          cbind(
            waves[[i+1]][,output_names],
            rep(i, nrow(waves[[i+1]]))), c(output_names, 'wave'))
      }
    }
  }
  total_data <- do.call('rbind', out_list)
  if (length(var_list) != 0) total_var <- do.call('rbind', var_list)
  total_data$wave <- factor(total_data$wave)
  output_ranges <- setNames(
    purrr::map(
      output_names,
      ~range(subset(total_data, wave == which_wave)[,.])), output_names)
  if (restrict) {
    targets_grid <- expand.grid(names(targets),
                                names(targets),
                                stringsAsFactors = FALSE)
    targets_grid <- targets_grid[targets_grid$Var1 != targets_grid$Var2,]
    targets_list <- c()
    for (i in seq_len(nrow(targets_grid))) {
      tg <- unlist(targets_grid[i,], use.names = F)
      dat_trunc <- total_data[,tg]
      any_match <- any(apply(dat_trunc, 1, function(x) {
        true1 <- x[1] <= targets[[tg[1]]][2] && x[1] >= targets[[tg[1]]][1]
        true2 <- x[2] <= targets[[tg[2]]][2] && x[2] >- targets[[tg[2]]][1]
        return(true1 && true2)
      }))
      if (!any_match)
        targets_list <- c(targets_list,
                          unlist(targets_grid[i,], use.names = FALSE))
    }
    if(is.null(targets_list))
      warning(paste("Expecting to restrict to failed output pairs but none exist.",
                    "Plotting all outputs."))
    else {
      targets <- targets[names(targets) %in% unique(targets_list)]
      output_names <- names(targets)
    }
  }
  pal <- viridis::viridis(max(length(waves), max(wave_numbers)),
                          option = "D", direction = -1)[wave_numbers+1]
  lfun <- function(data, mapping, targets, zoom = F) { #nocov start
    xname <- rlang::quo_get_expr(mapping$x)
    yname <- rlang::quo_get_expr(mapping$y)
    g <- ggplot(data = data, mapping = mapping) +
      geom_point(cex = p_size) +
      scale_colour_manual(values = pal) +
      theme_bw()
    if (surround) g <- g + geom_point(cex = p_size, pch = 1, colour = 'black')
    g <- g + geom_rect(xmin = targets[[xname]][1],
                       xmax = targets[[xname]][2],
                       ymin = targets[[yname]][1],
                       ymax = targets[[yname]][2],
                       colour = 'red', fill = NA, linewidth = l_wid)
    if (zoom) {
      xrange <- upper_scale*(targets[[xname]][2]-targets[[xname]][1])/2
      yrange <- upper_scale*(targets[[yname]][2]-targets[[yname]][1])/2
      g <- g + coord_cartesian(xlim = c(targets[[xname]][1] - xrange,
                                        targets[[xname]][2] + xrange),
                               ylim = c(targets[[yname]][1] - yrange,
                                        targets[[yname]][2] + yrange)) +
        theme_void()
    }
    else {
      g <- g + coord_cartesian(xlim = output_ranges[[xname]],
                               ylim = output_ranges[[yname]])
    }
    return(g)
  }
  lfun_var <- function(data, mapping, targets, var_data, zoom = F) {
    xname <- rlang::quo_get_expr(mapping$x)
    yname <- rlang::quo_get_expr(mapping$y)
    g <- ggplot(data = data, mapping = mapping) +
      geom_tile(aes(width = 6*var_data[,paste0(yname, 'V')],
                    height = 6*var_data[,paste0(yname, 'V')],
                    fill = data[,'wave']), colour = NA, alpha = 0.4) +
      scale_fill_manual(values = pal) +
      geom_point(cex = p_size/2) +
      scale_colour_manual(values = pal) +
      theme_bw()
    if (surround) g <- g + geom_point(cex = p_size/2, pch = 1, colour = 'black')
    g <- g + geom_rect(xmin = targets[[xname]][1],
                       xmax = targets[[xname]][2],
                       ymin = targets[[yname]][1],
                       ymax = targets[[yname]][2],
                       colour = 'red', fill = NA, linewidth = l_wid)
    if (zoom) {
      xrange <- upper_scale*(targets[[xname]][2]-targets[[xname]][1])/2
      yrange <- upper_scale*(targets[[yname]][2]-targets[[yname]][1])/2
      g <- g + coord_cartesian(xlim = c(targets[[xname]][1] - xrange,
                                        targets[[xname]][2] + xrange),
                               ylim = c(targets[[yname]][1] - yrange,
                                        targets[[yname]][2] + yrange))
    }
    else {
      g <- g + coord_cartesian(xlim = output_ranges[[xname]],
                               ylim = output_ranges[[yname]])
    }
    return(g)
  }
  dfun <- function(data, mapping, targets) {
    xname <- rlang::as_string(rlang::quo_get_expr(mapping$x))
    g <- ggally_densityDiag(data = data, mapping = mapping, alpha = 0.4) +
      geom_vline(xintercept = targets[[xname]], colour = 'red', linewidth = l_wid) +
      scale_fill_manual(values = pal) +
      theme_bw()
    return(g)
  } #nocov end
  if (length(var_list) != 0) {
    g <- ggpairs(total_data, columns = seq_along(output_names),
                 mapping = aes(colour = wave, group = wave),
                 lower = list(continuous = wrap(lfun_var, targets = targets,
                                                var_data = total_var)),
                 diag = list(continuous = wrap(dfun, targets = targets)),
                 upper = list(continuous = wrap(lfun_var, targets = targets,
                                                var_data = total_var, zoom = TRUE)),
                 progress = FALSE,
                 title = "Emulator plots with Targets", legend = 1)
  }
  else {
    g <- ggpairs(total_data, columns = seq_along(output_names),
                 mapping = aes(colour = wave, group = wave),
                 lower = list(continuous = wrap(lfun, targets = targets)),
                 diag = list(continuous = wrap(dfun, targets = targets)),
                 upper = list(continuous = wrap(lfun, targets = targets,
                                                zoom = TRUE)),
                 progress = FALSE,
                 title = "Output plots with targets", legend = 1)
  }
  return(g)
}

#' Multiple Wave Inputs vs Outputs
#'
#' Given multiple waves of points, produce input-output plots for each pair.
#'
#' It can be useful to consider what the dependencies between the input values and output
#' values are, to investigate the suitability of the chosen input ranges (i.e. if widening
#' an input range could result in the targets being matchable). This function provides those
#' plots.
#'
#' For each output-input pair, a points plot is produced with the input value on the x-axis
#' and the output value on the y-axis. The target bounds are superimposed as horizontal lines.
#' The points themselves are coloured by which wave of history matching they came from.
#'
#' These can show dependencies between specific outputs and inputs and, if points are clustering
#' at the far left or right edge of a plot, can give an indication that the input ranges are
#' unsuitable for matching the target.
#'
#' @importFrom GGally ggmatrix
#'
#' @param waves The list of data.frame objects, one for each set of outputs at that wave.
#' @param targets The target values of the outputs.
#' @param output_names The outputs to plot, if not all are wanted.
#' @param input_names The inputs to plot, if not all are wanted.
#' @param p_size Control for the point size on the plots: smaller is better for many plots.
#' @param l_wid Control for line width of superimposed targets.
#' @param normalize If true, plotting is done with target bounds equal size.
#' @param zero_in Is a wave 0 included in the waves list?
#' @param wave_numbers Which waves to plot
#' @param ... Optional parameters (not to be used directly)
#'
#' @return A grid of ggplot objects.
#'
#' @family visualisation tools
#' @export
#'
#' @examples
#'  wave_dependencies(SIRMultiWaveData, SIREmulators$targets, l_wid = 0.8, p_size = 0.8)
#'  wave_dependencies(SIRMultiWaveData, SIREmulators$targets, c('nS', 'nI'), c('aIR', 'aSI'))
#'  \donttest{
#'      # For many plots, it may be helpful to manually modify the font size
#'      wave_dependencies(SIRMultiWaveData, SIREmulators$targets) +
#'       ggplot2::theme(text = ggplot2::element_text(size = 5))
#'  }
wave_dependencies <- function(waves, targets, output_names = names(targets),
                              input_names = names(waves[[1]])[
                                !names(waves[[1]]) %in% names(targets)],
                              p_size = 1.5, l_wid = 1.5, normalize = FALSE,
                              zero_in = TRUE,
                              wave_numbers = ifelse(zero_in, 0, 1):
                                (length(waves)-ifelse(zero_in, 1, 0)), ...) {
  waves <- tryCatch(
    purrr::map(waves, data.frame),
    error = function(e) {
      stop(paste("Cannot convert wave objects to data.frame:", e))
    }
  )
  input_names <- input_names
  for (i in seq_along(targets)) {
    if (!is.atomic(targets[[i]]))
      targets[[i]] <- c(targets[[i]]$val - 3*targets[[i]]$sigma,
                        targets[[i]]$val + 3*targets[[i]]$sigma)
  }
  wave <- NULL
  for (i in ifelse(zero_in, 0, 1):(length(waves)-ifelse(zero_in, 1, 0))) {
    if (i %in% wave_numbers) waves[[i+1]]$wave <- i
  }
  total_data <- do.call('rbind', waves[wave_numbers+1])
  total_data$wave <- factor(total_data$wave)
  pal <- viridis::viridis(max(length(waves), max(wave_numbers)),
                          option = "D", direction = -1)[wave_numbers+1]
  plot_list <- list()
  for (i in seq_along(output_names)) {
    for (j in seq_along(input_names)) {
      plot_list[[(i-1)*length(input_names)+j]] <-
        local({
          i <- i
          j <- j
          line_limits <- targets[[output_names[i]]]
          y_limit <- c(min(line_limits[1], min(total_data[,output_names[i]])),
                       max(line_limits[2], max(total_data[,output_names[i]])))
          g <- ggplot(data = total_data, aes(x = total_data[,input_names[j]],
                                             y = total_data[,output_names[i]],
                                             colour = wave, group = wave)) +
            geom_point(cex = p_size) +
            labs(x = input_names[j], y = output_names[i]) +
            scale_colour_manual(values = pal) +
            xlim(range(total_data[,input_names[j]])) +
            ylim(y_limit) +
            geom_hline(yintercept = line_limits, colour = 'red', linewidth = l_wid) +
            theme_bw() +
            theme(axis.ticks = element_blank())
          if (normalize) {
            yrange <- diff(targets[[output_names[i]]])/2
            g <- g + coord_cartesian(
              ylim = c(targets[[output_names[i]]][1] - yrange,
                       targets[[output_names[i]]][2] + yrange))
          }
          return(suppressWarnings(g))
        })
    }
  }
  return(suppressWarnings(GGally::ggmatrix(plots = plot_list,
                                           ncol = length(input_names),
                                           nrow = length(output_names),
                                           xAxisLabels = input_names,
                                           yAxisLabels = output_names,
                                           progress = FALSE, legend = 1,
                                           title = "Outputs vs Inputs")))
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
#' The output can be plotted in a number of ways: raw; with outputs transformed to log scale;
#' or with targets normalised so that target bounds are all [-1, 1]. These two options may
#' be helpful in visualising behaviour when outputs have vastly different scales, but one
#' still wishes to see them all in the same plot: these options can be toggled by setting
#' \code{logscale = TRUE} or \code{normalize = TRUE} respectively. The data can be grouped in
#' two ways, either colouring by wave of emulation (default) or by the number of targets hit;
#' the latter option is enabled by setting \code{byhit = TRUE}.
#'
#' @import ggplot2
#' @importFrom dplyr mutate
#'
#' @param wave_points The set of wave points, as a list of data.frames
#' @param z The set of target values for each output
#' @param zero_in Is wave zero included? Default: TRUE
#' @param palette If a larger palette is required, it should be supplied here.
#' @param wave_numbers Which waves to plot. If not supplied, all waves are plotted.
#' @param normalize If true, plotting is done with rescaled target bounds.
#' @param logscale If true, targets are log-scaled before plotting.
#' @param byhit Should runs be grouped by number of targets hit, rather than wave?
#' @param barcol The colour of the target error bars/bounds
#' @param ... Optional parameters (not to be used directly)
#'
#' @return A ggplot object.
#'
#' @family visualisation tools
#' @export
#'
#' @examples
#'  simulator_plot(SIRMultiWaveData, SIREmulators$targets)
#'  simulator_plot(SIRMultiWaveData[2:4], SIREmulators$targets,
#'   zero_in = FALSE, wave_numbers = c(1,3))
#'  simulator_plot(SIRMultiWaveData, SIREmulators$targets, byhit = TRUE)
#'
simulator_plot <- function(wave_points, z, zero_in = TRUE, palette = NULL,
                           wave_numbers = seq(
                             ifelse(zero_in, 0, 1),
                             length(wave_points)-ifelse(zero_in, 1, 0)),
                           normalize = FALSE, logscale = FALSE, byhit = FALSE,
                           barcol = "#444444", ...) {
  wave_points <- tryCatch(
    purrr::map(wave_points, data.frame),
    error = function(e) {
      stop(paste("Cannot convert wave objects to data.frame:", e))
    }
  )
  if (normalize && logscale) {
    warning("Both normalize and logscale = TRUE; defaulting to logscale.")
    normalize <- FALSE
  }
  for (i in seq_along(z)) {
    if (!is.atomic(z[[i]]))
      z[[i]] <- c(z[[i]]$val - 3*z[[i]]$sigma, z[[i]]$val + 3*z[[i]]$sigma)
  }
  name <- value <- run <- wave <- val <- sigma <- NULL
  output_names <- names(z)
  target_hits <- function(result, targets, sum_func = "sum") {
    hits <- purrr::map_lgl(names(targets), function(t) {
        return(result[t] <= targets[[t]][2] && result[t] >= targets[[t]][1])
    })
    sum_function <- get(sum_func)
    return(sum_function(hits))
  }
  sim_runs <- do.call(
    'rbind',
    purrr::map(
      wave_numbers,
      ~data.frame(wave_points[[.+ifelse(zero_in, 1, 0)]][,output_names, drop = FALSE],
                  wave = .)))
  if (byhit) {
    t_hits <- apply(sim_runs, 1, target_hits, z)
    waves_by_hits <- purrr::map(0:length(z), ~sim_runs[t_hits == .,])
    w_numbers <- (0:length(z))[purrr::map_lgl(waves_by_hits, ~nrow(.) > 0)]
    return(simulator_plot(waves_by_hits, z, zero_in = TRUE, palette = palette,
                          wave_numbers = w_numbers, normalize = normalize,
                          logscale = logscale, byhit = FALSE, barcol = barcol,
                          change_legend = TRUE, ...))
  }
  sim_runs$run <- seq_along(sim_runs[,1])
  if (normalize) {
    for (i in names(z)) {
      sim_runs[[i]] <- (2*sim_runs[[i]] - z[[i]][1]-z[[i]][2])/(diff(z[[i]]))
      z[[i]] <- c(-1,1)
    }
  }
  if (logscale) {
    for (i in names(z)) {
      sim_runs[[i]] <- log(sim_runs[[i]])
      if(z[[i]][1] <= 0) z[[i]][1] <- 1e-4
      if(z[[i]][2] <= 0) z[[i]][2] <- 1e-4
      z[[i]] <- log(z[[i]])
    }
  }
  pivoted <- pivot_longer(sim_runs, cols = !c('run', 'wave'))
  pivoted$wave <- as.factor(pivoted$wave)
  pivoted$name <- factor(pivoted$name, levels = names(z))
  if (is.null(palette))
    pal <- viridis::viridis(length(wave_points),
                            option = 'plasma', direction = -1)
  else pal <- palette
  pal <- pal[seq_along(pal) %in% (wave_numbers+ifelse(zero_in, 1, 0))]
  obs <- data.frame(name = names(z), min = purrr::map_dbl(z, ~.[1]),
                    max = purrr::map_dbl(z, ~.[2]))
  optionals <- list(...)
  if (!is.null(optionals[['change_legend']]) && optionals[['change_legend']])
    legend_title = "# Hits"
  else
    legend_title = "Wave"
  if(length(output_names) == 1) {
      obs <- uncount(obs, length(unique(pivoted$wave))) %>%
        mutate(wave = factor(0:(length(unique(pivoted$wave)) - 1)))
      g <- ggplot(data = pivoted, aes(x = wave, y = value)) +
        ggbeeswarm::geom_beeswarm(aes(colour = wave, group = wave)) +
        scale_colour_manual(values = pal, name = legend_title) +
        geom_point(data = obs, aes(x = wave, y = (min+max)/2)) +
        geom_errorbar(data = obs,
                      aes(y = (min+max)/2, ymax = max, ymin = min),
                      width = 0.1, linewidth = 1.25, colour = barcol) +
        labs(title = paste0("Simulator evaluations (", output_names, ") at wave points",
                            (if (normalize) ": normalised" else (
                              if (logscale) ": log-scale" else "")))) +
        theme_minimal()
  } else {
      g <- ggplot(data = pivoted, aes(x = name, y = value)) +
        geom_line(aes(group = run, colour = wave)) +
        scale_colour_manual(values = pal, name = legend_title) +
        geom_point(data = obs, aes(x = name, y = (min+max)/2), colour = barcol) +
        geom_errorbar(data = obs,
                      aes(y = (min+max)/2, ymax = max, ymin = min),
                      width = 0.1, linewidth = 1.25, colour = barcol) +
        labs(title = paste0("Simulator evaluations at wave points",
                            (if (normalize) ": normalised" else (
                              if (logscale) ": log-scale" else "")))) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  if (normalize) g <- g + coord_cartesian(ylim = c(-3, 3))
  return(suppressWarnings(g))
}

#' Diagnostic plots for wave outputs
#'
#' A wrapper function for the set of diagnostic plots for multiple waves.
#'
#' The functions \code{\link{simulator_plot}}, \code{\link{wave_points}}, \code{\link{wave_points}},
#' and \code{\link{wave_dependencies}} are called, one after the other, to allow diagnosis of waves
#' of emulation.
#'
#' The \code{directory} option should be used as follows. If the desired location is in fact
#' a folder, it should end in "/"; if instead the structure requires each plot to be saved with a
#' prefix, then it should be provided. For example, \code{directory = "Plots/"} in the first event
#' or \code{directory = "Plots/unique-identifier"} in the second event.
#'
#' @importFrom grDevices dev.off png
#'
#' @param waves The wave points, as a list of data.frames.
#' @param targets The output targets.
#' @param output_names The outputs to plot.
#' @param input_names The inputs to plot.
#' @param directory The location of files to be saved (if required).
#' @param s.heights The heights of the saved pngs (if directory is not NULL).
#' @param s.widths The widths of the saved pngs (if directory is not NULL).
#' @param include.norm Should normalized versions of simulator_plot and wave_dependencies be made?
#' @param include.log Should the log-scale version of simulator_plot be made?
#' @param ... Optional parameters (eg \code{p_size}, \code{l_wid}, ...)
#'
#' @return The set of plots (either into console or saved).
#'
#' @family visualisation tools
#' @export
#'
#' @examples
#' \donttest{
#'  diagnostic_wrap(SIRMultiWaveData, SIREmulators$targets)
#'  diagnostic_wrap(SIRMultiWaveData, SIREmulators$targets,
#'   input_names = c('aSI', 'aIR'), output_names = c('nI', 'nR'),
#'   p_size = 0.8, l_wid = 0.8, wave_numbers = 1:3, zero_in = FALSE, surround = TRUE)
#'   }
diagnostic_wrap <- function(waves, targets, output_names = names(targets),
                            input_names = names(waves[[1]])[
                              !names(waves[[1]]) %in% names(targets)],
                            directory = NULL, s.heights = rep(1000, 4),
                            s.widths = s.heights, include.norm = TRUE,
                            include.log = TRUE, ...) {
  if (!is.null(directory) && !file.exists(sub("(.*)/[^/]*$", "\\1", directory))) #nocov
    stop("Specified directory does not exist. No plots saved.") #nocov
  s.widths[1] <- ceiling(length(output_names)/10)*1000
  get_dims <- function(initial, plots) { #nocov start
    if (length(initial) == 1)
      return(rep(initial, length(plots)))
    if (length(initial) < 4)
      return(rep(initial[1], length(plots)))
    if (length(initial) > 4 && length(initial) != length(plots))
      initial <- initial[1:4]
    if (length(initial) == 4 && length(plots) > 4) {
      replaced <- rep(0, length(plots))
      replaced[grepl("simulatorplot.*", names(plots))] <- initial[1]
      replaced[names(plots) == "posteriorplot"] <- initial[2]
      replaced[names(plots) == "outputsplot"] <- initial[3]
      replaced[grepl("dependencyplot.*", names(plots))] <- initial[4]
      return(replaced)
    }
  } #nocov end
  g <- list()
  g[["simulatorplot"]] <- simulator_plot(waves, targets, ...)
  if (include.norm)
    g[["simulatorplotnorm"]] <- simulator_plot(waves, targets, normalize = TRUE, ...)
  if (include.log)
    g[["simulatorplotlog"]] <- simulator_plot(waves, targets, logscale = TRUE, ...)
  g[["posteriorplot"]] <- wave_points(waves, input_names, ...)
  g[["outputsplot"]] <- wave_values(waves, targets, output_names, ...)
  g[["dependencyplot"]] <- wave_dependencies(waves, targets, output_names, input_names, ...)
  if (include.norm)
    g[["dependencyplotnorm"]] <- wave_dependencies(waves, targets, output_names,
                            input_names, normalize = TRUE, ...)
  if (!is.null(directory)) { #nocov start
    s.widths <- get_dims(s.widths, g)
    s.heights <- get_dims(s.heights, g)
    for (i in seq_along(g)) {
      png(filename = paste0(directory, names(g)[[i]], ".png"),
          width = s.widths[i], height = s.heights[i])
      print(g[[i]])
      dev.off()
    }
    return(NULL)
  } #nocov end
  return(g)
}

#' Output Hit Summary
#'
#' Provides a summary of numbers of points that hit n outputs
#'
#' Given a collection of wave points and the targets used in history matching,
#' it might be informative to consider the proportion of points whose model
#' output matches a given number of targets. This function provides by-wave
#' information about how many parameter sets are matches to 0,1,2,...,n outputs.
#'
#' The results of the analysis can be presented as a \code{data.frame} object
#' where each row is a wave and each column a number of outputs; if \code{plt = TRUE}
#' the results are instead presented visually, as a grid coloured by proportion of
#' total points (if \code{grid.plot = TRUE}, the default) or as a series of discrete
#' density lines, one per wave. The \code{as.per} argument determines whether the
#' output values are raw or if they are calculated percentages of the total number
#' of parameter sets for a given wave.
#'
#' When the data arise from a stochastic model, and therefore parameter sets have
#' multiple realisations, there are multiple ways to analyze the data (determined
#' by \code{measure}). The options are "mean" to compare the means of realisations
#' to the outputs; "real" to compare all individual realisations; and "stoch" to
#' consider an output matched to if the mean lies within 3 standard deviations of
#' the output, where the standard deviation is calculated over the realisations.
#'
#' @param waves The collection of waves, as a list of data.frames
#' @param targets The output targets
#' @param input_names The names of the input parameters
#' @param measure If stochastic, the measure to use to compare (see description)
#' @param plt If TRUE, results are plotted; else a data.frame is returned
#' @param as.per Should the data be percentages, or raw numbers?
#' @param grid.plot If \code{plt = TRUE}, determines the type of plot.
#'
#' @return Either a data.frame of results or a ggplot object plot
#'
#' @family visualisation tools
#' @export
#'
#' @examples
#'  # Default Usage
#'  hit_by_wave(SIRMultiWaveData, SIREmulators$targets, c('aSI', 'aIR', 'aSR'))
#'  # Plotting - line plot or raw figures
#'  hit_by_wave(SIRMultiWaveData, SIREmulators$targets, c('aSI', 'aIR', 'aSR'),
#'   plt = TRUE, as.per = FALSE, grid.plot = FALSE)
hit_by_wave <- function(waves, targets, input_names, measure = "mean",
                        plt = FALSE, as.per = TRUE, grid.plot = TRUE) {
  wave <- name <- value <- value_unsc <- NULL
  target_hits <- function(result, targets, sum_func = "sum") {
    hits <- purrr::map_lgl(names(targets), function(t) {
      if (is.atomic(targets[[t]]))
        return(result[t] <= targets[[t]][2] && result[t] >= targets[[t]][1])
      result[t] <= targets[[t]]$val + 3*targets[[t]]$sigma && result[t] >= targets[[t]]$val - 3*targets[[t]]$sigma
    })
    sum_function <- get(sum_func)
    return(sum_function(hits))
  }
  target_overlaps <- function(means, sds, reps, targets, sum_func = "sum", sd = 3) {
    check_overlap <- function(int1, int2) {
      if (int1[2] >= int2[1] && int1[2] <= int2[1]) return(TRUE)
      if (int1[1] <= int2[2] && int1[1] >= int2[1]) return(TRUE)
      if (int2[2] >= int1[1] && int2[2] <= int1[1]) return(TRUE)
      if (int2[1] <= int1[2] && int2[1] >= int1[1]) return(TRUE)
      return(FALSE)
    }
    intervals <- purrr::map(names(targets), ~as.numeric(c(means[.] - sd/sqrt(reps)*sds[.], means[.] + sd/sqrt(reps)*sds[.]))) |> setNames(names(targets))
    hits <- purrr::map_lgl(names(targets), function(t) {
      if (!is.atomic(targets[[t]]))
        targ_int <- c(targets[[t]]$val - 3*targets[[t]]$sigma, targets[[t]]$val + 3*targets[[t]]$sigma)
      else
        targ_int <- targets[[t]]
      return(check_overlap(targ_int, intervals[[t]]))
    })
    sum_function <- get(sum_func)
    return(sum_function(hits))
  }
  waves <- purrr::map(waves, ~.[,c(input_names, names(targets))])
  wave_uids <- purrr::map(waves, function(w) {
    apply(w[,input_names], 1, rlang::hash)
  })
  duplicated <- FALSE
  if (any(purrr::map_lgl(wave_uids, ~length(unique(.)) != length(.)))) duplicated <- TRUE
  if (duplicated == TRUE && measure != "real") {
    waves_grouped <- purrr::map(waves, ~. |> dplyr::group_by(across(all_of(input_names))))
    wave_means <- purrr::map(waves_grouped, ~. |> dplyr::summarise(.groups = "keep", across(everything(), mean)))
    if (measure == "stoch") {
      wave_reps <- purrr::map(waves_grouped, ~purrr::map_dbl(dplyr::group_rows(.), length))
      wave_sds <- purrr::map(waves_grouped, ~. |> dplyr::summarise(.groups = "keep", across(everything(), sd)))
    }
    else {
      wave_reps <- NULL
      wave_sds <- NULL
    }
    if (is.null(wave_reps)) {
      hit_data <- purrr::map(wave_means, function(w) {
        apply(w, 1, target_hits, targets)
      })
    }
    else {
      hit_data <- purrr::map(seq_along(wave_means), function(i) {
        purrr::map_dbl(seq_len(nrow(wave_means[[i]])), function(j) {
          target_overlaps(wave_means[[i]][j,], wave_sds[[i]][j,], wave_reps[[i]][j], targets)
        })
      })
    }
  }
  else {
    hit_data <- purrr::map(waves, function(w) {
      apply(w, 1, target_hits, targets)
    })
  }
  hit_data_df <- do.call('rbind.data.frame',
                         purrr::map(hit_data, function(h) purrr::map_dbl(0:length(targets), ~sum(h == .)))
  ) |> setNames(0:length(targets))
  if (as.per) hit_data_df <- sweep(hit_data_df, 1, apply(hit_data_df, 1, sum), "/")
  if (!plt) {
    row.names(hit_data_df) <- 0:(length(waves)-1)
    if (as.per) return(signif(hit_data_df, 2)*100)
    return(hit_data_df)
  }
  hit_data_df$wave <- seq_len(nrow(hit_data_df))-1
  plot.mat <- tidyr::pivot_longer(hit_data_df, cols = !wave)
  plot.mat$wave <- factor(plot.mat$wave, levels = rev(seq_along(hit_data_df)-1))
  plot.mat$name <- factor(plot.mat$name, levels = 0:length(targets))
  if (!as.per) {
    hit_data_df_max <- max(purrr::map_dbl(hit_data, length))
    hit_data_df_scale <- data.frame(t(apply(hit_data_df[,!names(hit_data_df) == 'wave'], 1, function(x) ceiling(x*hit_data_df_max/sum(x))))) |>
      setNames(0:length(targets))
    hit_data_df_scale$wave <- hit_data_df$wave
    scale_additional <- tidyr::pivot_longer(hit_data_df_scale, cols = !wave)
    holdout <- plot.mat$value
    plot.mat$value <- scale_additional$value
    plot.mat$value_unsc <- holdout
  }
  labr <- function(b) {
    if (as.per) as.numeric(b)*100
    else as.numeric(b)
  }
  if (duplicated) {
    if (measure == "mean") subtit <- "Mean of realisations"
    else if (measure == "real") subtit <- "Realisations"
    else if (measure == "stoch") subtit <- "Means with Stochasticity"
    else subtit <- ""
  }
  else subtit <- ""
  if (grid.plot)
    g <- ggplot(data = plot.mat, aes(x = name, y = wave)) +
    geom_tile(aes(fill = value)) +
    geom_text(data = subset(plot.mat, value != 0), aes(label = if (as.per) signif(value, 2)*100 else value_unsc), size = 2.5) +
    scale_fill_gradientn(
      colours = c("white", "yellow", "#77FF00", "green"),
      values = c(0, 1e-8, 0.2, 1),
      breaks = signif(seq(0, 1, length.out = 20), 2),
      name = ifelse(as.per, "%", "#"), labels = as_labeller(labr))
  else
    g <- ggplot(data = plot.mat, aes(x = name, y = value, group = wave, colour = wave)) +
    geom_line() +
    viridis::scale_colour_viridis(discrete = TRUE, name = "Wave")
  g <- g + theme_minimal() +
    labs(title = paste(ifelse(as.per, "Percentage", "Number"), "of Points Hitting # Outputs, by Wave"),
         subtitle = subtit, x = "# Outputs", y = ifelse(grid.plot, "Wave", ifelse(as.per, "Percentage", "Number")))
  if (grid.plot) g <- g + theme(legend.position = 'none')
  return(g)
}
