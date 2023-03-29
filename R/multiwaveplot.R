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
  for (i in ifelse(zero_in, 0, 1):(length(waves)-ifelse(zero_in, 1, 0))) {
    if (i %in% wave_numbers)
      out_list[[i+1]] <- setNames(
        cbind(waves[[i+1]][,input_names],
              rep(i, nrow(waves[[i+1]]))), c(input_names, 'wave'))
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
                          option = "D", direction = -1)[wave_numbers+1]
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
#'
simulator_plot <- function(wave_points, z, zero_in = TRUE, palette = NULL,
                           wave_numbers = seq(
                             ifelse(zero_in, 0, 1),
                             length(wave_points)-ifelse(zero_in, 1, 0)),
                           normalize = FALSE, logscale = FALSE,
                           barcol = "#444444", ...) {
  for (i in seq_along(z)) {
    if (!is.atomic(z[[i]]))
      z[[i]] <- c(z[[i]]$val - 3*z[[i]]$sigma, z[[i]]$val + 3*z[[i]]$sigma)
  }
  name <- value <- run <- wave <- val <- sigma <- NULL
  output_names <- names(z)
  sim_runs <- do.call(
    'rbind',
    purrr::map(
      wave_numbers,
      ~data.frame(wave_points[[.+ifelse(zero_in, 1, 0)]][,output_names, drop = FALSE],
                  wave = .)))
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
  if(length(output_names) == 1) {
      obs <- uncount(obs, length(unique(pivoted$wave))) %>%
        mutate(wave = factor(0:(length(unique(pivoted$wave)) - 1)))
      g <- ggplot(data = pivoted, aes(x = wave, y = value)) +
        ggbeeswarm::geom_beeswarm(aes(colour = wave, group = wave)) +
        scale_colour_manual(values = pal) +
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
        scale_colour_manual(values = pal) +
        geom_point(data = obs, aes(x = name, y = (min+max)/2)) +
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
