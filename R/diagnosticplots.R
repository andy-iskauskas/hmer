#' Output Plotting
#'
#' A simple diagnostic plot that compares the output values to input values, for
#' each possible combination. If emulators are provided, the emulator predictions
#' are plotted; otherwise the model outputs are plotted.
#'
#' If emulators are provided, then the \code{points} argument is optional: if given
#' then the emulator predictions will correspond to those at the points provided. If
#' no points are provided, 100*d (where d is the number of input parameters) are sampled
#' uniformly from the space and used to predict at.
#'
#' If no emulators are provided, then points must be provided, along with the names of
#' the outputs to plot; each named output must exist as a column in the points data.frame.
#'
#' @importFrom graphics plot par
#' @importFrom stats setNames
#'
#' @param ems A set of \code{\link{Emulator}} objects.
#' @param points A set of points at which to evaluate the emulator expectation
#' @param model If TRUE, use the model outputs; else use emulator expectation
#' @param out_names If no emulators are provided, use this argument to indicate outputs.
#' @param targets If targets are provided, these are added into the plots.
#'
#' @return The dependency plots.
#'
#' @family visualisation tools
#' @export
#'
#' @examples
#'  behaviour_plot(SIREmulators$ems, model = FALSE)
#'  behaviour_plot(points = SIRSample$training, out_names = names(SIREmulators$ems))
#'  #> Throws a warning
#'  behaviour_plot(SIRMultiWaveEmulators, model = TRUE, targets = SIREmulators$targets)
behaviour_plot <- function(ems, points, model = missing(ems), out_names = unique(names(collect_emulators(ems))), targets = NULL) {
  if (missing(ems) && missing(points))
    stop("One of 'ems' or 'points' must be supplied.")
  if (missing(points) && model) {
    warning("Cannot do model output plots with no data points. Using emulator expectation.")
    model <- FALSE
  }
  if (missing(ems) && !model)
    stop("Cannot perform emulator expectation comparison with no emulators (should model = TRUE?).")
  tryCatch(
    out_names,
    error = function(e) {
      stop("No output names provided and no emulators to determine them from.")
    }
  )
  if (!missing(ems))
    ems <- collect_emulators(ems)
  if (missing(points)) {
    ranges <- ems[[1]]$ranges
    points <- setNames(do.call('cbind.data.frame', purrr::map(ranges, ~runif(length(ranges)*100, .[1], .[2]))), names(ranges))
  }
  if (model && !all(out_names %in% names(points))) {
    warning("Not all outputs have a column in the data given by 'points'. Using emulator expectation.")
    model <- FALSE
  }
  if (!missing(ems))
    in_names <- names(ems[[1]]$ranges)
  else
    in_names <- names(points)[!names(points) %in% out_names]
  if (model) {
    data <- points[,c(in_names, out_names)]
  }
  else {
    data <- setNames(cbind.data.frame(points, purrr::map(ems[!duplicated(purrr::map_chr(ems, ~.$output_name))], ~.$get_exp(points))), c(in_names, out_names))
  }
  if (!is.null(targets))
    ylims <- setNames(purrr::map(out_names, function(i) {
      if (is.numeric(targets[[i]])) {
        ymin <- min(data[,i], targets[[i]][1]-0.05*diff(targets[[i]]))
        ymax <- max(data[,i], targets[[i]][2]+0.05*diff(targets[[i]]))
        return(c(ymin, ymax))
      }
      else {
        ymin <- min(data[,i], targets[[i]]$val - 3.5*targets[[i]]$sigma)
        ymax <- max(data[,i], targets[[i]]$val + 3.5*targets[[i]]$sigma)
        return(c(ymin, ymax))
      }
    }), out_names)
  else
    ylims <- setNames(purrr::map(out_names, ~c(min(data[,.]), max(data[,.]))), out_names)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  for (i in seq_along(out_names)) {
    op <- par(mfrow = c(min(4, ceiling(sqrt(length(in_names)))), min(4, ceiling(sqrt(length(in_names))))))
    for (j in seq_along(in_names)) {
      plot(data[,in_names[j]], data[,out_names[i]], pch = 16, xlab = in_names[[j]], ylab = out_names[[i]],
           ylim = ylims[[i]])
      if (!is.null(targets[[out_names[i]]])) {
        if(is.numeric(targets[[out_names[i]]]))
          abline(h = targets[[out_names[i]]], lty = 2)
        else
          abline(h = c(targets[[i]]$val-3*targets[[out_names[i]]]$sigma, targets[[out_names[i]]]$val+3*targets[[out_names[i]]]$sigma))
      }
    }
    par(op)
  }
}

#' Space Removal Diagnostics
#'
#' Finds the proportion of space removed as a function of implausibility cut-off and of one of
#' structural discrepancy, emulator variance, or correlation hyperparameter(s).
#'
#' The reduction in space is found by evaluating a p^d regular grid, where p is chosen by
#' \code{ppd} and d is the dimension of the input space. Larger values of p will give a more
#' accurate reflection of the space removed, at a corresponding computational cost. For the
#' purpose of quick-and-dirty diagnostics, \code{ppd = 5} is sufficient: the default is 10.
#'
#' The parameter \code{modified} can be one of three strings: \code{'obs'} corresponding
#' to observation uncertainty; \code{'disc'} corresponding to internal and external
#' discrepancy (as given in \code{Emulator$disc}); \code{'var'} corresponding to global
#' emulator variance (as given by \code{Emulator$u_sigma}), and \code{'hp'} corresponding to
#' the hyperparameters of the emulator correlation structure. In the first case, the
#' implausibilities are recalculated for each inflation value; in the other two cases the
#' emulators are retrained. For this reason, the \code{'var'} and \code{'hp'} options are
#' computationally more intensive. The default is \code{'obs'}.
#'
#' The inflationary/deflationary values are chosen by \code{u_mod}: the default is to take
#' 80\%, 90\%, 100\%, 110\%, and 120\% of the original value as the variation. The proportion of
#' points deemed non-implausible is checked at a set of implausibility cutoffs defined by
#' \code{intervals}, and a plot is returned showing the relevant data.
#'
#' @import ggplot2
#' @import tidyr
#' @importFrom stats setNames
#' @importFrom viridis scale_colour_viridis
#'
#' @param ems The \code{\link{Emulator}} objects.
#' @param targets The corresponding targets to match to.
#' @param ppd The number of points per input dimension to sample at.
#' @param u_mod The proportional values by which to inflate/deflate the relevant statistic.
#' @param intervals The interval values of the implausibility cutoff at which to evaluate.
#' @param modified The statistic to modify: obs, disc, var or hp (see above)
#' @param maxpoints The maximum number of points to evaluate at
#'
#' @return A ggplot object
#'
#' @family visualisation tools
#' @export
#'
#' @seealso \code{\link{space_removal}} for a numeric representation of space removed.
#'
#' @examples
#'  space_removed(SIREmulators$ems, SIREmulators$targets, ppd = 5)
#'  space_removed(SIREmulators$ems$nS, SIREmulators$targets,
#'   ppd = 5, u_mod = seq(0.75, 1.25, by = 0.25), intervals = seq(2, 6, by = 0.1))
space_removed <- function(ems, targets, ppd = 10, u_mod = seq(0.8, 1.2, by = 0.1), intervals = seq(0, 10, length.out = 200), modified = 'obs', maxpoints = 50000) {
  value <- name <- cutoff <- NULL
  if ("Emulator" %in% class(ems))
    ems <- setNames(list(ems), ems$output_name)
  if ("EmProto" %in% class(ems[[1]])) {
    if (modified == 'var' || modified == 'hp') {
      warning("Cannot consider space removed with respect to prior variance or hyperparameter for proto ems. Setting to observation error (obs).")
      modified <- 'obs'
    }
  }
  ranges <- ems[[1]]$ranges
  if (is.null(maxpoints)) maxpoints <- 100*length(ranges)
  if (ppd^length(ranges) > maxpoints) {
    ptgrid <- setNames(do.call('cbind.data.frame', purrr::map(ranges, ~runif(maxpoints, .[[1]], .[[2]]))), names(ranges))
  }
  else {
    ptgrid <- setNames(expand.grid(purrr::map(ranges, ~seq(.[[1]], .[[2]], length.out = ppd))), names(ranges))
  }
  imp_array <- array(0, dim = c(length(intervals), length(u_mod)))
  if (!modified %in% c('obs', 'disc', 'var', 'hp')) {
    warning("Unrecognised vary parameter. Setting to observation error (obs).")
    modified <- 'obs'
  }
  if (modified == "disc" && all(unlist(purrr::map(ems, ~all(.$discs == 0))))) {
    warning("'disc' chosen, but no emulators have any internal or external discrepancy. Setting to observation error (obs).")
    modified <- 'obs'
  }
  if (modified == 'obs') {
    for (i in u_mod) {
      targets_2 <- targets
      for (j in seq_along(targets)) {
        if (is.atomic(targets[[j]])) {
          targets_2[[j]] <- c(mean(targets[[j]]) - i * diff(targets[[j]])/2, mean(targets[[j]]) + i * diff(targets[[j]])/2)
        }
        else {
          targets_2[[j]] <- list(val = targets[[j]]$val, sigma = i * targets[[j]]$sigma)
        }
      }
      m_imps <- nth_implausible(ems, ptgrid, targets_2)
      coff <- purrr::map_dbl(intervals, ~1-length(m_imps[m_imps <= .])/length(m_imps))
      imp_array[, match(i, u_mod)] <- coff
    }
  }
  else {
    for (i in u_mod) {
      if (modified == 'var')
        ems_2 <- purrr::map(ems, ~.$mult_sigma(sqrt(i)))
      else if (modified == 'hp')
        ems_2 <- purrr::map(ems, ~.$set_hyperparams(purrr::map(.$corr$hyper_p, ~i*.)))
      else {
        ems_2 <- ems
        for (j in seq_along(ems_2)) {
          ems_2[[j]]$disc <- lapply(ems[[j]]$disc, function(a) a * i)
        }
      }
      m_imps <- nth_implausible(ems_2, ptgrid, targets)
      coff <- purrr::map_dbl(intervals, ~1-length(m_imps[m_imps <= .])/length(m_imps))
      imp_array[, match(i, u_mod)] <- coff
    }
  }
  df <- setNames(data.frame(imp_array), u_mod)
  df$cutoff <- intervals
  df_pivot <- pivot_longer(df, cols = !c('cutoff'))
  title <- switch(modified, 'obs' = "observational error", 'disc' = 'structural discrepancy', 'var' = 'variance inflation', 'hp' = 'hyperparameter inflation')
  subtitle <- switch(modified, 'obs' = "% Observational\nError", 'disc' = '% Structural\nDiscrepancy', 'var' = '% Variance', 'hp' = '% Hyperparameter')
  g <- ggplot(data = df_pivot, aes(x = cutoff, y = value, group = name, colour = name)) +
    geom_line(linewidth = 1.5) +
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
#' @family visualisation tools
#' @export
#'
#' @examples
#'  validation_pairs(SIREmulators$ems, SIRSample$validation, SIREmulators$targets)
#'  wider_ranges <- purrr::map(SIREmulators$ems[[1]]$ranges, ~.*c(-2, 2))
#'  validation_pairs(SIREmulators$ems, SIRSample$validation,
#'   SIREmulators$targets, ranges = wider_ranges, cb = TRUE)
validation_pairs <- function(ems, points, targets, ranges, nth = 1, cb = FALSE) {
  if ("Emulator" %in% class(ems))
    ems <- setNames(list(ems), ems$output_name)
  if (missing(ranges)) ranges <- ems[[1]]$ranges
  em_exp <- data.frame(purrr::map(ems, ~.$get_exp(points)))
  em_var <- data.frame(purrr::map(ems, ~.$get_cov(points) + .$disc$internal^2 + .$disc$external^2))
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
  g <- ggpairs(diag_vals, columns = seq_along(ranges), aes(colour = diag_vals['d']), legend = c(1,2),
               title = "Emulator Diagnostics (lower) and Emulator Implausibility (upper)",
               lower = list(continuous = wrap(plotfun), mapping = aes(colour = diag_vals[,'d'])),
               upper = list(continuous = wrap(plotfun), mapping = aes(colour = diag_vals[,'imp'])),
               diag = 'blank', progress = FALSE) +
    scale_colour_gradient2(low = cols$low, mid = cols$mid, high = cols$high, midpoint = 3, breaks = colour_breaks, name = "Scale", labels = colour_names) +
    theme(legend.position = 'right') +
    theme_minimal()
  return(g)
}

#' Find Effect Strength of Active Variables
#'
#' Collates the linear and quadratic contributions of the active variables to the global
#' emulators' behaviour
#'
#' For a set of emulators, it can be useful to see the relative contributions of various
#' parameters to the global part of the emulator (i.e. the regression surface). This
#' function extracts the relevant information from a list of emulator objects.
#'
#' The parameter \code{quadratic} controls whether quadratic effect strength is
#' calculated and plotted (an unnecessary plot if, say, linear emulators have been trained).
#' The remaining options control visual aspects of the plots: \code{line.plot} determines
#' whether a line or bar (default) plot should be produced, \code{grid.plot} determines
#' whether the results are plotted as a graph or a grid, and \code{labels} determines
#' if a legend should be provided with the plot (for large numbers of emulators, it is
#' advisable to set this to \code{FALSE}).
#'
#' @import ggplot2
#'
#' @param ems The Emulator object(s) to be analysed.
#' @param plt Should the results be plotted?
#' @param line.plot Should a line plot be produced?
#' @param grid.plot Should the effect strengths be plotted as a grid?
#' @param labels Whether or not the legend should be included.
#' @param quadratic Whether or not quadratic effect strength should be calculated.
#' @param xvar Should the inputs be used on the x-axis?
#'
#' @return A list of data.frames: the first is the linear strength, and the second quadratic.
#'
#' @family visualisation tools
#' @export
#'
#' @examples
#'  effect <- effect_strength(SIREmulators$ems)
#'  effect_line <- effect_strength(SIREmulators$ems, line.plot = TRUE)
#'  effect_grid <- effect_strength(SIREmulators$ems, grid.plot = TRUE)
effect_strength <- function(ems, plt = interactive(), line.plot = FALSE,
                            grid.plot = FALSE, labels = TRUE, quadratic = TRUE,
                            xvar = TRUE) {
  if ("EmProto" %in% unlist(purrr::map(ems, class), use.names = FALSE))
    stop("effect_strength not applicable for Proto_emulator objects.")
  get_effect_strength <- function(em, quad = FALSE) {
    es <- c()
    for (i in seq_along(em$ranges)) {
      plus_vec <- c(rep(0, i-1), 1, rep(0, length(em$ranges)-i))
      minus_vec <- c(rep(0, i-1), -1, rep(0, length(em$ranges)-i))
      zero_vec <- rep(0, length(em$ranges))
      p_f <- purrr::map_dbl(em$basis_f, purrr::exec, plus_vec) %*% em$beta_mu
      m_f <- purrr::map_dbl(em$basis_f, purrr::exec, minus_vec) %*% em$beta_mu
      z_f <- purrr::map_dbl(em$basis_f, purrr::exec, zero_vec) %*% em$beta_mu
      if (quad) es <- c(es, round((p_f+m_f-2*z_f)/2, 6))
      else es <- c(es, round((p_f - m_f)/2, 6))
    }
    return(es)
  }
  if ("Emulator" %in% class(ems)) ems <- setNames(list(ems), ems$output_name)
  ranges <- ems[[1]]$ranges
  linear_effect_strength <- setNames(
    data.frame(
      do.call('rbind', purrr::map(ems, get_effect_strength))), names(ranges))
  if (quadratic)
  {
    quadratic_effect_strength <- setNames(
      data.frame(
        do.call('rbind',
                purrr::map(ems, get_effect_strength, TRUE))), names(ranges))
    complete_set <- list(linear = linear_effect_strength,
                         quadratic = quadratic_effect_strength)
    quad.mat <- pivot_longer(quadratic_effect_strength,
                             cols = everything(), names_to = "Var2")
    quad.mat$Var1 <- rep(row.names(quadratic_effect_strength),
                         each = length(quadratic_effect_strength))
    quad.mat$Var1 <- factor(quad.mat$Var1,
                            levels = purrr::map_chr(ems, ~.$output_name))
    quad.mat$Var2 <- factor(quad.mat$Var2, levels = names(ranges))
  }
  else complete_set <- linear_effect_strength
  lin.mat <- pivot_longer(linear_effect_strength,
                          cols = everything(), names_to = "Var2")
  lin.mat$Var1 <- rep(row.names(linear_effect_strength),
                      each = length(linear_effect_strength))
  lin.mat$Var1 <- factor(lin.mat$Var1,
                         levels = purrr::map_chr(ems, ~.$output_name))
  lin.mat$Var2 <- factor(lin.mat$Var2, levels = names(ranges))
  Var1 <- Var2 <- value <- NULL
  if (plt) { #nocov start
    if (grid.plot) {
      g <- ggplot(data = lin.mat, aes(x = Var2, y = Var1, fill = value)) +
        geom_tile(colour = 'black') +
        scale_fill_gradient2(low = 'red', mid = 'white', high = 'blue',
                             midpoint = 0, name = "Strength") +
        labs(title = "Linear Effect Strength", x = "Parameter", y = "Output")
      print(g)
      if (quadratic) {
        g <- ggplot(data = quad.mat, aes(x = Var2, y = Var1, fill = value)) +
          geom_tile(colour = 'black') +
          scale_fill_gradient2(low = 'red', mid = 'white', high = 'blue',
                               midpoint = 0, name = "Strength") +
          labs(title = "Quadratic Effect Strength", x = "Parameter", y = "Output")
        print(g)
      }
    }
    else if (xvar) {
      if (line.plot) {
        g <- ggplot(data = lin.mat,
                    aes(x = Var2, y = value, group = Var1, colour = Var1)) +
            geom_line(linewidth = 1.2) +
            viridis::scale_color_viridis(discrete = TRUE, name = "Output") +
            theme_bw() +
            labs(title = "Linear Effect Strength",
                 x = "Parameter", y = "Strength") +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        if (!labels) g <- g + theme(legend.position = 'none')
        print(g)
        if(quadratic) {
          g <- ggplot(data = quad.mat,
                      aes(x = Var2, y = value, group = Var1, colour = Var1)) +
            geom_line(linewidth = 1.2) +
            viridis::scale_color_viridis(discrete = TRUE, name = "Output") +
            theme_bw() +
            labs(title = "Quadratic Effect Strength",
                 x = "Parameter", y = "Strength") +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
          if (!labels) g <- g + theme(legend.position = 'none')
          print(g)
        }
      }
      else {
        g <- ggplot(data = lin.mat,
                    aes(x = Var2, y = value, group = Var1, fill = Var1)) +
          geom_col(position = 'dodge', colour = 'black') +
          viridis::scale_fill_viridis(discrete = TRUE, name = "Output") +
          theme_bw() +
          labs(title = "Linear Effect Strength",
               x = "Parameter", y = "Strength") +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        if (!labels) g <- g + theme(legend.position = 'none')
        print(g)
        if (quadratic) {
          g <-ggplot(data = quad.mat,
                     aes(x = Var2, y = value, group = Var1, fill = Var1)) +
            geom_col(position = 'dodge', colour = 'black') +
            viridis::scale_fill_viridis(discrete = TRUE, name = "Output") +
            theme_bw() +
            labs(title = "Quadratic Effect Strength",
                 x = "Parameter", y = "Strength") +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
          if (!labels) g <- g + theme(legend.position = 'none')
          print(g)
        }
      }
    }
    else {
      if (line.plot) {
        g <- ggplot(data = lin.mat,
                    aes(x = Var1, y = value, group = Var2, colour = Var2)) +
          geom_line(linewidth = 1.2) +
          viridis::scale_color_viridis(discrete = TRUE, name = "Parameter") +
          theme_bw() +
          labs(title = "Linear Effect Strength", x = "Output", y = "Strength") +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        if (!labels) g <- g + theme(legend.position = 'none')
        print(g)
        if(quadratic) {
          g <- ggplot(data = quad.mat,
                      aes(x = Var1, y = value, group = Var2, colour = Var2)) +
            geom_line(linewidth = 1.2) +
            viridis::scale_color_viridis(discrete = TRUE, name = "Parameter") +
            theme_bw() +
            labs(title = "Quadratic Effect Strength",
                 x = "Output", y = "Strength") +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
          if (!labels) g <- g + theme(legend.position = 'none')
          print(g)
        }
      }
      else {
        g <- ggplot(data = lin.mat,
                    aes(x = Var1, y = value, group = Var2, fill = Var2)) +
          geom_col(position = 'dodge', colour = 'black') +
          viridis::scale_fill_viridis(discrete = TRUE, name = "Parameter") +
          theme_bw() +
          labs(title = "Linear Effect Strength", x = "Output", y = "Strength") +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        if (!labels) g <- g + theme(legend.position = 'none')
        print(g)
        if (quadratic) {
          g <- ggplot(data = quad.mat,
                     aes(x = Var1, y = value, group = Var2, fill = Var2)) +
            geom_col(position = 'dodge', colour = 'black') +
            viridis::scale_fill_viridis(discrete = TRUE, name = "Parameter") +
            theme_bw() +
            labs(title = "Quadratic Effect Strength",
                 x = "Output", y = "Strength") +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
          if (!labels) g <- g + theme(legend.position = 'none')
          print(g)
        }
      }
    }
  } #nocov end
  return(complete_set)
}
