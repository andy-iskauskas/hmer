% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{emulator_plot}
\alias{emulator_plot}
\title{Plot Emulator Outputs}
\usage{
emulator_plot(
  ems,
  plot_type = "exp",
  ppd = 30,
  targets = NULL,
  cb = FALSE,
  params = NULL,
  fixed_vals = NULL,
  nth = 1,
  imp_breaks = NULL,
  include_legend = TRUE
)
}
\arguments{
\item{ems}{An \code{\link{Emulator}} object, or a list thereof.}

\item{plot_type}{The statistic to plot (see description or examples).}

\item{ppd}{The number of points per plotting dimension}

\item{targets}{If required, the targets from which to calculate implausibility}

\item{cb}{A boolean representing whether a colourblind-friendly plot is produced.}

\item{params}{Which two input parameters should be plotted?}

\item{fixed_vals}{For fixed input parameters, the values they are held at.}

\item{nth}{If plotting nth maximum implausibility, which level maximum to plot.}

\item{imp_breaks}{If plotting nth maximum implausibility, defines the levels at
which to draw contours.}

\item{include_legend}{For multiple plots, should a combined legend be appended?}
}
\value{
A ggplot object, or collection thereof.
}
\description{
A function for plotting emulator expectations, variances, and implausibilities
}
\details{
Given a single emulator, or a set of emulators, the emulator statistics can be plotted
across a two-dimensional slice of the parameter space. Which statistic is plotted is
determined by \code{plot_type}: options are `exp', `var', `sd', `imp', and `nimp', which
correspond to expectation, variance, standard deviation, implausibility, and nth-max
implausibility.

By default, the slice varies in the first two parameters of the emulators, and all other
parameters are taken to be fixed at their mid-range values. This behaviour can be changed
with the \code{params} and \code{fixed_vals} parameters (see examples).

If the statistic is `exp', `var' or `sd', then the minimal set of parameters to pass to this
function are \code{ems} (which can be a list of emulators or a single one) and \code{plot_type}.
If the statistic is `imp' or `nimp', then the \code{targets} must be supplied - it is not
necessary to specify the individual target for a single emulator plot. If the statistic is
`nimp', then the level of maximum implausibility can be chosen with the parameter \code{nth}.

Implausibility plots are typically coloured from green (low implausibility) to red (high
implausibility): a colourblind-friendly option is available and can be turned on by setting
\code{cb = TRUE}.

The granularity of the plot is controlled by the \code{ppd} parameter, determining the number
of points per dimension in the grid. For higher detail, at the expense of longer computing
time, increase this value. The default is 30.
}
\examples{
 # Reducing ppd to 10 for speed.
 emulator_plot(SIREmulators$ems, ppd = 10)
 emulator_plot(SIREmulators$ems$nS, ppd = 10)
 emulator_plot(SIREmulators$ems, plot_type = 'var', ppd = 10, params = c('aIR', 'aSR'))
 \donttest{ # Excessive runtime
    emulator_plot(SIREmulators$ems, plot_type = 'imp', ppd = 10,
     targets = SIREmulators$targets,
     fixed_vals = list(aSR = 0.02))
    emulator_plot(SIREmulators$ems, plot_type = 'nimp', cb = TRUE,
     targets = SIREmulators$targets, nth = 2, ppd = 10)
 }

}
\seealso{
Other visualisation tools: 
\code{\link{behaviour_plot}()},
\code{\link{diagnostic_wrap}()},
\code{\link{effect_strength}()},
\code{\link{hit_by_wave}()},
\code{\link{output_plot}()},
\code{\link{plot_actives}()},
\code{\link{plot_lattice}()},
\code{\link{plot_wrap}()},
\code{\link{simulator_plot}()},
\code{\link{space_removed}()},
\code{\link{validation_pairs}()},
\code{\link{wave_dependencies}()},
\code{\link{wave_points}()},
\code{\link{wave_values}()}
}
\concept{visualisation tools}
