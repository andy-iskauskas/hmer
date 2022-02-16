#' Sample SIR data
#'
#' A small dataset containing points generated using the Gillespie algorithm.
#' The SIR model contains three input parameters, and generates three output
#' parameters. The initial populations are 950 susceptible (S), 50 infected (I),
#' and 0 recovered (R). The final values are taken at time t=20.
#'
#' @format A data frame with 30 rows and 6 variables:
#' \describe{
#'   \item{aSI}{Infection: transition rate from S to I}
#'   \item{aIR}{Recovery: transition rate from I to R}
#'   \item{aSR}{Immunisation: transition rate from S to R}
#'   \item{nS}{Final number of S}
#'   \item{nI}{Final number of I}
#'   \item{nR}{Final number of R}
#' }
"GillespieSIR"

#' Sample SIR validation data
#'
#' A small dataset containing points generated using the Gillespie algorithm.
#' Very similar to \code{\link{GillespieSIR}}, slightly larger in size.
#'
#' @format A data frame with 60 rows and 6 variables:
#' \describe{
#'   \item{aSI}{Infection: transition rate from S to I}
#'   \item{aIR}{Recovery: transition rate from I to R}
#'   \item{aSR}{Immunisation: transition rate from S to R}
#'   \item{nS}{Final number of S}
#'   \item{nI}{Final number of I}
#'   \item{nR}{Final number of R}
#' }
"GillespieValidation"

#' Sample Implausibility Data
#'
#' A dataset containing 1000 points from the region bounded by
#' [0.1, 0.8], [0, 0.5], [0, 0.05] for aSI, aIR and aSR respectively.
#' Implausibility has been calculated (for emulators trained on the
#' \code{\link{GillespieSIR}} dataset) for each of the outputs nS, nI, nR,
#' and the maximum implausibility is included.
#' The target values used in calculating implausibility were:
#' nS: 281 (sigma 10.43);
#' nI: 30 (sigma 11.16);
#' nR: 689 (sigma 14.32)
#'
#' @format A data frame with 1000 rows and 7 variables:
#' \describe{
#'   \item{aSI}{Infection: transition rate from S to I}
#'   \item{aIR}{Recovery: transition rate from I to R}
#'   \item{aSR}{Immunisation: transition rate from S to R}
#'   \item{nS}{Implausibility for nS}
#'   \item{nI}{Implausibility for nI}
#'   \item{nR}{Implausibility for nR}
#'   \item{I}{Maximum implausibility}
#' }
"GillespieImplausibility"

#' Sample Multi-wave Emulators
#'
#' An rda object containing three waves of emulation on the
#' Gillespie SIR model.
#'
#' @format A list containing \code{\link{Emulator}} objects:
#' \describe{
#'   \item{Wave 1}{Emulators trained on GillespieSIR to generate wave 2 points}
#'   \item{Wave 2}{Emulators trained on the results of the above wave 2 points}
#'   \item{Wave 3}{Emulators trained on the results of the wave 3 points}
#' }
"GillespieMultiWaveEmulators"

#' Sample Multi-wave Results
#'
#' An rda object containing four data.frames: an initial set of points
#' given by \code{GillespieSIR} and \code{GillespieValidation}, and
#' the 90 points generated at each of three subsequent waves. The trained
#' emulators are provided in \code{\link{GillespieMultiWaveEmulators}}.
#'
#' @format A list of data.frame objects:
#' \describe{
#'   \item{Wave 0}{The initial points used in other examples}
#'   \item{Wave 1}{Points generated from the wave 1 emulators}
#'   \item{Wave 2}{Points generated from the wave 2 emulators}
#'   \item{Wave 3}{Points generated from the wave 3 emulators}
#' }
"GillespieMultiWaveData"

#' Sample Emulators
#'
#' An RData object containing three trained emulators, and the associated
#' targets, for the Gillespie SIR example. The emulators have been trained
#' on the \code{\link{GillespieSIR}} dataset using methods documented in
#' this package.
#'
#' @format A list containing two objects:
#' \describe{
#'  \item{ems}{The trained \code{\link{Emulator}} objects.}
#'  \item{targets}{The targets to match to, as a set of val-sigma pairs.}
#' }
"sample_emulators"

#' Birth-Death Model Results
#'
#' An RData object containing two data.frames. The first consists of ten parameter
#' sets run through a simple, two-parameter, stochastic birth-death model; five of the
#' points have 500 replicates and the other five have only 5 replicates. The second
#' consists of ten further points, each with ten replicates. The objects are denoted
#' \code{training} and \code{validation}, representing their expected usage.
#'
#' The initial population for the simulations is 100 people; the model is run until
#' t = 15 to obtain the results to emulate.
#'
#' @format A list of two data.frames \code{training} and \code{validation}: each
#' data.frame has the following columns:
#' \describe{
#'  \item{lambda}{Birth rate}
#'  \item{mu}{Death rate}
#'  \item{Y}{The number of people at time t = 15}
#' }
"BirthDeath"

#' Stochastic SIR Data
#'
#' An RData object consisting of two data.frames (in a similar fashion to \code{BirthDeath}).
#' The first consists of 30 points in the parameter space (aSI, aIR, aSR), each of which has
#' been inputted into the Gillespie algorithm for the stochastic version of the model used
#' in \code{GillespieSIR} (but with changed starting conditions) 100 times. The second has
#' similar form but for 20 unique points, each with 100 repetitions.
#'
#' The outputs observed are the numbers of infected (I) and recovered (R) people at time points
#' t = 10, 25, 50. All outputs display some level of bimodality. The initial conditions to
#' generate the runs had S(0)=995, I(0)=5, R(0)=0.
#'
#' @format A list of two data.frames \code{training} and \code{validation}: each has the
#' following columns:
#' \describe{
#'  \item{aSI}{Infection rate}
#'  \item{aIR}{Recovery rate}
#'  \item{aSR}{Waning immunity rate}
#'  \item{I10 (25, 50)}{The number of infected people at t = 10 (25, 50)}
#'  \item{R10 (25, 50)}{The number of recovered people at t = 10 (25, 50)}
#' }
"SIR_stochastic"
