#' Sample SIR data
#'
#' A small dataset containing points generated from a simple deterministic SIR model.
#' The model contains three input parameters, and generates three output
#' parameters. The initial populations are 950 susceptible (S), 50 infected (I),
#' and 0 recovered (R). The final values are taken at time t=10.
#'
#' The model operates using simple differential equations, where
#'
#' S' = aSR*R - aSI*S*R/(S+I+R)
#'
#' I' = aSI*S*R/(S+I+R) - aIR*I
#'
#' R' = aIR*I - aSR*R.
#'
#' @format A list of two data frames. The first has 30 rows and 6 variables, the second
#' 60 rows and 6 variables. The structure is the same in both cases:
#' \describe{
#'   \item{aSI}{Infection: transition rate from S to I}
#'   \item{aIR}{Recovery: transition rate from I to R}
#'   \item{aSR}{Immunisation: transition rate from S to R}
#'   \item{nS}{Final number of S}
#'   \item{nI}{Final number of I}
#'   \item{nR}{Final number of R}
#' }
"SIRSample"

#' Sample Implausibility Data
#'
#' A dataset containing 1000 points from the region bounded by
#' [0.1, 0.8], [0, 0.5], [0, 0.05] for aSI, aIR and aSR respectively.
#' Implausibility has been calculated (for emulators trained on the
#' \code{\link{SIRSample}} training dataset) for each of the outputs
#' nS, nI, nR, and the maximum implausibility is included.
#' The target values used in calculating implausibility were:
#' \describe{
#'   \item{nS}{between 324 and 358}
#'   \item{nI}{mean 143 (sigma 7.15)}
#'   \item{nR}{between 490 and 542}
#' }
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
"SIRImplausibility"

#' Sample Multi-wave Emulators
#'
#' An rda object containing three waves of emulators applied to
#' SIR model (described in \code{\link{SIRSample}}). The corresponding points
#' (both training and validation) are stored in \code{\link{SIRMultiWaveData}}.
#'
#' @format A list containing \code{\link{Emulator}} objects:
#' \describe{
#'   \item{Wave 1}{Emulators trained on Wave 0, generating wave 1 points}
#'   \item{Wave 2}{Emulators trained on the results of the above wave 1 points}
#'   \item{Wave 3}{Emulators trained on the results of the above wave 2 points}
#' }
"SIRMultiWaveEmulators"

#' Sample Multi-wave Results
#'
#' An rda object containing four data.frames: an initial set of points
#' also provided in \code{\link{SIRSample}}, and
#' the 90 points generated at each of three subsequent waves. The trained
#' emulators are provided in \code{\link{SIRMultiWaveEmulators}}.
#'
#' @format A list of data.frame objects:
#' \describe{
#'   \item{Wave 0}{The initial points used in other examples}
#'   \item{Wave 1}{Points generated from the wave 1 emulators}
#'   \item{Wave 2}{Points generated from the wave 2 emulators}
#'   \item{Wave 3}{Points generated from the wave 3 emulators}
#' }
"SIRMultiWaveData"

#' Sample Emulators
#'
#' An RData object containing three trained emulators, and the associated
#' targets, for the SIR example. The emulators have been trained
#' on the \code{\link{SIRSample}} training dataset using methods documented in
#' this package.
#'
#' @format A list containing two objects:
#' \describe{
#'  \item{ems}{The trained \code{\link{Emulator}} objects.}
#'  \item{targets}{The targets to match to, as a named list.}
#' }
"SIREmulators"

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

#' Data for an interesting emulation problem
#'
#' An RData object consisting of four objects: a data.frame \code{data} of 208 points,
#' a set \code{targets} of 19 targets for outputs, a set \code{ranges} of 21 ranges for inputs,
#' and a data.frame \code{extra} of 26 additional points. This dataset is used to demonstrate
#' some of the subtleties of emulation in the vignettes, where data transformations can be
#' useful and careful attention should be paid to emulation at early waves.
#'
#' @format A list of objects:
#' \describe{
#'  \item{data}{The training data of 'space-filling' runs}
#'  \item{targets}{The output targets to match to}
#'  \item{ranges}{The input ranges over which the system is valid}
#'  \item{extra}{A set of 'extra' points, generated around a known point of best fit.}
#' }
"problem_data"
