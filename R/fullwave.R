#' Pre-flight Checks for Data
#'
#' Performs checks to see if there are redundancies or contradictions in output data.
#'
#' Given the set of data and the targets that will be used to train emulators, three
#' checks are performed. The most rudimentary is checking that each target is likely to
#' be matched to - that is, either there are output values within the target bounds
#' or at least there are both overestimates and underestimates.
#'
#' Checks are also performed as to whether targets are mutually contradictory: i.e. if
#' one target can only be hit if another is missed. The function also checks if one
#' target is hit every time a different target is, suggesting target redundancy. The
#' same analysis is performed for every triple of targets.
#'
#' @importFrom stats cor
#'
#' @param data The simulator data.
#' @param targets The targets to match to, either as bounds or (val, sigma) pairs.
#' @param coff The correlation value above which two outputs are assumed correlated.
#' @param verbose Should feedback be provided to console?
#'
#' @return NULL
#'
#' @keywords internal
#' @noRd
preflight <- function(data, targets, coff = 0.95, verbose = interactive()) {
  potential_problem <- FALSE
  applicable_targets <- intersect(names(data), names(targets))
  d_abridge <- data[,applicable_targets]
  targets <- targets[applicable_targets]
  for (i in seq_along(targets)) {
    if (!is.atomic(targets[[i]])) targets[[i]] <- c(targets[[i]]$val - 3*targets[[i]]$sigma, targets[[i]]$val + 3*targets[[i]]$sigma)
  }
  getHits <- function(data, targets) {
    hmm <- data.frame(do.call('cbind', purrr::map(seq_along(names(data)), ~data[,names(data)[.]] >= targets[[names(data)[[.]]]][1] & data[,names(data)[.]] <= targets[[names(data)[.]]][2])))
    hmm <- setNames(data.frame(t(apply(hmm, 1, as.numeric))), names(data))
    return(hmm)
  }
  hom <- getHits(d_abridge, targets)
  for (i in names(hom)) {
    if (all(hom[,i] == 0)) {
      if (all(data[,i] < targets[[i]][1])) {
        potential_problem <- TRUE
        if (verbose) cat("Target", i, "consistently underestimated.\n") #nocov
      }
      else {
        potential_problem <- TRUE
        if (verbose) cat("Target", i, "consistently overestimated.\n") #nocov
      }
      hom[,i] <- NULL
    }
  }
  missing_all <- (seq_len(nrow(hom)))[apply(hom, 1, function(x) all(x == 0))]
  hom <- hom[-missing_all,]
  checkPairs <- function(df) {
    for (i in 1:(length(df)-1)) {
      for (j in (i+1):length(df)) {
        cval <- cor(df[,c(i,j)])[1,2]
        if (cval < -coff) {
          if (verbose) cat("Strong negative correlation between points satisfying", #nocov
                            names(df)[i], "and", names(df)[j], "\n") #nocov
        }
        if (cval > coff) {
          if (verbose) cat("Strong positive correlation between points satisfying", #nocov
                            names(df)[i], "and", names(df)[j], "\n") #nocov
        }
      }
    }
  }

  checkTriples <- function(df) {
    for (i in 1:(length(df)-1)) {
      for (j in (i+1):length(df)) {
        log_and <- as.numeric(df[,i] & df[,j])
        if (all(log_and == 0) || all(log_and == 1)) next
        for (k in (seq_along(df))[-c(i,j)]) {
          cval <- cor(log_and, df[,k])
          if (cval < -coff) {
            if (verbose)
              cat("Strong negative correlation between points satisfying (", #nocov start
                    names(df)[i],", ",names(df)[j],
                    ") and ",
                    names(df)[k], "\n", sep = "") #nocov end
          }
          if (cval > coff) {
            if (verbose)
              cat("Strong positive correlation between points satisfying", #nocov
                    names(df)[i], "and", names(df)[j], "\n") #nocov
          }
        }
      }
    }
  }
  if (!is.null(nrow(hom)) && nrow(hom) >= 2) {
    checkPairs(hom)
    if (!is.null(nrow(hom)) && nrow(hom) >= 3) {
      checkTriples(hom)
    }
  }
  return(potential_problem)
}

#' Automatic Wave Calculation
#'
#' Performs a full wave of emulation and history matching, given data.
#'
#' This function uses all of the functionality from the package in a relatively conservative form.
#' The function performs the following steps:
#'
#' 1) Split the data into a training set and a validation set, where \code{prop_train} indicates
#' what proportion of the data is used to train.
#'
#' 2) Perform emulator training using \code{\link{emulator_from_data}}. If a more involved
#' specification is desired, optional arguments can be passed to \code{emulator_from_data} using
#' the \code{...} argument.
#'
#' 3) Perform diagnostics on the trained emulators, removing emulators that do not display
#' acceptable performance. Global emulator variance may also be modified to ensure that none of
#' the emulators demonstrate misclassification errors (from \code{\link{classification_diag}}).
#'
#' 4) Ordering the remaining emulators from most restrictive to least restrictive on the dataset
#' provided at this wave. Some point generation mechanisms terminate early if a point is ruled
#' out by a single emulator, so the ordering ensures this happens earlier rather than later.
#'
#' 5) Generate the new points using the default method of \code{\link{generate_new_runs}}, using
#' the normal procedure (for details, see the description for generate_new_runs).
#'
#' If the parameter \code{old_emulators} is provided, this should be a list of emulators used
#' at all previous waves - for example if \code{full_wave} is used to do a second wave of
#' history matching, then \code{old_emulators} would contain the list of first-wave emulators.
#'
#' The function returns a list of two objects: \code{emulators} corresponding to this wave's
#' emulators, and \code{points} corresponding to the new proposed points. The points can then
#' be put into the simulator to generate runs for a subsequent wave.
#'
#' @param data The data to train with.
#' @param ranges The ranges of the input parameters
#' @param targets The output targets to match to.
#' @param old_emulators Any emulators from previous waves.
#' @param prop_train What proportion of the data is used for training.
#' @param cutoff The implausibility cutoff for point generation and diagnostics.
#' @param nth The level of maximum implausibility to consider.
#' @param verbose Should progress be printed to console?
#' @param ... Any arguments to be passed to \code{\link{emulator_from_data}}.
#'
#' @return A list of two objects: \code{points} and \code{emulators}
#'
#' @export
#'
#' @examples
#' \donttest{
#'  ranges <- list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05))
#'  default <- full_wave(do.call('rbind.data.frame', SIRSample), ranges,
#'   SIREmulators$targets)
#'  non_quad <- full_wave(do.call('rbind.data.frame', SIRSample), ranges,
#'   SIREmulators$targets, quadratic = FALSE)
#'  second <- full_wave(SIRMultiWaveData[[2]], ranges, SIREmulators$targets,
#'   old_emulators = SIRMultiWaveEmulators[[1]])
#'  }
full_wave <- function(data, ranges, targets, old_emulators = NULL,
                      prop_train = 0.7, cutoff = 3, nth = 1,
                      verbose = interactive(), ...) {
  new_ranges <- setNames(
    purrr::map(
      names(ranges),
      ~c(max(ranges[[.]][1], min(data[,.]) - 0.05 * diff(range(data[,.]))),
         min(ranges[[.]][2], max(data[,.]) + 0.05 * diff(range(data[,.]))))),
    names(ranges))
  preflight(data, targets, verbose = verbose)
  samp <- sample(seq_len(nrow(data)), floor(prop_train*nrow(data)))
  train <- data[samp,]
  valid <- data[-samp,]
  if (verbose) cat("Training emulators...\n") #nocov
  tryCatch(
    ems <- emulator_from_data(train, names(targets), new_ranges, verbose = verbose, ...),
    error = function(e) {
      stop(paste("Problem with emulator training:", e))
    }
  )
  for (i in seq_along(ems)) {
    if (ems[[i]]$u_sigma^2 < 1e-8)
      ems[[i]]$set_sigma(targets[[ems[[i]]$output_name]]$sigma * 0.1)
  }
  if (verbose) cat("Performing diagnostics...\n") #nocov
  invalid_ems <- c()
  for (i in seq_along(ems)) {
    comp_fail <- nrow(comparison_diag(ems[[i]], targets, valid, plt = FALSE))
    if (comp_fail > nrow(valid)/4) {
      invalid_ems <- c(invalid_ems, i)
      if (verbose) cat("Emulator for output", #nocov start
                          ems[[i]]$output_name,
                          "fails comparison diagnostics.",
                          "It will not be matched to at this wave.\n") #nocov end
    }
  }
  if (length(invalid_ems) > 0)
    ems <- ems[-invalid_ems]
  if (!any(purrr::map_lgl(targets, is.atomic)))
    emulator_uncerts <- purrr::map_dbl(
      ems, ~(.$u_sigma^2 +
               targets[[.$output_name]]$sigma^2)/targets[[.$output_name]]$sigma)
  else emulator_uncerts <- NULL
  if (length(ems) == 0) stop("No emulator passed diagnostic checks.")
  for (i in seq_along(ems)) {
    misclass <- nrow(classification_diag(ems[[i]], targets, valid,
                                         cutoff = cutoff, plt = FALSE))
    while(misclass > 0) {
      ems[[i]] <- ems[[i]]$mult_sigma(1.1)
      misclass <- nrow(classification_diag(ems[[i]],
                                           targets, valid, plt = FALSE))
    }
  }
  emulator_order <- c(
    purrr::map_dbl(
      ems,
      ~sum(.$implausibility(data, targets[[.$output_name]], cutoff))),
    use.names = FALSE)
  ems <- ems[order(emulator_order)]
  if (!is.null(old_emulators))
    working_ems <- c(ems, old_emulators)
  else
    working_ems <- ems
  targets <- targets[purrr::map_chr(working_ems, ~.$output_name)]
  if (verbose) cat("Generating new points...\n") #nocov
  new_points <- generate_new_runs(working_ems, nrow(data), targets,
                                  cutoff = cutoff, verbose = FALSE, nth = nth)
  if (nrow(new_points) == 0)
    stop("Could not generate points in non-implausible space.")
  if (!any(purrr::map_lgl(targets, is.atomic)))
    if (!is.null(emulator_uncerts) &&
        all(emulator_uncerts < 1.02) &&
        length(emulator_uncerts) == length(targets))
      if (verbose) cat("All emulator uncertainties are comparable to target", #nocov
                  "uncertainties. More waves may be unnecessary.\n") #nocov
  return(list(points = new_points, emulators = ems))
}
