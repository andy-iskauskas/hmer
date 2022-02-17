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
#' @param logging If a logger exists, it is provided here.
#'
#' @return NULL
#'
#' @keywords internal
#' @noRd
preflight <- function(data, targets, coff = 0.95, logging = NULL) {
  d_abridge <- data[,names(targets)]
  for (i in 1:length(targets)) {
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
        if (!is.null(logging)) logging(level = "WARN", msg = paste("All output values for target",i,"are underestimates. This may suggest model inadequacy."))
        else print(paste("Target", i, "consistently underestimated."))
      }
      else {
        if (!is.null(logging)) logging(level = "WARN", msg = paste("All output values for target",i,"are overestimates. This may suggest model inadequacy."))
        else print(paste("Target", i, "consistently overestimated."))
      }
      hom[,i] <- NULL
    }
  }
  missing_all <- (1:nrow(hom))[apply(hom, 1, function(x) all(x == 0))]
  hom <- hom[-missing_all,]
  checkPairs <- function(df) {
    for (i in 1:(length(df)-1)) {
      for (j in (i+1):length(df)) {
        cval <- cor(df[,c(i,j)])[1,2]
        if (cval < -coff) {
          if (!is.null(logging)) logging(level = "WARN", msg = paste("Strong negative correlation between points satisfying", names(df)[i], "and", names(df)[j], "- this could mean that both targets cannot be hit simultaneously."))
          else print(paste("Strong negative correlation between points satisfying", names(df)[i], "and", names(df)[j]))
        }
        if (cval > coff) {
          if (!is.null(logging)) logging(level = "INFO", msg = paste("Strong positive correlation between points satisfying", names(df)[i], "and", names(df)[j], "- one of these targets may be sufficient."))
          else print(paste("Strong positive correlation between points satisfying", names(df)[i], "and", names(df)[j]))
        }
      }
    }
  }

  checkTriples <- function(df) {
    for (i in 1:(length(df)-1)) {
      for (j in (i+1):length(df)) {
        log_and <- as.numeric(df[,i] & df[,j])
        if (all(log_and == 0) || all(log_and == 1)) next
        for (k in (1:length(df))[-c(i,j)]) {
          cval <- cor(log_and, df[,k])
          if (cval < -coff) {
            if (!is.null(logging)) logging(level = "WARN", msg = paste("Strong negative correlation between points satisfying (",names(df)[i],", ",names(df)[j],") and ", names(df)[k], " - this could mean that the latter target cannot be hit if the other two are.", sep = ""))
            else print(paste("Strong negative correlation between points satisfying (",names(df)[i],", ",names(df)[j],") and ", names(df)[k], sep = ""))
          }
          if (cval > coff) {
            if (!is.null(logging)) logging(level = "INFO", msg = paste("Strong positive correlation between points satisfying (",names(df)[i],", ",names(df)[j],") and ", names(df)[k], " - perhaps not all of these targets are necessary.", sep = ""))
            else print(paste("Strong positive correlation between points satisfying", names(df)[i], "and", names(df)[j]))
          }
        }
      }
    }
  }
  if (length(hom) >= 2) {
    checkPairs(hom)
    if (length(hom) >= 3) {
      checkTriples(hom)
    }
  }
  return(NULL)
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
#' 5) Generate the new points using the defauly method of \code{\link{generate_new_runs}}. The
#' initial implausibility cutoff is set by \code{cutoff}; up to three inflations of the cutoff
#' can be performed. If points still cannot be generated, then the nth implausibility is
#' increased by one (so maximum implausibility becomes second-maximum implausibility).
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
#' @param ... Any arguments to be passed to \code{\link{emulator_from_data}}.
#'
#' @return A list of two objects: \code{points} and \code{emulators}
#'
#' @export
#'
#' @examples
#' \donttest{
#'  ranges <- list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05))
#'  default <- full_wave(rbind(GillespieSIR, GillespieValidation), ranges,
#'   sample_emulators$targets)
#'  non_quad <- full_wave(rbind(GillespieSIR, GillespieValidation), ranges,
#'   sample_emulators$targets, quadratic = FALSE)
#'  second <- full_wave(GillespieMultiWaveData[[2]], ranges, sample_emulators$targets,
#'   old_emulators = GillespieMultiWaveEmulators[[1]])
#'  }
full_wave <- function(data, ranges, targets, old_emulators = NULL, prop_train = 0.7, cutoff = 3, nth = 1, ...) {
  new_ranges <- setNames(purrr::map(names(ranges), ~c(max(ranges[[.]][1], min(data[,.]) - 0.05 * diff(range(data[,.]))), min(ranges[[.]][2], max(data[,.]) + 0.05 * diff(range(data[,.]))))), names(ranges))
  preflight(data, targets)
  samp <- sample(1:nrow(data), floor(prop_train*nrow(data)))
  train <- data[samp,]
  valid <- data[!seq_along(1:nrow(data)) %in% samp,]
  print("Training emulators...")
  tryCatch(
    ems <- emulator_from_data(train, names(targets), new_ranges, ...),
    error = function(e) {
      print(paste("Problem with emulator training:", e))
    }
  )
  for (i in 1:length(ems)) if (ems[[i]]$u_sigma^2 < 1e-8) ems[[i]]$set_sigma(targets[[ems[[i]]$output_name]]$sigma * 0.1)
  print("Performing diagnostics...")
  invalid_ems <- c()
  for (i in 1:length(ems)) {
    comp_fail <- nrow(comparison_diag(ems[[i]], targets, valid, plt = FALSE))
    if (comp_fail > nrow(valid)/4) {
      invalid_ems <- c(invalid_ems, i)
      print(paste("Emulator for output", ems[[i]]$output_name, "fails comparison diagnostics. It will not be matched to at this wave."))
    }
  }
  if (length(invalid_ems) > 0)
    ems <- ems[-invalid_ems]
  if (!any(purrr::map_lgl(targets, is.atomic)))
    emulator_uncerts <- purrr::map_dbl(ems, ~(.$u_sigma^2 + targets[[.$output_name]]$sigma^2)/targets[[.$output_name]]$sigma)
  else emulator_uncerts <- NULL
  if (length(ems) == 0) stop("No emulator passed diagnostic checks.")
  for (i in 1:length(ems)) {
    misclass <- nrow(classification_diag(ems[[i]], targets, valid, cutoff = cutoff, plt = FALSE))
    while(misclass > 0) {
      ems[[i]] <- ems[[i]]$mult_sigma(1.1)
      misclass <- nrow(classification_diag(ems[[i]], targets, valid, plt = FALSE))
    }
  }
  emulator_order <- c(purrr::map_dbl(ems, ~sum(.$implausibility(data, targets[[.$output_name]], cutoff))), use.names = FALSE)
  ems <- ems[order(emulator_order)]
  if (!is.null(old_emulators))
    working_ems <- c(ems, old_emulators)
  else
    working_ems <- ems
  targets <- targets[purrr::map_chr(working_ems, ~.$output_name)]
  print("Generating new points...")
  new_points <- generate_new_runs(working_ems, nrow(data), targets, cutoff = cutoff, verbose = FALSE, nth = nth)
  if (nrow(new_points) == 0) stop("Could not generate points in non-implausible space.")
  if (!any(purrr::map_lgl(targets, is.atomic)))
    if (!is.null(emulator_uncerts) && all(emulator_uncerts < 1.02) && length(emulator_uncerts) == length(targets)) print("All emulator uncertainties are comparable to target uncertainties. More waves may be unnecessary.")
  return(list(points = new_points, emulators = ems))
}
