# Implausibility colour palettes
redgreen <- c('#00FF00', '#18FF00', '#31FF00', '#49FF00', '#62FF00',
              '#7AFF00', '#93FF00', '#ABFF00', '#C4FF00', '#DDFF00',
              '#E0E200', '#E4C600', '#E8AA00', '#EC8D00', '#EF7100',
              '#F35500', '#F73800', '#FB1C00', '#FF0000', '#FF0000')
colourblind <- c('#1aff1a', '#2af219', '#3ae618', '#4ada17', '#5acd16',
                 '#6bc115', '#7bb514', '#8ba813', '#9b9c12', '#ac9011',
                 '#a2831d', '#98762a', '#8e6936', '#845d43', '#7a504f',
                 '#70435c', '#663768', '#5c2a75', '#521d81', '#48118e')
redgreencont <- list(low = '#00FF00', mid = '#DDFF00', high = '#FF0000')
colourblindcont <- list(low = '#1aff1a', mid = '#ac9011', high = '#48118e')

# Scales inputs: important since the emulators should take inputs purely in [-1,1]
scale_input <- function(x, r, forward = TRUE) {
  centers <- map(r, ~(.x[2]+.x[1])/2)
  scales <- map(r, ~(.x[2]-.x[1])/2)
  if (is.null(names(x))) {
    centers <- unlist(centers, use.names = F)
    scales <- unlist(scales, use.names = F)
  }
  if (forward) {
    if(is.null(names(x)) && !is.null(dim(x))) {
      output <- t(apply(
        t(apply(x, 1, function(y) y - centers)),
        1, function(z) z/scales))
    }
    else
      output <- (x-centers)/scales
  }
  else {
    if (is.null(names(x)) && !is.null(dim(x))) {
      if (is.null(dim(x))) x <- matrix(x, ncol = 1)
      output <- t(apply(
        t(apply(x, 1, function(y) y * scales)),
        1, function(z) z + centers
      ))
    }
    else
      output <- x * scales + centers
  }
  if (!is.data.frame(output)) {
    output <- data.frame(output)
    names(output) <- NULL
  }
  return(output)
}

# Helper to convert functions to names
function_to_names <- function(f, var_names, function_form = TRUE) {
  f_str <- deparse(body(f))
  first_sub <- gsub("x\\[{1,2}(\\d*)\\]{1,2}(^\\d)?",
                    "var_names[\\1]\\2", f_str)
  matching <- gregexec("var_names\\[\\d*\\]", first_sub)[[1]]
  if (length(matching) == 1 && attr(matching, "match.length") == -1) return(first_sub)
  match_locs <- matching[1,]
  match_lengths <- attr(matching, "match.length")
  subbed <- ""
  match_index <- 1
  i <- 1
  while (i <= nchar(first_sub)) {
    if (i %in% match_locs) {
      end_val <- i + match_lengths[match_index]
      subbed <- paste0(subbed, eval(parse(text = substr(first_sub, i, end_val-1))))
      i <- end_val
    }
    else {
      subbed <- paste0(subbed, substr(first_sub, i, i))
      i <- i + 1
    }
  }
  if (function_form) {
    subbed <- gsub("\\s\\*\\s", ":", subbed)
  }
  else {
    subbed <- gsub("\\s\\*\\s", "\\*", subbed)
  }
  return(subbed)
}
# function_to_names <- function(f, var_names, function_form = TRUE) {
#   f_str <- deparse(body(f))
#   subbed <- stringr::str_replace_all(
#     gsub("x\\[{1,2}(\\d*)\\]{1,2}(^\\d)?",
#          "var_names[\\1]\\2", f_str),
#     "var_names\\[\\d*\\]", function(x) eval(parse(text = x)))
#   if (function_form) {
#     subbed <- gsub("\\s\\*\\s", ":", subbed)
#   }
#   else {
#     subbed <- gsub("\\s\\*\\s", "\\*", subbed)
#   }
#   return(subbed)
# }

name_to_function <- function(str, var_names) {
  if (str == "(Intercept)") return(function(x) 1)
  str_power <- gsub("I\\((.*)\\^(\\d+)\\)", "\\1\\^\\2", str)
  str_times <- gsub(":", "*", str_power)
  str_par <- str_times
  var_names_ordered <- var_names[order(nchar(var_names), var_names,  decreasing = TRUE)]
  for (i in seq_along(var_names_ordered)) {
    str_par <- gsub(var_names_ordered[i], paste0("\u00a3\u00a3\u00a3\u00a3[[", which(var_names == var_names_ordered[i]), "]]\\1"), str_par)
  }
  str_par <- gsub("\u00a3\u00a3\u00a3\u00a3", "x", str_par)
  return(eval(parse(text = paste('function(x) return(', str_par, ')', sep = ''))))
}

# Evaluate multiple functions over points
eval_funcs <- function(funcs, points, ...) {
  pointsdim <- (length(dim(points)) != 0)
  manyfuncs <- (typeof(funcs) != "closure")
  if (manyfuncs && pointsdim)
    return(apply(points, 1,
                 function(x) map_dbl(funcs, exec, x, ...)))
  if (manyfuncs)
    return(map_dbl(funcs, exec, points, ...))
  if (pointsdim) {
    return(tryCatch(apply(points, 1, funcs, ...),
                    error = function(cond1) {
                      tryCatch(exec(funcs, points, ...),
                               error = function(cond2) {
                                 cat(cond1$message, "\n", cond2$message, "\n")
                                 stop()
                               })
                    }))
  }
  return(exec(funcs, points, ...))
}

# Split datasets by uniqueness of columns
split_dataset <- function(data, split_vars, method = "hash") {
  if (method != "dplyr") {
    uids <- apply(data[,split_vars, drop = FALSE], 1, hash)
    grouping <- map(unique(uids), ~data[which(uids == .),])
  }
  else {
    dplyr_group <- data |> group_by(across(all_of(split_vars)))
    grouping <- map(group_rows(dplyr_group), ~data[.,])
  }
  return(grouping)
}

# Inner modification of a function
multiply_function <- function(f, mult) {
  func_body <- body(f)
  if (length(func_body) == 2) {
    relev <- func_body[[2]]
    if (is.language(relev) && !(is.symbol(relev) || is.double(relev))) {
      innards <- relev[[2]]
      body(f)[[2]][[2]] <- substitute(mult * innards)
    }
    if (is.symbol(relev) || is.double(relev)) {
      body(f)[[2]] <- substitute(mult * relev)
    }
    return(f)
  }
  if (length(func_body) == 1) {
    body(f) <- substitute(mult * func_body)
    return(f)
  }
  relev <- func_body[[length(func_body)]]
  if (is.language(relev) && !(is.symbol(relev) || is.double(relev))) {
    innards <- relev[[2]]
    body(f)[[length(func_body)]][[2]] <- substitute(mult * innards)
    return(f)
  }
  if (is.symbol(relev) || is.double(relev)) {
    body(f)[[length(func_body)]] <- substitute(mult * relev)
    return(f)
  }
  warning("Function multiplication not successful. Returning original function.")
  return(f)
}

# Kurtosis estimator
kurtosis <- function(x, na.rm = FALSE) {
  if (is.matrix(x)) apply(x, 2, kurtosis, na.rm = na.rm)
  else if (is.vector(x)) {
    if(na.rm) x <- x[!is.na(x)]
    n <- length(x)
    n * sum((x-mean(x))^4)/sum((x-mean(x))^2)^2
  }
  else if (is.data.frame(x)) vapply(x, kurtosis, numeric(1), na.rm = na.rm)
  else kurtosis(as.vector(x), na.rm = na.rm)
}

#' Truncation of t-distribution
#'
#' Produces moments of a truncated t-distribution
#'
#' For stochastic disease models, it can be useful to truncate the available
#' output of an emulator. For example, if emulating the variance of a stochastic
#' output it is possible for an emulator to predict a negative value of the variance.
#' In this circumstance, we instead truncate by assuming a t-distribution over the
#' predicted expectation and variance and calculating the expectation and variance of
#' its truncated distribution to positive values.
#'
#' This function works for general truncations of data. The default behaviour is a
#' truncation to non-negative values, using a t-distribution with 6 degrees of freedom
#' (so that nu = 6, a = 0, b = Inf).
#'
#' @param e The original expectation
#' @param v The original variance
#' @param mu Should the modified expectation be returned? If FALSE, then variance is given.
#' @param nu The number of degrees of freedom of the underlying t-distribution
#' @param a The left truncation point
#' @param b The right truncation point
#'
#' @importFrom stats pt
#'
#' @keywords internal
#' @noRd
#'
#' @return Either the new expectation or the new variance.
get_truncation <- function(e, v, mu = TRUE, nu = 6, a = 0, b = Inf) {
  if (v <= 0) v <- 1e-10
  eff_v <- (nu-2)*v/nu
  new_a <- (a-e)/sqrt(eff_v)
  new_b <- (b-e)/sqrt(eff_v)
  a0 <- pt(new_b, nu) - pt(new_a, nu)
  kappa <- gamma((nu+1)/2)/(a0 * gamma(nu/2) * sqrt(nu*pi))
  tauj <- function(j) (nu-2*j)/nu
  m1 <- kappa*nu/(nu-1) *
    (((nu+new_a^2)/nu)^(-(nu-1)/2)-((nu+new_b^2)/nu)^(-(nu-1)/2))
  if (mu) return (sqrt(eff_v)*m1+e)
  m2 <- (nu-1)/tauj(1) *
    ((pt(new_b*sqrt(tauj(1)), nu-2)-pt(new_a*sqrt(tauj(1)), nu-2))/a0)-nu
  new_v <- m2-m1^2
  return(eff_v*new_v)
}

#' Subsetting for Bimodal/Variance Emulators
#'
#' Takes a collection of bimodal or stochastic emulators and subsets by output name.
#'
#' It can be useful to consider only a subset of outputs. In the normal case, this can be
#' easily achieved; however, when the emulators are in a nested structure such as that
#' provided by emulator_from_data with emulator_type = 'variance' or 'bimodal', it can
#' be more involved. This function allows the easy selecting of emulators by name, returning a
#' subset of them in the same form as the original object.
#'
#' This function is compatible with `standard' emulators; that is, those in a simple
#' list, equivalent to subsetting over the collection of output names of the emulators
#' that exist in \code{output_names}.
#'
#' @param emulators A set of emulators, often in nested form
#' @param output_names The names of the desired outputs
#'
#' @return An object of the same form as `emulators`.
#' @export
subset_emulators <- function(emulators, output_names) {
  if (!is.null(emulators$mode1)) {
    m1exp <- emulators$mode1$expectation[
      map_chr(emulators$mode1$expectation,
                     ~.$output_name) %in% output_names]
    m1var <- emulators$mode1$variance[
      map_chr(emulators$mode1$variance,
                     ~.$output_name) %in% output_names]
    m2exp <- emulators$mode2$expectation[
      map_chr(emulators$mode2$expectation,
                     ~.$output_name) %in% output_names]
    m2var <- emulators$mode2$variance[
      map_chr(emulators$mode2$variance,
                     ~.$output_name) %in% output_names]
    collated <- list(
      mode1 = list(
        variance = m1var,
        expectation = m1exp
      ),
      mode2 = list(
        variance = m2var,
        expectation = m2exp
      ),
      prop = emulators$prop
    )
  }
  else if (!is.null(emulators$variance)) {
    mexp <- emulators$expectation[
      map_chr(emulators$expectation,
                     ~.$output_name) %in% output_names]
    mvar <- emulators$variance[
      map_chr(emulators$variance,
                     ~.$output_name) %in% output_names]
    collated <- list(
      expectation = mexp,
      variance = mvar
    )
  }
  else {
    collated <- emulators[map(emulators, ~.$output_name) %in% output_names]
  }
  return(collated)
}

#' Collect and order emulators
#'
#' Manipulates lists (or lists of lists) of emulators into a useable form.
#'
#' Most often used as a pre-processing stage for \code{generate_new_design} or
#' \code{nth_implausible}, this takes a list of emulators in a variety of forms
#' coming from either multiple waves of history matching, hierarchical emulation
#' or bimodal emulation, and arrange them in a form suitable for sequential analysis.
#' Emulators are also ordered by a number of factors: number of parameters, size of
#' the minimum enclosing hyperrectangle, and implausibility (where applicable).
#'
#' If targets are provided, then the emulators can also be tested on the basis of
#' how restrictive they are: this is included in the ordering. The cutoff by which to
#' make a determination of implausibility for a point is governed by \code{cutoff}.
#' The number of points to sample to consider implausibility is chosen by the value
#' of \code{sample_size}: higher values are likely to be a more accurate reflection
#' but will take longer.
#'
#' The weighting of each of the three metrics can be chosen using the \code{ordering}
#' argument: metrics with higher weight are closer to the front of the character vector.
#' The metrics are denoted "params" for number of parameters, "imp" for restrictiveness,
#' and "volume" for volume of the hyperrectangle. For instance, a character vector
#' \code{c("volume", "imp")} would sort first by volume of the minimum enclosing hyperrectangle,
#' resolve ties by restrictiveness, and not consider the number of parameters.
#'
#' @param emulators The recursive list of emulators
#' @param targets If not NULL, uses implausibility to order the emulators.
#' @param cutoff The implausibility cutoff to use (if required)
#' @param ordering The order in which to apply the relevant metrics
#' @param sample_size The number of points to apply implausibility to (if required)
#' @param ... Any additional arguments to pass recursively to collect_emulators
#'
#' @return A list of emulators with the ordered property described above.
#'
#' @export
collect_emulators <- function(emulators, targets = NULL, cutoff = 3,
                              ordering = c("params", "imp", "volume"),
                              sample_size = 200, ...) {
  if (inherits(emulators, "Emulator"))
    return(setNames(list(emulators), emulators$output_name))
  if (all(map_lgl(emulators, ~inherits(., "Emulator")))) {
    em_names <- map_chr(emulators, ~.$output_name)
    em_range_lengths <- map_dbl(emulators, ~length(.$ranges))
    em_range_prods <- map_dbl(emulators,
                                     ~prod(map_dbl(.$ranges, diff)))
    if (is.null(targets)) {
      which_ordering <- match(ordering, c("params", "volume"))
    } else {
      which_ordering <- match(ordering, c("params", "imp", "volume"))
    }
    which_ordering <- which_ordering[!is.na(which_ordering)]
    if (length(which_ordering) == 0) {
      warning("No valid ordering parameters specified. Reverting to default.")
      if (is.null(targets)) which_ordering <- 1:2
      else which_ordering <- 1:3
    }
    if (!is.null(targets)) {
      if (!is.null(targets$expectation)) targets <- targets$expectation
      for (i in seq_along(targets)) {
        if (length(targets[[i]]) == 1) {
          targets[[i]] <- list(val = targets[[i]], sigma = 0.05*targets[[i]])
        }
      }
      maximal_ranges <- map(emulators, "ranges")[[which.max(map_dbl(emulators, ~length(.$ranges)))]]
      sample_points <- do.call('cbind.data.frame', map(maximal_ranges, ~runif(sample_size, min = .[[1]], max = .[[2]]))) |> setNames(names(maximal_ranges))
      if (length(cutoff) == 1) cutoff <- rep(cutoff, length(emulators))
      em_imps <- do.call('cbind', map(seq_along(emulators), ~emulators[[.]]$implausibility(sample_points, targets[[emulators[[.]]$output_name]], cutoff = cutoff[[.]])))
      imp_restriction <- apply(em_imps, 2, sum)
      em_ordering <- do.call(order, list(-em_range_lengths, imp_restriction, em_range_prods)[which_ordering])
    } else {
      em_ordering <- do.call(order, list(-em_range_lengths, em_range_prods)[which_ordering])
    }
    return(setNames(emulators[em_ordering], em_names[em_ordering]))
  }
  if ((!is.null(emulators$expectation) && sum(names(emulators) == "expectation") == 1) || (!is.null(emulators$mode1) && sum(names(emulators) == "mode1") == 1)) {
    return(emulators)
  }
  if ("expectation" %in% names(emulators)) {
    exp_ems <- c(emulators[names(emulators) == "expectation"], use.names = FALSE)
    var_ems <- c(emulators[names(emulators) == "variance"], use.names = FALSE)
    if (!is.null(targets$expectation)) {
      targs_exp <- targets$expectation
      targs_var <- targets$variance
    }
    else {
      targs_exp <- targets
      targs_var <- NULL
    }
    return(list(expectation = collect_emulators(exp_ems, targets = targs_exp),
                variance = collect_emulators(var_ems, targets = targs_var)))
  }
  if ("mode1" %in% names(emulators)) {
    m1ems <- c(emulators[names(emulators) == "mode1"], use.names = FALSE)
    m2ems <- c(emulators[names(emulators) == "mode2"], use.names = FALSE)
    prop_ems <- c(emulators[names(emulators) == "prop"], use.names = FALSE)
    return(list(mode1 = collect_emulators(m1ems, targets = NULL),
                mode2 = collect_emulators(m2ems, targets = NULL),
                prop = collect_emulators(prop_ems)))
  }
  if (!is.null(emulators[[1]]$mode1)) {
    m1ems <- map(emulators, ~.$mode1)
    m2ems <- map(emulators, ~.$mode2)
    prop_ems <- map(emulators, ~.$prop)
    return(list(mode1 = collect_emulators(m1ems, targets = NULL),
                mode2 = collect_emulators(m2ems, targets = NULL),
                prop = collect_emulators(prop_ems)))
  }
  if (!is.null(emulators[[1]]$expectation)) {
    exp_ems <- map(emulators, ~.$expectation)
    var_ems <- map(emulators, ~.$variance)
    if (!is.null(targets$expectation)) {
      targs_exp <- targets$expectation
      targs_var <- targets$variance
    }
    else {
      targs_exp <- targets
      targs_var <- NULL
    }
    return(list(expectation = collect_emulators(exp_ems, targets = targs_exp),
                variance = collect_emulators(var_ems, targets = targs_var)))
  }
  return(collect_emulators(unlist(emulators), targets = targets))
}

#' Obtain the parameter ranges from a collection of emulators
#'
#' This is a more complex version of the \code{em$ranges} command, which accommodates
#' the recursive structure of hierarchical, bimodal, and multiwave emulators. The minimal
#' argument determines whether we obtain the smallest (default) or largest set of ranges
#' in the collection of emulators.
#'
#' @param emulators The set of emulators (possibly as a recursive list)
#' @param minimal Whether to return the smallest (default) or largest ranges
#'
#' @return A list of paired numerics, corresponding to the parameter ranges
#'
#' @noRd
#' @keywords internal
getRanges <- function(emulators, minimal = TRUE) {
  emulators <- collect_emulators(emulators)
  if (!is.null(emulators$expectation)) emulators <- emulators$expectation
  if (!is.null(emulators$mode1))
    emulators <- c(emulators$mode1$expectation, emulators$mode2$expectation)
  range_lengths <- map_dbl(emulators, ~length(.$ranges))
  if (length(unique(range_lengths)) != 1) {
    emulators <- emulators[range_lengths == max(range_lengths)]
  }
  range_widths <- data.frame(
    do.call(
      'rbind', map(emulators, ~map(.$ranges, diff))))
  which_choose <- if (minimal)
    apply(range_widths, 2, which.min)
  else apply(range_widths, 2, which.max)
  return(setNames(map(names(range_widths),
                             ~emulators[[which_choose[[.]]]]$ranges[[.]]),
                  names(range_widths)))
}

### In progress: need to consider active variables in this function
exp_sq_helper <- function(points, data, first_points, theta) {
  orig_corr <- exp_sq(first_points, data, list(theta = theta))
  grid_indices <- which(map_dbl(seq_along(first_points),
                                       ~length(unique(first_points[,.]))) == 1)
  unique_grid_points <- unique(points[,grid_indices])
  starting_point <- unique_grid_points[1,]
  unique_grid_points <- unique_grid_points[-1,]
  diffs_list <- map(grid_indices, ~t(outer(first_points[,.], data[,.], "-")))
  all_corrs <- map(seq_len(nrow(unique_grid_points)), function(i) {
    trans_vector <- unlist(unique_grid_points[i,] - starting_point)
    modifier <- prod(map_dbl(trans_vector, ~exp(-.^2/theta^2)))
    cross_term <- Reduce("*", map(seq_along(diffs_list),
                                         ~exp(-2 * diff_list[[.]] * trans_vector[.]/theta^2)))
    return(modifier * orig_corr * cross_term)
  })
  calced_corrs <- do.call('cbind', all_corrs)
  return(cbind(orig_corr, calced_corrs))
}

#' Export Emulators
#'
#' Exports emulators to non-R format
#'
#' While having emulators saved in a native-R format (for example, as part of
#' an RData file) can be useful for fast loading, it might be prohibitive where
#' space is a consideration, or where emulators might be imported into other
#' languages. This function summarises the key features of deterministic
#' emulators, either returning the compressed details as an object or returning
#' a JSON file. If a filename is provided, then the JSON is automatically
#' saved to the relevant location.
#'
#' @importFrom jsonlite toJSON prettify minify
#'
#' @param ems The (list of) emulators to convert
#' @param inputs Any previous information from compressed emulators (mostly internal)
#' @param filename If provided, the location to save the JSON file created
#' @param output.type If filename = NULL, whether to return the JSON ("json") or object.
#'
#' @return Either the created object, or NULL if a filename is specified.
#'
#' @seealso import_emulator_from_json
#'
#' @examples
#' # Using the SIREmulators
#' ems <- SIREmulators$ems
#' single_em_json <- export_emulator_to_json(ems[[1]], output.type = "json")
#' single_em_json
#' all_em_json <- export_emulator_to_json(ems, output.type = "json")
#' # Checking with the corresponding import function
#' reconstruct_em <- import_emulator_from_json(single_em_json)
#' ems[[1]]
#' reconstruct_em
#'
#' @export
export_emulator_to_json <- function(ems, inputs = NULL, filename = NULL, output.type = "json") {
  if ("Emulator" %in% class(ems)) {
    input_ranges <- ems$ranges
    out_name <- ems$output_name
    input_uid <- NULL
    if (!is.null(ems$in_data)) {
      input_df <- cbind.data.frame(eval_funcs(scale_input, ems$in_data, input_ranges, forward = FALSE), ems$out_data) |>
        setNames(c(names(input_ranges), out_name))
      input_uid <- c(apply(input_df[,names(input_ranges)], 1, rlang::hash))
      if (!is.null(inputs$data)) {
        if (!out_name %in% names(inputs$data))
          inputs$data[,out_name] <- NA
        already_used <- input_uid %in% inputs$data$uid
        for (ind in seq_along(already_used)) {
          if (already_used[ind])
            inputs$data[inputs$data$uid == input_uid[ind], out_name] <- input_df[ind, out_name]
        }
        input_subset <- cbind.data.frame(input_uid, input_df)[!already_used,] |> setNames(c("uid", names(input_df)))
        if (nrow(input_subset) != 0) {
          for (nm in names(inputs$data)) {
            if (!nm %in% names(input_subset)) input_subset[,nm] <- NA
          }
          inputs$data <- rbind.data.frame(inputs$data, input_subset)
        }
      }
      else {
        inputs$data <- cbind.data.frame(input_uid, input_df) |> setNames(c("uid", names(input_df)))
      }
    }
    inputs <- list(em = list(in.ranges = input_ranges,
                             out.name = out_name,
                             input.uid = input_uid,
                             basis.f = ems$basis_f,
                             basis.beta = ems$beta_mu,
                             emulator.discrepancies = ems$disc,
                             corr.name = ems$corr$corr_name,
                             corr.sigma = ems$u_sigma,
                             corr.details = ems$corr$hyper_p,
                             nugget = ems$corr$nugget),
                   data = inputs$data)
  }
  else {
    if ("expectation" %in% names(ems) || "mode1" %in% names(ems))
      stop("Variance or multistate emulators are not currently supported.")
    ems <- collect_emulators(ems, ordering = "params")
    for (i in seq_along(ems)) {
      exported_em_details <- export_emulator_to_json(ems[[i]], inputs, output.type = "object")
      inputs$ems[[paste0("emulator", length(inputs$ems)+1)]] <- exported_em_details$em
      inputs$data <- exported_em_details$data
    }
  }
  if (is.null(filename)) {
    if (output.type == "json") return(toJSON(inputs, digits = NA) |> prettify())
    else return(inputs)
  }
  else {
    tryCatch(
      write(jsonlite::toJSON(inputs, digits = NA) |> minify(), file = filename),
      error = function(e) {
        cat("File could not be saved; returning JSON output.", e)
        return(toJSON(inputs, digits = NA) |> prettify())
      }
    )
    cat("Output JSON saved to", filename)
    return(NULL)
  }
}

#' Import JSON Emulator Data
#'
#' Given a file containing emulator details, reconstruct a collection of emulators.
#'
#' For data generated from \code{\link{export_emulator_to_json}} (for example),
#' emulators are recreated using the specifications therein. For each emulator,
#' a call is made to \code{\link{emulator_from_data}} with \code{specified_priors}
#' stipulated (so no retraining is performed, making the reconstruction fast).
#'
#' The structure of the JSON file used to import is relatively generic, and can
#' be created outside of this package. For examples of the structure, see the
#' code given in the companion export function.
#'
#' @importFrom jsonlite fromJSON
#'
#' @param filename Either a file location of a saved JSON file, or the string corresponding to it
#' @param details Mainly internal; any already reconstructed emulators and their input data
#'
#' @returns The emulator objects, as a list.
#'
#' @seealso export_emulator_to_json
#'
#' @export
import_emulator_from_json <- function(filename = NULL, details = NULL) {
  if (!is.null(filename)) details <- jsonlite::fromJSON(filename)
  input_data <- details$data
  if (!is.null(details[['em']])) {
    in_em_details <- details$em
    beta_funcs <- c(unlist(apply(in_em_details$basis.f, 1, function(y) {
      eval(parse(text = paste0("function(x)", sub("return\\((*.))\\)", "\\1", y[2]))))
    })))
    in_data <- input_data[input_data$uid %in% in_em_details$input.uid,]
    em <- emulator_from_data(in_data, in_em_details$out.name, in_em_details$in.ranges,
                             discrepancies = list(in_em_details$emulator.discrepancies),
                             specified_priors = list(func = list(beta_funcs), beta = list(list(mu = in_em_details$basis.beta)),
                                                     u = list(list(sigma = in_em_details$corr.sigma,
                                                                   corr = Correlator$new(in_em_details$corr.name,
                                                                                         in_em_details$corr.details,
                                                                                         nug = in_em_details$nugget)))),
                             verbose = FALSE, more_verbose = FALSE, check.ranges = FALSE)
  }
  else {
    em <- purrr::map(names(details$ems), ~import_emulator_from_json(details = list(data = input_data, em = details$ems[[.]]))[[1]]) |>
      setNames(purrr::map_chr(details$ems, "out.name"))
  }
  return(em)
}

# Pre-submission questions for CRAN submission
release_questions <- function() { # nocov start
  c(
    "Have you recompiled the vignettes using precompile.R?"
    )
} # nocov end
