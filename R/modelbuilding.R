## Helper function for converting ranges from data.frame or data.matrix to list
convertRanges <- function(object) {
  if (is.null(object) || missing(object)) return(object)
  if ("matrix" %in% class(object)) {
    if (ncol(object) == 2 &&
        length(row.names(object)[row.names(object)!=""]) == nrow(object))
      return(setNames(
        purrr::map(
          seq_len(nrow(object)),
          ~c(object[.,1], object[.,2], use.names = FALSE)),
        row.names(object)))
    else {
      warning(paste("Data.matrix of ranges is misspecified",
                    "(either row.names incomplete, or dimension of data.frame is not nx2)"))
      return(NULL)
    }
  }
  if (is.list(object) && !"data.frame" %in% class(object)) {
    if (length(names(object)) == length(object) &&
        all(purrr::map_dbl(object, length) == 2)) return(object)
    else {
      warning(paste("List of ranges is misspecified",
                    "(either not all named, or not all have a maximum and minimum)"))
      return(NULL)
    }
  }
  if (is.data.frame(object)) {
    if (length(object) == 2 &&
        !any(purrr::map_lgl(seq_len(nrow(object)), ~row.names(object)[.]==.)))
      return(setNames(
        purrr::map(seq_len(nrow(object)),
                   ~c(object[.,1], object[.,2])), row.names(object)))
    else {
      warning(paste("Data.frame of ranges is misspecified",
                    "(either row.names incomplete, or dimension of data.frame is not nx2)"))
      return(NULL)
    }
  }
}

#' Model Generation
#'
#' Creates a best fit of coefficients for a given data set.
#'
#' There are two ways to generate the model: we either start with all possible terms
#' (including cross-terms) up to order \code{n}, and then stepwise remove them; or we
#' start with an intercept and stepwise add terms up to order \code{n}, only retaining
#' a term if the information criterion is improved. Which method is chosen is dependent
#' on the value of \code{add} - if \code{add = FALSE} and there are not enough degrees
#' of freedom to accommodate all possible terms, a warning will be given.
#'
#' @importFrom stats lm step setNames as.formula anova
#'
#' @param data A \code{data.frame} containing the input and output values
#' @param ranges A named list consisting of the input parameter ranges
#' @param output_name A string corresponding to the output to fit to
#' @param add Should stepwise addition or deletion be performed?
#' @param order To what polynomial order should the model by fitted?
#' @param u_form An upper form for the model fit. Used internally.
#' @param printing Should the name of the output be printed?
#'
#' @keywords internal
#' @noRd
#'
#' @return The fitted \code{lm} model object
#'
get_coefficient_model <- function(data, ranges, output_name, add = FALSE,
                                  order = 2, u_form = NULL, printing = NULL) {
  if (!is.null(printing)) cat(printing, "\n") # nocov
  lower_form <- as.formula(paste(output_name, "1", sep = " ~ "))
  if (is.null(u_form)) {
    if (order == 1)
      upper_form <- as.formula(paste(output_name,
                                     " ~ ",
                                     paste0(c('1', names(ranges)),
                                            collapse = "+"),
                                     sep = ""))
    else {
      upper_form <- as.formula(
        paste(
          output_name,
          " ~ ",
          paste(paste0(c('1', names(ranges)), collapse = "+"),
                paste0("I(", names(ranges), "^2)", collapse = "+"),
                sep = "+"), sep = ""))
      start_model <- get_coefficient_model(data = data, ranges = ranges,
                                           output_name = output_name,
                                           add = add, order = order,
                                           u_form = upper_form)
      a_vars <- names(start_model$coefficients)[-1]
      in_model <- purrr::map_lgl(names(ranges), ~any(grepl(., a_vars)))
      a_vars <- names(ranges)[in_model]
      if (length(a_vars) == 0)
        upper_form <- as.formula(paste(output_name, "~ 1"))
      else
        upper_form <- as.formula(
          paste(
            output_name,
            " ~ ",
            paste(
              paste(
                "(",
                paste(c('1', a_vars), collapse = "+"), ")^", order, sep = ""),
              paste0("I(", a_vars, paste("^", order, ")", sep = ""),
                     collapse = "+"),
              paste0(c(outer(a_vars, names(ranges), paste, sep = ":")),
                     collapse = "+"),
              sep = "+"), sep = ""))
    }
  }
  else
    upper_form <- u_form
  if (!add & (choose(length(ranges) + order, length(ranges)) > nrow(data))) {
    warning(paste("Maximum number of regression terms is greater than",
                  "the available degrees of freedom. Changing to add = TRUE"))
    add <- TRUE
  }
  if (!"data.frame" %in% class(data))
    data <- setNames(data.frame(data), c(names(ranges), output_name))
  if (add) {
    model <- step(lm(formula = lower_form, data = data),
                  scope = list(lower = lower_form, upper = upper_form),
                  direction = "both", trace = 0, k = log(nrow(data)))
  }
  else {
    model <- step(lm(formula = upper_form, data = data),
                  scope = list(lower = lower_form, upper = upper_form),
                  direction = "backward", trace = 0, k = log(nrow(data)))
  }
  if (order == 2) {
    mod_coeffs <- summary(model)$coefficients[-1,]
    mod_anv <- anova(model)
    tot_sos <- sum(mod_anv$`Sum Sq`)
    quad_sos <- mod_anv[!row.names(mod_anv) %in% c(names(ranges), "Residuals"),
                        "Sum Sq"]/tot_sos
    quad_names <- row.names(mod_coeffs)[!row.names(mod_coeffs) %in%
                                          names(ranges)]
    quad_remove <- quad_names[quad_sos < 0.01]
    final_terms <- row.names(mod_coeffs)[!row.names(mod_coeffs) %in%
                                           quad_remove]
    model <- lm(data = data, formula = as.formula(
      paste(output_name, "~", paste0(c('1', final_terms), collapse = "+"))))
  }
  return(model)
}

#' Hyperparameter Estimation
#'
#' Performs hyperparameter fitting using MLE.
#'
#' The maximised regression coefficients beta and the overall SD sigma can be found in closed
#' form, given the hyperparameters of the correlation structure. Those hyperparameters are
#' estimated by first setting up a coarse grid over the main hyperparameters and finding the
#' value which maximises the likelihood. This initial guess is put into a Nelder-Mead
#' optimiser to finesse the estimate for the hyperparameters and the nugget term. Once this
#' is done, the final ML estimates of sigma and beta are calculated.
#'
#' @param inputs The input data
#' @param outputs The output values (usually as residuals from a fitted regression)
#' @param h The basis functions of the regression surface
#' @param corr_name The name of the correlation function to use.
#' @param hp_range The allowed range for the hyperparameters
#' @param beta If provided, the regression coefficients will be treated as known.
#' @param delta The value of the nugget term. If \code{NULL}, it is treated as a hyperparameter.
#' @param printing Should the output name be printed?
#'
#' @keywords internal
#' @noRd
#'
#' @return A list of hyperparameter values
likelihood_estimate <- function(inputs, outputs, h, corr_name = 'exp_sq',
                                hp_range, beta = NULL, delta = NULL,
                                nsteps = 30, printing = NULL) {
  if (!is.null(printing)) cat(printing, "\n")
  corr <- Correlator$new(corr_name,
                         hp = setNames(
                           purrr::map(names(hp_range),
                                      ~hp_range[[.]][[1]]), names(hp_range)))
  if (!"data.frame" %in% class(inputs)) inputs <- data.frame(inputs)
  H <- t(eval_funcs(h, inputs))
  if (!is.null(beta) && length(beta) != length(h))
    stop("Number of coefficients does not match number of regression functions.")
  av <- purrr::map_lgl(seq_along(names(inputs)), function(x) {
    point_vec <- c(rep(0, x-1), 1, rep(0, length(names(inputs))-x))
    func_vals <- purrr::map_dbl(h, purrr::exec, point_vec)
    sum(func_vals) > 1
  })
  if (all(av == FALSE)) av <- c(TRUE)
  corr_mat <- function(points, hp, delta) {
    this_corr <- corr$set_hyper_p(hp, delta)
    this_corr$get_corr(points, actives = av)
  }
  func_to_opt <- function(params, log_lik = TRUE,
                          b = beta, return_stats = FALSE) {
    hp <- params[seq_along(hp_range)]
    delta <- params[length(hp_range)+1]
    if (is.na(delta)) delta <- 0
    A <- corr_mat(inputs, hp, delta)
    if (suppressWarnings(is.nan(log(det(A)))) && !return_stats) return(-Inf)
    A_inv <- tryCatch(
      chol2inv(chol(A)),
      error = function(e) MASS::ginv(A)
    )
    if (is.null(b)) {
      b_ml <- tryCatch(
        chol2inv(chol(t(H) %*% A_inv %*% H)) %*% t(H) %*% A_inv %*% outputs,
        error = function(e) MASS::ginv(t(H) %*% A_inv %*% H) %*%
          t(H) %*% A_inv %*% outputs
      )
    }
    else b_ml <- b
    m_diff <- if (nrow(H) == 1)
      t(outputs - H * b_ml)
    else
      outputs - H %*% b_ml
    sigma_ml <- suppressWarnings(sqrt(t(m_diff) %*% A_inv %*% m_diff/length(outputs)))
    if (log_lik)
      lik <- -length(outputs) * log(2*pi*sigma_ml^2)/2 - log(det(A))/2 -
      (t(m_diff) %*% A_inv %*% m_diff)/(2*sigma_ml^2)
    else
      lik <-  1/((2*pi*sigma_ml^2)^(length(outputs)/2) * sqrt(det(A))) *
      exp(-(t(m_diff) %*% A_inv %*% m_diff)/(2 * sigma_ml^2))
    if (!is.finite(lik)) lik <- -Inf
    if (return_stats) return(list(beta = b_ml, sigma = sigma_ml))
    else return(lik)
  }
  if (length(hp_range[['theta']]) == 1) {
    if (is.null(delta)) delta <- 0.05
    best_point <- hp_range
    best_delta <- delta
    best_params <- c(best_point, best_delta)
  }
  else {
    if (corr_name == "matern") {
      grid_search <- expand.grid(
        theta = seq(hp_range$theta[[1]],
                    hp_range$theta[[2]],
                    length.out = nsteps),
        nu = c(0.5, 1.5, 2.5))
    }
    else {
      grid_search <- expand.grid(purrr::map(hp_range,
                                            ~seq(.[[1]], .[[2]],
                                                 length.out = nsteps)))
    }
    grid_liks <- apply(grid_search, 1, function(x) {
      if (is.null(delta)) return(func_to_opt(c(x, 0)))
      else return(func_to_opt(c(x, delta)))
    })
    maximin_distance <- max(
      apply(
        diag(Inf, nrow(inputs)) +
          as.matrix(dist(inputs, diag = TRUE, upper = TRUE)), 1, min))
    best_point <- setNames(
      as.list(c(grid_search[which.max(grid_liks),])), names(hp_range))
    if (maximin_distance < hp_range[['theta']][2])
      best_point$theta <- max(
        maximin_distance, grid_search[which.max(grid_liks), 'theta'])
    if (sum(av) == length(inputs)) delta <- 0
    if(is.null(delta))
      best_delta <- 0.05
    else
      best_delta <- ifelse (delta == 0, 1e-10, delta)
    initial_params <- unlist(c(best_point, best_delta), use.names = FALSE)
    # best_params <- tryCatch(
    #   nmkb(initial_params, func_to_opt, lower = c(purrr::map_dbl(hp_range,
    #   ~.[[1]]-1e-6), if(is.null(delta)) 0 else max(0, delta - 1e-6),
    #   use.names = F), upper = c(purrr::map_dbl(hp_range, ~.[[2]]+1e-6),
    #   if(is.null(delta)) 0.1 else delta + 1e-6, use.names = F),
    #   control = list(maximize = TRUE))$par,
    #   error = function(e) {
    #     return(initial_params)
    #   }
    # )
    #best_hp <- setNames(list(best_params[-length(best_params)]), names(hp_range))
    #best_delta <- min(best_params[length(best_params)], 0.05)
  }
  # if (all(unlist(best_hp, use.names = FALSE) - purrr::map_dbl(hp_range, ~.[[1]]) < 1e-5)) {
  #   best_hp <- purrr::map_dbl(hp_range, ~sum(.)/2)
  #   best_delta <- 0.05
  # }
  # if (sum(av) == length(inputs)) best_delta <- 0
  best_params <- unlist(c(best_point, best_delta), use.names = FALSE)
  other_pars <- func_to_opt(best_params, return_stats = TRUE)
  return(list(hp = best_point, delta = best_delta,
              sigma = other_pars$sigma, beta = other_pars$beta))
}

#' Generate Emulators from Data
#'
#' Given data from simulator runs, generates a set of univariate \code{\link{Emulator}} objects,
#' one for each output.
#'
#' Many of the parameters that can be passed to this function are optional: the minimal operating
#' example requires \code{input_data}, \code{output_names}, and one of \code{ranges} or
#' \code{input_names}. If \code{ranges} is supplied, the input names are generated from that
#' list, data.frame, or data.matrix; if only \code{input_names} is specified, then the ranges
#' are assumed to be [-1, 1] for every input.
#'
#' The ranges can be provided in alternative ways: either as a named list of length-2 numeric
#' vectors (corresponding to the maximum and minimum for each parameter); as a data.frame with
#' 2 columns where each row corresponds to a parameter; or as a data.matrix defined similarly
#' as the data.frame. In the cases where the ranges are provided as a data.frame or a data.matrix,
#' the \code{row.names} of the data object must be provided, corresponding to the names of the
#' parameters.
#'
#' If the minimum information is provided, then an emulator is fitted as follows. The basis
#' functions and associated regression coefficients are generated using \code{step} and \code{lm}
#' up to a desired order (default 2, determined by \code{quadratic}). These regression parameters
#' are assumed to be `known' unless \code{beta.var = TRUE}, in which case the derived parameter
#' variance is taken from the model fit too (and the regression coefficients themselves can
#' be modified by the maximum likelihood estimate performed below).
#'
#' The correlation function c(x, x') is assumed to be \code{\link{exp_sq}} and a corresponding
#' \code{\link{Correlator}} object is created. The hyperparameters of the correlation structure
#' are determined using a combination of maximum likelihood estimation and restriction to a
#' 'sensible' range of values, to avoid the correlation length tending to 0 or very large values.
#' This determines the variance \code{sigma^2}, correlation length \code{theta}, any other
#' hyperparameters (eg \code{nu} for the matern correlation function), and nugget
#' term \code{delta}. The hyperparameter priors can be overridden either by directly specifying
#' them using the \code{c_lengths} argument, or by supplying ranges to the \code{theta_ranges}
#' argument. Examples of this customisation can be found in the examples to this function.
#'
#' If \code{ev} is provided, then the ensemble variability is taken into account in the
#' determination of the nugget term via a two-stage training process.
#'
#' Some rudimentary data handling functionality is available but should be approached with
#' caution. The \code{na.rm} option will strip out rows of the training data that have NA values
#' in them; this of course may leave too few points to train to, and any consistent occurrence
#' of NAs in model data should be investigated. The \code{check.ranges} option allows a
#' redefinition of the ranges of the input parameters for emulator training; this is a common
#' practice in later waves in order to maximise the predictive power of the emulators, but should
#' only be used here if one is sure that the training set is representative of (and certainly
#' spanning) the full minimum enclosing hyperrectangle.
#'
#' @param input_data Required. A data.frame containing parameter and output values.
#' @param output_names Required. A character vector of output names.
#' @param ranges A named list of input parameter ranges.
#' @param input_names The names of the parameters (if \code{ranges} is not provided).
#' @param beta A list of regression coefficients for each output.
#' @param u A list of \code{\link{Correlator}} objects for each output.
#' @param c_lengths A list of correlation lengths for each output.
#' @param funcs A list of regression functions for each output.
#' @param deltas Nugget terms for each correlation structure.
#' @param ev Estimates of ensemble variability for each output.
#' @param quadratic Should a quadratic or linear fit be found?
#' @param beta.var Should regression coefficient uncertainty be included?
#' @param adjusted Are the raw emulators wanted, or Bayes Linear updated ones?
#' @param discrepancies Any internal or external discrepancies in the model.
#' @param has.hierarchy For hierarchical emulators, this will be TRUE.
#' @param verbose Should status updates be printed?
#' @param na.rm If NAs exist in the dataset, should those rows be removed?
#' @param check.ranges Should the ranges be modified in light of the data provided?
#' @param corr_name What correlation function to use. Defaults to exp_sq
#' @param targets If provided, outputs are checked for over/underestimation
#' @param ... Any additional parameters (eg for custom correlation functions)
#'
#' @return A list of \code{\link{Emulator}} objects.
#' @export
#'
#' @examples
#' # Use the SIRSample training dataset as an example.
#' ranges <- list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05))
#' out_vars <- c('nS','nI','nR')
#' ems_linear <- emulator_from_data(SIRSample$training, out_vars, ranges, quadratic = FALSE)
#' ems_linear # Printout of the key information.
#'
#' \donttest{
#'
#'   ems_quad <- emulator_from_data(SIRSample$training, out_vars, ranges)
#'   ems_quad # Now includes quadratic terms (but only where needed)
#'
#'   ems_unadjusted <- emulator_from_data(SIRSample$training, out_vars, ranges, adjusted = FALSE)
#'   ems_unadjusted # Looks the same as ems_quad, but the emulators are not BL adjusted
#'
#'   # Reproduce the linear case, but with slightly changed beta values
#'   basis_f <- list(
#'    c(function(x) 1, function(x) x[[1]], function(x) x[[2]], function(x) x[[3]]),
#'    c(function(x) 1, function(x) x[[1]], function(x) x[[2]]),
#'    c(function(x) 1, function(x) x[[1]], function(x) x[[2]])
#'   )
#'   beta_vals <- list(
#'    list(mu = c(300, -260, 220, -120)),
#'    list(mu = c(120, 110, -260)),
#'    list(mu = c(580, 160, 130))
#'   )
#'   ems2 <- emulator_from_data(SIRSample$training, out_vars, ranges,
#'                              funcs = basis_f, beta = beta_vals)
#'   ems2
#'   # Custom correlation functions
#'   corr_structs <- list(
#'     list(sigma = 83, corr = Correlator$new('exp_sq', list(theta = 0.5), nug = 0.1)),
#'     list(sigma = 95, corr = Correlator$new('exp_sq', list(theta = 0.4), nug = 0.25)),
#'     list(sigma = 164, corr = Correlator$new('exp_sq', list(theta = 0.2), nug = 0.45))
#'   )
#'   ems3 <- emulator_from_data(SIRSample$training, out_vars, ranges,
#'                              u = corr_structs)
#'   # Using alternative correlation functions and c_lengths
#'   # Allow code to choose hyperparameters
#'   ems_matern <- emulator_from_data(SIRSample$training, out_vars, ranges,
#'                                    corr_name = 'matern')
#'   # Providing hyperparameters to the function directly, via c_lengths
#'   matern_hp <- list(list(theta = 0.8, nu = 1.5), list(theta = 0.6, nu = 2.5),
#'    list(theta = 1.2, nu = 0.5))
#'   ems_matern2 <- emulator_from_data(SIRSample$training, out_vars, ranges,
#'                                     corr_name = 'matern', c_lengths = matern_hp)
#'   # If only one set of hyperparameters are provided to c_lengths, they are used for all
#'   ems_matern3 <- emulator_from_data(SIRSample$training, out_vars, ranges,
#'                                     corr_name = 'matern', c_lengths = matern_hp[[1]])
#'   # "Custom" correlation function with user-specified ranges: gamma exponential
#'   # 'gamma_exp' can be substituted for any correlation function - see Correlator documentation
#'   ems_gamma <- emulator_from_data(SIRSample$training, out_vars, ranges,
#'                                     corr_name = 'gamma_exp',
#'                                     theta_ranges = list(gamma = c(0.01, 2), theta = c(1/3, 2)))
#' }
#'
emulator_from_data <- function(input_data, output_names, ranges,
                               input_names = names(ranges), beta, u,
                               c_lengths, funcs, deltas, ev,
                               quadratic = TRUE, beta.var = FALSE,
                               adjusted = TRUE, discrepancies = NULL,
                               has.hierarchy = FALSE, verbose = interactive(),
                               na.rm = FALSE, check.ranges = FALSE,
                               corr_name = 'exp_sq', targets = NULL, ...) {
  if(!is.null(targets) &&
     length(intersect(names(targets), output_names) == length(output_names))) {
    do_preflight <- preflight(input_data, targets[output_names],
                              verbose = verbose, na.rm = na.rm)
    if (do_preflight && verbose) {
      cat("Some outputs may not be adequately emulated,", #nocov start
                "due to consistent over/underestimation in training data.\n")
      cat("Consider looking at the outputs (using, eg, behaviour_plot);",
                "some outputs may require extra runs and/or transformation.\n") #nocov end
    }
  }
  model_beta_mus <- model_u_sigmas <- model_u_corrs <- NULL
  if (missing(ranges)) {
    if (is.null(input_names))
      stop("Input ranges or names of inputs must be provided.")
    warning("No ranges provided. Inputs assumed to be in ranges [-1, 1].")
    ranges <- setNames(purrr::map(input_names, ~c(-1, 1)), input_names)
  }
  else {
    ranges <- convertRanges(ranges)
  }
  if (is.null(ranges)) stop("Ranges either not specified, or misspecified.")
  if (na.rm) input_data <- input_data[apply(
    input_data, 1, function(x) !any(is.na(x))),]
  if (check.ranges) {
    ranges <- setNames(
      purrr::map(
        names(ranges),
        ~c(
          max(
            ranges[[.]][1],
            min(input_data[,.]) - 0.05 * diff(range(input_data[,.]))),
          min(
            ranges[[.]][2],
            max(input_data[,.]) + 0.05 * diff(range(input_data[,.]))))),
      names(ranges))
  }
  if (nrow(input_data) < 10*length(ranges) && verbose)
    warning(paste("Fewer than", 10*length(ranges), #nocov start
                           "non-NA points in", length(ranges),
                           "dimensions - treat the emulated",
                           "outputs with caution, or include more",
                           "training points",
                           "(minimum 10 times the number of input parameters).")) #nocov end
  data <- setNames(
    cbind(
      eval_funcs(
        scale_input,
        input_data[,names(ranges)],
        ranges), input_data[,output_names]),
    c(names(ranges), output_names))
  if (!"data.frame" %in% class(data))
    data <- setNames(data.frame(data), c(names(ranges), output_names))
  if (is.null(list(...)[['more_verbose']]))
    more_verbose <- if (length(output_names) > 10) TRUE else FALSE
  else more_verbose <- list(...)[['more_verbose']]
  if (missing(funcs)) {
    if (verbose) cat("Fitting regression surfaces...\n") #nocov
    if (quadratic) {
      does_add <- (choose(length(input_names)+2,
                          length(input_names)) > nrow(data))
      models <- purrr::map(
        output_names,
        ~get_coefficient_model(data, ranges, ., add = does_add,
                               printing = (if(more_verbose) . else NULL)))
    }
    else {
      does_add <- (length(input_names)+1 > nrow(data))
      models <- purrr::map(
        output_names,
        ~get_coefficient_model(data, ranges, ., add = does_add, order = 1,
                               printing = (if(more_verbose) . else NULL)))
    }
    all_funcs <- c(function(x) 1,
                   purrr::map(seq_along(input_names), ~function(x) x[[.]]))
    all_coeffs <- c("(Intercept)", input_names)
    if (quadratic) {
      all_funcs <- c(all_funcs,
                     apply(
                       expand.grid(
                         seq_along(input_names),
                         seq_along(input_names)), 1,
                       function(y) function(x) x[[y[[1]]]] * x[[y[[2]]]]))
      all_coeffs <- c(all_coeffs,
                      apply(
                        expand.grid(
                          input_names,
                          input_names), 1, paste, collapse = ":"))
      all_coeffs <- sub("(.*):(\\1)$", "I(\\1^2)", all_coeffs)
    }
    model_basis_funcs <- purrr::map(
      models, ~all_funcs[match(variable.names(.), all_coeffs)])
    if (!beta.var) {
      model_beta_mus <- purrr::map(
        models,
        ~c(.$coefficients[all_coeffs[match(names(.$coefficients), all_coeffs)]],
           use.names = FALSE))
      model_beta_sigmas <- purrr::map(model_beta_mus, ~diag(0, nrow = length(.)))
    }
    else {
      model_beta_sigmas <- purrr::map(models, ~vcov(.))
    }
  }
  else {
    if (!(missing(beta) || is.null(beta[[1]]$mu))) {
      if (any(purrr::map_lgl(seq_along(beta), ~length(beta[[.]]$mu) != length(funcs[[.]]))))
        stop("Regression function and coefficient specifications do not match.")
      model_beta_mus <- purrr::map(beta, ~.$mu)
      if (is.null(beta[[1]]$sigma))
        model_beta_sigmas <- purrr::map(beta, ~diag(0, nrow = length(.$mu)))
      else model_beta_sigmas <- purrr::map(beta, ~.$sigma)
      model_basis_funcs <- funcs
    }
    else {
      model_basis_funcs <- funcs
      model_beta_sigmas <- purrr::map(model_basis_funcs, ~diag(0, nrow = length(.)))
    }
  }
  if (!missing(deltas)) model_deltas <- deltas
  else model_deltas <- NULL
  if (!(missing(u) || is.null(u[[1]]$sigma) || is.null(u[[1]]$corr))) {
    model_u_sigmas <- purrr::map(u, ~.$sigma)
    model_u_corrs <- purrr::map(u, ~.$corr)
  }
  if (verbose) cat("Building correlation structures...\n") #nocov
  if (any(is.null(model_beta_mus) ||
          is.null(model_u_sigmas) ||
          is.null(model_u_corrs))) {
    corr_func <- tryCatch(
      get(corr_name),
      error = function(e) {
        warning(paste("Can't find correlation function of type",
                    corr_name, "- reverting to exp_sq"))
        return(NULL)
      }
    )
    if (missing(c_lengths)) {
      if(is.null(corr_func)) corr_name <- 'exp_sq'
      else {
        if (!corr_name %in% c('exp_sq', 'orn_uhl', 'matern', 'rat_quad')) {
          th_ra <- list(...)[["theta_ranges"]]
          if (is.null(th_ra)) {
            warning(paste("User defined correlation function", corr_name,
                        "found but no corresponding hyperparameter ranges via",
                        "theta_ranges. Please provide -",
                        "reverting to exponential squared."))
            corr_name <- 'exp_sq'
          }
          else {
            theta_ranges <- list()
            for (i in seq_along(model_basis_funcs))
              theta_ranges[[length(theta_ranges)+1]] <- th_ra
          }
        }
      }
      if (corr_name == "exp_sq" || corr_name == "orn_uhl")
        theta_ranges <- purrr::map(
          model_basis_funcs,
          ~list(theta = c(1/3, ifelse(quadratic, 1, 2))))
      else if (corr_name == "matern")
        theta_ranges <- purrr::map(
          model_basis_funcs,
          ~list(theta = c(1/3, ifelse(quadratic, 1, 2)), nu = c(0.5, 2.5)))
      else if (corr_name == "rat_quad")
        theta_ranges <- purrr::map(
          model_basis_funcs,
          ~list(theta = c(1/3, ifelse(quadratic, 1, 2)), alpha = c(-1, 1)))
    }
    else {
      if (corr_name == "exp_sq" || corr_name == "orn_uhl") {
        if (length(c_lengths) == 1)
          theta_ranges <- purrr::map(
            seq_along(model_basis_funcs),
            ~list(theta = c_lengths))
        else theta_ranges <- purrr::map(
          seq_along(model_basis_funcs),
          ~list(theta = c_lengths[[.]]))
      }
      else {
        if (is.null(names(c_lengths)))
          theta_ranges <- purrr::map(
            seq_along(model_basis_funcs), ~c_lengths[[.]])
        else
          theta_ranges <- purrr::map(seq_along(model_basis_funcs), ~c_lengths)
      }
    }
    specs <- purrr::map(
      seq_along(model_basis_funcs),
      ~likelihood_estimate(data[,input_names],
                           data[,output_names[[.]]],
                           model_basis_funcs[[.]],
                           corr_name = corr_name,
                           hp_range = theta_ranges[[.]],
                           beta = model_beta_mus[[.]],
                           delta = model_deltas[[.]],
                           printing = (if (more_verbose) output_names[[.]] else NULL)))
    if (is.null(model_u_sigmas))
      model_u_sigmas <- purrr::map(specs, ~as.numeric(.$sigma))
    if (is.null(model_beta_mus))
      model_beta_mus <- purrr::map(specs, ~.$beta)
    if(is.null(model_u_corrs))
      model_u_corrs <- purrr::map(specs,
                                  ~Correlator$new(corr_name, hp = .$hp, nug = .$delta))
  }
  model_us <- purrr::map(
    seq_along(
      model_u_corrs),
    ~list(sigma = model_u_sigmas[[.]], corr = model_u_corrs[[.]]))
  model_betas <- purrr::map(
    seq_along(model_beta_mus),
    ~list(mu = model_beta_mus[[.]], sigma = model_beta_sigmas[[.]]))
  if (!is.null(discrepancies)) {
    if (is.numeric(discrepancies))
      discrepancies <- purrr::map(discrepancies,
                                  ~list(internal = ., external = 0))
  }
  if (verbose) cat("Creating emulators...\n") #nocov
  if (!has.hierarchy) {
    out_ems <- setNames(
      purrr::map(
        seq_along(model_us),
        ~Emulator$new(basis_f = model_basis_funcs[[.]], beta = model_betas[[.]],
                      u = model_us[[.]], ranges = ranges,
                      model = tryCatch(models[[.]], error = function(e) NULL),
                      discs = discrepancies[[.]])), output_names)
  } else {
    out_ems <- setNames(
      purrr::map(
        seq_along(model_us),
        ~HierarchicalEmulator$new(basis_f = model_basis_funcs[[.]],
                                  beta = model_betas[[.]],
                                  u = model_us[[.]], ranges = ranges,
                                  model = tryCatch(models[[.]], error = function(e) NULL),
                                  discs = discrepancies[[.]])), output_names)
  }
  if (!missing(ev)) {
    ev_deltas <- ev/purrr::map_dbl(
      out_ems, function(x) {
        if (is.numeric(x$u_sigma)) return(x$u_sigma)
        mean(apply(data[,names(ranges)], 1, x$u_sigma))
      })
    ev_deltas <- purrr::map_dbl(ev_deltas, ~min(1/3, .))
    return(emulator_from_data(input_data, output_names, ranges,
                              input_names, beta, u,
                              purrr::map_dbl(out_ems, ~.$corr$hyper_p$theta),
                              funcs, ev_deltas, quadratic = quadratic,
                              beta.var = beta.var, verbose = verbose,
                              discrepancies = discrepancies,
                              has.hierarchy = has.hierarchy))
  }
  if (!is.null(model_deltas)) {
    for (i in seq_along(model_deltas))
      out_ems[[i]]$corr$nugget <- model_deltas[[i]]
  }
  for (i in seq_along(out_ems))
    out_ems[[i]]$output_name <- output_names[[i]]
  if (adjusted) {
    if (verbose) cat("Performing Bayes linear adjustment...\n") #nocov
    out_ems <- purrr::map(out_ems, ~.$adjust(input_data, .$output_name))
  }
  return(out_ems)
}

#' Variance Emulator Creation
#'
#' Trains hierarchical emulators to stochastic systems
#'
#' For stochastic systems, one may emulate the variance as well as the function itself.
#' This is particularly true if one expects the variance to be very different in different
#' areas of the parameter space (for example, in an epidemic model). This function performs
#' the requisite two-stage Bayes Linear update.
#'
#' All observations are required (including replicates at points) - this function collects
#' them into the required chunks and calculates the summary statistics as required.
#'
#' All other parameters passed to this function are equivalent to those in
#' emulators are the Bayes Linear adjusted forms.
#'
#' @importFrom stats var
#'
#' @param input_data All model runs at all points.
#' @param output_names The observation names.
#' @param ranges A named list of parameter ranges
#' @param input_names The names of the parameters (if \code{ranges} is not provided).
#' @param verbose Should status updates be printed to console?
#' @param na.rm Should NA values be removed before training?
#' @param ... Optional parameters that can be passed to \code{link{emulator_from_data}}.
#'
#' @return A list of lists: one for the variance emulators and one for the function emulators.
#'
#' @references Goldstein & Vernon (2016) in preparation
#'
#' @examples
#'  # A simple example using the BirthDeath dataset
#'  v_ems <- variance_emulator_from_data(BirthDeath$training, c("Y"),
#'   list(lambda = c(0, 0.08), mu = c(0.04, 0.13)), c_lengths = c(0.75))
#'
#' @export
variance_emulator_from_data <- function(input_data, output_names, ranges,
                                        input_names = names(ranges),
                                        verbose = interactive(), na.rm = FALSE, ...) {
  if (na.rm) input_data <- input_data[apply(
    input_data, 1, function(x) !any(is.na(x))),]
  unique_points <- unique(input_data[, input_names])
  uids <- apply(unique_points, 1, hash)
  data_by_point <- purrr::map(uids, function(x) {
    input_data[apply(input_data[,names(ranges)], 1, hash) == x,]
  })
  data_by_point <- data_by_point[purrr::map_lgl(data_by_point, ~nrow(.)>1)]
  if (verbose) cat("Separated dataset by unique points...\n") #nocov
  collected_stats <- do.call('rbind', lapply(data_by_point, function(x) {
    n_points <- nrow(x)
    if (length(output_names) == 1) {
      means <- mean(x[,output_names])
      vars <- var(x[,output_names])
      kurts <- kurtosis(x[,output_names])
    }
    else {
      means <- apply(x[,output_names], 2, mean)
      vars <- apply(x[,output_names], 2, var)
      kurts <- apply(x[,output_names], 2, kurtosis)
    }
    vofv <- vars^2 * (kurts - 1 + 2/(n_points-1))/n_points
    kurt_relevant <- purrr::map_dbl(seq_along(kurts), function(x) {
        if(!is.nan(kurts[x]) && (vofv/vars^2)[x] <= 1) kurts[x] else NA
      })
    return(c(x[1, input_names], means, vars, kurt_relevant,
             n_points, use.names = FALSE))
  }))
  collected_df <- setNames(data.frame(collected_stats),
                           c(input_names, paste0(output_names, "mean"),
                             paste0(output_names, "var"),
                             paste0(output_names, "kurt"), "n"))
  collected_df <- data.frame(apply(collected_df, 2, unlist))
  if (length(output_names) == 1)
    collected_df_var <- collected_df[
      !is.na(collected_df[,paste0(output_names, "var")]),]
  else
    collected_df_var <- collected_df[
      apply(collected_df[,paste0(output_names, "var")], 1,
            function(a) !any(is.na(a))),]
  if (verbose) cat("Computed summary statistics...\n") #nocov
  variance_emulators <- list()
  for (i in output_names) {
    is_high_rep <- !is.na(collected_df_var[,paste0(i,"kurt")])
    all_var <- setNames(
      collected_df_var[,c(input_names, paste0(i, 'var'))], c(input_names, i))
    all_n <- collected_df_var$n
    if (all(is_high_rep)) kurt_ave <- mean(collected_df_var[,paste0(i,'kurt')])
    else if (!any(is_high_rep)) kurt_ave <- 3
    else kurt_ave <- mean(collected_df_var[is_high_rep, paste0(i, 'kurt')])
    if (all(is_high_rep) || any(is_high_rep)) {
      var_df <- setNames(
        collected_df_var[is_high_rep, c(input_names, paste0(i, "var"))],
        c(input_names, i))
      npoints <- collected_df_var[is_high_rep, 'n']
      if (sum(is_high_rep) == 1) {
        point <- all_var[is_high_rep,]
        temp_corr <- Correlator$new(hp = list(theta = 0.5))
        variance_em <- HierarchicalEmulator$new(
          basis_f = c(function(x) 1),
          beta = list(mu = c(point[1,i]), sigma = matrix(0, nrow = 1, ncol = 1)),
          u = list(sigma = point[1,i]^2, corr = temp_corr),
          ranges = ranges, out_name = i, verbose = FALSE)
      }
      else {
        variance_em <- emulator_from_data(
          var_df, i, ranges, quadratic = FALSE,
          adjusted = FALSE, has.hierarchy = TRUE, verbose = FALSE, ...)[[1]]
      }
    }
    else {
      variance_em <- emulator_from_data(
        all_var, i, ranges, quadratic = FALSE,
        adjusted = FALSE, has.hierarchy = TRUE, verbose = FALSE, ...)[[1]]
    }
    if (round(variance_em$u_sigma, 10) <= 0) {
      s_vars <- all_var[is_high_rep, i]
      s_n <- all_n[is_high_rep]
      sig_est <- s_vars^2 * (kurt_ave-1+2/(s_n-1))/s_n
      variance_em$u_sigma <- sqrt(mean(sig_est))
    }
    var_mod <- function(x, n) {
      if (n > 1)
        return(variance_em$get_exp(x)^2 +
                 variance_em$get_cov(x))/n * (kurt_ave - 1 + 2/(n-1))
      return(0)
    }
    variance_em$s_diag <- var_mod
    variance_em$em_type <- "variance"
    if (all(is_high_rep) || !any(is_high_rep) || sum(!is_high_rep) == 1) {
      variance_em$samples <- all_n
      v_em <- variance_em$adjust(all_var, i)
    }
    else {
      variance_em$samples <- all_n[!is_high_rep]
      v_em <- variance_em$adjust(
        setNames(
          collected_df[!is_high_rep, c(input_names, paste0(i, 'var'))],
          c(input_names, i)), i)
    }
    variance_emulators <- c(variance_emulators, v_em)
  }
  variance_emulators <- setNames(variance_emulators, output_names)
  if (verbose) cat("Completed variance emulators. Training mean emulators...\n") #nocov
  exp_mods <- purrr::map(variance_emulators, ~function(x, n) .$get_exp(x)/n)
  exp_data <- setNames(
    collected_df[,c(input_names, paste0(output_names, 'mean'))],
    c(input_names, output_names))
  exp_em <- emulator_from_data(
    exp_data, output_names, ranges, input_names,
    adjusted = FALSE, has.hierarchy = TRUE,
    verbose = verbose, more_verbose = FALSE, ...)
  for (i in seq_along(exp_em)) {
    exp_em[[i]]$s_diag <- exp_mods[[i]]
    exp_em[[i]]$samples <- collected_df$n
  }
  expectation_emulators <- setNames(
    purrr::map(
      seq_along(exp_em),
      ~exp_em[[.]]$adjust(exp_data, output_names[[.]])),
    output_names)
  return(list(variance = variance_emulators,
              expectation = expectation_emulators))
}

#' Bimodal Emulation
#'
#' Performs emulation of bimodal outputs and/or systems.
#'
#' In many stochastic systems, particularly disease models, the outputs exhibit bimodality - a
#' familiar example is where a disease either takes off or dies out. In these cases, it is not
#' sensible to emulate the outputs based on all realisations, and instead we should emulate each
#' mode separately.
#'
#' This function first tries to identify bimodality. If detected, it determines which of the
#' outputs in the data exhibits the bimodality: to these two separate emulators are trained, one
#' to each mode. The emulators are provided with any data that is relevant to their training; for
#' example, bimodality can exist in some regions of parameter space but not others. Points where
#' bimodality is present have their realisations allocated between the two modes while points
#' where no bimodality exists have their realisations provided to both modes. Targets that do not
#' exhibit bimodality are trained as a normal stochastic output: that is, using the default of
#' \code{\link{variance_emulator_from_data}}.
#'
#' The function also estimates the proportion of realisations in each mode for the set of outputs.
#' This value is also emulated as a deterministic emulator and included in the output.
#'
#' The output of the function is a list, containing three objects: \code{mode1}, \code{mode2}, and
#' \code{prop}. The first two objects have the form produced by \code{variance_emulator_from_data}
#' while \code{prop} has the form of an \code{emulator_from_data} output.
#'
#' @importFrom rlang hash
#' @importFrom mclust Mclust mclustBIC emControl mclust.options
#'
#' @param data The data to train emulators on (as in variance_emulator_from_data)
#' @param output_names The names of the outputs to emulate
#' @param ranges The parameter ranges
#' @param input_names The names of the parameters (by default inferred from \code{ranges})
#' @param verbose Should status updates be provided?
#' @param na.rm Should NA values be removed before training?
#' @param ... Any other parameters to pass to emulator training
#'
#' @return A list \code{(mode1, mode2, prop)} of emulator lists and objects.
#' @export
#'
#' @examples
#'  \donttest{
#'   # Use the stochastic SIR dataset
#'   SIR_ranges <- list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05))
#'   SIR_names <- c("I10", "I25", "I50", "R10", "R25", "R50")
#'   b_ems <- bimodal_emulator_from_data(SIR_stochastic$training, SIR_names, SIR_ranges)
#'  }
#'
bimodal_emulator_from_data <- function(data, output_names, ranges,
                                       input_names = names(ranges),
                                       verbose = interactive(), na.rm = FALSE, ...) {
  if (na.rm) input_data <- input_data[apply(
    input_data, 1, function(x) !any(is.na(x))),]
  unique_points <- unique(data[,input_names])
  uids <- apply(unique_points, 1, hash)
  param_sets <- purrr::map(uids, function(x) {
    data[apply(data[,input_names], 1, hash) == x,]
  })
  if(verbose) cat("Separated dataset by unique points.\n") #nocov
  modNames <- mclust.options("emModelNames")[
    !mclust.options("emModelNames") == "EEE"]
  proportion <- purrr::map_dbl(param_sets, function(x) {
    prop_clust <- Mclust(x[, output_names], G = 1:2, verbose = FALSE,
                         modelNames = modNames,
                         control = emControl(tol = 1e-3))$classification
    return(sum(prop_clust ==1)/length(prop_clust))
  })
  prop_df <- setNames(
    data.frame(cbind(unique_points, proportion)),
    c(names(unique_points), 'prop'))
  if (verbose) cat("Training emulator to proportion in modes.\n") #nocov
  prop_em <- emulator_from_data(prop_df,
                                c('prop'), ranges, verbose = FALSE, ...)
  if (verbose) cat("Performing clustering to identify modes.\n") #nocov
  has_bimodality <- setNames(
    do.call('rbind.data.frame', purrr::map(param_sets, function(x) {
    purrr::map_lgl(output_names, function(y) {
      if (length(unique(x[,y])) == 1) return(FALSE)
      return(Mclust(x[,y], G = 1:2, verbose = FALSE)$G == 2)
    })
  })), output_names)
  is_bimodal_target <- apply(
    has_bimodality, 2,
    function(x) sum(x)/length(x) >= 0.1)
  if (!any(is_bimodal_target))
    return(variance_emulator_from_data(
      data, output_names, ranges, verbose = FALSE, ...))
  if (!all(is_bimodal_target)) {
    if (verbose) cat("Training to unimodal targets.\n") #nocov
    non_bimodal <- variance_emulator_from_data(
      data, output_names[!is_bimodal_target], ranges, verbose = FALSE, ...)
  }
  else {
    if (verbose) cat("No targets appear to be unimodal.\n") #nocov
    non_bimodal <- NULL
  }
  if (verbose) cat("Training to bimodal targets.\n") #nocov
  bimodal <- purrr::map(output_names[is_bimodal_target], function(x) {
    c1_data <- list()
    c2_data <- list()
    param_bimodal <- has_bimodality[,x]
    for (i in seq_along(param_bimodal)) {
      if (!param_bimodal[[i]]) {
        c1_data[[length(c1_data)+1]] <- param_sets[[i]]
        c2_data[[length(c2_data)+1]] <- param_sets[[i]]
      }
      else {
        this_clust <- Mclust(param_sets[[i]][,x], G = 1:2,
                             verbose = FALSE)$classification
        c1_data[[length(c1_data)+1]] <- param_sets[[i]][this_clust == 1,]
        c2_data[[length(c2_data)+1]] <- param_sets[[i]][this_clust == 2,]
      }
    }
    mode1_dat <- do.call('rbind', c1_data)
    mode2_dat <- do.call('rbind', c2_data)
    m1em <- tryCatch(
      variance_emulator_from_data(mode1_dat, x, ranges, verbose = FALSE, ...),
      error = function(e) {
        cat("Problem training mode 1 emulator for target ", x, "\n", e, sep = "")
        return(list(expectation = NA, variance = NA))
      }
    )
    m2em <- tryCatch(
      variance_emulator_from_data(mode2_dat, x, ranges, verbose = FALSE, ...),
      error = function(e) {
        cat("Problem training mode 2 emulator for target ", x, "\n", e, sep = "")
        return(list(expectation = NA, variance = NA))
      }
    )
    return(list(m1 = m1em, m2 = m2em))
  })
  bimodals <- setNames(bimodal, output_names[is_bimodal_target])
  if (verbose) cat("Trained emulators. Collating.\n") #nocov
  m1exps <- m2exps <- m1vars <- m2vars <- list()
  for (i in output_names) {
    if (i %in% names(non_bimodal$expectation)) {
      m1exps <- c(m1exps, non_bimodal$expectation[[i]])
      m2exps <- c(m2exps, non_bimodal$expectation[[i]])
      m1vars <- c(m1vars, non_bimodal$variance[[i]])
      m2vars <- c(m2vars, non_bimodal$variance[[i]])
    }
    else {
      m1exps <- c(m1exps, bimodals[[i]]$m1$expectation)
      m2exps <- c(m2exps, bimodals[[i]]$m2$expectation)
      m1vars <- c(m1vars, bimodals[[i]]$m1$variance)
      m2vars <- c(m2vars, bimodals[[i]]$m2$variance)
    }
  }
  names(m1exps) <- names(m1vars) <- output_names
  names(m2exps) <- names(m2vars) <- output_names
  return(list(
    mode1 = list(expectation = m1exps, variance = m1vars),
    mode2 = list(expectation = m2exps, variance = m2vars),
    prop = prop_em[[1]]
  ))
}
