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
#' @importFrom stats lm step setNames as.formula
#'
#' @param data A \code{data.frame} containing the input and output values
#' @param ranges A named list consisting of the input parameter ranges
#' @param output_name A string corresponding to the output to fit to
#' @param add Should stepwise addition or deletion be performed?
#' @param order To what polynomial order should the model by fitted?
#' @param u_form An upper form for the model fit. Used internally.
#'
#' @keywords internal
#' @noRd
#'
#' @return The fitted \code{lm} model object
#'
get_coefficient_model <- function(data, ranges, output_name, add = FALSE, order = 2, u_form = NULL) {
  lower_form <- as.formula(paste(output_name, "1", sep = " ~ "))
  if (is.null(u_form)) {
    if (order == 1)
      upper_form <- as.formula(paste(output_name, " ~ ", paste0(c('1', names(ranges)), collapse = "+"), sep = ""))
    else {
      upper_form <- as.formula(paste(output_name, " ~ ", paste(paste0(c('1', names(ranges)), collapse = "+"), paste0("I(", names(ranges), "^2)", collapse = "+"), sep = "+"), sep = ""))
      start_model <- get_coefficient_model(data = data, ranges = ranges, output_name = output_name, add = add, order = order, u_form = upper_form)
      a_vars <- names(start_model$coefficients)[-1]
      in_model <- purrr::map_lgl(names(ranges), ~any(grepl(., a_vars)))
      a_vars <- names(ranges)[in_model]
      if (length(a_vars) == 0)
        upper_form <- as.formula(paste(output_name, "~ 1"))
      else
        upper_form <- as.formula(paste(output_name, " ~ ", paste(paste("(", paste(c('1', a_vars), collapse = "+"), ")^", order, sep = ""), paste0("I(", a_vars, paste("^", order, ")", sep = ""), collapse = "+"), paste0(c(outer(a_vars, names(ranges), paste, sep = ":")), collapse = "+"), sep = "+"), sep = ""))
    }
  }
  else
    upper_form <- u_form
  if (!add & (choose(length(ranges) + order, length(ranges)) > nrow(data))) {
    warning("Maximum number of regression terms is greater than the available degrees of freedom. Changing to add = TRUE")
    add <- TRUE
  }
  #scaled_input_data <- scale_input(data[, names(ranges)], ranges)
  #full_scaled_data <- setNames(cbind(scaled_input_data, data[,output_name]), c(names(ranges), output_name))
  if (!"data.frame" %in% class(data)) data <- setNames(data.frame(data), c(names(ranges), output_name))
  if (add) {
    model <- step(lm(formula = lower_form, data = data),
                  scope = list(lower = lower_form, upper = upper_form),
                  direction = "forward", trace = 0, k = log(nrow(data)))
  }
  else {
    model <- step(lm(formula = upper_form, data = data),
                  scope = list(lower = lower_form, upper = upper_form),
                  direction = "backward", trace = 0, k = log(nrow(data)))
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
#' @importFrom dfoptim nmkb
#'
#' @param inputs The input data
#' @param outputs The output values (usually as residuals from a fitted regression)
#' @param h The basis functions of the regression surface
#' @param corr The \code{\link{Correlator}} object to use
#' @param hp_range The allowed range for the hyperparameters
#' @param beta If provided, the regression coefficients will be treated as known.
#' @param delta The value of the nugget term. If \code{NULL}, it is treated as a hyperparameter.
#'
#'
#' @keywords internal
#' @noRd
#'
#' @return A list of hyperparameter values
likelihood_estimate <- function(inputs, outputs, h, corr = Correlator$new(), hp_range, beta = NULL, delta = NULL) {
  if (!"data.frame" %in% class(inputs)) inputs <- data.frame(inputs)
  H <- t(eval_funcs(h, inputs))
  if (!is.null(beta) && length(beta) != length(h)) stop("Number of coefficients does not match number of regression functions.")
  av <- purrr::map_lgl(seq_along(names(inputs)), function(x) {
    point_vec <- c(rep(0, x-1), 1, rep(0, length(names(inputs))-x))
    func_vals <- purrr::map_dbl(h, purrr::exec, point_vec)
    sum(func_vals) > 1
  })
  if (all(av == FALSE)) av <- c(TRUE)
  corr_mat <- function(points, hp, delta) {
    this_corr <- corr$set_hyper_p(hp, delta)
    apply(points, 1, function(a) apply(points, 1, this_corr$get_corr, a, av))
  }
  func_to_opt <- function(params, log_lik = TRUE, b = beta, return_stats = FALSE) {
    hp <- params[1:length(hp_range)]
    delta <- params[length(hp_range)+1]
    if (is.na(delta)) delta <- 0
    A <- corr_mat(inputs, hp, delta)
    A_inv <- tryCatch(
      chol2inv(chol(A)),
      error = function(e) MASS::ginv(A)
    )
    if (is.null(b)) {
      b_ml <- tryCatch(
        chol2inv(chol(t(H) %*% A_inv %*% H)) %*% t(H) %*% A_inv %*% outputs,
        error = function(e) MASS::ginv(t(H) %*% A_inv %*% H) %*% t(H) %*% A_inv %*% outputs
      )
    }
    else b_ml <- b
    m_diff <- if (nrow(H) == 1) t(outputs - H * b_ml) else outputs - H %*% b_ml
    sigma_ml <- sqrt(t(m_diff) %*% A_inv %*% m_diff/length(outputs))
    if (log_lik) lik <- -length(outputs) * log(2*pi*sigma_ml^2)/2 - det(A)/2 - (t(m_diff) %*% A_inv %*% m_diff)/(2*sigma_ml^2)
    else lik <-  1/((2*pi*sigma_ml^2)^(length(outputs)/2) * sqrt(det(A))) * exp(-(t(m_diff) %*% A_inv %*% m_diff)/(2 * sigma_ml^2))
    if (return_stats) return(list(beta = b_ml, sigma = sigma_ml))
    else return(lik)
  }
  if (length(hp_range[[1]]) == 1) {
    if (is.null(delta)) delta <- 0.1
    best_hp <- hp_range
    best_delta <- delta
    best_params <- c(best_hp, best_delta)
  }
  else {
    grid_search <- expand.grid(purrr::map(hp_range, ~seq(.[[1]], .[[2]], length.out = 20)))
    grid_liks <- apply(grid_search, 1, function(x) {
      if (is.null(delta)) return(func_to_opt(c(x, 0)))
      else return(func_to_opt(c(x, delta)))
    })
    best_initial <- grid_search[which.max(grid_liks),]
    if (sum(av) == length(inputs)) delta <- 0
    if(is.null(delta)) {
      initial_params <- c(best_initial, 0.01, use.names = F)
    }
    else
      initial_params <- c(best_initial, if (delta == 0) 1e-10 else delta, use.names = F)
    best_params <- nmkb(initial_params, func_to_opt, lower = c(purrr::map_dbl(hp_range, ~.[[1]]-1e-6), if(is.null(delta)) 0 else max(0, delta - 1e-6), use.names = F), upper = c(purrr::map_dbl(hp_range, ~.[[2]]+1e-6), if(is.null(delta)) 0.25 else delta + 1e-6, use.names = F), control = list(maximize = TRUE))$par
    best_hp <- best_params[-length(best_params)]
    best_delta <- best_params[length(best_params)]
  }
  # if (all(unlist(best_hp, use.names = FALSE) - purrr::map_dbl(hp_range, ~.[[1]]) < 1e-5)) {
  #   best_hp <- purrr::map_dbl(hp_range, ~sum(.)/2)
  #   best_delta <- 0.05
  # }
  # if (sum(av) == length(inputs)) best_delta <- 0
  return(c(hp = best_hp, delta = best_delta, func_to_opt(best_params, return_stats = TRUE)))
}

#' Generate Emulators from Data
#'
#' Given data from simulator runs, generates a set of univariate \code{\link{Emulator}} objects,
#' one for each output.
#'
#' Many of the parameters that can be passed to this function are optional: the minimal operating
#' example requires \code{input_data}, \code{output_names}, and one of \code{ranges} or
#' \code{input_names}. If \code{ranges} is supplied, the input names are generated from that
#' list; if only \code{input_names} is specified, then the ranges are assumed to be [-1, 1]
#' for every input.
#'
#' If the minimum information is provided, then an emulator is fitted as follows. The basis
#' functions and associated regression coefficients are generated using \code{step} and \code{lm}
#' up to a desired order (default 2, determined by \code{quadratic}). These regression parameters
#' are assumed to be 'known' unless \code{beta.var = TRUE}, in which case the derived parameter
#' variance is taken from the model fit too (and the regression coefficients themselves can
#' be modified by the maximum likelihood estimate performed below).
#'
#' The correlation function c(x, x') is assumed to be \code{\link{exp_sq}} and a corresponding
#' \code{\link{Correlator}} object is created. The hyperparameters of the correlation structure
#' is determined using a combination of the Durham heuristic and maximum likelihood estimation.
#' This determines the variance \code{sigma^2}, correlation length \code{theta}, and nugget
#' term \code{delta}.
#'
#' If \code{ev} is provided, then the ensemble variability is taken into account in the
#' determination of the nugget term via a two-stage training process.
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
#'
#' @return A list of \code{\link{Emulator}} objects.
#' @export
#'
#' @examples
#' # Use the \code{\link{GillespieSIR}} dataset as an example.
#' ranges <- list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05))
#' out_vars <- c('nS','nI','nR')
#' ems_linear <- emulator_from_data(GillespieSIR, out_vars, ranges, quadratic = FALSE)
#' ems_linear # Printout of the key information.
#'
#' \donttest{
#'
#'   ems_quad <- emulator_from_data(GillespieSIR, out_vars, ranges)
#'   ems_quad # Now includes quadratic terms (but only where needed)
#'
#'   ems_unadjusted <- emulator_from_data(GillespieSIR, out_vars, ranges, adjusted = FALSE)
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
#'   ems2 <- emulator_from_data(GillespieSIR, out_vars, ranges,
#'                              funcs = basis_f, beta = beta_vals)
#'   ems2
#'   # Custom correlation functions
#'   corr_structs <- list(
#'     list(sigma = 83, corr = Correlator$new('exp_sq', list(theta = 0.5), nug = 0.1)),
#'     list(sigma = 95, corr = Correlator$new('exp_sq', list(theta = 0.4), nug = 0.25)),
#'     list(sigma = 164, corr = Correlator$new('exp_sq', list(theta = 0.2), nug = 0.45))
#'   )
#'   ems3 <- emulator_from_data(GillespieSIR, out_vars, ranges,
#'                              u = corr_structs)
#' }
#'
emulator_from_data <- function(input_data, output_names, ranges,
                               input_names = names(ranges), beta, u,
                               c_lengths, funcs, deltas, ev,
                               quadratic = TRUE, beta.var = FALSE,
                               adjusted = TRUE) {
  model_beta_mus <- model_u_sigmas <- model_u_corrs <- NULL
  if (missing(ranges)) {
    if (is.null(input_names)) stop("Input ranges or names of inputs must be provided.")
    warning("No ranges provided. Inputs assumed to be in ranges [-1, 1].")
    ranges <- setNames(purrr::map(input_names, ~c(-1, 1)), input_names)
  }
  data <- setNames(cbind(eval_funcs(scale_input, input_data[,names(ranges)], ranges), input_data[,output_names]), c(names(ranges), output_names))
  if (!"data.frame" %in% class(data)) data <- setNames(data.frame(data), c(names(ranges), output_names))
  if (missing(funcs)) {
    if (quadratic) {
      does_add <- (choose(length(input_names)+2, length(input_names)) > nrow(data))
      models <- purrr::map(output_names, ~get_coefficient_model(data, ranges, ., add = does_add))
    }
    else {
      does_add <- (length(input_names)+1 > nrow(data))
      models <- purrr::map(output_names, ~get_coefficient_model(data, ranges, ., add = does_add, order = 1))
    }
    all_funcs <- c(function(x) 1, purrr::map(seq_along(input_names), ~function(x) x[[.]]))
    all_coeffs <- c("(Intercept)", input_names)
    if (quadratic) {
      all_funcs <- c(all_funcs, apply(expand.grid(1:length(input_names), 1:length(input_names)), 1, function(y) function(x) x[[y[[1]]]] * x[[y[[2]]]]))
      all_coeffs <- c(all_coeffs, apply(expand.grid(input_names, input_names), 1, paste, collapse = ":"))
      all_coeffs <- sub("(.*):(\\1)$", "I(\\1^2)", all_coeffs)
    }
    model_basis_funcs <- purrr::map(models, ~all_funcs[all_coeffs %in% variable.names(.)])
    if (!beta.var) {
      model_beta_mus <- purrr::map(models, ~c(.$coefficients[all_coeffs[all_coeffs %in% names(.$coefficients)]], use.names = FALSE))
      model_beta_sigmas <- purrr::map(model_beta_mus, ~diag(0, nrow = length(.)))
    }
    else {
      model_beta_sigmas <- purrr::map(models, ~vcov(.))
    }
  }
  else {
    if (!(missing(beta) || is.null(beta[[1]]$mu))) {
      if (any(purrr::map_lgl(seq_along(beta), ~length(beta[[.]]$mu) != length(funcs[[.]])))) stop("Regression function and coefficient specifications do not match.")
      model_beta_mus <- purrr::map(beta, ~.$mu)
      if (is.null(beta[[1]]$sigma)) model_beta_sigmas <- purrr::map(beta, ~diag(0, nrow = length(.$mu)))
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
  if (any(is.null(model_beta_mus) || is.null(model_u_sigmas) || is.null(model_u_corrs))) {
    if (missing(c_lengths)) theta_ranges <- purrr::map(model_basis_funcs, ~list(theta = c(0.3, ifelse(quadratic, 1, 2))))
    else theta_ranges <- c_lengths
    specs <- purrr::map(seq_along(model_basis_funcs), ~likelihood_estimate(data[,input_names], data[,output_names[[.]]], model_basis_funcs[[.]], hp_range = theta_ranges[[.]], beta = model_beta_mus[[.]], delta = model_deltas[[.]]))
    if (is.null(model_u_sigmas)) model_u_sigmas <- purrr::map(specs, ~as.numeric(.$sigma))
    if (is.null(model_beta_mus)) model_beta_mus <- purrr::map(specs, ~.$beta)
    if(is.null(model_u_corrs)) model_u_corrs <- purrr::map(specs, ~Correlator$new()$set_hyper_p(.$hp, .$delta))
  }
  model_us <- purrr::map(seq_along(model_u_corrs), ~list(sigma = model_u_sigmas[[.]], corr = model_u_corrs[[.]]))
  model_betas <- purrr::map(seq_along(model_beta_mus), ~list(mu = model_beta_mus[[.]], sigma = model_beta_sigmas[[.]]))
  out_ems <- setNames(purrr::map(seq_along(model_us),
                                 ~Emulator$new(basis_f = model_basis_funcs[[.]], beta = model_betas[[.]], u = model_us[[.]], ranges = ranges, model = tryCatch(models[[.]], error = function(e) NULL))),
  output_names)
  if (!missing(ev)) {
    ev_deltas <- ev/purrr::map_dbl(out_ems, ~.$u_sigma)
    ev_deltas <- purrr::map_dbl(ev_deltas, ~min(1/3, .))
    return(emulator_from_data(input_data, output_names, ranges, input_names, beta, u, purrr::map_dbl(out_ems, ~.$corr$hyper_p$theta), funcs, ev_deltas, quadratic = quadratic, beta.var = beta.var))
  }
  if (!is.null(model_deltas)) {
    for (i in 1:length(model_deltas)) out_ems[[i]]$corr$nugget <- model_deltas[[i]]
  }
  for (i in 1:length(out_ems)) out_ems[[i]]$output_name <- output_names[[i]]
  if (adjusted) {
    out_ems <- purrr::map(out_ems, ~.$adjust(input_data, .$output_name))
  }
  return(out_ems)
}

#' Variance Emulator Creation
#'
#' For stochastic systems, it can be helpful to emulate the variance as well as the function.
#' This is particularly true if one expects the variance to be very different in different
#' areas of the parameter space (for example, in an epidemic model). This function performs
#' the requisite two-stage Bayes Linear update.
#'
#' Two sets of data are required: the observed variance of the stochastic runs at each point,
#' and their corresponding mean output values. The data.frames should have the same named
#' inputs and outputs. A parameter \code{npoints} is also required: this is a single numeric
#' or vector of numerics representing the number of replicates done at each input point.
#'
#' All other parameters passed to this function are equivalent to those in
#' \code{\link{emulator_from_data}} - one notable absence is that by default, the returned
#' emulators are the Bayes Linear adjusted forms.
#'
#' @param input_data_var The variance data.
#' @param input_data_exp The function data.
#' @param npoints The number of replicates per observation.
#' @param output_names The observation names.
#' @param ranges A named list of parameter ranges
#' @param input_names The names of the parameters (if \code{ranges} is not provided).
#' @param kurt The expected kurtosis of the data.
#' @param beta A list of regression coefficients for each output.
#' @param u A list of \code{\link{Correlator}} objects for each output.
#' @param c_lengths A list of correlation lengths for each output.
#' @param funcs A list of regression functions for each output.
#' @param deltas Nugget terms for each correlation structure.
#' @param ev Estimates of ensemble variability for each output.
#' @param quadratic Should a quadratic or linear fit be found?
#' @param beta.var Should regression coefficient uncertainty be included?
#'
#' @return A list of lists: one for the variance emulators and one for the function emulators.
#' @export
#'
#' @examples
#'  # Use the BirthDeath dataset
#'  ranges <- list(lambda = c(0, 0.08), mu = c(0.04, 0.13))
#'  v_ems <- variance_emulator_from_data(BirthDeath$var, BirthDeath$mean, BirthDeath$reps,
#'   paste0('t', c(1, 7, 15)), ranges)
#'  v_ems$varance
#'  v_ems$expectation
variance_emulator_from_data <- function(input_data_var, input_data_exp, npoints,
                                        output_names, ranges, input_names = names(ranges),
                                        kurt = 3, beta, u,
                                        c_lengths, funcs, deltas, ev,
                                        quadratic = TRUE, beta.var = FALSE) {
  p_v_e <- emulator_from_data(input_data_var, output_names, ranges, input_names, beta, u, c_lengths, funcs, deltas, ev, quadratic, beta.var, adjusted = FALSE)
  var_mods <- purrr::map(p_v_e, ~(.$get_exp(input_data_var)^2 + .$get_cov(input_data_var))/c(npoints) * (kurt - 1 + 2/(c(npoints)-1)))
  for (i in 1:length(p_v_e)) p_v_e[[i]]$s_diag <- c(var_mods[[i]])
  t_v_e <- setNames(purrr::map(seq_along(p_v_e), ~p_v_e[[.]]$adjust(input_data_var, output_names[[.]])), output_names)
  exp_mods <- purrr::map(t_v_e, ~.$get_exp(input_data_exp)/c(npoints))
  p_e_e <- emulator_from_data(input_data_exp, output_names, ranges, input_names, beta, u, c_lengths, funcs, deltas, ev, quadratic, beta.var, adjusted = FALSE)
  for (i in 1:length(p_e_e)) p_e_e[[i]]$s_diag <- c(exp_mods[[i]])
  t_e_e <- setNames(purrr::map(seq_along(p_e_e), ~p_e_e[[.]]$adjust(input_data_exp, output_names[[.]])), output_names)
  return(list(variance = t_v_e, expectation = t_e_e))
}
