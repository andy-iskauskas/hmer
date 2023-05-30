## Helper function for converting ranges from data.frame or data.matrix to list
convertRanges <- function(object) {
  if (is.null(object) || missing(object)) return(object)
  if ("matrix" %in% class(object)) {
    if (ncol(object) == 2 &&
        length(row.names(object)[row.names(object) != ""]) == nrow(object))
      return(setNames(
        purrr::map(seq_len(nrow(object)),
                   ~c(min(object[.,]), max(object[.,]), use.names = FALSE)
        ),
        row.names(object)
      ))
    else {
      warning(paste("Data.matrix or ranges is misspecified",
                    "(either row.names incomplete, or dimension of data.frame is not nx2)."))
      return(NULL)
    }
  }
  if (is.list(object) && !"data.frame" %in% class(object)) {
    if (length(names(object)) == length(object) &&
        all(purrr::map_dbl(object, length) == 2))
      return(purrr::map(object, sort))
    else {
      warning(paste("List of ranges is misspecified,",
                    "(either not all named, or not all have maximum and minimum value)."))
      return(NULL)
    }
  }
  if (is.data.frame(object)) {
    if (length(object) == 2 &&
        !any(purrr::map_lgl(seq_len(nrow(object)), ~row.names(object)[.]==.)))
      return(purrr::map(seq_len(nrow(object)),
                        ~c(min(object[.,]), max(object[.,]))) |>
               setNames(row.names(object)))
    else {
      warning(paste("Data.frame of ranges is misspecified",
                    "(either row.names incomplete, or dimension of data.frame is not n by 2."))
      return(NULL)
    }
  }
}

#' Model Generation
#'
#' Creates a best fit of coefficients for a given data set.
#'
#' There are two ways to generate the model: either start with all possible terms
#' (including cross-terms) up to order \code{n}, and then stepwise-remove; or start
#' with an intercept and stepwise-add terms up to order \code{n}, only retaining
#' a term is the information critereon is improved. Which method is chosen is
#' dependent on the value of \code{add} - if \code{add = FALSE} and there are not
#' enough degrees of freedom to accommodate all possible terms, a warning will
#' be given.
#'
#' @importFrom stats lm step setNames as.formula anova
#'
#' @param data A \code{data.frame} containing the input and output values
#' @param ranges A named list consisting of the input parameter ranges
#' @param output_name A string corresponding to the output name to fit to
#' @param add Should stepwise-add or stepwise-delete be performed?
#' @param order The order to which a polynomial should be fitted
#' @param u_form An 'upper form' for the model fit (used internally)
#' @param verbose Should the name of the output be printed?
#'
#' @keywords internal
#' @noRd
#'
#' @return The fitted \code{lm} model object.
get_coefficient_model <- function(data, ranges, output_name, add = FALSE,
                                  order = 2, u_form = NULL, verbose = FALSE) {
  if (verbose) cat(output_name, "\n") #nocov
  lower_form <- as.formula(paste(output_name, "1", sep = "~"))
  if (is.null(u_form)) {
    if (order == 1)
      upper_form <- as.formula(paste(output_name,
                                     "~",
                                     paste0(c("1", names(ranges)), collapse = "+"),
                                            sep = ""))
    else {
      start_string <- paste(output_name, "~",
                              paste0(c("1", names(ranges)), collapse = "+"), sep = "")
      string_vec <- c(start_string)
      for (i in 2:order) {
        string_vec <-c(string_vec,
                        paste0("I(", names(ranges), "^", i, ")", collapse = "+"))
      }
      upper_form <- as.formula(paste0(string_vec, collapse = "+"))
      start_model <- get_coefficient_model(data = data, ranges = ranges,
                                           output_name = output_name, add = add,
                                           order = order, u_form = upper_form)
      a_vars <- names(start_model$coefficients)[-1]
      a_vars <- names(ranges)[purrr::map_lgl(names(ranges), ~any(grepl(., a_vars)))]
      if (length(a_vars) == 0)
        upper_form <- as.formula(paste(output_name, "~ 1"))
      else {
        start_string <- paste(c(output_name, "~",
                                paste0(c("1", a_vars), collapse = "+")), sep = "")
        for (i in 2:order) {
          string_vec <- c(string_vec,
                          paste0("I(", a_vars, "^", i, ")", collapse = "+"),
                          paste("(", paste0(c("1", a_vars), collapse = "+"), ")^", order, sep = "")
                          )
          string_vec <- c(string_vec,
                          paste0(c(outer(a_vars, names(ranges), paste, sep = ":"))))
        }
        upper_form <- as.formula(paste0(string_vec, collapse = "+"))
      }
    }
  }
  else
    upper_form <- u_form
  if (!add && verbose && choose(length(ranges) + order, length(ranges)) > nrow(data)) {
    warning(paste("Maximum number of regression terms is greater than",
                  "the available degrees of freedom. Changing to add = TRUE."))
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
  if (order != 1) {
    mod_coeffs <- summary(model)$coefficients[-1,]
    mod_anv <- anova(model)
    tot_sos <- sum(mod_anv$`Sum Sq`)
    pow_sos <- mod_anv[!row.names(mod_anv) %in% c(names(ranges), "Residuals"), "Sum Sq"]/tot_sos
    pow_names <- row.names(mod_coeffs)[!row.names(mod_coeffs) %in% names(ranges)]
    pow_remove <- pow_names[pow_sos < 0.01]
    final_terms <- row.names(mod_coeffs)[!row.names(mod_coeffs) %in% pow_remove]
    model <- lm(data = data, formula = as.formula(
      paste(output_name, "~", paste0(c("1", final_terms), collapse = "+"))
    ))
  }
return(model)
}

#' Hyperparameter Estimation
#'
#' Performs hyperparameter fitting using MLE.
#'
#' The maximised regression coefficients \code{beta} and the overall standard deviation
#' \code{sigma} can be found in closed form, given the hyperparameters of the correlatio
#' structure. Those hyperparameters are estimated by first setting up a coarse grid over
#' the main hyperparameters and finding the value which maximises the likelihood. This
#' initial guess is used as a seed for a Nelder-Mead optimiser to finesse the estimate
#' for the hyperparameters and the nugget term. Once done, the final ML estimates of beta
#' and sigma are calculated.
#'
#' @param inputs The input data
#' @param outputs The output values (usually as residuals from a fitted regression)
#' @param model The basis functions of the regression surface
#' @param corr_name The name of the type of correlation function to use
#' @param hp_range The allowed range of the hyperparameter(s)
#' @param beta If provided, the regression coefficients will be treated as known
#' @param delta The value of the nugget term. If \code{delta = NULL}, it will be treated
#' as a hyperparameter
#' @param verbose Should the output name be printed?
#'
#' @importFrom stats optim
#'
#' @keywords internal
#' @noRd
#'
#' @return A list of hyperparameter values
hyperparameter_estimate <- function(inputs, outputs, model, corr_name = "exp_sq",
                                    hp_range, beta = NULL, delta = NULL,
                                    nsteps = 30, verbose = FALSE,
                                    log_likelihood = TRUE, perform_opt = FALSE) {
  if (verbose) cat(names(outputs)[[1]], "\n")
  corr <- Correlator$new(corr_name,
                         hp = setNames(
                           purrr::map(names(hp_range),
                                      ~hp_range[[.]][[1]]), names(hp_range))
  )
  if (!"data.frame" %in% class(inputs)) inputs <- data.frame(inputs)
  if (!is.null(beta) &&
      (("lm" %in% class(model) && length(coef(model)) != length(beta)) ||
       (!"lm" %in% class(model) && length(beta) != length(model))))
    stop("Number of coefficients does not match number of regression functions")
  if ("lm" %in% class(model)) H <- model.matrix(model) #%*% diag(coef(model), nrow = length(coef(model)))
  else H <- t(eval_funcs(model, inputs))
  if (dim(H)[1] == 1) H <- t(H)
  av <- purrr::map_lgl(seq_along(names(inputs)), function(x) {
    point_vec <- c(rep(0, x-1), 1, rep(0, length(names(inputs))-x))
    if ("lm" %in% class(model))
      func_vals <- model.matrix(model$terms,
                                setNames(data.frame(matrix(c(point_vec, 0), nrow = 1)), c(names(inputs), names(outputs))))
    else
      func_vals <- purrr::map_dbl(model, purrr::exec, point_vec)
    sum(func_vals) > 1
  })
  if (all(av == FALSE)) av <- c(TRUE)
  corr_mat <- function(points, hp, delta) {
    this_corr <- corr$set_hyper_p(hp, delta)
    this_corr$get_corr(points, actives = av)
  }
  func_to_opt <- function(params, log_lik = log_likelihood, return_stats = FALSE) {
    hp <- params[seq_along(hp_range)]
    delta <- params[length(params)]
    if (is.na(delta)) delta <- 0
    A <- corr_mat(inputs, hp, delta)
    a_det <- det(A)
    if (a_det < 0) a_det <- 1e-20
    A_inv <- tryCatch(chol2inv(chol(A)), error = function(e) MASS::ginv(A))
    outputs <- purrr::map_dbl(seq_len(nrow(outputs)), ~as.numeric(outputs[.,]))
    if (is.null(beta)) {
      inv_mat <- tryCatch(chol2inv(chol(t(H) %*% A_inv %*% H)), error = function(e) MASS::ginv(t(H) %*% A_inv %*% H))
      b_ml <- inv_mat %*% t(H) %*% A_inv %*% outputs
    }
    else b_ml <- beta
    mod_diff <- if (nrow(H) == 1) t(outputs - H * b_ml) else outputs - H %*% b_ml
    sigmasq_ml <- t(mod_diff) %*% A_inv %*% mod_diff/length(outputs)
    if (log_lik)
      lik <- -length(outputs) * log(sigmasq_ml)/2 - log(a_det)/2
    else
      lik <- 1/sqrt(sigmasq_ml^length(outputs)) * 1/sqrt(a_det)
    if (is.infinite(lik)) lik <- -Inf
    if (return_stats) return(list(beta = b_ml, sigma = as.numeric(sqrt(sigmasq_ml))))
    else return(lik)
  }
  if (length(hp_range[['theta']]) == 1) {
    if (is.null(delta)) delta <- 0.05
    best_point <- hp_range
    best_delta <- delta
    best_params <- c(best_point, best_delta)
  }
  func_grad <- function(params, log_lik = log_likelihood) {
    hp <- params[seq_along(hp_range)]
    delta <- params[length(params)]
    if (is.na(delta)) delta <- 0
    A <- corr_mat(inputs, hp, delta)
    a_det <- det(A)
    if (a_det < 0) a_det <- 1e-20
    A_inv <- tryCatch(chol2inv(chol(A)), error = function(e) MASS::ginv(A))
    outputs <- purrr::map_dbl(seq_len(nrow(outputs)), ~as.numeric(outputs[.,]))
    A_diff_theta <- -2*(1-delta) * as.matrix(dist(inputs, diag = TRUE, upper = TRUE))^2/params[1]^3 * A
    A_diff_delta <- -A + diag(1, nrow = length(outputs))
    if (is.null(beta)) {
      inv_mat <- tryCatch(chol2inv(chol(t(H) %*% A_inv %*% H)), error = function(e) MASS::ginv(t(H) %*% A_inv %*% H))
      b_ml <- inv_mat %*% t(H) %*% A_inv %*% outputs
    }
    else b_ml <- beta
    mod_diff <- if (nrow(H) == 1) t(outputs - H * b_ml) else outputs - H %*% b_ml
    sigmasq_ml <- t(mod_diff) %*% A_inv %*% mod_diff/length(outputs)
    if (log_lik)
      lik <- purrr::map_dbl(list(A_diff_theta, A_diff_delta),
                            ~-(length(outputs) * (t(mod_diff) %*% A_inv %*% . %*% A_inv %*% mod_diff)/(t(mod_diff) %*% A_inv %*% mod_diff) +
                                sum(diag(A_inv %*% .)))/2
                            )
    else
      lik <- purrr::map_dbl(list(A_diff_theta, A_diff_delta),
                            ~-(t(mod_diff) %*% A_inv %*% . %*% A_inv %*% mod_diff)/(2*length(outputs)*sigmasq_ml^(length(outputs)/2 + 1) * sqrt(a_det)) -
                              sum(diag(A_inv %*% .))/(2 * sigmasq_ml^(length(outputs)/2) * sqrt(a_det)))
    if (any(is.infinite(lik))) lik <- c(0,0)
    return(lik)
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
        theta = seq(hp_range$theta[[1]], hp_range$theta[[2]], length.out = nsteps*4),
        nu = c(0.5, 1.5, 2.5)
      )
    }
    else {
      grid_search <- expand.grid(purrr::map(hp_range,
                                            ~seq(.[[1]], .[[2]], length.out = nsteps)))
    }
    grid_liks <- apply(grid_search, 1, function(x) {
      if (is.null(delta)) return(func_to_opt(c(x, 0.01)))
      else return(func_to_opt(c(x, delta)))
    })
    dists <- diag(Inf, nrow(inputs)) +
      as.matrix(dist(inputs, diag = TRUE, upper = TRUE))
    maximin_distance <- max(apply(dists, 1, min))
    best_point <- setNames(
      as.list(c(grid_search[which.max(grid_liks),])),
      names(hp_range)
    )
    if (maximin_distance < hp_range[['theta']][2])
      best_point$theta <- max(
        maximin_distance, grid_search[which.max(grid_liks), 'theta']
      )
    if (sum(av) == length(inputs)) delta <- 0
    if (is.null(delta)) best_delta <- 0.05
    else best_delta <- ifelse(delta == 0, 1e-10, delta)
    initial_params <- unlist(c(best_point, best_delta), use.names = FALSE)
  }
  if (perform_opt) {
    if (corr$corr_name == "exp_sq")
      optimise <- optim(initial_params, fn = func_to_opt,
                        gr = func_grad, method = "L-BFGS-B",
                        lower = c(purrr::map_dbl(hp_range, ~.[[1]]), 0),
                        upper = c(purrr::map_dbl(hp_range, ~.[[2]]), 0.2),
                        control = list(fnscale = -1))
    else
      optimise <- optim(initial_params, fn = func_to_opt,
                        method = "L-BFGS-B",
                        lower = c(purrr::map_dbl(hp_range, ~.[[1]]), 0),
                        upper = c(purrr::map_dbl(hp_range, ~.[[2]]), 0.2),
                        control = list(fnscale = -1))
    best_params <- optimise$par
  }
  else {
    best_params <- unlist(c(best_point, best_delta), use.names = FALSE)
  }
  other_pars <- func_to_opt(best_params, return_stats = TRUE)
  return(list(hp = setNames(as.list(best_params[-length(best_params)]), names(hp_range)),
              delta = best_params[length(best_params)],
              sigma = other_pars$sigma, beta = other_pars$beta))
}

#' Generate Emulators from Data
#'
#' Given data from simulator runs, generates a set of \code{\link{Emulator}} objects,
#' one for each output.
#'
#' Many of the parameters that can be passed to this function are optional: the minimal operating
#' example requires \code{input_data}, \code{output_names}, and one of \code{ranges} or
#' \code{input_names}. If \code{ranges} is supplied, the input names are intuited from that list,
#' data.frame, or data.matrix; if only \code{input_names} is supplied, then ranges are
#' assumed to be [-1, 1] for each input.
#'
#' The ranges can be provided in a few different ways: either as a named list of length-2
#' numeric vectors (corresponding to upper and lower bounds for each parameter); as a
#' data.frame with 2 columns and each row corresponding to a parameter; or as a data.matrix
#' defined similarly as the data.frame. In the cases where the ranges are provided as a
#' data.frame or data.matrix, the \code{row.names} of the data object must be provided, and
#' a warning will be given if not.
#'
#' If the set \code{(input_data, output_names, ranges)} is provided and nothing else,
#' then emulators are fitted as follows. The basis functions and associated regression
#' coefficients are generated using linear regression up to quadratic order, allowing for
#' cross-terms. These regression parameters are assumed 'known'.
#'
#' The correlation function c(x, x') is assumed to be \code{\link{exp_sq}} and a corresponding
#' \code{\link{Correlator}} object is created. The hyperparameters of the correlation
#' structure are determined using a constrained maximum likelihood argument. This determines
#' the variance, correlation length, and nugget term.
#'
#' The maximum allowed order of the regression coefficients is controlled by \code{order};
#' the regression coefficients themselves can be deemed uncertain by setting
#' \code{beta.var = TRUE} (in which case their values can change in the hyperparameter
#' estimation); the hyperparameter search can be overridden by specifying ranges for
#' each using \code{hp_range}.
#'
#' In the presence of expert beliefs about the structure of the emulators, information
#' can be supplied directly using the \code{specified_priors} argument. This can contain
#' specific regression coefficient values \code{beta} and regression functions \code{func},
#' correlation structures \code{u}, hyperparameter values \code{hyper_p} and nugget term
#' values \code{delta}.
#'
#' Some rudimentary data handling functionality exists, but is not a substitute for
#' sense-checking input data directly. The \code{na.rm} option will remove rows of
#' training data that contain NA values if true; the \code{check.ranges} option allows
#' a redefinition of the ranges of input parameters for emulator training if true. The
#' latter is a common practice in later waves of emulation in order to maximise the
#' predictive power of the emulators, but should only be used if it is believed that
#' the training set provided is truly representative of and spans the full space of
#' interest.
#'
#' Various different classes of emulator can be created using this function, depending
#' on the nature of the model. The \code{emulator_type} argument accepts a few different
#' options:
#'
#' \describe{
#' \item{"variance"}{Create emulators for the mean and variance surfaces, for each stochastic output}
#' \item{"covariance}{Create emulators for the mean surface, and a covariance matrix for the variance surface}
#' \item{"multistate"}{Create sets of emulators per output for multistate stochastic systems}
#' \item{"default"}{Deterministic emulators with no covariance structure}
#' }
#'
#' The "default" behaviour will apply if the \code{emulator_type} argument is not supplied, or
#' does not match any of the above options. If the data provided looks to display stochasticity,
#' but default behaviour is used, a warning will be generated and only the first model result
#' for each individual parameter set will be used in training.
#'
#' For examples of this function's usage (including optinal argument behaviour), see the examples.
#'
#' @param input_data Required. A data.frame containing parameter and output values
#' @param output_names Required. A character vector of output names
#' @param ranges Required if input_names is not given. A named list of input parameter ranges
#' @param input_names Required if ranges is not given. The names of the parameters
#' @param emulator_type Selects between deterministic, variance, covariance, and multistate emulation
#' @param specified_priors A collection of user-determined priors (see description)
#' @param order To what polynomial order should regression surfaces be fitted?
#' @param beta.var Should uncertainty in the regression coefficients be included?
#' @param corr_name If not exp_sq, the name of the correlation structures to fit
#' @param adjusted Should the return emulators be Bayes linear adjusted?
#' @param discrepancies Any known internal or external discrepancies of the model
#' @param verbose Should status updates be provided?
#' @param na.rm If TRUE, removes output values that are NA
#' @param check.ranges If TRUE, modifies ranges to a conservative minimum enclosing hyperrectangle
#' @param targets If provided, outputs are checked for consistent over/underestimation
#' @param has.hierarchy Internal - distinguishes deterministic from hierarchical emulators
#' @param covariance_opts User-specified options for emulating covariance matrices
#' @param ... Any additional parameters for custom correlators or additional verbosity options
#'
#' @importFrom rlang hash
#' @importFrom cluster daisy fanny
#' @importFrom stats coef formula model.matrix var
#'
#' @return An appropriately structured list of \code{\link{Emulator}} objects
#' @export
#'
#' @examples
#' # Deterministic: use the SIRSample training dataset as an example.
#' ranges <- list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05))
#' out_vars <- c('nS', 'nI', 'nR')
#' ems_linear <- emulator_from_data(SIRSample$training, out_vars, ranges, order = 1)
#' ems_linear # Printout of the key information.
#'
#' # Stochastic: use the BirthDeath training dataset
#' v_ems <- emulator_from_data(BirthDeath$training, c("Y"),
#'  list(lambda = c(0, 0.08), mu = c(0.04, 0.13)), emulator_type = 'variance')
#'
#' # If different specifications are wanted for variance/expectation ems, then
#' # enter a list with entries 'variance', 'expectation'. Eg corr_names
#' v_ems_corr <- emulator_from_data(BirthDeath$training, c("Y"),
#'  list(lambda = c(0, 0.08), mu = c(0.4, 0.13)), emulator_type = 'variance',
#'  corr_name = list(variance = "matern", expectation = "exp_sq")
#' )
#'
#' \donttest{
#'   ems_quad <- emulator_from_data(SIRSample$training, out_vars, ranges)
#'   ems_quad # Now includes quadratic terms
#'   ems_cub <- emulator_from_data(SIRSample$training, out_vars, ranges, order = 3)
#'   ems_cub # Up to cubic order in the parameters
#'
#'   ems_unadjusted <- emulator_from_data(SIRSample$training, out_vars, ranges, adjusted = FALSE)
#'   ems_unadjusted # Looks the same as ems_quad, but the emulators are not Bayes Linear adjusted
#'
#'   # Reproduce the linear case, but with slightly adjusted beta values
#'   basis_f <- list(
#'    c(function(x) 1, function(x) x[[1]], function(x) x[[2]]),
#'    c(function(x) 1, function(x) x[[1]], function(x) x[[2]]),
#'    c(function(x) 1, function(x) x[[1]], function(x) x[[3]])
#'   )
#'   beta_val <- list(
#'    list(mu = c(550, -400, 250)),
#'    list(mu = c(200, 200, -300)),
#'    list(mu = c(200, 200, -50))
#'   )
#'   ems_custom_beta <- emulator_from_data(SIRSample$training, out_vars, ranges,
#'    specified_priors = list(func = basis_f, beta = beta_val)
#'   )
#'   # Custom correlation functions
#'   corr_structs <- list(
#'    list(sigma = 83, corr = Correlator$new('exp_sq', list(theta = 0.5), nug = 0.1)),
#'    list(sigma = 95, corr = Correlator$new('exp_sq', list(theta = 0.4), nug = 0.25)),
#'    list(sigma = 164, corr = Correlator$new('matern', list(theta = 0.2, nu = 1.5), nug = 0.45))
#'   )
#'   ems_custom_u <- emulator_from_data(SIRSample$training, out_vars, ranges,
#'   specified_priors = list(u = corr_structs))
#'   # Allowing the function to choose hyperparameters for 'non-standard' correlation functions
#'   ems_matern <- emulator_from_data(SIRSample$training, out_vars, ranges, corr_name = 'matern')
#'   # Providing hyperparameters directly
#'   matern_hp <- list(
#'    list(theta = 0.8, nu = 1.5),
#'    list(theta = 0.6, nu = 2.5),
#'    list(theta = 1.2, nu = 0.5)
#'   )
#'   ems_matern2 <- emulator_from_data(SIRSample$training, out_vars, ranges, corr_name = 'matern',
#'    specified_priors = list(hyper_p = matern_hp))
#'   # "Custom" correaltion function with user-specified ranges: gamma exponential
#'   # Any named, defined, correlation function can be passed. See Correlator documentation
#'   ems_gamma <- emulator_from_data(SIRSample$training, out_vars, ranges, corr_name = 'gamma_exp',
#'    specified_priors = list(hyper_p = list(gamma = c(0.01, 2), theta = c(1/3, 2))))
#'
#'   # Multistate emulation: use the stochastic SIR dataset
#'   SIR_names <- c("I10", "I25", "I50", "R10", "R25", "R50")
#'   b_ems <- emulator_from_data(SIR_stochastic$training, SIR_names,
#'    ranges, emulator_type = 'multistate')
#'
#'   # Covariance emulation, with specified non-zero matrix elements
#'   which_cov <- matrix(rep(TRUE, 16), nrow = 4)
#'   which_cov[2,3] <- which_cov[3,2] <- which_cov[1,4] <- which_cov[4,1] <- FALSE
#'   c_ems <- emulator_from_data(SIR_stochastic$training, SIR_names[-c(3,6)], ranges,
#'    emulator_type = 'covariance', covariance_opts = list(matrix = which_cov))
#' }
#'
emulator_from_data <- function(input_data, output_names, ranges,
                               input_names = names(ranges), emulator_type = NULL,
                               specified_priors = NULL, order = 2, beta.var = FALSE,
                               corr_name = "exp_sq", adjusted = TRUE, discrepancies = NULL,
                               verbose = interactive(), na.rm = FALSE, check.ranges = FALSE,
                               targets = NULL, has.hierarchy = FALSE, covariance_opts = NULL, ...) {
  if (is.data.frame(input_data)) {
    ## Data cleaning and checking
    if (!all(output_names %in% names(input_data)))
      stop("output_names do not match data. Check data.frame.")
    if (!is.null(targets) && length(intersect(names(targets), output_names) == length(output_names))) {
      do_preflight <- preflight(input_data, targets[output_names], verbose = verbose, na.rm = na.rm)
      if (do_preflight && verbose) {
        cat("Some outputs may not be adequately emulated", #nocov start
            "due to consistent over/underestimation of outputs in training data.\n")
        cat("Consider looking at the outputs (e.g. using behaviour_plot);",
            "some outputs may require extra runs and/or transformation applied.\n") #nocov end
      }
    }
  }
  if (missing(ranges)) {
    if (is.null(input_names))
      stop("One of input_names or ranges must be provided.")
    warning("No ranges provided: inputs assumed to be in ranges [-1,1].")
    ranges <- setNames(purrr::map(input_names, ~c(-1,1)), input_names)
  }
  else {
    ranges <- convertRanges(ranges)
  }
  if (is.null(ranges)) stop("Ranges either not specified, or misspecified.")
  if (any(grepl("\u00a3", names(ranges)))) {
    stop("Character \u00a3 not permitted in input names - please rename.")
  }
  if (is.null(emulator_type) || !is.character(emulator_type))
    emulator_type <- "default"
  if (!emulator_type %in% c("default", "variance", "covariance", "multistate"))
    emulator_type <- "default"
  if (emulator_type == "multistate")
    input_data_backup <- input_data
  if (is.data.frame(input_data)) {
    unique_hash <- apply(unique(input_data[, input_names, drop = FALSE]), 1, hash)
    data_by_point <- purrr::map(unique_hash, function(x) {
      input_data[apply(input_data[,names(ranges),drop = FALSE], 1, hash) == x,]
    })
    if (any(purrr::map_dbl(data_by_point, nrow) > 1) && emulator_type == "default") {
      if (FALSE) {
        ## Hidden for now. Might change this at some point!
        # if (verbose) {
        want_variance <- readline(paste("emulator_type is default but multiple runs provided for some parameter sets.\n",
                                        "If variance or covariance was wanted, enter 'variance' or 'covariance'.\n",
                                        "Else hit RETURN.\n "))
        if (want_variance == "variance" || want_variance == "covariance") {
          cat("Changed\n")
          emulator_type <- want_variance
        }
        else
          input_data <- do.call('rbind.data.frame', purrr::map(data_by_point, ~.[1,,drop=FALSE]))
      }
      else {
        warning("emulator_type is default but multiple runs provided for some parameter sets.
            Training to the mean of realisations.")
        input_data <- do.call('rbind.data.frame', purrr::map(data_by_point, ~c(.[1,names(ranges)], apply(.[,output_names,drop=FALSE], 2, mean)))) |>
          setNames(names(input_data))
        variabilities <- do.call('rbind.data.frame', purrr::map(data_by_point, ~apply(.[,output_names,drop=FALSE], 2, var))) |> setNames(output_names)
        ave_vars <- apply(variabilities, 2, mean)
        if (is.null(discrepancies))
          discrepancies <- purrr::map(ave_vars, ~list(external = 0, internal = sqrt(.))) |> setNames(output_names)
        else {
          for (i in output_names) {
            if (is.null(discrepancies[[i]])) discrepancies[[i]] <- list(external = 0, internal = sqrt(ave_vars[[i]]))
            else if (is.null(discrepancies[[i]]$internal)) discrepancies[[i]]$internal <- sqrt(ave_vars[[i]])
            else discrepancies[[i]]$internal <- sqrt(discrepancies[[i]]$internal^2 + ave_vars[[i]])
          }
        }
      }
    }
    if (all(purrr::map_dbl(data_by_point, nrow) == 1) && emulator_type != "default") {
      warning("emulator_type is not default but only one model run per point.
            Changing to emulator_type = 'default'.")
      emulator_type <- "default"
    }
    if (na.rm) {
      input_data <- purrr::map(output_names, function(o) {
        new_df <- input_data[,c(names(ranges), o)]
        return(new_df[apply(new_df, 1, function(x) !any(is.na(x))),])
      }) |> setNames(output_names)
    }
    else
      input_data <- purrr::map(output_names, ~input_data[,c(names(ranges), .)])
  }
  if (check.ranges) {
    ranges <- setNames(
      purrr::map(
        names(ranges), function(nm) {
          all_pts <- do.call('c', purrr::map(input_data, ~.[,nm]))
          c(
            max(ranges[[nm]][1], min(all_pts)-0.05*diff(range(all_pts))),
            min(ranges[[nm]][2], max(all_pts)+0.05*diff(range(all_pts)))
          )
        }
      ),
    names(ranges))
  }
  not_enough_points <- purrr::map_lgl(input_data, ~nrow(.) < 10*length(ranges))
  if (verbose && any(not_enough_points)) {
    not_enough_output <- output_names[which(not_enough_points)]
    print_message <- paste("Fewer than", 10*length(ranges), #nocov start
                           "valid points in", length(ranges),
                           "dimensions for outputs", paste0(not_enough_output, collapse = ", "),
                           "- treat these emulated outputs with caution or include more",
                           "training points (min 10 times number of input parameters)." #nocov end
                           )
  }
  data <- purrr::map(input_data, function(dat) {
    temp_names <- names(dat)
    cbind.data.frame(eval_funcs(scale_input, dat[,names(ranges)], ranges), dat[,length(dat)]) |>
      setNames(temp_names)
  })
  if (length(ranges) != length(input_names)) {
    input_names <- intersect(names(ranges), input_names)
    ranges <- ranges[input_names]
  }
  if (is.null(list(...)[['more_verbose']]))
    more_verbose <- (length(output_names) > 10)
  else
    more_verbose <- list(...)[['more_verbose']]
  model_beta_mus <- model_u_sigmas <- model_u_corrs <- NULL
  if (emulator_type == "default") {
    if (is.null(specified_priors$func)) {
      if (verbose) cat("Fitting regression surfaces...\n") #nocov
      models <- purrr::map(data, ~get_coefficient_model(., ranges, names(.)[length(.)],
                                                        order = order, verbose = more_verbose))
      model_beta_mus <- purrr::map(models, coef)
      model_basis_funcs <- purrr::map(models, function(m) {
        purrr::map(names(m$coefficients), name_to_function, names(ranges))
      })
      if (!beta.var)
        model_beta_sigmas <- purrr::map(models, ~diag(0, nrow = length(coef(.))))
      else
        model_beta_sigmas <- purrr::map(models, ~vcov(.))
    }
    else {
      models <- NULL
      if (!(is.null(specified_priors$beta) ||
            any(purrr::map_lgl(specified_priors$beta, ~is.null(.$mu))))) {
        beta_priors <- specified_priors$beta
        if (any(purrr::map_lgl(seq_along(beta_priors), ~length(beta_priors[[.]]$mu) != length(specified_priors$func[[.]]))))
          stop("Provided regression function and coefficient specifications do not match.")
        model_beta_mus <- purrr::map(beta_priors, "mu")
        if (!any(purrr::map_lgl(beta_priors, ~is.null(.$sigma))))
          model_beta_sigmas <- purrr::map(beta_priors, "sigma")
        else
          model_beta_sigmas <- purrr::map(beta_priors, function(bp) {
            if (is.null(bp$sigma) || nrow(bp$sigma) != length(bp$mu) || ncol(bp$sigma) != length(bp$mu))
              return(diag(0, nrow = length(bp$mu)))
            return(bp$sigma)
          })
        model_basis_funcs <- specified_priors$func
      }
      else {
        model_basis_funcs <- specified_priors$func
        model_beta_sigmas <- purrr::map(model_basis_funcs, ~diag(0, nrow = length(.)))
      }
    }
    if (!is.null(specified_priors$delta)) model_deltas <- specified_priors$delta
    else model_deltas <- NULL
    if (!(is.null(specified_priors$u) || any(purrr::map_lgl(specified_priors$u, ~is.null(.$sigma) || is.null(.$corr))))) {
      model_u_sigmas <- purrr::map(specified_priors$u, "sigma")
      model_u_corrs <- purrr::map(specified_priors$u, "corr")
    }
    if (verbose) cat("Building correlation structures...\n") #nocov
    if (is.null(model_beta_mus) || is.null(model_u_sigmas) || is.null(model_u_corrs)) {
      corr_func <- tryCatch(
        get(corr_name),
        error = function(e) {
          warning(paste("Can't find correlation function of type",
                        corr_name, "- reverting to 'exp_sq'"))
          return(NULL)
        }
      )
      if (is.null(specified_priors$hyper_p)) {
        if (is.null(corr_func)) corr_name <- "exp_sq"
        else {
          if (!corr_name %in% c("exp_sq", "orn_uhl", "matern", "rat_quad")) {
            th_ra <- list(...)[["theta_ranges"]]
            if (is.null(th_ra)) {
              warning(paste("User-defined correlation function", corr_name,
                            "found but no corresponding hyperparameter ranges
                            via theta_ranges. Reverting to 'exp_sq'"))
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
          theta_ranges <- purrr::map(model_basis_funcs,
                                     ~list(theta = c(min(2/(order+1), 1/3), 2/order)))
        else if (corr_name == "matern")
          theta_ranges <- purrr::map(model_basis_funcs,
                                     ~list(theta = c(min(2/(order+1), 1/3), 2/order),
                                           nu = c(0.5, 2.5)))
        else if (corr_name == "rat_quad")
          theta_ranges <- purrr::map(model_basis_funcs,
                                     ~list(theta = c(min(2/(order+1), 1/3), 2/order),
                                           alpha = c(-1, 1)))
      }
      else {
        if (corr_name == "exp_sq" || corr_name == "orn_uhl") {
          if (length(specified_priors$hyper_p) == 1)
            theta_ranges <- purrr::map(model_basis_funcs,
                                       ~list(theta = specified_priors$hyper_p))
          else
            theta_ranges <- purrr::map(seq_along(model_basis_funcs),
                                       ~list(theta = specified_priors$hyper_p[[.]]))
        }
        else {
          if (is.null(names(specified_priors$hyper_p)))
            theta_ranges <- purrr::map(seq_along(model_basis_funcs),
                                       ~specified_priors$hyper_p[[.]])
          else
            theta_ranges <- purrr::map(seq_along(model_basis_funcs),
                                       ~specified_priors$hyper_p)
        }
      }
      if (is.null(models))
        specs <- purrr::map(seq_along(model_basis_funcs),
                            ~hyperparameter_estimate(
                              data[[.]][,input_names],
                              data[[.]][,length(data),drop=FALSE],
                              model_basis_funcs[[.]],
                              corr_name = corr_name,
                              hp_range = theta_ranges[[.]],
                              beta = model_beta_mus[[.]],
                              delta = model_deltas,
                              verbose = more_verbose
                            ))
      else
        specs <- purrr::map(seq_along(model_basis_funcs),
                            ~hyperparameter_estimate(
                              data[[.]][,input_names,drop=FALSE],
                              data[[.]][,length(data[[.]]),drop=FALSE],
                              models[[.]],
                              corr_name = corr_name,
                              hp_range = theta_ranges[[.]],
                              beta = model_beta_mus[[.]],
                              delta = model_deltas,
                              verbose = more_verbose
                            ))
      if (is.null(model_u_sigmas)) model_u_sigmas <- purrr::map(specs, ~as.numeric(.$sigma))
      if (is.null(model_beta_mus)) model_beta_mus <- purrr::map(specs, ~.$beta)
      if (is.null(model_u_corrs)) model_u_corrs <- purrr::map(specs, ~Correlator$new(corr_name, hp = .$hp, nug = .$delta))
    }
    model_us <- purrr::map(seq_along(model_u_corrs),
                           ~list(sigma = model_u_sigmas[[.]], corr = model_u_corrs[[.]]))
    model_betas <- purrr::map(seq_along(model_beta_mus),
                              ~list(mu = model_beta_mus[[.]], sigma = model_beta_sigmas[[.]]))
    if (!is.null(discrepancies)) {
      if (is.numeric(discrepancies)) discrepancies <- purrr::map(discrepancies,
                                                                 ~list(internal = ., external = 0))
    }
    else {
      discrepancies <- purrr::map(model_betas, ~list(internal = 0, external = 0))
    }
    if (verbose) cat("Creating emulators...\n") #nocov
    if (!has.hierarchy) {
      out_ems <- setNames(purrr::map(seq_along(model_us),
                                     ~Emulator$new(basis_f = model_basis_funcs[[.]],
                                                   beta = model_betas[[.]],
                                                   u = model_us[[.]],
                                                   ranges = ranges,
                                                   model = tryCatch(models[[.]], error = function(e) NULL),
                                                   discs = discrepancies[[.]],
                                                   )),
                          output_names)
    } else {
      out_ems <- setNames(purrr::map(seq_along(model_us),
                                     ~HierarchicalEmulator$new(basis_f = model_basis_funcs[[.]],
                                                               beta = model_betas[[.]],
                                                               u = model_us[[.]],
                                                               ranges = ranges,
                                                               model = tryCatch(models[[.]], error = function(e) NULL),
                                                               discs = discrepancies[[.]])),
                          output_names)
    }
    if (!is.null(model_deltas))
      for (i in seq_along(model_deltas)) out_ems[[i]]$corr$nugget <- model_deltas[[i]]
    for (i in seq_along(out_ems)) out_ems[[i]]$output_name <- output_names[[i]]
    if (adjusted) {
      if (verbose) cat("Performing Bayes linear adjustment...\n") #nocov
      out_ems <- purrr::map(seq_along(out_ems), ~out_ems[[.]]$adjust(input_data[[.]], out_ems[[.]]$output_name))
    }
    return(setNames(out_ems, output_names))
  }
  ## Variance emulation starts here
  if (verbose) cat("Multiple model runs per point detected. Splitting by input parameters...\n") #nocov
  unique_uids <- Reduce(union, purrr::map(input_data, function(dat) {
    c(apply(unique(dat[,input_names]), 1, hash), use.names = FALSE)
  }))
  data_by_point <- purrr::map(input_data, function(dat) {
    unique_hash <- apply(unique(dat[,input_names]), 1, hash)
    purrr::map(unique_hash, function(x) {
      dat[apply(dat[,names(ranges)], 1, hash) == x,]
    })
  })
  data_by_point <- purrr::map(data_by_point, function(dat) {
    dat[purrr::map_lgl(dat, ~nrow(.) > 1)]
  })
  param_sets <- purrr::map(unique_uids, function(uid) {
    output_vals <- purrr::map(data_by_point, function(d) {
      hashes <- purrr::map_chr(d, ~unique(apply(.[,input_names], 1, hash)))
      which_matches <- which(hashes == uid)
      if (length(which_matches) == 0) return(NULL)
      return(c(d[[which_matches[1]]][,length(d[[which_matches[1]]])]))
    })
    largest_length <- max(purrr::map_dbl(output_vals, length))
    output_vals_padded <- purrr::map(output_vals, ~c(., rep(NA, largest_length-length(.))))
    which_matches <- which(apply(do.call('rbind.data.frame', purrr::map(input_data, ~.[,input_names])), 1, hash) == uid)[1]
    which_point <- do.call('rbind.data.frame', purrr::map(input_data, ~.[,input_names]))[which_matches, input_names]
    return(cbind.data.frame(which_point[rep(1,largest_length),], do.call('cbind.data.frame', output_vals_padded)) |>
             setNames(c(input_names, output_names)))
  })
  if (verbose) cat("Separated dataset by unique points...\n") #nocov
  ## Multistate emulation goes here, so that if it isn't multistate we can go to variance emulation
  if (emulator_type == "multistate") {
    proportion <- purrr::map_dbl(param_sets, function(x) {
      p_clust <- suppressWarnings(fanny(suppressWarnings(daisy(x[,output_names, drop = FALSE])), k = 2)$clustering)
      return(sum(p_clust == 1)/length(p_clust))
    })
    unique_points <- unique(input_data_backup[,names(ranges)])
    prop_df <- cbind.data.frame(unique_points, proportion) |> setNames(c(names(unique_points), "prop"))
    has_bimodality <- do.call('rbind.data.frame', purrr::map(param_sets, function(x) {
      purrr::map_lgl(output_names, function(y) {
        if (length(unique(x[,y])) == 1) return(FALSE)
        clust1 <- suppressWarnings(fanny(suppressWarnings(daisy(x[,y,drop=FALSE])), k = 1))
        clust2 <- suppressWarnings(fanny(suppressWarnings(daisy(x[,y,drop=FALSE])), k = 2))
        if (clust1$objective[["objective"]] < clust2$objective[["objective"]]) return(FALSE)
        return(TRUE)
      })
    })) |> setNames(output_names)
    is_bimodal_target <- apply(has_bimodality, 2, function(x) sum(x)/length(x) >= 0.1)
    if (!any(is_bimodal_target)) {
      if (verbose) cat("No targets appear to be multistate. Reverting to variance emulation.\n") #nocov
      emulator_type <- "variance"
    }
    if (emulator_type == "multistate") {
      if (verbose) cat("Training an emulator to proportion in each mode.\n") #nocov
      prop_em <- emulator_from_data(prop_df, c('prop'), ranges, verbose = FALSE,
                                    specified_priors = specified_priors$prop, order = order,
                                    beta.var = beta.var, corr_name = corr_name, adjusted = adjusted,
                                    discrepancies = discrepancies$prop, na.rm = na.rm, check.ranges = check.ranges,
                                    targets = targets$prop, has.hierarchy = FALSE, covariance_opts = NULL, ...)
      if (any(!is_bimodal_target)) {
        if (verbose) cat("Training to single-state targets.\n") #nocov
        non_bimodal <- emulator_from_data(input_data_backup, output_names[!is_bimodal_target], ranges, verbose = FALSE,
                                          specified_priors = specified_priors, order = order, emulator_type = "variance",
                                          beta.var = beta.var, corr_name = corr_name, adjusted = adjusted,
                                          discrepancies = discrepancies, na.rm = na.rm, check.ranges = check.ranges,
                                          targets = targets, has.hierarchy = TRUE, covariance_opts = covariance_opts, ...)

      }
      else {
        if (verbose) {
          cat("No targets appear to be single-state.\n") #nocov
        }
        non_bimodal <- NULL
      }
      if (verbose) cat("Training to multistate targets.\n") #nocov
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
            this_clust <- fanny(suppressWarnings(daisy(param_sets[[i]][,x,drop=FALSE])), k=2)$clustering
            c1_data[[length(c1_data)+1]] <- param_sets[[i]][this_clust == 1,]
            c2_data[[length(c2_data)+1]] <- param_sets[[i]][this_clust == 2,]
          }
        }
        mode1_dat <- do.call('rbind.data.frame', c1_data)[,c(names(ranges), x)]
        mode2_dat <- do.call('rbind.data.frame', c2_data)[,c(names(ranges), x)]
        return(list(mode1 = mode1_dat, mode2 = mode2_dat))
      }) |> setNames(output_names[is_bimodal_target])
      mode1_dats <- purrr::map(bimodal, "mode1") |> setNames(output_names[is_bimodal_target])
      mode2_dats <- purrr::map(bimodal, "mode2") |> setNames(output_names[is_bimodal_target])
      mode1_ems <- emulator_from_data(mode1_dats, names(mode1_dats), ranges, verbose = FALSE, emulator_type = "variance",
                                      specified_priors = specified_priors, order = order, beta.var = beta.var,
                                      corr_name = corr_name, adjusted = adjusted, discrepancies = discrepancies,
                                      na.rm = na.rm, check.ranges = check.ranges, targets = targets, has.hierarchy = TRUE,
                                      covariance_opts = covariance_opts, more_verbose = FALSE, ...)
      mode2_ems <- emulator_from_data(mode2_dats, names(mode2_dats), ranges, verbose = FALSE, emulator_type = "variance",
                                      specified_priors = specified_priors, order = order, beta.var = beta.var,
                                      corr_name = corr_name, adjusted = adjusted, discrepancies = discrepancies,
                                      na.rm = na.rm, check.ranges = check.ranges, targets = targets, has.hierarchy = TRUE,
                                      covariance_opts = covariance_opts, more_verbose = FALSE, ...)
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
          m1exps <- c(m1exps, mode1_ems$expectation[[i]])
          m2exps <- c(m2exps, mode2_ems$expectation[[i]])
          m1vars <- c(m1vars, mode1_ems$variance[[i]])
          m2vars <- c(m2vars, mode2_ems$variance[[i]])
        }
      }
      names(m1exps) <- names(m1vars) <-
        names(m2exps) <- names(m2vars) <- output_names
      return(list(
        mode1 = list(expectation = m1exps, variance = m1vars),
        mode2 = list(expectation = m2exps, variance = m2vars),
        prop = prop_em[[1]]
      ))
    }
  }
  if (length(output_names) == 1 && emulator_type == "covariance") {
    warning("Covariance emulation selected, but only one output. Changing to emulator_type = 'variance")
    emulator_type <- "variance"
  }
  if (emulator_type == "variance") {
    ## This is a bit of a cludge to deal with bimodal emulators where
    ## one parameter doesn't contribute to an output
    data_by_point <- purrr::map(data_by_point, function(dat) {
      dat_out_name <- names(dat)[length(names(dat))]
      unique_points_in_dat <- purrr::map_chr(dat, function(subdat) {
        apply(subdat[,names(ranges)], 1, hash)[1]
      })
      which_missing <- which(!unique_uids %in% unique_points_in_dat)
      if (length(which_missing) == 0) return(dat)
      purrr::map(seq_along(unique_uids), function(i) {
        if (!i %in% which_missing) return(dat[[which(unique_points_in_dat == unique_uids[i])]])
        else {
          ## This ALMOST works. But I think it's having downstream effects.
          fake_df <- data.frame(matrix(
            c(rep(0, length(ranges)), NA, rep(0, length(ranges)), NA),
            nrow = 2, byrow = TRUE
          )) |> setNames(c(names(ranges), dat_out_name))
          return(fake_df)
        }
      })
    })
    collected_stats <- purrr::map(data_by_point, function(dat) {
      lapply(dat, function(x) {
        n_points <- nrow(x)
        out_vals <- x[,length(x)]
        mu <- mean(out_vals, na.rm = TRUE)
        sigsq <- var(out_vals, na.rm = TRUE)
        kurt <- kurtosis(out_vals, na.rm = TRUE)
        vofv <- sigsq^2 * (kurt - 1 + 2/(n_points-1))/n_points
        kurt_rel <- ifelse(!is.nan(kurt) && !is.na(kurt) && (vofv/sigsq^2) <= 1, kurt, NA)
        return(list(
          point = as.numeric(x[1,input_names], use.names = FALSE),
          mean = mu,
          var = sigsq,
          kurt = kurt_rel,
          np = n_points
        ))
      })
    })
    collected_df_pts <- do.call('rbind.data.frame', purrr::map(collected_stats[[1]], "point")) |>
      setNames(input_names)
    collected_df_var <- cbind.data.frame(
      collected_df_pts,
      do.call('cbind.data.frame', purrr::map(collected_stats, ~purrr::map_dbl(., "var")))
    ) |>
      setNames(c(input_names, output_names))
    collected_df_mean <- cbind.data.frame(
      collected_df_pts,
      do.call('cbind.data.frame', purrr::map(collected_stats, ~purrr::map_dbl(., "mean")))
    ) |>
      setNames(c(input_names, output_names))
    collected_df_kurt <- cbind.data.frame(
      collected_df_pts,
      do.call('cbind.data.frame', purrr::map(collected_stats, ~purrr::map_dbl(., "kurt")))
    ) |>
      setNames(c(input_names, output_names))
    collected_df_kurt <- cbind.data.frame(collected_df_kurt,
                                          do.call('cbind.data.frame', purrr::map(collected_stats, ~purrr::map_dbl(., "np")))) |>
      setNames(c(input_names, output_names, paste0(output_names, "_n")))
  }
  if (emulator_type == "covariance") {
    cov_out_names <- outer(output_names, output_names, paste0)
    data_by_point <- purrr::map(seq_along(data_by_point[[1]]), function(i) {
      cbind.data.frame(data_by_point[[1]][[i]][,names(ranges)],
                       do.call('cbind.data.frame', purrr::map(data_by_point, ~.[[i]][,length(.[[i]])]))) |>
        setNames(c(names(ranges), output_names))
    })
    collected_stats <- lapply(data_by_point, function(x) {
      n_points <- nrow(x)
      means <- apply(x[,output_names], 2, mean)
      covs <- cov(x[,output_names])
      kurts <- apply(x[,output_names], 2, kurtosis)
      vofv <- diag(covs)^2 * (kurts - 1 + 2/(n_points-1))/n_points
      kurt_relev <- purrr::map_dbl(seq_along(kurts), function(y) {
        if (!is.nan(kurts[y]) && (vofv/diag(covs)^2)[y] <= 1) kurts[y] else NA
      })
      if (!is.null(covariance_opts$logged) && covariance_opts$logged) {
        c_est <- eigen(covs)
        c_vals <- c_est$values
        c_vals[c_vals <= 0] <- 1e-8
        c_vals <- log(c_vals)
        covs <- c_est$vectors %*% diag(c_vals) %*% t(c_est$vectors)
      }
      return(
        list(
          point = as.numeric(x[1,input_names], use.names = FALSE),
          means = means,
          covs = covs,
          kurt = kurt_relev,
          np = n_points
        )
      )
    })
    collected_df_var <- do.call('rbind.data.frame',
                                purrr::map(collected_stats, function(x) {
                                  c(
                                    x$point, diag(x$covs)
                                  )
                                })) |>
      setNames(c(input_names, output_names))
    collected_df_cov <- do.call('rbind.data.frame',
                                purrr::map(collected_stats, function(x) {
                                  c(
                                    x$point, c(x$covs[upper.tri(x$covs)])
                                  )
                                })) |>
      setNames(c(input_names, cov_out_names[upper.tri(cov_out_names)]))
    collected_df_mean <- do.call('rbind.data.frame',
                                 purrr::map(collected_stats, function(x) {
                                   c(
                                     x$point, x$means
                                   )
                                 })) |>
      setNames(c(input_names, output_names))
    collected_df_kurt <- do.call('rbind.data.frame',
                                  purrr::map(collected_stats, function(x) {
                                    c(
                                      x$point, x$kurt, x$np
                                    )
                                  })) |>
      setNames(c(input_names, output_names, "n"))
  }
  if (verbose) cat("Computed summary statistics...\n") #nocov
  if (verbose) cat("Building variance emulator priors...\n") #nocov
  which_high_rep <- do.call('cbind.data.frame',
                            purrr::map(output_names,
                                       ~!is.na(collected_df_kurt[,.]))) |>
    setNames(output_names)
  which_high_rep[!which_high_rep] <- NA
  collected_df_var_model <- collected_df_var
  collected_df_var_model[,output_names] <- collected_df_var_model[,output_names] * which_high_rep
  if (!is.character(corr_name) && !is.null(corr_name$variance)) corr_name_var <- corr_name$variance
  else corr_name_var <- corr_name
  variance_emulators <- emulator_from_data(collected_df_var_model, output_names, ranges,
                                           input_names = input_names, emulator_type = "default",
                                           specified_priors = specified_priors$variance,
                                           order = max(order-1, 1), beta.var = beta.var,
                                           corr_name = corr_name_var, adjusted = FALSE,
                                           discrepancies = discrepancies$variance,
                                           verbose = FALSE, na.rm = TRUE, check.ranges = FALSE,
                                           targets = targets$variance, has.hierarchy = TRUE, ...)
  kurt_aves <- apply(collected_df_kurt[,output_names,drop=FALSE], 2, function(x) {
    if (sum(x[!is.na(x)]) < 2) return(3)
    return(mean(x[!is.na(x)]))
  })
  trained_var_ems <- purrr::map(names(variance_emulators), function(v_name) {
    if (round(variance_emulators[[v_name]]$u_sigma, 10) <= 0) {
      s_vars <- collected_df_var[!is.na(collected_df_kurt[,v_name]), v_name]
      if (emulator_type == "covariance")
        s_n <- collected_df_kurt[!is.na(collected_df_kurt[,v_name]), "n"]
      else
        s_n <- collected_df_kurt[!is.na(collected_df_kurt[,v_name]), paste0(v_name, "_n")]
      sig_est <- s_vars^2 * (kurt_aves[v_name]-1+2/(s_n-1))/s_n
      variance_emulators[[v_name]]$u_sigma <- sqrt(mean(sig_est))
    }
    var_mod <- function(x, n) {
      if (n > 1)
        return((variance_emulators[[v_name]]$get_exp(x)^2 + variance_emulators[[v_name]]$get_cov(x))/n *
                 (kurt_aves[v_name] - 1 + 2/(n-1)))
      return(0)
    }
    variance_emulators[[v_name]]$s_diag <- var_mod
    variance_emulators[[v_name]]$em_type <- "variance"
    valid_cov_indices <- which(!is.na(collected_df_var[,v_name]))
    if (emulator_type == "covariance")
      variance_emulators[[v_name]]$samples <- collected_df_kurt$n
    else
      variance_emulators[[v_name]]$samples <- collected_df_kurt[valid_cov_indices,paste0(v_name, "_n")]
    adjust_df <- collected_df_var[valid_cov_indices ,c(input_names, v_name)]
    return(variance_emulators[[v_name]]$adjust(adjust_df, v_name))
  }) |> setNames(output_names)
  if (verbose && emulator_type != "covariance") cat("Completed variance emulators...\n") #nocov
  if (verbose) cat("Creating prior (untrained) mean emulators...\n") #nocov
  if (!is.character(corr_name) && !is.null(corr_name$expectation)) corr_name_exp <- corr_name$expectation
  else corr_name_exp <- corr_name
  mean_emulators <- emulator_from_data(collected_df_mean, output_names, ranges,
                                       input_names = input_names, emulator_type = "default",
                                       specified_priors = specified_priors$expectation,
                                       order = order, beta.var = beta.var,
                                       corr_name = corr_name_exp, adjusted = FALSE,
                                       discrepancies = discrepancies$expectation,
                                       verbose = FALSE, na.rm = TRUE, check.ranges = FALSE,
                                       targets = targets$expectation, has.hierarchy = TRUE, ...)
  if (verbose && emulator_type == "covariance") cat("Variance emulators trained. Creating covariance emulators...\n") #nocov
  ### Covariance emulation here...
  if (emulator_type == "covariance") {
    name_to_time <- function(names) {
      t_names <- sub("\\.*[^\\d](\\d+)$", "\\1", names)
      return(suppressWarnings(as.numeric(t_names)))
    }
    multi_point_cov <- function(i, j, x, theta.t, rho, v_ems) {
      part_indices <- partition_by_output(data, output_names, FALSE, TRUE)
      i_name <- v_ems[[i,i]]$output_name
      j_name <- v_ems[[j,j]]$output_name
      if (nrow(rho) == 1) group1 <- group2 <- 1
      else {
        group1 <- which(purrr::map_lgl(part_indices, ~i_name %in% .))
        group2 <- which(purrr::map_lgl(part_indices, ~j_name %in% .))
      }
      t_nam <- name_to_time(c(i_name, j_name))
      if (any(is.na(t_nam))) t_nam <- rep(0, 2)
      rho[group1, group2] * exp(-diff(t_nam)^2/theta.t^2) * sqrt(
        v_ems[[i,i]]$get_cov(x) * v_ems[[j,j]]$get_cov(x)
      )
    }
    comb_rv <- function(i, j, x, v_ems, k_est = 3) {
      if (length(k_est) == 1) k_est <- c(k_est, k_est)
      rvi <- (k_est[1]-1)*(v_ems[[i,i]]$get_cov(x) + v_ems[[i,i]]$get_exp(x)^2)
      rvj <- (k_est[2]-1)*(v_ems[[j,j]]$get_cov(x) + v_ems[[j,j]]$get_exp(x)^2)
      return(purrr::map2_dbl(rvi, rvj, min))
    }
    cov_vt <- function(x, i, j, n, cov_ems, theta.t, rho, k_est = 3) {
      if (i > j) {
        temp <- j
        j <- i
        i <- temp
      }
      1/n * comb_rv(i, j, x, cov_ems, k_est) + 1/(n*(n-1)) * (
        multi_point_cov(i, j, x, theta.t, rho, cov_ems) +
          cov_ems[[i,i]]$get_exp(x) * cov_ems[[j,j]]$get_exp(x) +
          cov_ems[[i,j]]$get_exp(x)^2
      )
    }
    if (is.null(covariance_opts$matrix))
      covariance_opts$matrix <- matrix(TRUE, nrow = length(variance_emulators), ncol = length(variance_emulators))
    which_outputs <- cov_out_names[upper.tri(cov_out_names) & covariance_opts$matrix]
    init_cov_ems <- emulator_from_data(collected_df_cov, which_outputs, ranges, input_names = input_names,
                                       emulator_type = "default", specified_priors = specified_priors$covariance,
                                       order = max(1, order - 1), beta.var = beta.var, corr_name = corr_name,
                                       adjusted = FALSE, discrepancies = discrepancies$covariance,
                                       verbose = FALSE, more_verbose = FALSE, na.rm = TRUE,
                                       check.ranges = FALSE, targets = targets$covariance, has.hierarchy = TRUE,
                                       ...)
    zero_em <- Proto_emulator$new(
      output_name = "zero_emulator",
      ranges = init_cov_ems[[1]]$ranges,
      predict_func <- function(x) return(rep(0, max(1,nrow(x)))),
      variance_func <- function(x) return(rep(1, max(1,nrow(x))))
    )
    init_cov_mat <- matrix(nrow = length(variance_emulators), ncol = length(variance_emulators))
    init_cov_mat[upper.tri(init_cov_mat) & covariance_opts$matrix] <- init_cov_ems
    init_cov_mat <- matrix(init_cov_mat, nrow = length(variance_emulators), ncol = length(variance_emulators))
    how_many_zero <- sum(upper.tri(init_cov_mat) & !covariance_opts$matrix)
    if (how_many_zero != 0)
      init_cov_mat[upper.tri(init_cov_mat) & !covariance_opts$matrix] <- list(zero_em)
    init_cov_mat[col(init_cov_mat) == row(init_cov_mat)] <- variance_emulators
    recomb_input_data <- cbind.data.frame(
      input_data[[1]][,input_names],
      do.call('cbind.data.frame', purrr::map(input_data, ~.[,length(.)]))
    ) |> setNames(c(input_names, output_names))
    if (is.null(covariance_opts$rho)) rho_mat <- get_mpc_rho_est(recomb_input_data, output_names)
    else rho_mat <- covariance_opts$rho
    if (is.null(covariance_opts$theta)) theta_val <- get_mpc_theta_est(recomb_input_data, output_names, variance_emulators, rho_mat)
    else theta_val <- covariance_opts$theta
    if (theta_val == 0 || is.nan(theta_val)) theta_val <- 1
    trained_cov_ems <- purrr::map(1:length(init_cov_ems), function(i) {
      indices <- as.numeric(which(cov_out_names == init_cov_ems[[i]]$output_name, arr.ind = TRUE))
      kurts <- kurt_aves[indices]
      init_cov_ems[[i]]$s_diag <- function(x, n) {
        cov_vt(x, indices[1], indices[2], n, init_cov_mat, theta_val, rho_mat, kurts)
      }
      init_cov_ems[[i]]$samples <- collected_df_kurt$n
      init_cov_ems[[i]]$em_type <- "covariance"
      init_cov_ems[[i]]$adjust(collected_df_cov, init_cov_ems[[i]]$output_name)
    })
    trained_cov_mat <- matrix(nrow = length(trained_var_ems), ncol = length(trained_var_ems))
    trained_cov_mat[upper.tri(trained_cov_mat) & covariance_opts$matrix] <- trained_cov_ems
    trained_cov_mat <- matrix(trained_cov_mat, nrow = length(trained_var_ems), ncol = length(trained_var_ems))
    trained_cov_mat[upper.tri(trained_cov_mat) & !covariance_opts$matrix] <- list(zero_em)
    trained_cov_mat[col(trained_cov_mat) == row(trained_cov_mat)] <- trained_var_ems
  }
  ## Mean emulators
  if (verbose) {
    if (emulator_type == "covariance")
      cat("Completed covariance emulators. Training mean emulators...\n")
    else
      cat("Training mean emulators...\n")
  }
  trained_mean_ems <- purrr::map(output_names, function(m_name) {
    if (!is.null(covariance_opts$logged) && covariance_opts$logged)
      mean_emulators[[m_name]]$s_diag <- function(x, n) exp(trained_var_ems[[m_name]]$get_exp(x))/n
    else
      mean_emulators[[m_name]]$s_diag <- function(x, n) trained_var_ems[[m_name]]$get_exp(x)/n
    valid_mean_indices <- which(!is.nan(collected_df_mean[,m_name]))
    if (emulator_type == "covariance")
      mean_emulators[[m_name]]$samples <- collected_df_kurt$n
    else
      mean_emulators[[m_name]]$samples <- collected_df_kurt[valid_mean_indices, paste0(m_name, "_n")]
    valid_mean <- collected_df_mean[valid_mean_indices, c(input_names, m_name)]
    return(mean_emulators[[m_name]]$adjust(valid_mean, m_name))
  }) |> setNames(output_names)
  if (emulator_type == "covariance")
    return(list(variance = EmulatedMatrix$new(trained_cov_mat, theta_val, rho_mat, logged = (!is.null(covariance_opts$logged) && covariance_opts$logged)),
                expectation = trained_mean_ems))
  return(list(variance = trained_var_ems, expectation = trained_mean_ems))
}

#' Variance Emulator Creation (Deprecated)
#'
#' Trains hierarchical emulators to stochastic systems
#'
#' This function is deprecated in favour of using \code{\link{emulator_from_data}}
#' with argument \code{emulator_type = "variance"}. See the associated help file.
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
#' \donttest{
#'  # A simple example using the BirthDeath dataset
#'  v_ems <- variance_emulator_from_data(BirthDeath$training, c("Y"),
#'   list(lambda = c(0, 0.08), mu = c(0.04, 0.13)), c_lengths = c(0.75))
#' }
#'
#' @export
variance_emulator_from_data <- function(input_data, output_names, ranges,
                                        input_names = names(ranges),
                                        verbose = interactive(), na.rm = FALSE, ...) {
  .Deprecated(new = "emulator_from_data", msg = "variance_emulator_from_data(...) is deprecated in favour of emulator_from_data(..., emulator_type = 'variance').")
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
  variance_emulators <- purrr::map(output_names, function(i) {
    is_high_rep <- !is.na(collected_df_var[,paste0(i,"kurt")])
    all_var <- setNames(
      collected_df_var[,c(input_names, paste0(i, 'var'))], c(input_names, i))
    all_n <- collected_df_var$n
    if (all(is_high_rep)) kurt_ave <- mean(collected_df_var[,paste0(i,'kurt')])
    else if (!any(is_high_rep)) kurt_ave <- 3
    else kurt_ave <- mean(collected_df_var[is_high_rep, paste0(i, 'kurt')])
    # if (verbose) print(paste0(i, " kurtosis:", kurt_ave))
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
          var_df, i, ranges,
          adjusted = FALSE, has.hierarchy = TRUE, verbose = FALSE, ...)[[1]]
      }
    }
    else {
      variance_em <- emulator_from_data(
        all_var, i, ranges,
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
        return((variance_em$get_exp(x)^2 +
                  variance_em$get_cov(x))/n * (kurt_ave - 1 + 2/(n-1)))
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
    return(v_em)
  })
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
#' This function is deprecated in favour of using \code{\link{emulator_from_data}}
#' with argument \code{emulator_type = "multistate"}. See the associated help file.
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
#' @importFrom cluster daisy fanny
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
  .Deprecated(new = "emulator_from_data", msg = "bimodal_emulator_from_data(...) is deprecated in favour of emulator_from_data(..., emulator_type = 'multistate').")
  if (na.rm) input_data <- input_data[apply(
    input_data, 1, function(x) !any(is.na(x))),]
  unique_points <- unique(data[,input_names])
  uids <- apply(unique_points, 1, hash)
  param_sets <- purrr::map(uids, function(x) {
    data[apply(data[,input_names], 1, hash) == x,]
  })
  if(verbose) cat("Separated dataset by unique points.\n") #nocov
  proportion <- purrr::map_dbl(param_sets, function(x) {
    p_clust <- suppressWarnings(fanny(suppressWarnings(daisy(x[, output_names, drop = FALSE])), k = 2)$clustering)
    return(sum(p_clust == 1)/length(p_clust))
  })
  prop_df <- setNames(
    data.frame(cbind(unique_points, proportion)),
    c(names(unique_points), 'prop'))
  if (verbose) cat("Training emulator to proportion in modes.\n") #nocov
  prop_em <- emulator_from_data(prop_df,
                                c('prop'), ranges, verbose = FALSE, ...)
  if (verbose) cat("Performing clustering to identify modes.\n") #nocov
  has_bimodality <- do.call('rbind.data.frame', purrr::map(param_sets, function(x) {
    purrr::map_lgl(output_names, function(y) {
      if (length(unique(x[,y])) == 1) return(FALSE)
      clust1 <- suppressWarnings(fanny(suppressWarnings(daisy(x[,y, drop = FALSE])), k = 1))
      clust2 <- suppressWarnings(fanny(suppressWarnings(daisy(x[,y, drop = FALSE])), k = 2))
      if (clust1$objective[["objective"]] < clust2$objective[["objective"]]) return(FALSE)
      return(TRUE)
    })
  })) |> setNames(output_names)
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
        this_clust <- suppressWarnings(fanny(suppressWarnings(daisy(param_sets[[i]][,x, drop = FALSE])), k = 2)$clustering)
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
