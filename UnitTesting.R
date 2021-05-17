unit_test <- function() {
  ro <- TRUE
  ## Check the correlation structure stuff
  test_u <- Correlator$new(nug = 0.1)
  ro <- test_u$get_corr(c(a = 0.1, b = 0.2, c = 0.3)) == 1
  ro <- round(test_u$get_corr(c(a = 0.1, b = 0.2, c = 0.3), c(a = 0.15, b = 0.18, c = 0.295)), 7) == 0.6717557
  ro <- round(test_u$get_corr(c(a = 0.1, b = 0.2, c = 0.3), c(a = 0.15, b = 0.18, c = 0.295), actives = c(T,T,F)), 7) == 0.6734372

  new_u <- test_u$set_hyper_p(list(theta = 0.5), 0.01)
  ro <- round(new_u$get_corr(c(a = 0.1, b = 0.2, c = 0.3), c(a = 0.15, b = 0.18, c = 0.295)), 7) == 0.9784845
  ro <- round(new_u$get_corr(c(a = 0.1, b = 0.2, c = 0.3), c(a = 0.15, b = 0.18, c = 0.295), actives = c(T,T,F)), 7) == 0.9785824

  if (!ro) stop("Problem with Correlator construction.")
  print("Correlator construction passes...")

  ## Check the Emulator stuff
  basis_f <- c(
    function(x) 1,
    function(x) x[[1]],
    function(x) x[[2]],
    function(x) x[[1]]*x[[2]]
  )
  beta <- list()
  beta$mu <- c(0.7, 1.2, 0.9, -0.01)
  beta$sigma <- diag(0, nrow = 4)
  u <- list()
  u$sigma <- 2
  u$corr <- test_u
  ranges <- list(a = c(0, 1), b = c(-2, 2))

  ## Boring, non-trained, emulator
  test_e <- Emulator$new(basis_f, beta, u, ranges)
  ro <- all(test_e$beta_mu == c(0.7, 1.2, 0.9, -0.01)) && test_e$u_sigma == 2 && test_e$corr$hyper_p$theta == 0.1
  if (!ro) stop("Problem with Emulator construction.")
  print("Emulator construction passes...")

  ## Checking expectation and variance
  test_lhs <- lhs::randomLHS(20, 2)
  test_data <- setNames(data.frame(t(apply(test_lhs, 1, function(x) (x - c(0, 0.5))*c(1, 4)))), c('a', 'b'))
  exp_out <- test_e$get_exp(test_data)
  test_data_prescaled <- data.frame(t(apply(test_data, 1, function(x) {
    (2*x - purrr::map_dbl(ranges, sum))/purrr::map_dbl(ranges, diff)
  })))
  other_way <- apply(test_data_prescaled, 1, function(x) purrr::map_dbl(test_e$basis_f, purrr::exec, x) %*% test_e$beta_mu)
  ro <- all(other_way == exp_out)

  var_out <- test_e$get_cov(test_data)
  ro <- all(var_out == 4)

  if (!ro) stop("Problem with non-trained emulator predictions.")
  print("Untrained emulator predictions as expected...")

  ## Checking adjust...
  train_data <- test_data
  train_data$f <- train_data$a + 2*train_data$b
  test_t_e <- test_e$adjust(train_data, 'f')

  ## Expectations should match at training points, and variance should be 0
  ro <- all(test_t_e$get_exp(train_data[,names(ranges)]) - train_data$f < 1e-10)
  ro <- all(test_t_e$get_cov(train_data[,names(ranges)]) < 1e-10)

  ## Look at some non-training points
  grid_data <- expand.grid(a = seq(0, 1, length.out = 20), b = seq(-2, 2, length.out = 20))
  test_data <- grid_data[sample(1:nrow(grid_data), 10),]

  ## The different ways of doing get_cov should all match up.
  ro <- all(test_t_e$get_cov(test_data) == diag(test_t_e$get_cov(test_data, full = TRUE)))
  ro <- all(test_t_e$get_cov(test_data) == diag(test_t_e$get_cov(test_data, test_data, full = TRUE)))
  ro <- all(test_t_e$get_cov(test_data) == test_t_e$get_cov(test_data, test_data))
  ro <- all(test_t_e$get_cov(test_data, full = TRUE) == test_t_e$get_cov(test_data, test_data, full = TRUE))
  ro <- all(test_t_e$get_cov(test_data, full = TRUE)[1:5, 6:10] == test_t_e$get_cov(test_data[1:5,], test_data[6:10,], full = TRUE))
  ro <- all(test_t_e$get_cov(test_data, full = TRUE)[3:7, 5:9] == test_t_e$get_cov(test_data[3:7,], test_data[5:9,], full = TRUE))

  if (!ro) stop("Emulator adjustments not behaving as expected.")
  print("Adjusted emulators behaving correctly...")

  ## Visualisation testing
  tryCatch({
    grid_exp <- test_t_e$get_exp(grid_data)
    grid_var <- test_t_e$get_cov(grid_data)
    full_data <- setNames(data.frame(cbind(grid_data, grid_exp, grid_var)), c(names(grid_data), "E", "V"))
    print(ggplot(data = full_data, aes(x = a, y = b)) +
      geom_contour_filled(aes(z = E), bins = 30))
    print(ggplot(data = full_data, aes(x = a, y = b)) +
      geom_contour_filled(aes(z = V), breaks = seq(0, 4.1, length.out = 30)))

    ## Modifying the correlation structure
    mod_t_e <- test_t_e$set_hyperparams(list(theta = 0.5), nugget = 0.05)
    mod_grid_exp <- mod_t_e$get_exp(grid_data)
    mod_grid_var <- mod_t_e$get_cov(grid_data)
    mod_full_data <- setNames(data.frame(cbind(grid_data, mod_grid_exp, mod_grid_var)), c(names(grid_data), "E", "V"))
    print(ggplot(data = mod_full_data, aes(x = a, y = b)) +
      geom_contour_filled(aes(z = E), bins = 30))
    print(ggplot(data = mod_full_data, aes(x = a, y = b)) +
      geom_contour_filled(aes(z = V), bins = 30))

    target <- list(val = 0.643 - 2*1.057, sigma = 0.01)
    imps <- mod_t_e$implausibility(grid_data, target)
    imp_lgl <- mod_t_e$implausibility(grid_data, target, 3)
    length(imp_lgl[imp_lgl]) == length(imps[imps < 3])
    mod_full_data$I <- imps
    print(ggplot(data = mod_full_data, aes(x = a, y = b)) +
      geom_contour_filled(aes(z = I), bins = 30))
  },
  error = function(e) {
    stop(paste("Plots of expectation, variance and implausibility not working as expected.", e))
  }
  )
  print("Emulator working with ggplot...")

  ### Testing model building and emulator_from_data
  ranges <- list(
    a = c(0, 1),
    b = c(-2, 2),
    c = c(3, 4),
    d = c(-2, 1),
    e = c(7, 10)
  )
  test_lhs <- lhs::randomLHS(50, 5)
  test_points <- data.frame(t(apply(test_lhs, 1, function(x) x*purrr::map_dbl(ranges, diff) + purrr::map_dbl(ranges, ~.[[1]]))))
  test_func <- function(x) {
    10*x[[1]] - x[[2]]/20 + x[[3]]*x[[4]]/5 - x[[5]]^2/10 + rnorm(1)
  }
  test_points$f <- apply(test_points, 1, test_func)
  scaled_points <- data.frame(t(apply(test_points[,c('a','b','c','d','e')], 1, function(x) (x - c(0.5, 0, 3.5, -0.5, 8.5))*c(2, 0.5, 2, 2/3, 2/3))))
  scaled_points$f <- test_points$f

  ### Different options for specifying emulator parts in emulator_from_data
  tryCatch({
    tdm0 <- emulator_from_data(test_points, c('f'), ranges = ranges)
    tdm1 <- emulator_from_data(test_points, c('f'), ranges = ranges, adjusted = FALSE)
    tdm2 <- emulator_from_data(test_points, c('f'), ranges = ranges, adjusted = TRUE)
    tdm3 <- emulator_from_data(test_points, c('f'), ranges = ranges, adjusted = TRUE, beta.var = TRUE)
    tdm4 <- suppressWarnings(emulator_from_data(scaled_points, c('f'), input_names = c('a','b','c','d','e')))
    tdm5 <- emulator_from_data(test_points, c('f'), ranges, quadratic = FALSE)
    tdm6 <- emulator_from_data(test_points, c('f'), ranges, funcs = list(c(function(x) 1, function(x) x[[1]], function(x) x[[5]], function(x) x[[1]]*x[[4]])))
    tdm7 <- emulator_from_data(test_points, c('f'), ranges, funcs = list(c(function(x) 1, function(x) x[[1]], function(x) x[[5]], function(x) x[[1]]*x[[4]])), beta = list(list(mu = c(-2.845, 5.13, -2.6, -1.43))))
    tdm8 <- emulator_from_data(test_points, c('f'), ranges, c_lengths = c(0.6))
    tdm9 <- emulator_from_data(test_points, c('f'), ranges, u = list(list(sigma = 0.9, corr = Correlator$new('exp_sq', list(theta = 0.25), nug = 0.02))))
    tdm10 <- emulator_from_data(test_points, c('f'), ranges, c_lengths = c(0.4), deltas = c(0.09))
    tdm11 <- emulator_from_data(test_points, c('f'), ranges, ev = c(0.5))
    tdm12 <- emulator_from_data(test_points, c('f'), ranges, deltas = c(0.143))
    tdm_whole <- emulator_from_data(test_points, c('f'), ranges,
                                    funcs = list(c(function(x) 1, function(x) x[[1]], function(x) x[[5]], function(x) x[[1]]*x[[4]])),
                                    beta = list(list(mu = c(-2.845, 5.13, -2.16, -1.43))),
                                    c_lengths = c(0.45),
                                    deltas = c(0.09))
    tdm_whole2 <- emulator_from_data(test_points, c('f'), ranges,
                                     funcs = list(c(function(x) 1, function(x) x[[1]], function(x) x[[5]], function(x) x[[1]]*x[[4]])),
                                     beta = list(list(mu = c(-2.845, 5.13, -2.16, -1.43))),
                                     u = list(list(sigma = 6, corr = Correlator$new('exp_sq', list(theta = 0.4), nug = 0.02))))

    test_func_2 <- function(x) {
      -5*x[[1]] + x[[2]]*x[[3]]/10 + x[[4]]^2 + rnorm(1)
    }
    test_points$g <- apply(test_points, 1, test_func_2)

    ### Multiple outputs
    test_derived_ems <- emulator_from_data(test_points, c('f','g'), ranges = ranges, c_lengths = c(0.4, 0.35))
  },
  error = function(e) stop(paste("Something went wrong with emulator_from_data specifications", e))
  )
  print("emulator_from_data working as expected...")

  ### Diagnostics: Use GillespieSIR dataset for simplicity.
  ranges <- list(
    aSI = c(0.1, 0.8),
    aIR = c(0, 0.5),
    aSR = c(0, 0.05)
  )
  out_vars <- c('nS', 'nI', 'nR')
  ems <- emulator_from_data(GillespieSIR, out_vars, ranges)
  targets <- list(
    nS = list(val = 281, sigma = 10.43),
    nI = list(val = 30, sigma = 11.16),
    nR = list(val = 689, sigma = 14.32)
  )
  tryCatch({
    for (i in names(targets)) {
      standard_errors(ems[[i]], GillespieValidation)
    }
    for (i in names(targets)) {
      comparison_diag(ems[[i]], GillespieValidation, targets, sd = 0.5)
    }
    for (i in names(targets)) {
      classification_diag(ems[[i]], GillespieValidation, targets, cutoff = 1.8)
    }
  }, error = function(e) stop(paste("Validation plotting failed", e))
  )
  ro <- nrow(validation_diagnostics(ems, GillespieValidation, targets)) == 0
  ro <- nrow(validation_diagnostics(ems, GillespieValidation, targets, c('ce','cd'))) == 0
  ro <- nrow(validation_diagnostics(ems, GillespieValidation, targets, cutoff = 2, sd = 2)) == 1
  if (!ro) stop("Combined validation diagnostics not working as expected")
  print("Validation diagnostic functionality working...")

  ## Nth Implausibility Checks

  ro <- round(nth_implausible(ems, GillespieValidation, targets)[1], 6) == 3.30738
  ro <- round(nth_implausible(ems, GillespieValidation, targets, n = 2)[1], 6) == 1.30417
  ro <- round(nth_implausible(ems, GillespieValidation, targets, n = 6)[1], 6) == 0.121097
  ro <- sum(nth_implausible(ems, GillespieValidation, targets, cutoff = c(4, 1, 5))) == 25
  ro <- round(nth_implausible(ems, data.frame(aSI = 0.4, aIR = 0.25, aSR = 0.025), targets), 6) == 1.038521

  if (!ro) stop("Nth Implausibility not working as expected.")
  print("Nth implausibility working...")

  ## Point Proposal Checks
  tryCatch(
    {
      all_methods <- generate_new_runs(ems, 200, targets, method = c('lhs', 'line', 'importance', 'slice', 'optical'), cutoff = 0.5)
      default_method <- generate_new_runs(ems, 200, targets, cutoff = 0.5)
      sample_points <- default_method[sample(nrow(default_method), 100),]
      with_plausible <- generate_new_runs(ems, 200, targets, method = c('line', 'importance'), cutoff = 0.5, plausible_set = sample_points, distro = 'normal')
      no_points <- generate_new_runs(ems, 100, targets, cutoff = 0)

      plot(rbind(all_methods, default_method, with_plausible), pch = 16, cex = 0.5, col = c(rep('black', nrow(all_methods)), rep('blue', nrow(default_method)), rep('red', nrow(with_plausible))))
    },
    error = function(e) {
      stop("Point proposal failed.")
    }
  )
  print("Point proposals working...")

  tryCatch(
    {
      emulator_plot(ems, ppd = 10)
      emulator_plot(ems$nS, ppd = 10)
      emulator_plot(ems, var_name = 'var', ppd = 10, params = c('aIR', 'aSR'))
      emulator_plot(ems, var_name = 'imp', ppd = 10, targets = targets, fixed_vals = list(aSR = 0.02))
      emulator_plot(ems, var_name = 'nimp', cb = TRUE, targets = targets, nth = 2, ppd = 10)
      print("Standard plots successful...")
      plot_lattice(ems, targets, ppd = 10)
      plot_lattice(ems[c('nS', 'nI')], targets[c('nS', 'nI')], ppd = 10)
      print("Lattice plots successful...")
    },
    error = function(e) {
      stop("Emulator plotting failed.")
    }
  )

  tryCatch(
    {
      behaviour_plot(ems, GillespieValidation)
      behaviour_plot(ems$nS, GillespieValidation)
      output_plot(ems, targets)
      output_plot(ems, targets, points = GillespieSIR)
      output_plot(ems[c('nS', 'nI')], targets = targets[c('nS', 'nI')])
      print("General behaviour plots successful...")
      space_removed(ems, targets, ppd = 10)
      space_removed(ems$nS, targets$nS, ppd = 10, modified = 'hp')
      space_removed(ems, targets, ppd = 10, modified = 'var')
      validation_pairs(ems, GillespieValidation, targets)
      validation_pairs(ems$nS, GillespieValidation, targets)
      print("Diagnostic plotting successful...")
    },
    error = function(e) {
      stop("Validation plots failed.")
    }
  )

  tryCatch(
    {
      simulator_plot(GillespieMultiWaveData, targets)
      simulator_plot(GillespieMultiWaveData[2:4], targets, zero_in = FALSE, wave_numbers = c(1,3))
      wave_points(GillespieMultiWaveData, names(ranges))
      wave_points(GillespieMultiWaveData, names(ranges), TRUE, 0.8)
      wave_values(GillespieMultiWaveData, targets)
      wave_values(GillespieMultiWaveData, targets[c('nS', 'nI')], surround = TRUE, p_size = 1, l_wid = 0.5)
      print("Wave plots working...")
      wave_variance(GillespieMultiWaveEmulators, names(targets), ppd = 10)
      wave_variance(GillespieMultiWaveEmulators, names(targets), plot_dirs = c('aIR', 'aSR'), wave_numbers = c(2,3), ppd = 20, sd = TRUE)
    },
    error = function(e) {
      stop("Wave plots failed.")
    }
  )

  tryCatch(
    {
      default <- full_wave(rbind(GillespieSIR, GillespieValidation), ranges, targets)
      non_quad <- full_wave(rbind(GillespieSIR, GillespieValidation), ranges, targets, quadratic = FALSE)
      second <- full_wave(GillespieMultiWaveData[[2]], ranges, targets, old_emulators = GillespieMultiWaveEmulators[[1]])
    },
    error = function(e) {
      stop("Automated wave generation failed.")
    }
  )
  print("Automatic wave generation working.")

  tryCatch(
    {
      deriv_exp <- purrr::map(ems, ~purrr::map(names(ranges), function(x) .$get_exp_d(GillespieValidation, x)))
      deriv_var <- purrr::map(ems, ~purrr::map(names(ranges), function(x) .$get_cov_d(GillespieValidation, x)))
    },
    error = function(e) {
      stop("Derivative testing failed.")
    }
  )
  print("Derivative calculations successful.")

  tryCatch(
    {
      ve <- variance_emulator_from_data(BirthDeath$var, BirthDeath$mean, BirthDeath$reps, paste0('t', c(1, 7, 15)), list(lambda = c(0, 0.08), mu = c(0.04, 0.13)))
      emulator_plot(ve$variance)
      emulator_plot(ve$expectation)
    },
    error = function(e) {
      stop("Variance emulator construction failed.")
    }
  )
  print("Variance emulators working.")

  if (ro) return("All unit tests passed. Hooray!")
}

library(hmer)
library(ggplot2)
unit_test()


