test_that("Basic emulator construction", {
  test_em <- Emulator$new(
    basis_f <- c(function(x) 1, function(x) x[[1]]),
    beta = list(mu = c(1, 1),
                sigma = diag(0, nrow = 2)),
    u <- list(sigma = 2, corr = Correlator$new()),
    ranges <- list(x = c(0, 2))
  )
  expect_equal(
    test_em$active_vars,
    c(TRUE)
  )
  expect_equal(
    test_em$beta_sigma,
    diag(0, nrow = 2)
  )
  expect_equal(
    purrr::map_dbl(
      test_em$basis_f, purrr::exec, data.frame(x = 2)
    ),
    c(1, 2)
  )
  expect_equal(
    test_em$u_sigma,
    2
  )
  expect_equal(
    test_em$corr$corr_name,
    "exp_sq"
  )
  expect_equal(
    test_em$corr$hyper_p$theta,
    0.1
  )
})

data <- data.frame(x = seq(-1, 1, by = 0.2), y = seq(0, 2, by = 0.2),
                   f = 0.8*sin((seq(-1, 1, by = 0.2)-0.2)*pi/0.35))
data_em <- Emulator$new(
  basis_f = c(function(x) 1),
  beta = list(mu = c(1), sigma = diag(0, nrow = 1)),
  u = list(
    corr = Correlator$new('matern', hp = list(theta = 0.2, nu = 1.5)),
    sigma = 1
  ),
  ranges = list(x = c(-1, 1), y = c(0, 2)),
)
data_em$output_name <- 'f'

test_that("Emulator with data", {
  test_data <- data.frame(x = c(-0.2, 0, 0.2), y = c(0.8, 1, 1.2))
  expect_equal(
    data_em$get_exp(test_data),
    c(1, 1, 1)
  )
  expect_equal(
    data_em$get_cov(test_data),
    c(1, 1, 1)
  )
  data_em_adj <- data_em$adjust(data, 'f')
  expect_equal(
    data_em_adj$get_exp(test_data),
    data$f[5:7],
    tolerance = 1e-5
  )
  expect_equal(
    c(data_em_adj$get_cov(test_data), use.names = FALSE),
    c(0, 0, 0),
    tolerance = 1e-5
  )
})

em <- emulator_from_data(SIRSample$training,
                         c('nI'),
                         list(aSI = c(0.1, 0.8),
                              aIR = c(0, 0.5),
                              aSR = c(0, 0.05)),
                         verbose = FALSE)$nI

test_that("Trained emulator covariance", {
  expect_equal(
    length(
      em$get_cov(SIRSample$validation[1:5,],
                 SIRSample$validation[5:9,])
    ),
    5
  )
  expect_equal(
    dim(
      em$get_cov(SIRSample$validation[1:5,],
                 SIRSample$validation[5:10,],
                 full = TRUE)
    ),
    c(5, 6)
  )
})

test_that("Modifying priors and functional sigma", {
  em_2 <- em$set_sigma(2)
  expect_equal(
    em_2$u_sigma,
    2
  )
  em_3 <- em_2$mult_sigma(2)
  expect_equal(
    em_3$u_sigma,
    4
  )
  em_4 <- em_2$set_hyperparams(
    hp = list(theta = 0.75),
    nugget = 0.1
  )
  expect_equal(
    em_4$corr$hyper_p$theta,
    0.75
  )
  expect_equal(
    em_4$corr$nugget,
    0.1
  )
  em_sigma <- em$set_sigma(function(x) x[[1]]*5)
  expect_false(
    all(em_sigma$get_cov(SIRSample$training[1:3,]) == 0)
  )
  expect_equal(
    dim(em_sigma$get_cov(SIRSample$training[1:3,],
                         SIRSample$training[2:5,],
                         full = TRUE)),
    c(3, 4)
  )
  expect_equal(
    c(em_sigma$get_exp(SIRSample$validation[1:3,],
                     include_c = FALSE), use.names = FALSE),
    c(85.11743, 59.98822, 338.93812),
    tolerance = 1e-4
  )
  em_sigma_2 <- em_sigma$mult_sigma(2)
  expect_equal(
    em_sigma_2$u_sigma(c(1, 0, 0)),
    10
  )
})

test_that("Modifying priors and functional sigma - untrained", {
  em_o <- em$o_em
  em_o2 <- em_o$set_sigma(2)
  expect_equal(
    em_o2$u_sigma,
    2
  )
  em_o3 <- em_o2$mult_sigma(2)
  expect_equal(
    em_o3$get_cov(SIRSample$validation[1,,drop=FALSE]),
    4
  )
  expect_equal(
    em_o3$u_sigma,
    2
  )
  em_o4 <- em_o2$set_hyperparams(
    hp = list(theta = 0.7),
    nugget = 0.3
  )
  expect_equal(
    em_o4$corr$hyper_p$theta,
    0.7
  )
  expect_equal(
    em_o4$corr$nugget,
    0.3
  )
})

test_that("Derivative functions", {
  expect_equal(
    nrow(
      em$get_exp_d(SIRSample$training[1:5,], 'aSI')
    ),
    5
  )
  expect_equal(
    c(em$get_exp_d(SIRSample$training[1:5,], 'aSR'), use.names = FALSE),
    rep(0, 5)
  )
  expect_equal(
    length(em$get_cov_d(SIRSample$training[1:5,], 'aSI')),
    5
  )
  expect_equal(
    dim(
      em$get_cov_d(SIRSample$training[1:5,], 'aSI',
                   SIRSample$training[2:7,], 'aIR',
                   full = TRUE)
    ),
    c(5, 6)
  )
  oem <- em$o_em
  expect_equal(
    nrow(
      oem$get_exp_d(SIRSample$training[1:5,], 'aSI')
    ),
    5
  )
  expect_equal(
    c(oem$get_exp_d(SIRSample$training[1:5,], 'aSR'), use.names = FALSE),
    rep(0, 5)
  )
  expect_equal(
    length(unique(oem$get_cov_d(SIRSample$training[1:5,], 'aSI'))),
    1
  )
  expect_equal(
    dim(
      oem$get_cov_d(SIRSample$training[1:5,], 'aSI',
                   SIRSample$training[2:7,], 'aIR',
                   full = TRUE)
    ),
    c(5, 6)
  )
})

test_that("Batch processing is called for >1000 points", {
  many_points <- data.frame(
    aSI = runif(2400, 0.1, 0.8),
    aIR = runif(2400, 0, 0.5),
    aSR = runif(2400, 0, 0.05)
  )
  expect_equal(
    length(c(
      em$get_exp(many_points))),
    2400
  )
  expect_equal(
    length(c(
      em$get_cov(many_points))),
    2400
  )
  expect_equal(
    length(em$implausibility(many_points, SIREmulators$targets$nS)),
    2400
  )
})

test_that("Emulator with no variable dependence", {
  fake_grid <- expand.grid(x = seq(1, 10, by = 1), y = seq(1, 10, by = 1))
  fake_output <- runif(100, -5, 5)
  fake_data <- cbind.data.frame(fake_grid, fake_output) |> setNames(letters[24:26])
  fake_em <- emulator_from_data(fake_data, c('z'), list(x = c(1, 10), y = c(1, 10)),
                                beta_var = TRUE, verbose = FALSE)$z
  expect_equal(
    c(fake_em$get_exp(fake_data), use.names = FALSE),
    c(fake_data$z),
    tolerance = 1e-6
  )
  expect_true(
    all(fake_em$get_cov(fake_data) < 1e-6)
  )
  expect_equal(
    length(c(fake_em$o_em$get_exp(fake_data))),
    100
  )
  expect_equal(
    length(c(fake_em$o_em$get_cov(fake_data))),
    100
  )
  expect_equal(
    length(c(fake_em$get_exp_d(fake_data, 'x'))),
    100
  )
  expect_equal(
    length(c(fake_em$get_cov_d(fake_data, 'x'))),
    100
  )
})

test_that("Printing works", {
  expect_output(
    print(em),
    "Parameters and ranges"
  )
  expect_output(
    print(em),
    "Regression surface Variance"
  )
  expect_output(
    print(em),
    "Bayes-adjusted emulator - prior specifications listed"
  )
})
