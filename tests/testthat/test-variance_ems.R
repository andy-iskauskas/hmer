v_em <- variance_emulator_from_data(BirthDeath$training,
                                    c('Y'),
                                    list(lambda = c(0, 0.08), mu = c(0.04, 0.13)),
                                    verbose = FALSE, beta.var = TRUE)

test_that("Variance Emulators", {
  expect_equal(
    class(v_em$expectation$Y),
    c("Hierarchical", "Emulator", "R6")
  )
  expect_equal(
    v_em$expectation$Y$em_type,
    "mean"
  )
  expect_equal(
    v_em$variance$Y$em_type,
    "variance"
  )
})

test_that("Batch runs", {
  many_points <- data.frame(lambda = runif(1400, 0, 0.08),
                            mu = runif(1400, 0.04, 0.13))
  expect_equal(
    length(c(v_em$variance$Y$get_exp(many_points))),
    1400
  )
  expect_equal(
    length(c(v_em$expectation$Y$get_cov(many_points))),
    1400
  )
  expect_equal(
    length(c(v_em$expectation$Y$implausibility(many_points,
                                               list(Y = c(90, 110))$Y))),
    1400
  )
})

em <- v_em$variance$Y
test_train <- unique(BirthDeath$training[,c('lambda', 'mu')])[1:5,]
test_points <- unique(BirthDeath$validation[,c('lambda', 'mu')])[1:5,]
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
    all(em_sigma$get_cov(test_train) == 0)
  )
  expect_equal(
    dim(em_sigma$get_cov(test_train[1:3,],
                         test_train[2:5,],
                         full = TRUE, check_neg = FALSE)),
    c(3, 4)
  )
  em_sigma_2 <- em_sigma$mult_sigma(2)
  expect_equal(
    em_sigma_2$u_sigma(c(0.01, 0)),
    0.1
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
    unname(em_o3$get_cov(test_points[1,,drop=FALSE])),
    359.2923,
    tolerance = 1e-4
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
