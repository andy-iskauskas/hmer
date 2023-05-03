em <- emulator_from_data(SIRSample$training, c('nS', 'nI', 'nR'),
                         list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05)),
                         verbose = FALSE)

test_that("Obtaining derivative information", {
  only_exp <- get_deriv_info(
    SIREmulators$ems$nI,
    SIRSample$validation[1,]
  )
  exp_and_var <- get_deriv_info(
    SIREmulators$ems$nR,
    SIRSample$validation[5,],
    var = TRUE
  )
  expect_equal(
    length(only_exp),
    1
  )
  expect_equal(
    length(only_exp$exp),
    length(SIREmulators$ems[[1]]$ranges)
  )
  expect_equal(
    c(only_exp$exp, use.names = FALSE),
    c(
      SIREmulators$ems$nI$get_exp_d(SIRSample$validation[1,], 'aSI'),
      SIREmulators$ems$nI$get_exp_d(SIRSample$validation[1,], 'aIR'),
      SIREmulators$ems$nI$get_exp_d(SIRSample$validation[1,], 'aSR')
    )
  )
  expect_true(
    !is.null(exp_and_var$exp) &&
      !is.null(exp_and_var$var)
  )
  expect_equal(
    length(exp_and_var$exp),
    length(SIREmulators$ems[[1]]$ranges)
  )
  expect_equal(
    dim(exp_and_var$var),
    rep(length(SIREmulators$ems[[1]]$ranges), 2)
  )
  expect_equal(
    exp_and_var$var,
    t(exp_and_var$var)
  )
})

test_that("Directional derivative function works", {
  direction <- c(get_deriv_info(
    SIREmulators$ems$nS,
    SIRSample$training[1,]
  )$exp, use.names = FALSE)
  d_deriv_exp <- directional_deriv(
    SIREmulators$ems$nS,
    SIRSample$training[1,],
    c(1, 1, 1)
  )
  expect_equal(
    d_deriv_exp,
    c(1,1,1) %*% direction/sqrt(3*sum(direction^2))
  )
  var_mod_direct <- directional_deriv(
    SIREmulators$ems$nS,
    SIRSample$training[1,],
    c(1,1,1), sd = 1
  )
  expect_equal(
    length(var_mod_direct),
    2
  )
  expect_true(
    d_deriv_exp[[1]] >= var_mod_direct[[1]] &&
      d_deriv_exp[[1]] <= var_mod_direct[[2]]
  )
})

test_that("Derivative proposals" ,{
  exp_measure <- directional_proposal(
    SIRMultiWaveEmulators[[3]],
    SIRSample$validation[2,],
    SIREmulators$targets
  )
  imp_measure <- directional_proposal(
    SIRMultiWaveEmulators[[3]],
    SIRSample$validation[2,],
    SIREmulators$targets,
    iteration.measure = "imp"
  )
  expect_true(
    nth_implausible(SIRMultiWaveEmulators[[3]],
                    SIRSample$validation[2,],
                    SIREmulators$targets) >=
      nth_implausible(SIRMultiWaveEmulators[[3]],
                      exp_measure,
                      SIREmulators$targets)
  )
  expect_true(
    nth_implausible(SIRMultiWaveEmulators[[3]],
                    SIRSample$validation[2,],
                    SIREmulators$targets) >=
      nth_implausible(SIRMultiWaveEmulators[[3]],
                      imp_measure,
                      SIREmulators$targets)
  )
})
