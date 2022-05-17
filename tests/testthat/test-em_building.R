test_that("Emulator training: default", {
  ems <- emulator_from_data(
    SIRSample$training,
    c('nS', 'nI', 'nR'),
    list(aSI = c(0.1, 0.8),
         aIR = c(0, 0.5),
         aSR = c(0, 0.05)),
    verbose = FALSE
  )
  expect_equal(
    ems[[1]]$beta_mu,
    SIREmulators$ems[[1]]$beta_mu
  )
  expect_equal(
    c(ems[[2]]$get_exp(SIRSample$training[1:5,]),
      use.names = FALSE),
    SIRSample$training[1:5, 'nI']
  )
  expect_equal(
    c(ems[[3]]$get_cov(SIRSample$training[1:5,]),
      use.names = FALSE),
    rep(0, 5)
  )
  expect_equal(
    matrix(ems[[2]]$get_cov(SIRSample$validation[1:3,], full = TRUE),
           dimnames = NULL, nrow = 3),
    matrix(c(237.1851790, 0.9288659, 5.660917,
             0.9288659, 234.3065391, 0.000001,
             5.6609166, 0.000001, 216.613111), nrow = 3),
    tolerance = 1e-6
  )
})

test_that("Emulator training: matern and preflight", {
  ems_extra <- emulator_from_data(
    SIRSample$training,
    c('nS', 'nI', 'nR'),
    list(aSI = c(0.1, 0.8),
         aIR = c(0, 0.5),
         aSR = c(0, 0.05)),
    corr_name = "matern",
    na.rm = TRUE,
    check.ranges = TRUE,
    targets = SIREmulators$targets,
    verbose = FALSE,
  )
  expect_equal(
    ems_extra[[1]]$corr$corr_name,
    "matern"
  )
})

test_that("Emulator training: orn_uhl, no ranges, few points", {
  expect_warning(
    ems_extra_2 <- emulator_from_data(
      SIRSample$training[1:20,],
      c('nS', 'nI', 'nR'),
      input_names = c('aSI', 'aIR', 'aSR'),
      corr_name = "orn_uhl",
      verbose = FALSE
    )
  )
  expect_equal(
    ems_extra_2[[1]]$ranges,
    list(aSI = c(-1, 1),
         aIR = c(-1, 1),
         aSR = c(-1, 1))
  )
})

test_that("Emulator training: no correlation", {
  expect_warning(
    ems_extra_3 <- emulator_from_data(
      SIRSample$training,
      c('nS', 'nI', 'nR'),
      list(
        aSI = c(0.1, 0.8),
        aIR = c(0, 0.5),
        aSR = c(0, 0.05)
      ),
      check.ranges = TRUE,
      corr_name = "not_a_correlation",
      verbose = FALSE
    )
  )
  expect_equal(
    ems_extra_3[[2]]$corr$corr_name,
    "exp_sq"
  )
})

test_that("Emulator training: bad correlation", {
  expect_warning(
    ems_extra_4 <- emulator_from_data(
      SIRSample$training,
      c('nS', 'nI', 'nR'),
      list(aSI = c(0.1, 0.8),
           aIR = c(0, 0.5),
           aSR = c(0, 0.05)),
      corr_name = "gamma_exp",
      verbose = FALSE
    )
  )
  expect_equal(
    ems_extra_4[[3]]$corr$corr_name,
    "exp_sq"
  )
})

test_that("Emulator training: provided hyperparams", {
  ems_5 <- emulator_from_data(
    SIRSample$training,
    c('nS', 'nI', 'nR'),
    list(
      aSI = c(0.1, 0.8),
      aIR = c(0, 0.5),
      aSR = c(0, 0.05)
    ),
    c_lengths = rep(0.75, 3),
    ev = c(30, 20, 10),
    verbose = FALSE
  )
  expect_equal(
    ems_5[[1]]$corr$hyper_p$theta,
    0.75
  )
  expect_equal(
    ems_5[[3]]$corr$nugget,
    0.2587429,
    tolerance = 1e-6
  )
})

test_that("Emulator training: uncertain beta", {
  discreps <- c(2, 4, 3.5)
  ems <- emulator_from_data(
    SIRSample$training,
    c('nS', 'nI', 'nR'),
    list(aSI = c(0.1, 0.8),
         aIR = c(0, 0.5),
         aSR = c(0, 0.05)),
    beta.var = TRUE,
    discrepancies = discreps,
    verbose = FALSE)
  expect_equal(
    c(
      ems$nS$get_exp(SIRSample$training[1:5,]),
      use.names = FALSE),
    SIRSample$training[1:5, 'nS'],
    tolerance = 1e-6
  )
  expect_equal(
    c(
      ems$nI$get_cov(SIRSample$training[1:5,]),
      use.names = FALSE),
    rep(0, 5),
    tolerance = 1e-6
  )
  expect_equal(
    ems$nR$disc,
    list(internal = 3.5, external = 0)
  )
})

test_that("Range handling", {
  standard_ranges <- list(
    aSI = c(0.1, 0.8),
    aIR = c(0, 0.5),
    aSR = c(0, 0.05)
  )
  ranges1 <- matrix(c(0.1, 0, 0, 0.8, 0.5, 0.05),
                    nrow = 3)
  row.names(ranges1) <- c('aSI', 'aIR', 'aSR')
  expect_equal(
    convertRanges(ranges1),
    standard_ranges
  )
  expect_warning(
    convertRanges(t(ranges1))
  )
  ranges2 <- list(aSI = c(0.1), aIR = c(0, 0.5), aSR = c(0, 0.05))
  expect_warning(
    convertRanges(ranges2)
  )
  ranges3 <- data.frame(min = c(0.1, 0, 0), max = c(0.8, 0.5, 0.05))
  row.names(ranges3) <- c("aSI", "aIR", "aSR")
  expect_equal(
    convertRanges(ranges3),
    standard_ranges
  )
})

test_that("Full wave behaves", {
  skip_on_cran()
  fw <- full_wave(rbind(
    SIRSample$training,
    SIRSample$validation),
                  list(aSI = c(0.1, 0.8),
                       aIR = c(0, 0.5),
                       aSR = c(0, 0.05)),
                  targets = SIREmulators$targets,
    verbose = FALSE)
  expect_equal(
    length(fw$emulators),
    3
  )
  expect_equal(
    nrow(fw$points),
    90
  )
})

test_that("Full wave with all atomic targets", {
  skip_on_cran()
  all_atomic <- SIREmulators$targets
  all_atomic$nI <- c(169-3*8.45, 169+3*8.45)
  fw <- full_wave(rbind(
    SIRSample$training,
    SIRSample$validation),
    list(aSI = c(0.1, 0.8),
         aIR = c(0, 0.5),
         aSR = c(0, 0.05)),
    targets = all_atomic,
    verbose = FALSE,
    old_emulators <- SIREmulators$ems
  )
  expect_equal(
    length(fw$emulators),
    3
  )
  expect_equal(
    nrow(fw$points),
    90
  )
})
