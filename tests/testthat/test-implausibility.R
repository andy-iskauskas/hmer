## Implausibility calculations

## Deterministic emulators
test_that("Implausibility with no cutoff" ,{
  n_imps <- nth_implausible(SIREmulators$ems,
                            SIRSample$validation,
                            SIREmulators$targets)
  expect_equal(
    n_imps[[1]],
    c(nth_implausible(SIREmulators$ems,
                      SIRSample$validation[1,],
                      SIREmulators$targets), use.names = FALSE)
  )
  expect_equal(
    length(n_imps),
    nrow(SIRSample$validation)
  )
})

test_that("Implausibility with cutoff", {
  n_imps_lgl <- nth_implausible(SIREmulators$ems,
                                SIRSample$validation,
                                SIREmulators$targets,
                                cutoff = 3)
  expect_true(
    isTRUE(all(is.logical(n_imps_lgl)))
  )
  expect_equal(
    sum(n_imps_lgl),
    3
  )
})

test_that("Nth maximum implausibility", {
  ems <- SIREmulators$ems
  pts <- SIRSample$validation
  targ <- SIREmulators$targets
  imps1 <- nth_implausible(ems, pts, targ)
  imps2 <- nth_implausible(ems, pts, targ, n = 2)
  imps3 <- nth_implausible(ems, pts, targ, n = 3)
  expect_true(
    sum(imps1 <= 3) <= sum(imps2 <= 3)
  )
  expect_true(
    sum(imps2 <= 3) <= sum(imps3 <= 3)
  )
  expect_equal(
    imps2 <= 3,
    nth_implausible(ems, pts, targ, cutoff = 3, n = 2)
  )
})

test_that("Warn if nth too large", {
  expect_warning(
    nth_implausible(
      SIREmulators$ems,
      SIRSample$validation,
      SIREmulators$targets,
      n = 4
    ),
    paste("n cannot be greater than the number of targets to match to.",
                    "Switching to minimum implausibility.")
  )
})

test_that("Warning if target is misspecified", {
  fake_targets <- SIREmulators$targets
  fake_targets$nI <- 169
  expect_warning(
    nth_implausible(SIREmulators$ems, SIRSample$validation, fake_targets),
    "Target nI is a single value. Assuming it's a val with sigma = 5%."
  )
})

test_that("Automatic 2nd max works", {
  fake_targets <- setNames(rep(SIREmulators$targets, 4), paste0("X", 1:12))
  fake_ems <- c()
  for (i in 1:12) {
    em_to_clone <- SIREmulators$ems[[(i-1) %% 3 + 1]]
    fake_ems <- c(fake_ems, em_to_clone$clone())
    fake_ems[[i]]$output_name <- paste0("X", i)
  }
  expect_equal(
    nth_implausible(fake_ems, SIRSample$validation, fake_targets),
    nth_implausible(fake_ems, SIRSample$validation, fake_targets, n = 2)
  )
  expect_equal(
    nth_implausible(fake_ems, SIRSample$validation, fake_targets),
    nth_implausible(SIREmulators$ems,
                    SIRSample$validation,
                    SIREmulators$targets, n = 1)
  )
})

test_that("Sequential implausibility works as expected", {
  expect_equal(
    nth_implausible(SIREmulators$ems,
                    SIRSample$validation,
                    SIREmulators$targets,
                    cutoff = 3, n = 2,
                    sequential = TRUE),
    c(nth_implausible(SIREmulators$ems,
                    SIRSample$validation,
                    SIREmulators$targets,
                    cutoff = 3, n = 2), use.names = FALSE)
  )
})

test_that("Processed raw = TRUE is same as raw = FALSE", {
  raw_imps <- nth_implausible(SIREmulators$ems,
                              SIRSample$validation,
                              SIREmulators$targets,
                              get_raw = TRUE)
  expect_equal(
    apply(raw_imps, 1, max),
    nth_implausible(SIREmulators$ems,
                    SIRSample$validation,
                    SIREmulators$targets)
  )
})

test_that("Fails if an emulator is without a target", {
  expect_error(
    nth_implausible(
      SIREmulators$ems[[1]],
      SIRSample$validation,
      SIREmulators$targets[-1]
    ),
  )
})

## Variance Emulation

var_em <- emulator_from_data(BirthDeath$training, c('Y'),
                             list(lambda = c(0, 0.08), mu = c(0.04, 0.13)),
                             emulator_type = "variance", verbose = FALSE)

test_that("Variance emulation implausibility basics", {
  var_target <- list(Y = c(90, 110))
  var_points <- unique(BirthDeath$validation[,c('lambda', 'mu')])
  var_imps <- nth_implausible(var_em, var_points, var_target)
  expect_equal(
    length(var_imps),
    nrow(var_points)
  )
  expect_equal(
    sum(var_imps <= 3),
    sum(nth_implausible(var_em, var_points, var_target, cutoff = 3))
  )
})

test_that("Different combinations of variance targets and emulators", {
  variance_targets <- list(
    expectation = list(Y = c(90, 110)),
    variance = list(Y = c(30, 35))
  )
  var_points <- unique(BirthDeath$validation[,c('lambda', 'mu')])
  expect_equal(
    ncol(nth_implausible(var_em, var_points, variance_targets, get_raw = TRUE)),
    2
  )
  expect_equal(
    apply(nth_implausible(var_em, var_points, variance_targets, get_raw = TRUE), 1, max),
    nth_implausible(var_em, var_points, variance_targets)
  )
  expect_equal(
    nth_implausible(var_em$expectation, var_points, variance_targets),
    nth_implausible(var_em$expectation, var_points, variance_targets$expectation)
  )
})


## Bimodal Emulation

bim_targets <- list(
  I10 = list(val = 35, sigma = 3.5),
  I25 = list(val = 147, sigma = 14.7),
  I50 = list(val = 55, sigma = 5.5),
  R10 = list(val = 29, sigma = 2.9),
  R25 = list(val = 276, sigma = 27.6),
  R50 = list(val = 579, sigma = 57.9)
)
bim_em <- emulator_from_data(SIR_stochastic$training, names(bim_targets),
                             SIREmulators$ems[[1]]$ranges, verbose = FALSE,
                             emulator_type = "multistate")

test_that("Bimodal emulation basics", {
  bim_points <- unique(SIR_stochastic$validation[,names(SIREmulators$ems[[1]]$ranges)])
  expect_equal(
    sort(names(nth_implausible(bim_em, bim_points, bim_targets, get_raw = TRUE))),
    sort(names(bim_targets))
  )
})

test_that("Separated modes: same result", {
  bim_points <- unique(SIR_stochastic$validation[,names(SIREmulators$ems[[1]]$ranges)])
  first_imps <- nth_implausible(bim_em$mode1, bim_points, bim_targets, get_raw = TRUE)
  second_imps <- nth_implausible(bim_em$mode2, bim_points, bim_targets, get_raw = TRUE)
  expect_equal(
    apply(as.data.frame(Map(pmin, first_imps, second_imps[,names(first_imps)])), 1, max),
    nth_implausible(bim_em, bim_points, bim_targets)
  )
})

### Point proposal - makes more sense here since emulators are trained
test_that("Bimodal point proposal", {
  skip_on_cran()
  bim_points <- generate_new_design(
    bim_em,
    100,
    bim_targets,
    resample = 0, verbose = FALSE
  )
  expect_equal(
    nrow(bim_points),
    100
  )
})
