test_that("Diagnostic summary measures", {
  expect_equal(
    summary_diag(
      SIREmulators$ems[[1]],
      SIRSample$validation,
      verbose = FALSE
    ),
    c(TRUE, TRUE)
  )
})

test_that("Residuals", {
  expect_equal(
    nrow(
      residual_diag(
        SIREmulators$ems[[2]],
        histogram = FALSE
      )
    ),
    0
  )
})

test_that("Space removal", {
  expect_equal(
    space_removal(
      SIREmulators$ems,
      SIREmulators$targets,
    ),
    c(nS = 29/30, nI = 9/10, nR = 1)
  )
  expect_equal(
    space_removal(
      SIREmulators$ems,
      SIREmulators$targets,
      cutoff = 4.5
    ),
    c(nS = 29/30, nI = 26/30, nR = 1)
  )
  expect_equal(
    space_removal(
      SIREmulators$ems,
      SIREmulators$targets,
      individual = FALSE
    ),
    1
  )
})

test_that("Validation diagnostics - validation set", {
  v1 <- validation_diagnostics(
    SIREmulators$ems,
    SIREmulators$targets,
    SIRSample$validation,
    plt = FALSE
  )
  expect_equal(
    nrow(v1),
    0
  )
  expect_warning(
    v2 <- validation_diagnostics(
      SIREmulators$ems,
      validation = SIRSample$validation,
      plt = FALSE
    )
  )
  expect_equal(
    nrow(v2),
    8
  )
  expect_warning(
    v3 <- validation_diagnostics(
      SIREmulators$ems,
      SIREmulators$targets,
      SIRSample$validation,
      which_diag = c('cd', 'ce', 'se', 'ft'),
      target_viz = "hatched"
    )
  )
  expect_equal(
    nrow(v3),
    0
  )
})

test_that("Validation diagnostics - no validation set", {
  v1 <- validation_diagnostics(
    SIREmulators$ems,
    SIREmulators$targets,
    k = 4,
    target_viz = "interval"
  )
  expect_equal(
    nrow(v1),
    0
  )
  v2 <- validation_diagnostics(
    SIREmulators$ems,
    which_diag = c('cd', 'se'),
    k = 15,
    target_viz = "solid"
  )
  expect_true(
    nrow(v2) >= 0
  )
})

v_em <- emulator_from_data(BirthDeath$training, c('Y'),
                                    list(lambda = c(0, 0.08), mu = c(0.04, 0.13)),
                                    verbose = FALSE, emulator_type = "variance")
v_targs <- list(Y = c(90, 110))

test_that("Variance emulator validation", {
  vv1 <- validation_diagnostics(
    v_em,
    v_targs,
    BirthDeath$validation,
    plt = FALSE
  )
  vv2 <- validation_diagnostics(
    v_em,
    v_targs,
    k = 10,
    plt = FALSE
  )
  expect_equal(
    nrow(vv1),
    nrow(vv2)
  )
})

bim_em <- emulator_from_data(SIR_stochastic$training,
                                     c('I10', 'I25', 'I50',
                                       'R10', 'R25', 'R50'),
                                     list(aSI = c(0.1, 0.8),
                                          aIR = c(0, 0.5),
                                          aSR = c(0, 0.05)),
                                     verbose = FALSE, emulator_type = "multistate")
bim_targets <- list(
  I10 = list(val = 35, sigma = 3.5),
  I25 = list(val = 147, sigma = 14.7),
  I50 = list(val = 55, sigma = 5.5),
  R10 = list(val = 29, sigma = 2.9),
  R25 = list(val = 276, sigma = 27.6),
  R50 = list(val = 579, sigma = 57.9)
)

test_that("Bimodal emulation validation", {
  skip_on_cran()
  vb1 <- validation_diagnostics(
    bim_em,
    bim_targets,
    SIR_stochastic$validation
  )
  expect_true(
    nrow(vb1) > 0
  )
  vb2 <- validation_diagnostics(
    bim_em,
    bim_targets,
    k = 10
  )
  expect_true(
    nrow(vb2) > 0
  )
})

test_that("Individual errors", {
  em <- SIREmulators$ems[[2]]
  i1 <- individual_errors(em, SIRSample$validation)
  i2 <- individual_errors(em, SIRSample$validation, "chol", "em")
  i3 <- individual_errors(em, SIRSample$validation, "eigen", plottype = "qq")
  i4 <- individual_errors(em, SIRSample$validation, "cholpivot", xtype = "aSI")
  expect_equal(
    nrow(i1),
    60
  )
  expect_equal(
    nrow(i2),
    nrow(i1)
  )
  expect_equal(
    nrow(i3),
    nrow(i2)
  )
  expect_equal(
    nrow(i4),
    nrow(i3)
  )
  expect_warning(
    individual_errors(
      em, SIRSample$validation,
      plottype = 'qq', errtype = 'normal'
    )
  )
  expect_warning(
    individual_errors(
      em, SIRSample$validation,
      errtype = 'eigen', xtype = 'em'
    )
  )
})

test_that("Alias functions", {
  expect_equal(
    nrow(
      classification_diag(
        SIREmulators$ems[[1]],
        SIREmulators$targets,
        SIRSample$validation,
        plt = FALSE
      )
    ),
    0
  )
  expect_equal(
    nrow(
      comparison_diag(
        SIREmulators$ems[[2]],
        SIREmulators$targets,
        SIRSample$validation,
        plt = FALSE
      )
    ),
    0
  )
  expect_equal(
    nrow(
      standard_errors(
        SIREmulators$ems[[1]],
        SIREmulators$targets,
        SIRSample$validation,
        plt = FALSE
      )
    ),
    0
  )
})

test_that("Automated Diagnostics - all pass", {
  new_ems <- diagnostic_pass(SIREmulators$ems,
                             SIREmulators$targets,
                             SIRSample$validation)
  expect_equal(new_ems[[1]]$beta, SIREmulators$ems[[1]]$beta)
  expect_equal(new_ems[[2]]$beta, SIREmulators$ems[[2]]$beta)
  expect_equal(new_ems[[3]]$beta, SIREmulators$ems[[3]]$beta)
})

test_that("Automated Diagnostics - modify sigma", {
  smaller_sigma_ems <- purrr::map(SIREmulators$ems,
                                  ~.$mult_sigma(0.2))
  new_ems <- diagnostic_pass(smaller_sigma_ems,
                             SIREmulators$targets,
                             SIRSample$validation)
  expect_true(all(purrr::map_dbl(smaller_sigma_ems, "u_sigma") <= purrr::map_dbl(new_ems, "u_sigma")))
})

all_pts <- do.call('rbind.data.frame', SIRSample)
all_pts_by_input <- all_pts[order(all_pts$aSI),]
new_ems_by_input <- emulator_from_data(
  all_pts_by_input[1:30,],
  names(SIREmulators$targets),
  SIREmulators$ems[[1]]$ranges
)
all_pts_by_output <- all_pts[order(all_pts$nR, decreasing = TRUE),]
new_ems_by_output <- emulator_from_data(
  all_pts_by_output[1:30,],
  names(SIREmulators$targets),
  SIREmulators$ems[[1]]$ranges
)

test_that("Automated Diagnostics: trained only on subset of input", {
  fixed_input_ems <- diagnostic_pass(new_ems_by_input, SIREmulators$targets, all_pts_by_input[31:90,], threshhold = 0.3)
  expect_equal(length(fixed_input_ems), 1)
  expect_equal(nrow(validation_diagnostics(fixed_input_ems, SIREmulators$targets, all_pts_by_input[31:90,], plt = FALSE)), 0)
})

test_that("Automated Diagnostics: trained only on subset of output", {
  fixed_output_ems <- diagnostic_pass(new_ems_by_output, SIREmulators$targets, all_pts_by_output[31:90,], threshhold = 0.25)
  expect_equal(length(fixed_output_ems), 0)
})

test_that("Automated Diagnostics: checking output suitability", {
  new_ems <- diagnostic_pass(SIREmulators$ems,
                             SIREmulators$targets,
                             SIRSample$validation, check_output = TRUE)
  expect_equal(length(new_ems), 3)
})
