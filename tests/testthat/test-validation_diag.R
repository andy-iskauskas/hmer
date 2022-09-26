test_that("Diagnostic summary measures", {
  expect_equal(
    summary_diag(
      SIREmulators$ems[[1]],
      SIRSample$validation,
      verbose = FALSE
    ),
    c(FALSE, TRUE)
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
    c(nS = 1/30, nI = 1/10, nR = 0)
  )
  expect_equal(
    space_removal(
      SIREmulators$ems,
      SIREmulators$targets,
      cutoff = 4.5
    ),
    c(nS = 1/30, nI = 4/30, nR = 0)
  )
  expect_equal(
    space_removal(
      SIREmulators$ems,
      SIREmulators$targets,
      individual = FALSE
    ),
    0
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
    6
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

v_em <- variance_emulator_from_data(BirthDeath$training, c('Y'),
                                    list(lambda = c(0, 0.08), mu = c(0.04, 0.13)),
                                    verbose = FALSE)
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

bim_em <- bimodal_emulator_from_data(SIR_stochastic$training,
                                     c('I10', 'I25', 'I50',
                                       'R10', 'R25', 'R50'),
                                     list(aSI = c(0.1, 0.8),
                                          aIR = c(0, 0.5),
                                          aSR = c(0, 0.05)),
                                     verbose = FALSE)
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
  em <- SIREmulators$ems[[1]]
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
