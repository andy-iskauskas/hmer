## From plotting.R

test_that("Expectation plot", {
  skip_on_cran()
  g <- emulator_plot(SIREmulators$ems,
                     ppd = 10)
  g_other_pars <- emulator_plot(SIREmulators$ems,
                                params = c('aIR', 'aSR'))
  g_single_em <- emulator_plot(SIREmulators$ems[[1]])
  expect_equal(
    length(g$plots),
    8
  )
})

test_that("Variance plot", {
  skip_on_cran()
  g <- emulator_plot(SIREmulators$ems,
                     ppd = 10, plot_type = 'var')
  g_other_fixed <- emulator_plot(SIREmulators$ems,
                                plot_type = 'var', ppd = 10,
                                fixed_vals = c(aSR = 0.04))
  g_single_em <- emulator_plot(SIREmulators$ems[[1]],
                               plot_type = 'var')
  g_sd <- emulator_plot(SIREmulators$ems, plot_type = 'sd',
                        ppd = 10)
  expect_equal(
    length(g$plots),
    8
  )
})

test_that("Implausibility plot", {
  skip_on_cran()
  g <- emulator_plot(SIREmulators$ems,
                     ppd = 10, plot_type = "imp",
                     targets = SIREmulators$targets)
  g_nimp <- emulator_plot(SIREmulators$ems,
                          ppd = 10, plot_type = 'nimp',
                          targets = SIREmulators$targets,
                          cb = TRUE)
  expect_equal(
    length(g$plots),
    8
  )
  expect_warning(
    g_4imp <- emulator_plot(SIREmulators$ems,
                          ppd = 10, plot_type = 'nimp',
                          targets = SIREmulators$targets,
                          nth = 4)
  )
})

test_that("Plots from direct call", {
  skip_on_cran()
  em <- SIREmulators$ems$nS
  g_exp <- exp_plot(em, ppd = 10)
  g_var <- var_plot(em, ppd = 10)
  g_var_const <- var_plot(em$o_em, ppd = 10)
  expect_true(
    length(unique(g_var_const$data$V)) == 1
  )
})

test_that("Call from Emulator object", {
  skip_on_cran()
  g_orig <- emulator_plot(SIREmulators$ems$nS, ppd = 10)
  g <- plot(SIREmulators$ems$nS, ppd = 10)
  g2 <- SIREmulators$ems$nS$plot(ppd = 10)
  expect_equal(
    g_orig$data,
    g$data
  )
  expect_equal(
    g$data,
    g2$data
  )
})

test_that("Variance and Expectation plotting", {
  skip_on_cran()
  v_ems <- emulator_from_data(BirthDeath$training, c('Y'),
                                       list(
                                         lambda = c(0, 0.08),
                                         mu = c(0.04, 0.13)),
                                       verbose = FALSE,
                              emulator_type = "variance")
  g_v <- emulator_plot(v_ems, ppd = 10)
  expect_match(
    g_v$plots[[2]]$labels$title,
    "Y Mean Emulator Expectation"
  )
  g_v2 <- emulator_plot(v_ems$expectation, ppd = 10)
  expect_equal(
    g_v$plots[[2]]$data,
    g_v2$plots[[2]]$data
  )
  g_vv <- v_ems$variance$Y$plot(ppd = 10)
  expect_equal(
    g_vv$labels$title,
    "Y Variance Emulator Expectation"
  )
})

test_that("Output plotting", {
  skip_on_cran()
  g <- output_plot(
    SIREmulators$ems,
    SIREmulators$targets,
    npoints = 500
  )
  expect_equal(
    nrow(g$data),
    500 * length(SIREmulators$ems)
  )
  g_pts <- output_plot(
    SIREmulators$ems,
    SIREmulators$targets,
    points = SIRSample$validation
  )
  expect_equal(
    nrow(g_pts$data),
    60 * length(SIREmulators$ems)
  )
})

test_that("Plot lattice", {
  skip_on_cran()
  pl <- plot_lattice(
    SIREmulators$ems,
    SIREmulators$targets,
    ppd = 10
  )
  expect_equal(
    length(pl$plots),
    length(SIREmulators$ems)^2
  )
  pl_alt <- plot_lattice(
    SIREmulators$ems,
    SIREmulators$targets,
    ppd = 20,
    maxpoints = 4000,
    cb = TRUE
  )
  expect_true(
    is.list(pl_alt$plots[[1]]$data)
  )
})

test_that("Active Variable plots", {
  skip_on_cran()
  p1 <- plot_actives(
    SIREmulators$ems
  )
  expect_equal(
    nrow(p1$data),
    9
  )
  expect_equal(
    levels(p1$data$Var1),
    c("nS", "nI", "nR")
  )
  expect_equal(
    levels(p1$data$Var2),
    c('aSI', 'aIR', 'aSR')
  )
  expect_true(
    all(p1$data$value %in% c("TRUE", "FALSE"))
  )
  p2 <- plot_actives(
    SIREmulators$ems,
    c('nS'),
    c('aSI', 'aIR')
  )
  expect_equal(
    nrow(p2$data),
    2
  )
})

## Diagnostic plots
test_that("Behaviour plot error and warning", {
  skip_on_cran()
  expect_error(
    behaviour_plot(out_names = c('nS', 'nI', 'nR'),
                   targets = SIREmulators$targets),
    "One of 'ems' or 'points' must be supplied."
  )
  expect_error(
    behaviour_plot(points = SIRSample$training, model = FALSE),
    "Cannot perform emulator expectation"
  )
  expect_error(
    behaviour_plot(points = SIRSample$validation, targets = SIREmulators$targets),
    "No output names"
  )
  expect_warning(
    behaviour_plot(SIREmulators$ems, model = TRUE),
    "Cannot do model output"
  )
  expect_warning(
    behaviour_plot(SIREmulators$ems, SIRSample$validation[,-c(5)], model = TRUE),
    "Not all outputs"
  )
})

test_that("Behaviour plot behaviour", {
  skip_on_cran()
  p1 <- behaviour_plot(SIREmulators$ems, targets = SIREmulators$targets)
  p2 <- behaviour_plot(points = SIRSample$validation,
                       targets = SIREmulators$targets,
                       out_names = names(SIREmulators$targets))
  expect_true(
    is.null(p1) && is.null(p2)
  )
})

test_that("Space removed plot", {
  skip_on_cran()
  g1 <- space_removed(
    SIREmulators$ems, SIREmulators$targets, ppd = 5
  )
  expect_equal(
    nrow(g1$data),
    1000
  )
  expect_warning(
    g2 <- space_removed(
      SIREmulators$ems, SIREmulators$targets, ppd = 5,
      u_mod = c(0.8, 1, 1.2), intervals = seq(0, 4, length.out = 100),
      modified = 'disc'
    ),
    "'disc' chosen"
  )
  expect_equal(
    nrow(g2$data),
    300
  )
  g_reduced <- space_removed(
    SIREmulators$ems, SIREmulators$targets, ppd = 10,
    modified = 'var', maxpoints = 990, u_mod = c(1, 1.1, 1.2),
    intervals = seq(0, 3, length.out = 80)
  )
  expect_equal(
    nrow(g_reduced$data),
    240
  )
  g_hp <- space_removed(
    SIREmulators$ems, SIREmulators$targets, ppd = 5,
    modified = 'hp'
  )
  expect_equal(
    nrow(g_hp$data),
    1000
  )
})

test_that("Validation pairs plot",  {
  skip_on_cran()
  v1 <- validation_pairs(SIREmulators$ems, SIRSample$validation,
                         SIREmulators$targets, nth = 2)
  v2 <- validation_pairs(SIRMultiWaveEmulators[[3]],
                         SIRMultiWaveData[[4]],
                         SIREmulators$targets,
                         SIREmulators$ems[[1]]$ranges,
                         cb = TRUE)
  expect_equal(
    nrow(v1$data),
    60
  )
  expect_equal(
    nrow(v2$data),
    90
  )
})

test_that("Emulator effect strength", {
  skip_on_cran()
  e1 <- effect_strength(SIREmulators$ems, plt = FALSE, quadratic = FALSE)
  e2 <- effect_strength(SIREmulators$ems, plt = FALSE)
  e_plot <- effect_strength(SIREmulators$ems,
                            line.plot = TRUE,
                            xvar = TRUE)
  e_plot2 <- effect_strength(SIREmulators$ems,
                             grid.plot = TRUE)
  expect_true(
    !is.null(e2$linear) && !is.null(e2$quadratic)
  )
  expect_equal(
    dim(e1),
    c(3,3)
  )
  expect_equal(
    e_plot$linear,
    e_plot2$linear
  )
})

test_that("Multi-wave: diagnostic_wrap", {
  skip_on_cran()
  mw <- diagnostic_wrap(
    SIRMultiWaveData,
    SIREmulators$targets,
    input_names = c('aSI', 'aIR'),
    output_names = c('nS', 'nR'),
    palette = c('red', 'blue', 'green', 'yellow'),
    wave_numbers = 1:3,
    p_size = 0.8,
    l_wid = 0.8,
    upper_scale = 1.2,
  )
  expect_true(
    all(names(mw) %in%
    c("simulatorplot", "simulatorplotnorm",
      "simulatorplotlog", "posteriorplot",
      "outputsplot", "dependencyplot",
      "dependencyplotnorm")
    )
  )
  surr <- wave_points(
    SIRMultiWaveData,
    c('aSI', 'aIR', 'aSR'),
    surround = TRUE,
    wave_numbers = c(0, 2, 3)
  )
  expect_equal(
    nrow(surr$data),
    nrow(do.call('rbind.data.frame', SIRMultiWaveData[c(0,2,3)+1]))
  )
  expect_warning(
    wv <- wave_values(
      SIRMultiWaveData, SIREmulators$targets,
      ems = SIRMultiWaveEmulators[[3]], restrict = TRUE
    ),
    "Expecting to restrict"
  )
})
