ems <- SIREmulators$ems
targs <- SIREmulators$targets

test_that("Point generation methods", {
  skip_on_cran()
  points <- generate_new_design(
    ems, 100, targs, verbose = FALSE
  )
  expect_equal(
    nrow(points),
    100
  )
  expect_warning(
    points_no_clust <- generate_new_design(
      ems, 100, targs, verbose = FALSE,
      resample = 0, cluster = TRUE, method = 'lhs'
    ),
    "Cannot distinguish"
  )
  points_line <- generate_new_design(
    ems, 100, targs, method = c('line'),
    plausible_set = points, cutoff = 2.5,
    resample = 0, verbose = FALSE
  )
  expect_equal(
    nrow(points_line),
    100
  )
  points_slice <- generate_new_design(
    ems, 100, targs, method = c('slice'),
    plausible_set = points_line[1:50,], cutoff = 3,
    resample = 0, verbose = FALSE
  )
  expect_equal(
    nrow(points_slice),
    100
  )
  points_importance <- generate_new_design(
    ems, 100, targs, method = c('importance'),
    plausible_set = points_line[1:50,], cutoff = 4,
    resample = 0, verbose = FALSE
  )
  expect_equal(
    nrow(points_importance),
    100
  )
  points_importance_sphere <- generate_new_design(
    ems, 100, targs, method = c("importance"),
    plausible_set = points_line[1:50,], cutoff = 4,
    resample = 0, verbose = FALSE, imp_distro = "normal"
  )
  expect_equal(
    nrow(points_importance_sphere),
    100
  )
  points_seek_good <- generate_new_design(
    ems, 100, targs, cutoff = 4, seek = 5,
    verbose = FALSE
  )
  expect_equal(
    nrow(points_seek_good),
    100
  )
})

test_that("Forced laddering of implausibility", {
  skip_on_cran()
  bad_targets <- list(
    nS = c(680, 751),
    nI = list(val = 229, sigma = 8.45),
    nR = c(199, 221)
  )
  g1 <- generate_new_design(ems, 100, bad_targets, verbose = FALSE)
  expect_equal(
    nrow(g1),
    100
  )
  really_bad_targets <- list(
    nS = c(780, 851),
    nI = list(val = 229, sigma = 8.45),
    nR = c(199, 221)
  )
  g2 <- generate_new_design(ems, 100, really_bad_targets,
                          resample = 0, verbose = FALSE)
  expect_true(
    nrow(g2) < 100
  )
})

test_that("Optional pca arguments", {
  skip_on_cran()
  lhd_pca_prop <- generate_new_design(ems, 100, targs, verbose = FALSE,
                                    opts = list(pca_lhs = TRUE))
  expect_equal(nrow(lhd_pca_prop), 100)
  slice_pca_prop <- generate_new_design(ems, 100, targs, method = c('lhs', 'slice'),
                                      verbose = FALSE, opts = list(pca_slice = TRUE))
  expect_equal(nrow(slice_pca_prop), 100)
})

test_that("Handling options via opts", {
  skip_on_cran()
  all_options <- generate_new_design(ems, 100, targs, verbose = FALSE, opts = list(
    accept_measure = "default", cluster = FALSE, cutoff_tolerance = 0.01,
    ladder_tolerance = 0.1, nth = 1, resample = 0, seek = 0
  ))
  expect_true(nrow(all_options) == 100)
})
