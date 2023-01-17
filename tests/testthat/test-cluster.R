all_ranges <- list(
  alpha1 = c(0, 0.3),
  beta1 = c(0, 0.05),
  gamma1 = c(0, 0.005),
  alpha2 = c(0, 0.4),
  beta2 = c(0, 0.025),
  gamma2 = c(0, 0.01),
  birth = c(0, 0.001),
  death = c(0, 0.001)
)
ems_split <- emulator_from_data(split_data, names(targets_split), all_ranges[1:6],
                                verbose = FALSE, more_verbose = FALSE)
ems_unsplit <- emulator_from_data(unsplit_data, names(targets_unsplit), all_ranges,
                                  verbose = FALSE, more_verbose = FALSE)
fake_custom <- function(ems, x, z, cutoff, n = 1, ...) {
  nth_implausible(ems, x, z, n = n, cutoff = cutoff, ...)
}

test_that("Cluster generation disjoint", {
  skip_on_cran()
  test_cluster <- generate_new_runs(ems_split, 100, targets_split, verbose = FALSE,
                                    opts = list(nth = 1, cluster = TRUE))
  expect_true(
    max(nth_implausible(ems_split, test_cluster, targets_split)) <= 3
  )
})

test_that("Cluster generation disjoint - custom implausibility", {
  skip_on_cran()
  test_cluster_custom <- generate_new_runs(ems_split, 100, targets_split,
                                           verbose = FALSE, opts = list(nth = 1,
                                                                        cluster = TRUE,
                                                                        accept_measure = fake_custom))
  expect_true(
    max(nth_implausible(ems_split, test_cluster_custom, targets_split)) <= 3
  )
})

test_that("Cluster generation shared actives", {
  skip_on_cran()
  test_cluster_shared <- generate_new_runs(ems_unsplit, 100, targets_unsplit,
                                           verbose = FALSE, opts = list(nth = 1,
                                                                        cluster = TRUE))
  expect_true(
    max(nth_implausible(ems_unsplit, test_cluster_shared, targets_unsplit)) <= 3
  )
})

test_that("Cluster generation shared - custom implausibility", {
  skip_on_cran()
  rest_ems <- subset_emulators(ems_unsplit, c("I1100", "R110", "R1100",
                                                   "I2100", "R210", "R2100"))
  test_cluster_custom_shared <- generate_new_runs(rest_ems, 100, targets_unsplit,
                                                  verbose = FALSE, opts = list(nth = 1,
                                                                               cluster = TRUE,
                                                                               accept_measure = fake_custom))
  expect_true(
    max(nth_implausible(rest_ems, test_cluster_custom_shared, targets_unsplit)) <= 3
  )
})

test_that("Cluster generation - multi wave", {
  skip_on_cran()
  multi_ems <- list(ems_split, ems_split)
  test_multi <- generate_new_runs(multi_ems, 100, targets_split, verbose = FALSE,
                                  opts = list(nth = 1, cluster = TRUE))
  expect_true(
    max(nth_implausible(ems_split, test_multi, targets_split)) <= 3
  )
})

test_that("No cluster determined", {
  skip_on_cran()
  rest_ems2 <- subset_emulators(ems_unsplit, c("R110", "R1100",
                                                    "R210", "R2100"))
  expect_warning(
    tcnc <- generate_new_runs(rest_ems2, 100, targets_unsplit,
                              verbose = FALSE, opts = list(nth = 1, cluster = TRUE))
  )
})

