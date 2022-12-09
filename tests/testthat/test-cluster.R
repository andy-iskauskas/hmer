fake_custom <- function(ems, x, z, cutoff, ...) {
  nth_implausible(ems, x, z, n = 1, cutoff = cutoff, ...)
}

test_that("Cluster generation disjoint", {
  skip_on_cran()
  test_cluster <- generate_new_runs(data$ems_split, 200, data$targets_split, verbose = FALSE,
                                    opts = list(nth = 1, cluster = TRUE))
  expect_true(
    max(nth_implausible(data$ems_split, test_cluster, data$targets_split)) <= 3
  )
})

test_that("Cluster generation disjoint - custom implausibility", {
  skip_on_cran()
  test_cluster_custom <- generate_new_runs(data$ems_split, 200, data$targets_split,
                                           verbose = FALSE, opts = list(nth = 1,
                                                                        cluster = TRUE,
                                                                        accept_measure = fake_custom))
  expect_true(
    max(nth_implausible(data$ems_split, test_cluster_custom, data$targets_split)) <= 3
  )
})

test_that("Cluster generation shared actives", {
  skip_on_cran()
  test_cluster_shared <- generate_new_runs(data$ems_unsplit, 200, data$targets_unsplit,
                                           verbose = FALSE, opts = list(nth = 1,
                                                                        cluster = TRUE))
  expect_true(
    max(nth_implausible(data$ems_unsplit, test_cluster_shared, data$targets_unsplit)) <= 3
  )
})

test_that("Cluster generation shared - custom implausibility", {
  skip_on_cran()
  rest_ems <- subset_emulators(data$ems_unsplit, c("I1100", "R110", "R1100",
                                                   "I2100", "R210", "R2100"))
  test_cluster_custom_shared <- generate_new_runs(rest_ems, 200, data$targets_unsplit,
                                                  verbose = FALSE, opts = list(nth = 1,
                                                                               cluster = TRUE,
                                                                               accept_measure = fake_custom))
  expect_true(
    max(nth_implausible(rest_ems, test_cluster_custom_shared, data$targets_unsplit)) <= 3
  )
})

test_that("No cluster determined", {
  skip_on_cran()
  rest_ems2 <- subset_emulators(data$ems_unsplit, c("R110", "R1100",
                                                    "R210", "R2100"))
  expect_warning(
    tcnc <- generate_new_runs(rest_ems2, 200, data$targets_unsplit,
                              verbose = FALSE, opts = list(nth = 1, cluster = TRUE))
  )
})

