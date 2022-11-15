ems <- SIREmulators$ems
targs <- SIREmulators$targets
bad_targs <- targs
bad_targs$nS <- c(380, 451)

custom_measure <- function(ems, x, z, cutoff, ...) {
  imps_df <- nth_implausible(ems, x, z, get_raw = TRUE)
  sorted_imps <- t(apply(imps_df, 1, sort, decreasing = TRUE))
  imps1 <- sorted_imps[,1] <= cutoff
  imps2 <- sorted_imps[,2] <= cutoff - 0.5
  constraint <- apply(x, 1, function(y) y[[1]] <= 0.4)
  return(imps1 & imps2 & constraint)
}

test_that("Custom generation behaves", {
  skip_on_cran()
  points <- generate_new_runs(
    ems, 100, targs, verbose = FALSE,
    opts = list(accept_measure = custom_measure)
  )
  expect_equal(
    nrow(points),
    100
  )
  bad_points <- generate_new_runs(
    ems, 200, bad_targs, verbose = FALSE,
    opts = list(accept_measure = custom_measure)
  )
  expect_equal(
    nrow(bad_points),
    200
  )
})
