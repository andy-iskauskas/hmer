# Kurtosis check

points <- seq(0, 4, by = 0.01)

test_that("works on vectors", {
  expect_equal(
    kurtosis(points),
    9/5,
    tolerance = 1e-3)
})

test_that("works on data.frames", {
  expect_equal(
    kurtosis(points),
    c(kurtosis(data.frame(x = points)),
      use.names = FALSE)
  )
})

test_that("works on matrices", {
  expect_equal(
    kurtosis(points),
    kurtosis(matrix(points, ncol = 1))
  )
})

test_that("na gives na", {
  points[30] <- NA
  expect_true(
    is.na(kurtosis(points))
  )
})

test_that("removing na works", {
  npoints <- points
  npoints[30] <- NA
  expect_true(
    !is.na(kurtosis(npoints, na.rm = TRUE))
  )
  expect_false(
    kurtosis(points) == kurtosis(npoints,na.rm = TRUE)
  )
})
