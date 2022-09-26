### Matern correlation function

test_that("correlation with self is 1", {
  expect_equal(
    c(matern(
      data.frame(a = 1),
      data.frame(a = 1),
      list(theta = 0.1, nu = 1.5)
    ), use.names = FALSE),
    1
  )
  expect_equal(
    c(diag(matern(
      data.frame(a = c(1, 2, 3), b = c(0.1, 0.4, 0.3)),
      data.frame(a = c(1, 2, 3), b = c(0.1, 0.4, 0.3)),
      list(theta = 0.2, nu = 1.5)
    )), use.names = FALSE),
    rep(1, 3)
  )
})

test_that("one-dimensional matern; single point", {
  expect_equal(
    c(matern(
      data.frame(a = 1),
      data.frame(a = 2),
      list(theta = 0.1, nu = 2.5)
    ), use.names = FALSE),
    3.6956962e-08)
})

test_that("one-dimensional matern; multi point", {
  expect_equal(
    matern(
      data.frame(a = c(1, 2)),
      data.frame(a = c(1.1, 2.9)),
      list(theta = 0.4, nu = 2.5)
    ),
    matrix(c(0.950959922, 0.09449877, 0.001200627, 0.09449877),
           nrow = 2, byrow = TRUE),
    tolerance = 1e-7
    )
})

test_that("multi-dimensional matern; single point", {
  expect_equal(
    c(matern(
      data.frame(a = 1, b = 2, c = -1),
      data.frame(a = 1.5, b = 2.9, c = -0.7),
      list(theta = 0.2, nu = 0.5)
    ), use.names = FALSE),
    0.0046919704,
    tolerance = 1e-6
  )
})

test_that("multi-dimensional matern; multi point", {
  expect_equal(
    matern(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.8, 2.4, 3.2), b = c(0.5, 0, -0.5)),
      list(theta = 1, nu = 1.5)
    ),
    matrix(
      c(0.8392642, 0.6764411, 0.2350773,
        0.7786323, 0.8949942, 0.4436402,
        0.2914432, 0.3986634, 0.5259420),
      nrow = 3, byrow = TRUE
    ),
    tolerance = 1e-7
  )
})

test_that("dimensionality checks", {
  expect_equal(
    dim(matern(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.8, 2.4), b = c(0.5, 0)),
      list(theta = 1, nu = 1.5)
    )),
    c(2,3)
  )
  expect_equal(
    dim(matern(
      data.frame(a = c(1.8), b = c(0.5)),
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      list(theta = 1, nu = 1.5)
    )),
    c(3,1)
  )
  expect_equal(
    dim(matern(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.8), b = c(0.5)),
      list(theta = 1, nu = 2.5)
    )),
    c(1,3)
  )
})

test_that("same points gives symmetric matrix", {
  corr_out <- matern(
    data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
    data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
    list(theta = 0.2, nu = 1.5)
  )
  expect_equal(
    corr_out,
    t(corr_out)
  )
})

test_that("fails with no theta", {
  expect_error(
    matern(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      list(nu = 1.5)
    )
  )
})

test_that("fails with no nu", {
  expect_error(
    matern(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      list(theta = 0.8)
    )
  )
})

test_that("fails with missing data.frame", {
  expect_error(
    matern(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      hp = list(theta = 0.1, nu = 1.5)
    )
  )
})

test_that("fails if nu not half-integer", {
  df <- df <- data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4))
  expect_error(
    matern(df, df, list(theta = 0.8, nu = 0.6))
  )
  expect_error(
    matern(df, df, list(theta = 0.8, nu = 2))
  )
})

test_that("works with data.matrix or data.frame", {
  df <- data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4))
  dm <- data.matrix(df)
  expect_equal(
    matern(df, df, list(theta = 0.2, nu = 1.5)),
    matern(dm, dm, list(theta = 0.2, nu = 1.5))
  )
})

## Derivative of matern correlation function

test_that("correlation with self is 0", {
  expect_equal(
    unname(matern_d(
      data.matrix(data.frame(a = 1)),
      data.matrix(data.frame(a = 1)),
      list(theta = 0.1, nu = 1.5),
      1
    )),
    matrix(0, nrow = 1)
  )
  expect_equal(
    c(diag(matern_d(
      data.matrix(data.frame(a = c(1, 2, 3), b = c(0.1, 0.4, 0.3))),
      data.matrix(data.frame(a = c(1, 2, 3), b = c(0.1, 0.4, 0.3))),
      list(theta = 0.2, nu = 1.5),
      2
    )), use.names = FALSE),
    rep(0, 3)
  )
})

test_that("one-dimensional matern derivative; single point", {
  expect_equal(
    matern_d(
      data.matrix(data.frame(a = 1)),
      data.matrix(data.frame(a = 2)),
      list(theta = 0.1, nu = 1.5),
      1
    ),
    matrix(-9.0140544e-06, nrow = 1))
})

test_that("one-dimensional matern derivative; multi point", {
  expect_equal(
    matern_d(
      data.matrix(data.frame(a = c(1, 2))),
      data.matrix(data.frame(a = c(1.1, 2.9))),
      list(theta = 0.4, nu = 2.5),
      1
    ),
    matrix(c(-0.928542145, 0.3692918, -0.005609912, -0.3692918),
           nrow = 2, byrow = TRUE),
    tolerance = 1e-6
  )
})

test_that("multi-dimensional matern derivative; single point", {
  expect_equal(
    matern_d(
      data.matrix(data.frame(a = 1, b = 2, c = -1)),
      data.matrix(data.frame(a = 1.5, b = 2.9, c = -0.7)),
      list(theta = 0.2, nu = 1.5),
      3
    ),
    matrix(-0.0020837784, nrow = 1),
    tolerance = 1e-5
  )
})

test_that("multi-dimensional matern derivative; multi point and multi deriv", {
  expect_equal(
    matern_d(
      data.matrix(data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4))),
      data.matrix(data.frame(a = c(1.8, 2.4, 3.2), b = c(0.5, 0, -0.5))),
      list(theta = 1, nu = 2.5),
      1,
      2
    ),
    matrix(
      c(0.132580306, 0.334695240, 0.036993702,
        0.133234551, -0.123267173, -0.299888029,
        0.264540759,  0.267678812, -0.190884317),
      nrow = 3, byrow = TRUE
    ),
    tolerance = 1e-6
  )
})

test_that("fails with no derivative direction", {
  expect_error(
    matern_d(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4))
    )
  )
})

test_that("fails if matern not differentiable", {
  df <- data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4))
  expect_error(
    matern_d(df, df, list(theta = 0.8, nu = 0.5), 1)
  )
  expect_error(
    matern_d(df, df, list(theta = 0.8, nu = 1.5), 1, 2)
  )
})
