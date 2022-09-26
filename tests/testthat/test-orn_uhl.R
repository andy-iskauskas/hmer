### Ornstein-Uhlenbeck correlation function

test_that("correlation with self is 1", {
  expect_equal(
    c(orn_uhl(
      data.frame(a = 1),
      data.frame(a = 1),
      list(theta = 0.1)
    ), use.names = FALSE),
    1
  )
  expect_equal(
    unname(diag(orn_uhl(
      data.frame(a = c(1, 2, 3), b = c(0.1, 0.4, 0.3)),
      data.frame(a = c(1, 2, 3), b = c(0.1, 0.4, 0.3)),
      list(theta = 0.2)
    ))),
    rep(1, 3)
  )
})

test_that("one-dimensional ornstein-uhlenbeck; single point", {
  expect_equal(
    c(orn_uhl(
      data.frame(a = 1),
      data.frame(a = 2),
      list(theta = 0.1)
    ), use.names = FALSE),
    4.539993e-05,
    tolerance = 1e-7)
})

test_that("one-dimensional ornstein-uhlenbeck; multi point", {
  expect_equal(
    orn_uhl(
      data.frame(a = c(1, 2)),
      data.frame(a = c(1.1, 2.9)),
      list(theta = 0.4)
    ),
    matrix(c(0.778800783, 0.1053992, 0.008651695, 0.1053992),
           nrow = 2, byrow = TRUE),
    tolerance = 1e-7
    )
})

test_that("multi-dimensional ornstein-uhlenbeck; single point", {
  expect_equal(
    c(orn_uhl(
      data.frame(a = 1, b = 2, c = -1),
      data.frame(a = 1.5, b = 2.9, c = -0.7),
      list(theta = 0.2)
    ), use.names = FALSE),
    0.00469197,
    tolerance = 1e-6
  )
})

test_that("multi-dimensional ornstein-uhlenbeck; multi point", {
  expect_equal(
    orn_uhl(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.8, 2.4, 3.2), b = c(0.5, 0, -0.5)),
      list(theta = 1)
    ),
    matrix(
      c(0.6621186, 0.5112889, 0.2012672,
        0.6005545, 0.7288934, 0.3406046,
        0.2388828, 0.3102211, 0.3977409),
      nrow = 3, byrow = TRUE
    ),
    tolerance = 1e-6
  )
})

test_that("dimensionality checks", {
  expect_equal(
    dim(orn_uhl(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.8, 2.4), b = c(0.5, 0)),
      list(theta = 1)
    )),
    c(2,3)
  )
  expect_equal(
    dim(orn_uhl(
      data.frame(a = c(1.8), b = c(0.5)),
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      list(theta = 1)
    )),
    c(3,1)
  )
  expect_equal(
    dim(orn_uhl(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.8), b = c(0.5)),
      list(theta = 1)
    )),
    c(1,3)
  )
})

test_that("same points gives symmetric matrix", {
  corr_out <- orn_uhl(
    data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
    data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
    list(theta = 0.2)
  )
  expect_equal(
    corr_out,
    t(corr_out)
  )
})

test_that("fails with no theta", {
  expect_error(
    orn_uhl(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4))
    )
  )
})

test_that("fails with missing data.frame", {
  expect_error(
    orn_uhl(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      hp = list(theta = 0.1)
    )
  )
})

test_that("works with data.matrix or data.frame", {
  df <- data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4))
  dm <- data.matrix(df)
  expect_equal(
    orn_uhl(df, df, list(theta = 0.2)),
    orn_uhl(dm, dm, list(theta = 0.2))
  )
})

test_that("orh-uhl is the same as matern for nu = 0.5", {
  df1 <- data.frame(a = runif(10), b = runif(10, -1, 1))
  df2 <- data.frame(a = runif(5), b = runif(5, -1, 1))
  expect_equal(
    matern(df1, df2, list(theta = 0.2, nu = 0.5)),
    orn_uhl(df1, df2, list(theta = 0.2))
  )
})
