### Gamma Exponential correlation function

test_that("correlation with self is 1", {
  expect_equal(
    c(gamma_exp(
      data.frame(a = 1),
      data.frame(a = 1),
      list(theta = 0.1, gamma = 1)
    ), use.names = FALSE),
    1
  )
  expect_equal(
    c(diag(gamma_exp(
      data.frame(a = c(1, 2, 3), b = c(0.1, 0.4, 0.3)),
      data.frame(a = c(1, 2, 3), b = c(0.1, 0.4, 0.3)),
      list(theta = 0.2, gamma = 1.5)
    )), use.names = FALSE),
    rep(1, 3)
  )
})

test_that("one-dimensional gamma-exponential; single point", {
  expect_equal(
    c(gamma_exp(
      data.frame(a = 1),
      data.frame(a = 2),
      list(theta = 0.1, gamma = 1.2)
    )),
    1.3088694e-07,
    tolerance = 1e-7)
})

test_that("one-dimensional gamma-exponential; multi point", {
  expect_equal(
    gamma_exp(
      data.frame(a = c(1, 2)),
      data.frame(a = c(1.1, 2.9)),
      list(theta = 0.4, gamma = 0.6)
    ),
    matrix(c(6.470865e-01, 0.196575705, 7.832213e-02, 0.196575705),
           nrow = 2, byrow = TRUE),
    tolerance = 1e-7
    )
})

test_that("multi-dimensional gamma-exponential; single point", {
  expect_equal(
    c(gamma_exp(
      data.frame(a = 1, b = 2, c = -1),
      data.frame(a = 1.5, b = 2.9, c = -0.7),
      list(theta = 0.2, gamma = 1.1)
    )),
    0.0017601453,
    tolerance = 1e-6
  )
})

test_that("multi-dimensional gamma-exponential; multi point", {
  expect_equal(
    gamma_exp(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.8, 2.4, 3.2), b = c(0.5, 0, -0.5)),
      list(theta = 1, gamma = 0.7)
    ),
    matrix(
      c(0.5840054, 0.4694570, 0.24870724,
        0.5357538, 0.6397463, 0.34877791,
        0.2764781, 0.3274293, 0.38879390),
      nrow = 3, byrow = TRUE
    ),
    tolerance = 1e-7
  )
})

test_that("dimensionality checks", {
  expect_equal(
    dim(gamma_exp(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.8, 2.4), b = c(0.5, 0)),
      list(theta = 1, gamma = 2)
    )),
    c(2,3)
  )
  expect_equal(
    dim(gamma_exp(
      data.frame(a = c(1.8), b = c(0.5)),
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      list(theta = 1, gamma = 1.4)
    )),
    c(3,1)
  )
  expect_equal(
    dim(gamma_exp(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.8), b = c(0.5)),
      list(theta = 1, gamma = 1.2)
    )),
    c(1,3)
  )
})

test_that("same points gives symmetric matrix", {
  corr_out <- gamma_exp(
    data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
    data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
    list(theta = 0.2, gamma = 0.9)
  )
  expect_equal(
    corr_out,
    t(corr_out)
  )
})

test_that("fails with no theta", {
  expect_error(
    gamma_exp(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      list(gamma = 2)
    )
  )
})

test_that("fails with no gamma", {
  expect_error(
    gamma_exp(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      list(theta = 0.5)
    )
  )
})

test_that("fails with missing data.frame", {
  expect_error(
    gamma_exp(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      hp = list(theta = 0.1, gamma = 0.7)
    )
  )
})

test_that("fails if gamma not correctly specified", {
  expect_error(
    gamma_exp(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      hp = list(theta = 0.1, gamma = 2.1)
    )
  )
  expect_error(
    gamma_exp(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      hp = list(theta = 0.1, gamma = 0)
    )
  )
})

test_that("works with data.matrix or data.frame", {
  df <- data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4))
  dm <- data.matrix(df)
  expect_equal(
    gamma_exp(df, df, list(theta = 0.2, gamma = 1.5)),
    gamma_exp(dm, dm, list(theta = 0.2, gamma = 1.5))
  )
})
