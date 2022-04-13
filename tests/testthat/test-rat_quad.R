### Rational Quadratic correlation function

test_that("correlation with self is 1", {
  expect_equal(
    rat_quad(
      data.frame(a = 1),
      data.frame(a = 1),
      list(theta = 0.1, alpha = 1.5)
    ),
    1
  )
  expect_equal(
    diag(rat_quad(
      data.frame(a = c(1, 2, 3), b = c(0.1, 0.4, 0.3)),
      data.frame(a = c(1, 2, 3), b = c(0.1, 0.4, 0.3)),
      list(theta = 0.2, alpha = 1.5)
    )),
    rep(1, 3)
  )
})

test_that("one-dimensional rational quadratic; single point", {
  expect_equal(
    rat_quad(
      data.frame(a = 1),
      data.frame(a = 2),
      list(theta = 0.1, alpha = 2.5)
    ),
    0.00049482515)
})

test_that("one-dimensional rational quadratic; multi point", {
  expect_equal(
    rat_quad(
      data.frame(a = c(1, 2)),
      data.frame(a = c(1.1, 2.9)),
      list(theta = 0.4, alpha = 2.5)
    ),
    matrix(c(0.96942099, 0.1740445, 0.01401614, 0.1740445),
           nrow = 2, byrow = TRUE),
    tolerance = 1e-7
    )
})

test_that("multi-dimensional rational quadratic; single point", {
  expect_equal(
    rat_quad(
      data.frame(a = 1, b = 2, c = -1),
      data.frame(a = 1.5, b = 2.9, c = -0.7),
      list(theta = 0.2, alpha = 0.5)
    ),
    0.1833397
  )
})

test_that("multi-dimensional rational quadratic; multi point", {
  expect_equal(
    rat_quad(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.8, 2.4, 3.2), b = c(0.5, 0, -0.5)),
      list(theta = 1, alpha = 1.5)
    ),
    matrix(
      c(0.9206467, 0.8108737, 0.3952748,
        0.8827861, 0.9520052, 0.6124095,
        0.4578728, 0.5688002, 0.6878453),
      nrow = 3, byrow = TRUE
    ),
    tolerance = 1e-7
  )
})

test_that("dimensionality checks", {
  expect_equal(
    dim(rat_quad(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.8, 2.4), b = c(0.5, 0)),
      list(theta = 1, alpha = 1.5)
    )),
    c(2,3)
  )
  expect_equal(
    dim(rat_quad(
      data.frame(a = c(1.8), b = c(0.5)),
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      list(theta = 1, alpha = 1.5)
    )),
    c(3,1)
  )
  expect_equal(
    dim(rat_quad(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.8), b = c(0.5)),
      list(theta = 1, alpha = 2.5)
    )),
    NULL
  )
})

test_that("same points gives symmetric matrix", {
  corr_out <- rat_quad(
    data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
    data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
    list(theta = 0.2, alpha = 1.5)
  )
  expect_equal(
    corr_out,
    t(corr_out)
  )
})

test_that("fails with no theta", {
  expect_error(
    rat_quad(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      list(alpha = 1.5)
    )
  )
})

test_that("fails with no alpha", {
  expect_error(
    rat_quad(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      list(theta = 0.8)
    )
  )
})

test_that("fails with missing data.frame", {
  expect_error(
    rat_quad(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      hp = list(theta = 0.1, alpha = 1.5)
    )
  )
})

test_that("works with data.matrix or data.frame", {
  df <- data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4))
  dm <- data.matrix(df)
  expect_equal(
    rat_quad(df, df, list(theta = 0.2, alpha = 1.5)),
    rat_quad(dm, dm, list(theta = 0.2, alpha = 1.5))
  )
})

## Derivative of rat_quad correlation function

test_that("correlation with self is 0", {
  expect_equal(
    rat_quad_d(
      data.matrix(data.frame(a = 1)),
      data.matrix(data.frame(a = 1)),
      list(theta = 0.1, alpha = 1.5),
      1
    ),
    matrix(0, nrow = 1)
  )
  expect_equal(
    diag(rat_quad_d(
      data.matrix(data.frame(a = c(1, 2, 3), b = c(0.1, 0.4, 0.3))),
      data.matrix(data.frame(a = c(1, 2, 3), b = c(0.1, 0.4, 0.3))),
      list(theta = 0.2, alpha = 1.5),
      2
    )),
    rep(0, 3)
  )
})

test_that("one-dimensional rational quadratic derivative; single point", {
  expect_equal(
    rat_quad_d(
      data.matrix(data.frame(a = 1)),
      data.matrix(data.frame(a = 2)),
      list(theta = 0.1, alpha = 1.5),
      1
    ),
    matrix(-0.01447805, nrow = 1))
})

test_that("one-dimensional rational quadratic derivative; multi point", {
  expect_equal(
    rat_quad_d(
      data.matrix(data.frame(a = c(1, 2))),
      data.matrix(data.frame(a = c(1.1, 2.9))),
      list(theta = 0.4, alpha = 2.5),
      1
    ),
    matrix(c(-0.5984080, 0.4864598, -0.0301935, -0.4864598),
           nrow = 2, byrow = TRUE),
    tolerance = 1e-6
  )
})

test_that("multi-dimensional rational quadratic derivative; single point", {
  expect_equal(
    rat_quad_d(
      data.matrix(data.frame(a = 1, b = 2, c = -1)),
      data.matrix(data.frame(a = 1.5, b = 2.9, c = -0.7)),
      list(theta = 0.2, alpha = 1.5),
      3
    ),
    matrix(-0.0205828295, nrow = 1)
  )
})

test_that("multi-dimensional rational quadratic derivative; multi point and multi deriv", {
  expect_equal(
    rat_quad_d(
      data.matrix(data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4))),
      data.matrix(data.frame(a = c(1.8, 2.4, 3.2), b = c(0.5, 0, -0.5))),
      list(theta = 1, alpha = 2.5),
      1,
      2
    ),
    matrix(
      c(0.04817765, 0.17099417, 0.03464827,
        0.05572207, -0.03841922, -0.21899785,
        0.23266799,  0.20716591, -0.12432663),
      nrow = 3, byrow = TRUE
    ),
    tolerance = 1e-7
  )
})

test_that("fails with no derivative direction", {
  expect_error(
    rat_quad_d(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4))
    )
  )
})
