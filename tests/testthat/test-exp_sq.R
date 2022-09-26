### Exponential-squared correlation function

test_that("correlation with self is 1", {
  expect_equal(
    c(exp_sq(
      data.frame(a = 1),
      data.frame(a = 1),
      list(theta = 0.1)
    ), use.names = FALSE),
    1
  )
  expect_equal(
    c(diag(exp_sq(
      data.frame(a = c(1, 2, 3), b = c(0.1, 0.4, 0.3)),
      data.frame(a = c(1, 2, 3), b = c(0.1, 0.4, 0.3)),
      list(theta = 0.2)
    )), use.names = FALSE),
    rep(1, 3)
  )
})

test_that("one-dimensional exp-squared; single point", {
  expect_equal(
    c(exp_sq(
      data.frame(a = 1),
      data.frame(a = 2),
      list(theta = 0.1)
    )),
    3.720076e-44)
})

test_that("one-dimensional exp-squared; multi point", {
  expect_equal(
    exp_sq(
      data.frame(a = c(1, 2)),
      data.frame(a = c(1.1, 2.9)),
      list(theta = 0.4)
    ),
    matrix(c(9.394131e-01, 0.006329715, 1.589391e-10, 0.006329715),
           nrow = 2, byrow = TRUE),
    tolerance = 1e-7
    )
})

test_that("multi-dimensional exp-squared; single point", {
  expect_equal(
    c(exp_sq(
      data.frame(a = 1, b = 2, c = -1),
      data.frame(a = 1.5, b = 2.9, c = -0.7),
      list(theta = 0.2)
    )),
    3.266131e-13
  )
})

test_that("multi-dimensional exp=squared; multi point", {
  expect_equal(
    exp_sq(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.8, 2.4, 3.2), b = c(0.5, 0, -0.5)),
      list(theta = 1)
    ),
    matrix(
      c(0.8436648, 0.6376282, 0.07653555,
        0.7710516, 0.9048374, 0.31348618,
        0.1287349, 0.2541070, 0.42741493),
      nrow = 3, byrow = TRUE
    ),
    tolerance = 1e-7
  )
})

test_that("dimensionality checks", {
  expect_equal(
    dim(exp_sq(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.8, 2.4), b = c(0.5, 0)),
      list(theta = 1)
    )),
    c(2,3)
  )
  expect_equal(
    dim(exp_sq(
      data.frame(a = c(1.8), b = c(0.5)),
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      list(theta = 1)
    )),
    c(3,1)
  )
  expect_equal(
    dim(exp_sq(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.8), b = c(0.5)),
      list(theta = 1)
    )),
    c(1,3)
  )
})

test_that("different theta per dimension", {
  expect_equal(
    exp_sq(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.8, 2.4, 3.2), b = c(0.5, 0, -0.5)),
      list(theta = c(1, 0.5))
    ),
    matrix(
      c(0.5220457768, 0.0009396529, 0.07427358,
        0.0002139004, 0.6907343306, 0.07427358,
        0.0437177973, 0.1636541368, 0.03762826),
      nrow = 3, byrow = TRUE
    ),
    tolerance = 1e-6
  )
})

test_that("same points gives symmetric matrix", {
  corr_out <- exp_sq(
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
    exp_sq(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4))
    )
  )
})

test_that("fails with missing data.frame", {
  expect_error(
    exp_sq(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      hp = list(theta = 0.1)
    )
  )
})

test_that("works with data.matrix or data.frame", {
  df <- data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4))
  dm <- data.matrix(df)
  expect_equal(
    exp_sq(df, df, list(theta = 0.2)),
    exp_sq(dm, dm, list(theta = 0.2))
  )
})

## Derivative of exponential-squared correlation function

test_that("correlation with self is 0", {
  expect_equal(
    unname(exp_sq_d(
      data.matrix(data.frame(a = 1)),
      data.matrix(data.frame(a = 1)),
      list(theta = 0.1),
      1
    )),
    matrix(0, nrow = 1)
  )
  expect_equal(
    c(diag(exp_sq_d(
      data.matrix(data.frame(a = c(1, 2, 3), b = c(0.1, 0.4, 0.3))),
      data.matrix(data.frame(a = c(1, 2, 3), b = c(0.1, 0.4, 0.3))),
      list(theta = 0.2),
      2
    )), use.names = FALSE),
    rep(0, 3)
  )
})

test_that("one-dimensional exp-squared derivative; single point", {
  expect_equal(
    exp_sq_d(
      data.matrix(data.frame(a = 1)),
      data.matrix(data.frame(a = 2)),
      list(theta = 0.1),
      1
    ),
    matrix(7.440152e-42, nrow = 1))
})

test_that("one-dimensional exp-squared derivative; multi point", {
  expect_equal(
    exp_sq_d(
      data.matrix(data.frame(a = c(1, 2))),
      data.matrix(data.frame(a = c(1.1, 2.9))),
      list(theta = 0.4),
      1
    ),
    matrix(c(1.174266, -0.0712093, 3.774804e-09, 0.0712093),
           nrow = 2, byrow = TRUE),
    tolerance = 1e-6
  )
})

test_that("multi-dimensional exp-squared derivative; single point", {
  expect_equal(
    exp_sq_d(
      data.matrix(data.frame(a = 1, b = 2, c = -1)),
      data.matrix(data.frame(a = 1.5, b = 2.9, c = -0.7)),
      list(theta = 0.2),
      3
    ),
    matrix(4.899197e-12, nrow = 1)
  )
})

test_that("multi-dimensional exp-squared derivative; single point and multi-deriv", {
  expect_equal(
    exp_sq_d(
      data.matrix(data.frame(a = 1, b = 2, c = -1)),
      data.matrix(data.frame(a = 1.5, b = 2.9, c = -0.7)),
      list(theta = 0.2),
      3,
      2
    ),
    matrix(-2.204639e-10, nrow = 1)
  )
})

test_that("multi-dimensional exp-squared derivative; multi point and multi deriv", {
  expect_equal(
    exp_sq_d(
      data.matrix(data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4))),
      data.matrix(data.frame(a = c(1.8, 2.4, 3.2), b = c(0.5, 0, -0.5))),
      list(theta = 1),
      1,
      2
    ),
    matrix(
      c(0.1349864, 0.4590923, 0.04898275,
        0.1542103, -0.1085805, -0.50157789,
        0.4016529, 0.4472282, -0.30773875),
      nrow = 3, byrow = TRUE
    ),
    tolerance = 1e-7
  )
})

test_that("fails with no derivative direction", {
  expect_error(
    exp_sq_d(
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4)),
      data.frame(a = c(1.9, 2.1, 3.4), b = c(0.1, -0.1, 0.4))
    )
  )
})
