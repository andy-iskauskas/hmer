test_that("simple distances work", {
  expect_equal(
    get_dist(
      data.frame(a = c(1, 2)),
      data.frame(a = c(0.9, 2.1))
    ),
    matrix(c(0.1, 1.1, 1.1, 0.1), nrow = 2, byrow = TRUE)
  )
})

test_that("multi-dimensional, different sized works", {
  dist_res <- get_dist(
    data.frame(a = c(1, 2, 4), b = c(0.1, 0.4, -0.5)),
    data.frame(a = c(0, 1.5), b = c(0, 0.76))
  )
  expect_equal(
    dim(dist_res),
    c(2, 3)
  )
  expect_equal(
    dist_res[1,1],
    1.0049876,
    tolerance = 1e-7
  )
})

test_that("large point number using apply method", {
  df1 <- data.frame(a = runif(20),
                    b = runif(20, 1, 2),
                    c = runif(20, -1, 1),
                    d = runif(20, -5, -2),
                    e = runif(20, -2, 2)
                    )
  df2 <- data.frame(a = runif(1000),
                    b = runif(1000, 1, 2),
                    c = runif(1000, -1, 1),
                    d = runif(1000, -5, -2),
                    e = runif(1000, -2, 2)
                    )
  dists <- get_dist(df1, df2)
  expect_equal(
    dim(dists),
    c(1000, 20)
  )
})

test_that("large dimension number: dist method", {
  df1 <- setNames(data.frame(
    matrix(runif(1000), ncol = 20)
  ), paste0("X", 1:20))
  df2 <- setNames(data.frame(
    matrix(runif(800), ncol = 20)
  ), paste0("X", 1:20))
  dists <- get_dist(df1, df2)
  expect_equal(
    dim(dists),
    c(40, 50)
  )
})

test_that("fails if data.frames have different column dimensions", {
  df1 <- data.frame(a = runif(20),
                    b = runif(20, 1, 2),
                    c = runif(20, -1, 1),
                    d = runif(20, -5, -2),
                    e = runif(20, -2, 2)
  )
  df2 <- data.frame(a = runif(1000),
                    b = runif(1000, 1, 2),
                    c = runif(1000, -1, 1),
                    d = runif(1000, -5, -2)
  )
  expect_error(
    get_dist(df1, df2)
  )
})
