## Correlator Object

## Initialisation

test_that("Correlator default initialises", {
  test_corr <- Correlator$new()
  expect_match(
    test_corr$corr_name,
    "exp_sq"
  )
  expect_match(
    names(test_corr$hyper_p),
    c("theta")
  )
  expect_equal(
    test_corr$hyper_p$theta,
    0.1
  )
  expect_equal(
    test_corr$nugget,
    0
  )
})

test_that("Custom Correlator initialises", {
  test_corr <- Correlator$new("matern", list(theta = 0.8, nu = 1.5), nug = 0.15)
  expect_match(
    test_corr$corr_name,
    "matern"
  )
  expect_equal(
    names(test_corr$hyper_p),
    c("theta", "nu")
  )
  expect_equal(
    test_corr$hyper_p$theta,
    0.8
  )
  expect_equal(
    test_corr$hyper_p$nu,
    1.5
  )
  expect_equal(
    test_corr$nugget,
    0.15
  )
})

test_that("Print statement behaves,", {
  test_corr <- Correlator$new("matern", list(theta = 0.8, nu = 1.5), nug = 0.15)
  expect_output(
    print(test_corr),
    "Correlation type"
  )
  expect_output(
    print(test_corr),
    "Hyperparameters"
  )
  expect_output(
    print(test_corr),
    "Nugget term"
  )
})

test_that("Querying hyper-parameters", {
  expect_equal(
    Correlator$new()$get_hyper_p(),
    list(theta = 0.1)
  )
})

test_that("Changing hyper-parameters", {
  c1 <- Correlator$new()
  c2 <- c1$set_hyper_p(list(theta = 0.6))
  c2alt <- c1$set_hyper_p(list(0.6))
  c3 <- c1$set_hyper_p(new_hp = c1$hyper_p, nug = 0.64)
  expect_equal(
    c2$hyper_p$theta,
    0.6
  )
  expect_equal(
    c2alt$hyper_p$theta,
    0.6
  )
  expect_equal(
    c3$nugget,
    0.64
  )
})

## Calculating correlations

test_corr <- Correlator$new("matern", list(theta = 0.8, nu = 1.5), nug = 0.15)

test_that("single point correlation", {
  expect_equal(
    test_corr$get_corr(data.frame(a = runif(1), b = runif(1, -1, 1))),
    1
  )
})

test_that("repeated data.frame same as single data.frame argument", {
  pts <- data.frame(a = runif(10), b = runif(10))
  expect_equal(
    test_corr$get_corr(pts),
    test_corr$get_corr(pts, pts)
  )
})

test_that("correlation function is symmetric", {
  pts <- data.frame(a = runif(10), b = runif(10))
  expect_equal(
    c(test_corr$get_corr(pts[1,], pts[2,]), use.names = FALSE),
    c(test_corr$get_corr(pts[2,], pts[1,]), use.names = FALSE)
  )
})

test_that("Correlation with and without nugget", {
  pts <- data.frame(a = runif(2), b = runif(2, -1, 1))
  corrs1 <- test_corr$get_corr(pts)
  corrs2 <- test_corr$get_corr(pts, use.nugget = FALSE)
  expect_equal(
    corrs1[1,1],
    corrs2[1,1]
  )
  expect_false(
    isTRUE(all.equal(corrs1, corrs2))
    )
})

test_that("Active variable changes result", {
  pts <- data.frame(a = runif(2), b = runif(2, -1, 1))
  corrs1 <- test_corr$get_corr(pts)
  corrs2 <- test_corr$get_corr(pts, actives = c(TRUE, FALSE))
  expect_equal(
    corrs1[1,1],
    corrs2[1,1]
  )
  expect_true(
    isTRUE(all(corrs2 >= corrs1))
    )
})

## Calculating derivatives

test_corr <- Correlator$new()

test_that("1d derivative works", {
  pts <- data.matrix(data.frame(a = runif(2), b = runif(2, -1, 1)))
  expect_equal(
    test_corr$get_corr_d(pts, p1 = 2)[1,1],
    0
  )
})

test_that("2d derivative works", {
  pts <- data.matrix(data.frame(a = runif(2), b = runif(2, -1, 1)))
  pts2 <- data.matrix(data.frame(a = runif(4), b = runif(4, -1, 1)))
  expect_equal(
    dim(test_corr$get_corr_d(pts, pts2, p1 = 1, p2 = 1)),
    c(4, 2)
  )
})

test_that("One data.frame = repeated data.frame", {
  pts <- data.matrix(data.frame(a = runif(3), b = runif(3, -1, 1)))
  expect_equal(
    test_corr$get_corr_d(pts, p1 = 1),
    test_corr$get_corr_d(pts, pts, p1 = 1)
  )
})

test_that("Swapped data.frames = minus sign transposition", {
  pts1 <- data.matrix(data.frame(a = runif(3), b = runif(3, -1, 1)))
  pts2 <- data.matrix(data.frame(a = runif(2), b = runif(2, -1, 1)))
  expect_equal(
    test_corr$get_corr_d(pts1, pts2, p1 = 1),
    -t(test_corr$get_corr_d(pts2, pts1, p1 = 1))
  )
})

test_that("Not active = derivative 0", {
  pts <- data.frame(a = runif(2), b = runif(2, -1, 1))
  expect_equal(
    test_corr$get_corr_d(pts, p1 = 1, actives = c(FALSE, TRUE)),
    matrix(rep(0, 4), nrow = 2)
  )
  expect_equal(
    test_corr$get_corr_d(pts, p1 = 1, p2 = 2, actives = c(TRUE, FALSE)),
    matrix(rep(0, 4), nrow = 2)
  )
})

test_that("Dimensions work", {
  pts1 <- data.frame(a = runif(2), b = runif(2, -1, 1))
  pts2 <- data.frame(a = runif(3), b = runif(3, -1, 1))
  expect_equal(
    dim(test_corr$get_corr_d(pts1, pts2, 1, 2, actives = c(FALSE, TRUE))),
    c(3, 2)
  )
})

test_that("Error if no derivative function", {
  t_corr <- Correlator$new("orn_uhl", list(theta = 0.3), nug = 0.1)
  expect_error(
    t_Corr$get_corr_d(
      data.matrix(data.frame(a = runif(2), b = runif(2, -1, 1))),
      p1 = 1)
  )
})
