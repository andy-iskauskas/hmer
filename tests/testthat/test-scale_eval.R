## Utility functions: scale_input, eval_funcs, function_to_names, multiply_function

## scale_input
test_that("Data.frame scaling works", {
  initial_df <- data.frame(a = c(-2, 0, 2), b = c(0, 2, 4))
  ranges <- list(a = c(-2, 2), b = c(0, 4))
  expect_equal(
    scale_input(initial_df, ranges),
    data.frame(a = c(-1, 0, 1), b = c(-1, 0, 1))
  )
  expect_equal(
    scale_input(
      scale_input(initial_df, ranges),
      ranges, forward = FALSE
    ),
    initial_df
  )
})

test_that("Scaling with eval_funcs", {
  initial_df <- data.frame(
    a = c(-3, -2, -1, -2.5),
    b = c(0.75, 1, 1.5, 1)
  )
  ranges <- list(a = c(-4, -1), b = c(0, 2))
  expect_equal(
    eval_funcs(scale_input,
               initial_df,
               ranges),
    data.frame(
      a = c(-1/3, 1/3, 1, 0),
      b = c(-1/4, 0, 1/2, 0)
    )
  )
  expect_equal(
    eval_funcs(
      scale_input,
      eval_funcs(
        scale_input, initial_df, ranges, forward = FALSE),
      ranges
      ),
    initial_df
  )
})

test_that("no names", {
  initial_df <- data.frame(a = c(-2, 0, 2), b = c(0, 2, 4))
  ranges <- list(a = c(-2, 2), b = c(0, 4))
  names(initial_df) <- NULL
  expect_equal(
    scale_input(initial_df, ranges),
    setNames(
      data.frame(a = c(-1, 0, 1), b = c(-1, 0, 1)),
      NULL)
  )
  expect_equal(
    scale_input(
      scale_input(initial_df, ranges),
      ranges, forward = FALSE
    ),
    initial_df
  )
})

# eval_funcs
test_that("eval_funcs with single function", {
  func <- function(x) x[[1]] + x[[2]]^2
  points <- data.frame(a = c(1, 2), b = c(0.5, 3))
  expect_equal(
    c(eval_funcs(func, points[1,]), use.names = FALSE),
    1.25
  )
  expect_equal(
    c(eval_funcs(func, points), use.names = FALSE),
    c(1.25, 11)
  )
  expect_equal(
    length(eval_funcs(func, rbind(points, points))),
    4
  )
})

test_that("eval_funcs with multiple functions", {
  funcs <- c(
    function(x) x[[1]],
    function(x) x[[2]]-x[[1]],
    function(x) 1)
  points <- data.frame(a = c(1, 2), b = c(0.5, 3))
  expect_equal(
    dim(eval_funcs(funcs, points)),
    c(3, 2)
  )
  expect_equal(
    eval_funcs(funcs, points),
    matrix(c(1, -0.5, 1, 2, 1, 1), nrow = 3)
  )
  points_extra <- data.frame(
    a = c(1, 2),
    b = c(0.5, 3),
    c = c(7, pi),
    d = c(-1, -1)
  )
  expect_equal(
    eval_funcs(funcs, points),
    eval_funcs(funcs, points_extra)
  )
  funcs_list <- list(function(x) x[[1]],
                     function(x) x[[2]]-x[[1]],
                     function(x) 1)
  expect_equal(
    eval_funcs(funcs, points),
    eval_funcs(funcs_list, points)
  )
})

# function_to_names
test_that("simple function - linear in one var", {
  func <- function(x) x[[1]]
  variable_names = "x"
  expect_match(
    function_to_names(func, variable_names),
    "x"
  )
  expect_match(
    function_to_names(func, variable_names),
    function_to_names(func, variable_names, FALSE)
  )
})

test_that("cubic in multiple variables", {
  func <- function(x) x[[1]] + x[[2]]^3 - x[[1]]^2*x[[2]]
  var_names <- c('x', 'y')
  expect_match(
    function_to_names(func, var_names),
    "x + y^3 - x^2:y",
    fixed = TRUE
  )
  expect_match(
    function_to_names(func, var_names, FALSE),
    "x + y^3 - x^2*y",
    fixed = TRUE
  )
})

# multiply_function

test_that("simple return function, variety of settings", {
  func <- function(x) {
    return(x)
  }
  expect_equal(
    multiply_function(func, 2)(6),
    12
  )
  func <- function(x) x
  expect_equal(
    multiply_function(func, -1)(-1),
    1
  )
  func <- function(x) {
    x <- x
    return(x)
  }
  expect_equal(
    multiply_function(func, 0.5)(1/4),
    1/8
  )
})

test_that("constant function", {
  func <- function(x) 7
  expect_equal(
    multiply_function(func, 2)(1),
    14
  )
  func <- function(x) return(6)
  expect_equal(
    multiply_function(func, 3)(1),
    18
  )
  func <- function(x) {
    -1
  }
  expect_equal(
    multiply_function(func, 0.1)(1),
    -0.1
  )
})

test_that("complex functions", {
  func <- function(x) {
    x <- x^2
    return(x/10)
  }
  expect_equal(
    multiply_function(func, 5)(-2),
    2
  )
  func <- function(x) {
    y <- sin(x)
    return(ifelse(y < 0, 1, 0))
  }
  expect_equal(
    multiply_function(func, 2)(3*pi/2),
    2
  )
})
