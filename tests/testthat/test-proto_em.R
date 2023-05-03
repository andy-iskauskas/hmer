ranges <- list(aSI = c(0.1, 0.8),
               aIR = c(0, 0.5),
               aSR = c(0, 0.05))
targets <- SIREmulators$targets
lms <- purrr::map(names(targets),
                  ~step(lm(data = SIRSample$training,
                           formula = as.formula(paste0(., "~(", paste0(names(ranges), collapse = "+"), ")^2"))),
                        trace = 0))
relev_lm <- lms[[1]]
pf <- function(x) predict(relev_lm, x)
vf <- function(x) {
  pred <- predict(relev_lm, x, se.fit = TRUE)
  return(pred$se.fit^2 + pred$residual.scale^2)
}
ifunc <- function(x, z, cutoff) return(x-z)

test_that("Basic proto_emulator behaviour works", {
  pe1 <- Proto_emulator$new(
    ranges, targets, pf, vf
  )
  expect_true("EmProto" %in% class(pe1))
  pe2 <- Proto_emulator$new(
    ranges, targets, pf, vf, ifunc
  )
  expect_true("EmProto" %in% class(pe2))
})

test_that("Proto em failure states",  {
  expect_error(
    Proto_emulator$new(
      ranges, targets, 7, vf
    ),
    "does not appear to be a function"
  )
  expect_error(
    Proto_emulator$new(
      ranges, targets, function() pf(x), vf
    ),
    "requires at least one argument"
  )
  expect_error(
    Proto_emulator$new(
      ranges, targets, pf, 7
    ),
    "does not appear to be a function"
  )
  expect_error(
    Proto_emulator$new(
      ranges, targets, pf, function() vf(x)
    ),
    "requires at least one argument"
  )
  expect_error(
    Proto_emulator$new(
      ranges, targets, pf, vf, 7
    ),
    "does not appear to be a function"
  )
  expect_error(
    Proto_emulator$new(
      ranges, targets, pf, vf, function(x) ifunc(x, z, cutoff)
    ),
    "at least two arguments"
  )
})

vpe <- purrr::map(seq_along(lms), function(l) {
  Proto_emulator$new(
    ranges,
    names(targets)[l],
    function(x) predict(lms[[l]], x),
    function(x) {
      pred <- predict(lms[[l]], x, se.fit = TRUE)
      return(pred$se.fit^2 + pred$residual.scale^2)
    },
    print_func = function() print(summary(lms[[l]]))
  )
}) |> setNames(names(targets))

test_that("Unavailable functions fail gracefully", {
  expect_error(
    effect_strength(vpe),
    "not applicable"
  )
  expect_error(
    individual_errors(vpe[[2]], SIRSample$validation),
    "not applicable"
  )
  expect_error(
    residual_diag(vpe[[3]]),
    "not applicable"
  )
  expect_error(
    space_removal(vpe, targets),
    "Cannot access"
  )
  expect_error(
    summary_diag(vpe[[1]], SIRSample$validation),
    "not applicable"
  )
  expect_error(
    validation_diagnostics(vpe, targets),
    "require a validation set"
  )
  expect_error(
    get_diagnostic(vpe[[2]], targets),
    "requires validation set"
  )
})

test_that("Proto behaviour - validation", {
  expect_equal(
    nrow(classification_diag(vpe[[1]], targets, SIRSample$validation, plt = FALSE)),
    0
  )
  expect_equal(
    nrow(comparison_diag(vpe[[2]], targets, SIRSample$validation, plt = FALSE)),
    0
  )
  expect_equal(
    nrow(standard_errors(vpe[[3]], targets, SIRSample$validation, plt = FALSE)),
    0
  )
  expect_equal(
    nrow(validation_diagnostics(vpe, targets, SIRSample$validation, plt = FALSE)),
    0
  )
})

test_that("Proto behaviour - point generation", {
  g1 <- generate_new_design(vpe, 100, targets, verbose = FALSE)
  expect_equal(nrow(g1), 100)
  g2 <- generate_new_design(vpe, 100, targets,
                          plausible_set = g1[1:20,],
                          method = c('line', 'importance'), verbose = FALSE)
  expect_equal(nrow(g2), 100)
  cimp <- function(ems, x, z, n = n, cutoff = cutoff, ...) {
    imps <- nth_implausible(ems, x, z, n = n, cutoff = cutoff, ...)
    constraint <- apply(x, 1, function(y) y[[1]] < y[[2]])
    return(imps & constraint)
  }
  g3 <- generate_new_design(vpe, 100, targets, accept_measure = cimp, verbose = FALSE)
  expect_equal(nrow(g3), 100)
})

test_that("Proto behaviour - other output functions", {
  expect_equal(
    as.numeric(space_removal(vpe, targets, SIRSample$validation)),
    c(27, 14, 4)/60
  )
  expect_true(
    all(is.logical(nth_implausible(vpe, SIRSample$validation, targets, n=2, cutoff = 1.5)))
  )
  expect_equal(
    sum(nth_implausible(vpe, SIRSample$validation, targets, cutoff = 1.5)),
    11
  )
  expect_equal(
    max(nth_implausible(vpe, SIRSample$validation, targets)),
    8.854471,
    tolerance = 1e-4
  )
})

test_that("Proto derivative fails gracefully", {
  expect_error(
    directional_deriv(vpe[[3]], SIRSample$validation[1,], c(1,1,1)),
    "No derivative expression"
  )
  expect_error(
    directional_proposal(vpe, SIRSample$validation[1,], targets),
    "No derivative expression"
  )
})
