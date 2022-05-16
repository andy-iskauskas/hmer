# Testing subsetting and collating emulators

emulators <- SIREmulators$ems

## Subsetting

test_that("subsetting variance ems", {
  fake_var <- list(
    expectation = emulators,
    variance = emulators
  )
  subset_names <- c('nS', 'nR')
  subsetted <- subset_emulators(fake_var, subset_names)
  expect_equal(
    names(subsetted),
    c("expectation", "variance")
  )
  expect_equal(
    names(
      subsetted$expectation
      ),
    subset_names
  )
})

test_that("subsetting bimodal ems", {
  fake_bim <- list(
    mode1 = list(
      expectation = emulators,
      variance = emulators),
    mode2 = list(
      expectation = emulators,
      variance = emulators),
    prop = emulators[[1]]
  )
  subset_names <- c('nI', 'nS')
  subsetted <- subset_emulators(fake_bim, subset_names)
  expect_equal(
    names(subsetted),
    c('mode1', 'mode2', 'prop')
  )
  expect_equal(
    names(subsetted$mode1),
    c('variance', 'expectation')
  )
  expect_equal(
    names(subsetted$mode1$expectation),
    rev(subset_names)
  )
  expect_true(
    "Emulator" %in% class(subsetted$mode1$expectation[[1]])
  )
})

test_that("subsetting a normal collection", {
  expect_equal(
    names(
      subset_emulators(emulators, c('nI', 'nR'))
    ),
    c('nI', 'nR')
  )
})

## Getting ranges

ems2 <- setNames(purrr::map(emulators, ~.$clone()), names(emulators))
for (i in seq_along(ems2)) ems2[[i]]$ranges <- list(
  aSI = c(0, 0.3),
  aIR = c(0, 0.2),
  aSR = c(0, 0.5)
)

test_that("Collecting emulators: same ranges", {
  combine_ems <- list(emulators, emulators)
  expect_equal(
    getRanges(combine_ems)$aSI,
    c(0.1, 0.8)
  )
})

test_that("Collecting emulators: different ranges", {
  combine_ems <- list(emulators, ems2)
  expect_equal(
    getRanges(combine_ems),
    list(aSI = c(0, 0.3),
         aIR = c(0, 0.2),
         aSR = c(0, 0.05))
  )
  expect_equal(
    getRanges(combine_ems, minimal = FALSE),
    list(aSI = c(0.1, 0.8),
         aIR = c(0, 0.5),
         aSR = c(0, 0.5))
  )
})

test_that("Variance and bimodal emulators", {
  v_ems <- list(
    expectation = emulators,
    variance = emulators
  )
  expect_equal(
    getRanges(v_ems),
    emulators[[1]]$ranges
  )
  b_ems <- list(
    mode1 = v_ems,
    mode_2 = v_ems,
    prop = emulators$nI
  )
  expect_equal(
    getRanges(b_ems),
    emulators[[1]]$ranges
  )
})

## Collecting emulators

test_that("Single emulator", {
  expect_true(
    is.list(collect_emulators(emulators$nS))
  )
  expect_match(
    names(collect_emulators(emulators$nS))[1],
    "nS"
  )
})

test_that("List of simple emulators", {
  expect_equal(
    names(collect_emulators(emulators)),
    names(emulators)
  )
  expect_true(
    all(
      class(collect_emulators(emulators)) ==
        class(emulators)
      )
  )
})

test_that("list of lists", {
  l_em <- list(emulators, emulators, emulators$nS)
  expect_equal(
    length(collect_emulators(l_em)),
    7
  )
  expect_equal(
    names(collect_emulators(l_em)),
    c(rep(names(emulators), 2), 'nS')
  )
})

test_that("variance ems", {
  v_ems <- list(
    expectation = emulators,
    variance = emulators
  )
  expect_equal(
    names(collect_emulators(v_ems)),
    c("expectation", "variance")
  )
  vv_ems <- list(v_ems, v_ems)
  expect_equal(
    length(collect_emulators(vv_ems)),
    2
  )
  expect_equal(
    length(collect_emulators(vv_ems)$expectation),
    6
  )
})

test_that("bimodal ems", {
  bim_ems <- list(
    mode1 = list(
      expectation = emulators,
      variance = emulators
    ),
    mode2 = list(
      expectation = emulators,
      variance = emulators
    ),
    prop = emulators$nS
  )
  expect_equal(
    names(collect_emulators(bim_ems)),
    c("mode1", "mode2", "prop")
  )
  bbim_ems <- list(bim_ems, bim_ems)
  bbim_collect <- collect_emulators(bbim_ems)
  expect_equal(
    names(bbim_collect),
    c("mode1", "mode2", "prop")
  )
  expect_equal(
    names(bbim_collect$mode1),
    c("expectation", "variance")
  )
  expect_equal(
    names(bbim_collect$mode1$expectation),
    rep(names(emulators), 2)
  )
  bbbim_ems <- c(bim_ems, bim_ems)
  bbbim_collect <- collect_emulators(bbbim_ems)
  expect_equal(
    names(bbbim_collect),
    c("mode1", "mode2", "prop")
  )
  expect_equal(
    names(bbbim_collect$mode1),
    c("expectation", "variance")
  )
  expect_equal(
    names(bbbim_collect$mode1$expectation),
    rep(names(emulators), 2)
  )
})
