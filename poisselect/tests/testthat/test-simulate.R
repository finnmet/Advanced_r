# Tests for the data generator: output structure, NA pattern, and input
# validation.

test_that("simulate_poisselect returns the documented data frame", {
  set.seed(3)
  d <- simulate_poisselect(100,
    beta = c(1, 0.5), gamma = c(0.5, 1, -0.8),
    sigma = 0.5, rho = 0.6
  )
  expect_s3_class(d, "data.frame")
  expect_identical(names(d), c("y", "s", "x1", "w1"))
  expect_identical(nrow(d), 100L)
  expect_true(all(d$s %in% c(0L, 1L)))
  # y must be NA exactly where s == 0.
  expect_identical(is.na(d$y), d$s == 0L)
  y_observed <- d$y[d$s == 1L]
  expect_true(all(y_observed >= 0))
  expect_true(all(y_observed == round(y_observed)))
})

test_that("simulate_poisselect rejects invalid arguments", {
  valid_beta <- c(1, 0.5)
  valid_gamma <- c(0.5, 1, -0.8)
  expect_error(simulate_poisselect(10, valid_beta, valid_gamma, 0.5, 0.5))
  expect_error(
    simulate_poisselect(
      100, c(1, 0.5, 0.3, 0.2), valid_gamma,
      0.5, 0.5
    ),
    "gamma"
  )
  expect_error(simulate_poisselect(100, valid_beta, valid_gamma, 0, 0.5))
  expect_error(simulate_poisselect(100, valid_beta, valid_gamma, -1, 0.5))
  expect_error(
    simulate_poisselect(100, valid_beta, valid_gamma, 0.5, 1),
    "between -1 and 1"
  )
  expect_error(
    simulate_poisselect(100, valid_beta, valid_gamma, 0.5, -1),
    "between -1 and 1"
  )
  expect_error(
    simulate_poisselect(100, valid_beta, valid_gamma, 0.5, 1.5),
    "between -1 and 1"
  )
  expect_error(
    simulate_poisselect(100, "a", valid_gamma, 0.5, 0.5),
    "beta"
  )
})
