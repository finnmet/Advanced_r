# One expect_error() per faulty input from the specification; the regexp
# fragments pin the informative part of each message.

test_that("argument-level checks catch malformed arguments", {
  d <- test_data
  expect_error(
    poisselect(y ~ x1, s ~ x1 + w1, data = "nope"),
    "data.frame"
  )
  expect_error(
    poisselect(y ~ x1, s ~ x1 + w1, data = data.frame()),
    "at least 2"
  )
  expect_error(poisselect("y ~ x1", s ~ x1 + w1, data = d), "formula")
  expect_error(poisselect(~x1, s ~ x1 + w1, data = d), "two-sided")
  expect_error(poisselect(y ~ x1, "s ~ x1 + w1", data = d), "formula")
  expect_error(
    poisselect(y ~ x1 + nonexistent, s ~ x1 + w1, data = d),
    "nonexistent"
  )
})

test_that("the quadrature node argument k is validated", {
  d <- test_data
  expect_error(poisselect(y ~ x1, s ~ x1 + w1, data = d, k = 1), "k")
  expect_error(poisselect(y ~ x1, s ~ x1 + w1, data = d, k = 0), "k")
  expect_error(poisselect(y ~ x1, s ~ x1 + w1, data = d, k = 2.5), "k")
  expect_error(poisselect(y ~ x1, s ~ x1 + w1, data = d, k = NA), "k")
  expect_error(
    poisselect(y ~ x1, s ~ x1 + w1, data = d, k = c(10, 20)),
    "k"
  )
})

test_that("user-supplied starting values are validated", {
  d <- test_data
  expect_error(
    poisselect(y ~ x1, s ~ x1 + w1, data = d, start = c(1, 2, 3)),
    "length"
  )
  expect_error(
    poisselect(y ~ x1, s ~ x1 + w1,
      data = d,
      start = c(1, 0.5, 0.5, 1, -0.8, -1, 0)
    ),
    "positive"
  )
  expect_error(
    poisselect(y ~ x1, s ~ x1 + w1,
      data = d,
      start = c(1, 0.5, 0.5, 1, -0.8, 0.5, 1.2)
    ),
    "between -1 and 1"
  )
  expect_error(
    poisselect(y ~ x1, s ~ x1 + w1,
      data = d,
      start = c(1, 0.5, 0.5, 1, -0.8, 0.5, NA)
    ),
    "start"
  )
})

test_that("the optim control list is validated", {
  d <- test_data
  expect_error(
    poisselect(y ~ x1, s ~ x1 + w1, data = d, control = "x"),
    "list"
  )
  expect_error(
    poisselect(y ~ x1, s ~ x1 + w1, data = d, control = list(bogus = 1)),
    "subset"
  )
})

test_that("the selection indicator is validated", {
  d <- test_data
  d_na <- d
  d_na$s[1] <- NA
  expect_error(
    poisselect(y ~ x1, s ~ x1 + w1, data = d_na),
    "selection indicator"
  )
  d_two <- d
  d_two$s[1] <- 2
  expect_error(
    poisselect(y ~ x1, s ~ x1 + w1, data = d_two),
    "only 0 and 1"
  )
  d_frac <- d
  d_frac$s <- d_frac$s + 0.5
  expect_error(
    poisselect(y ~ x1, s ~ x1 + w1, data = d_frac),
    "only 0 and 1"
  )
  d_none <- d
  d_none$s <- 0
  d_none$y <- NA
  expect_error(
    poisselect(y ~ x1, s ~ x1 + w1, data = d_none),
    "No selected observations"
  )
  d_all <- d
  d_all$s <- 1
  d_all$y[is.na(d_all$y)] <- 0
  expect_error(
    poisselect(y ~ x1, s ~ x1 + w1, data = d_all),
    "All observations are selected"
  )
})

test_that("a logical selection indicator is accepted", {
  d <- test_data
  d$s <- d$s == 1
  fit <- poisselect(y ~ x1, s ~ x1 + w1, data = d, k = 5)
  expect_s3_class(fit, "poisselect")
})

test_that("the outcome counts are validated on the selected rows", {
  d <- test_data
  selected_row <- which(d$s == 1)[1]
  d_na <- d
  d_na$y[selected_row] <- NA
  expect_error(
    poisselect(y ~ x1, s ~ x1 + w1, data = d_na),
    "NA on selected rows"
  )
  d_neg <- d
  d_neg$y[selected_row] <- -1
  expect_error(
    poisselect(y ~ x1, s ~ x1 + w1, data = d_neg),
    "non-negative integer"
  )
  d_frac <- d
  d_frac$y[selected_row] <- 2.5
  expect_error(
    poisselect(y ~ x1, s ~ x1 + w1, data = d_frac),
    "non-negative integer"
  )
  d_inf <- d
  d_inf$y[selected_row] <- Inf
  expect_error(poisselect(y ~ x1, s ~ x1 + w1, data = d_inf), "finite")
})

test_that("missing values in the covariates are caught and located", {
  d <- test_data
  selected_row <- which(d$s == 1)[1]
  d_z <- d
  d_z$w1[3] <- NA
  expect_error(
    poisselect(y ~ x1, s ~ x1 + w1, data = d_z),
    "Selection covariates.*w1"
  )
  # x1 appears in both equations there, so use a model where x1 is
  # outcome-only to hit the outcome-covariate check specifically.
  d_x <- d
  d_x$x1[selected_row] <- NA
  expect_error(
    poisselect(y ~ x1, s ~ w1, data = d_x),
    "Outcome covariates.*x1"
  )
})

test_that("NA in an outcome covariate on an unselected row is allowed", {
  d <- test_data
  unselected_row <- which(d$s == 0)[1]
  d$x1[unselected_row] <- NA
  fit <- poisselect(y ~ x1, s ~ w1, data = d, k = 5)
  expect_s3_class(fit, "poisselect")
})

test_that("collinear model matrices are rejected", {
  d <- test_data
  d$x1_copy <- d$x1
  expect_error(
    poisselect(y ~ x1 + x1_copy, s ~ x1 + w1, data = d),
    "collinear"
  )
  d$w1_copy <- d$w1
  expect_error(
    poisselect(y ~ x1, s ~ x1 + w1 + w1_copy, data = d),
    "collinear"
  )
})

test_that("too few selected observations are rejected", {
  d <- test_data
  d$s <- 0L
  keep <- which(!is.na(d$y))[1:3]
  d$s[keep] <- 1L
  expect_error(
    poisselect(y ~ x1, s ~ x1 + w1, data = d),
    "Too few selected"
  )
})

test_that("a missing exclusion restriction triggers a warning", {
  d <- test_data
  expect_warning(
    poisselect(y ~ x1, s ~ x1, data = d, k = 5),
    "exclusion restriction"
  )
})

test_that("non-convergence triggers a warning", {
  d <- test_data
  expect_warning(
    poisselect(y ~ x1, s ~ x1 + w1,
      data = d, k = 5,
      control = list(maxit = 1)
    ),
    "converge"
  )
})
