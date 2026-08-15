# Reproducibility: identical seeds give identical data, and refitting the
# same data gives the same estimates.

test_that("simulation is reproducible under a fixed seed", {
  set.seed(99)
  d_first <- simulate_poisselect(200,
    beta = c(1, 0.5),
    gamma = c(0.5, 1, -0.8), sigma = 0.5,
    rho = 0.6
  )
  set.seed(99)
  d_second <- simulate_poisselect(200,
    beta = c(1, 0.5),
    gamma = c(0.5, 1, -0.8), sigma = 0.5,
    rho = 0.6
  )
  expect_identical(d_first, d_second)
})

test_that("fitting the same data twice gives equal results", {
  d <- make_test_data(n = 400, seed = 5)
  fit_first <- poisselect(y ~ x1, s ~ x1 + w1, data = d, k = 8)
  fit_second <- poisselect(y ~ x1, s ~ x1 + w1, data = d, k = 8)
  expect_equal(fit_first$coefficients, fit_second$coefficients)
  expect_equal(fit_first$loglik, fit_second$loglik)
})
