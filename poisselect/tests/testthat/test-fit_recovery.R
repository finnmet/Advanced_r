# The assignment's headline question: does the estimator recover known
# parameters on simulated data? Plus quadrature stability and start-value
# robustness. The large dataset and the default k = 20 fit are shared
# fixtures for this file.

set.seed(2024)
recovery_data <- simulate_poisselect(4000,
  beta = c(1, 0.5),
  gamma = c(0.5, 1, -0.8), sigma = 0.5,
  rho = 0.6
)
recovery_fit <- poisselect(y ~ x1, s ~ x1 + w1,
  data = recovery_data,
  k = 20
)

test_that("the estimator recovers the known parameters", {
  expect_identical(recovery_fit$convergence, 0L)
  beta_hat <- recovery_fit$coefficients$outcome
  gamma_hat <- recovery_fit$coefficients$selection
  expect_lt(max(abs(beta_hat - c(1, 0.5))), 0.15)
  expect_lt(max(abs(gamma_hat - c(0.5, 1, -0.8))), 0.15)
  expect_lt(abs(recovery_fit$coefficients$sigma - 0.5), 0.15)
  expect_lt(abs(recovery_fit$coefficients$rho - 0.6), 0.2)
  standard_errors <- unlist(recovery_fit$standard_errors)
  expect_true(all(is.finite(standard_errors)))
  expect_true(all(standard_errors > 0))
})

test_that("with ignorable selection the estimator matches the plain GLM", {
  set.seed(7)
  d <- simulate_poisselect(4000,
    beta = c(1, 0.5), gamma = c(0.5, 1, -0.8),
    sigma = 0.5, rho = 0
  )
  fit <- poisselect(y ~ x1, s ~ x1 + w1, data = d, k = 20)
  expect_lt(max(abs(fit$coefficients$outcome - c(1, 0.5))), 0.1)
  glm_fit <- glm(y ~ x1, family = poisson, data = d[d$s == 1, ])
  # The GLM slope is directly comparable, but its intercept estimates
  # beta_0 + sigma^2 / 2, not beta_0: with log-normal heterogeneity eps
  # the observed mean is exp(x'beta + sigma^2 / 2), and the plain GLM
  # absorbs the constant sigma^2 / 2 into its intercept. So compare the
  # slope against poisselect and the GLM intercept against the shifted
  # target - never intercept against intercept.
  expect_lt(
    abs(coef(glm_fit)[["x1"]] - fit$coefficients$outcome[["x1"]]),
    0.1
  )
  expect_lt(abs(coef(glm_fit)[["(Intercept)"]] - (1 + 0.5^2 / 2)), 0.1)
})

test_that("the fit is stable in the number of quadrature nodes", {
  fit_low <- poisselect(y ~ x1, s ~ x1 + w1, data = recovery_data, k = 10)
  fit_high <- poisselect(y ~ x1, s ~ x1 + w1, data = recovery_data, k = 30)
  difference <- unlist(fit_low$coefficients) -
    unlist(fit_high$coefficients)
  expect_lt(max(abs(difference)), 0.05)
})

test_that("a user start at the truth reaches the same optimum", {
  fit_true_start <- poisselect(y ~ x1, s ~ x1 + w1,
    data = recovery_data,
    k = 20,
    start = c(1, 0.5, 0.5, 1, -0.8, 0.5, 0.6)
  )
  expect_equal(unlist(fit_true_start$coefficients),
    unlist(recovery_fit$coefficients),
    tolerance = 1e-2
  )
  expect_equal(fit_true_start$loglik, recovery_fit$loglik,
    tolerance = 1e-6
  )
})
