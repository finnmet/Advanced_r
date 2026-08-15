# Tests for the approximated log-likelihood: finiteness, log-sum-exp
# robustness for large counts, and the analytic factorisation limit.

test_that("the log-likelihood is finite at the starting values", {
  model <- build_model_data(y ~ x1, s ~ x1 + w1, test_data)
  start <- compute_starting_values(model$y, model$x, model$z, model$s)
  gh <- compute_gh_rule(20)
  loglik <- compute_poisselect_loglik(
    start$beta, start$gamma, start$sigma,
    start$rho, model$y, model$x, model$z,
    model$s, gh
  )
  expect_true(is.finite(loglik))
})

test_that("log-sum-exp keeps the likelihood finite for large counts", {
  # A large intercept produces counts of several hundred; forming the
  # inner sum on the probability scale would underflow to 0 (log -> -Inf).
  set.seed(11)
  d <- simulate_poisselect(800,
    beta = c(4.5, 0.4),
    gamma = c(0.5, 1, -0.8), sigma = 0.3, rho = 0.5
  )
  expect_gt(max(d$y, na.rm = TRUE), 150)
  model <- build_model_data(y ~ x1, s ~ x1 + w1, d)
  gh <- compute_gh_rule(10)
  loglik <- compute_poisselect_loglik(
    c(4.5, 0.4), c(0.5, 1, -0.8), 0.3,
    0.5, model$y, model$x, model$z,
    model$s, gh
  )
  expect_true(is.finite(loglik))
  # The full estimation must also run through without any NaN/Inf issues.
  fit <- expect_no_warning(poisselect(y ~ x1, s ~ x1 + w1,
    data = d,
    k = 10
  ))
  expect_s3_class(fit, "poisselect")
  expect_true(is.finite(fit$loglik))
})

test_that("the quadrature collapses to the analytic limit", {
  # With rho = 0 the Phi term no longer depends on the node and factors
  # out of the quadrature sum (sum_k w_k / sqrt(pi) = 1), and with
  # sigma -> 0 the Poisson mixture collapses to a plain Poisson at
  # exp(x'beta): the likelihood factorises into Poisson x probit.
  model <- build_model_data(y ~ x1, s ~ x1 + w1, test_data)
  beta <- c(1, 0.5)
  gamma <- c(0.5, 1, -0.8)
  gh <- compute_gh_rule(20)
  approx_loglik <- compute_poisselect_loglik(
    beta, gamma, 1e-8, 0,
    model$y, model$x, model$z,
    model$s, gh
  )
  selected <- model$s == 1L
  lin_out <- drop(model$x %*% beta)
  lin_sel <- drop(model$z %*% gamma)
  exact_loglik <-
    sum(dpois(model$y[selected], exp(lin_out[selected]), log = TRUE)) +
    sum(pnorm(lin_sel[selected], log.p = TRUE)) +
    sum(pnorm(lin_sel[!selected], lower.tail = FALSE, log.p = TRUE))
  expect_equal(approx_loglik, exact_loglik, tolerance = 1e-4)
})
