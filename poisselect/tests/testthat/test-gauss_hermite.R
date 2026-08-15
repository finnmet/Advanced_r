# Tests for the Golub-Welsch Gauss-Hermite rule: known small rules, exact
# quadrature identities, and input validation.

test_that("compute_gh_rule reproduces the known K = 1 and K = 2 rules", {
  rule_one <- compute_gh_rule(1)
  expect_equal(rule_one$nodes, 0, tolerance = 1e-10)
  expect_equal(rule_one$weights, sqrt(pi), tolerance = 1e-10)
  rule_two <- compute_gh_rule(2)
  expect_equal(rule_two$nodes, c(-1, 1) / sqrt(2), tolerance = 1e-10)
  expect_equal(rule_two$weights, rep(sqrt(pi) / 2, 2), tolerance = 1e-10)
})

test_that("compute_gh_rule satisfies the exact GH identities for K = 20", {
  rule <- compute_gh_rule(20)
  # Integrating 1 and t^2 against exp(-t^2) is exact for any K:
  # sum(w) = sqrt(pi) and sum(w * t^2) = sqrt(pi) / 2.
  expect_equal(sum(rule$weights), sqrt(pi), tolerance = 1e-10)
  expect_equal(sum(rule$weights * rule$nodes^2), sqrt(pi) / 2,
    tolerance = 1e-10
  )
  expect_equal(rule$nodes, sort(rule$nodes))
  expect_equal(rule$nodes, -rev(rule$nodes), tolerance = 1e-10)
  expect_equal(rule$log_weights, log(rule$weights))
})

test_that("compute_gh_rule rejects invalid node counts", {
  expect_error(compute_gh_rule(0))
  expect_error(compute_gh_rule(-3))
  expect_error(compute_gh_rule(2.5))
  expect_error(compute_gh_rule("a"))
  expect_error(compute_gh_rule(c(2, 3)))
  expect_error(compute_gh_rule(NA))
})
