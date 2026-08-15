# Shared fixtures: one moderate simulated dataset and one fitted model,
# reused across test files to keep the whole suite fast. Every simulated
# dataset is seeded (reproducibility rule).

make_test_data <- function(n = 600, seed = 42) {
  set.seed(seed)
  simulate_poisselect(n,
    beta = c(1, 0.5), gamma = c(0.5, 1, -0.8),
    sigma = 0.5, rho = 0.6
  )
}

test_data <- make_test_data()
test_fit <- poisselect(y ~ x1, s ~ x1 + w1, data = test_data, k = 10)
