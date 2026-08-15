# Data generator matching the model exactly; used by the recovery tests,
# the examples, and for trying out the package.

#' Simulate data from the Poisson selection model
#'
#' Generates a dataset that follows the poisselect model exactly. Shared
#' covariates `x1, ..., x_{p-1}` (iid standard normal) enter both
#' equations; additional covariates `w1, ..., w_{q-p}` (iid standard
#' normal) enter only the selection equation and act as exclusion
#' restrictions. The errors are built as `eps = sigma * v` and
#' `u = rho * v + sqrt(1 - rho^2) * e` with independent standard normal
#' `v`, `e`, which yields `Var(u) = 1` and `Corr(eps, u) = rho` exactly.
#' The outcome is drawn as `y ~ Poisson(exp(x'beta + eps))`, selection as
#' `s = 1{z'gamma + u > 0}`, and `y` is set to `NA` wherever `s = 0` -
#' exactly like real data for this model.
#'
#' @param n Single integerish value `>= 50`, the number of units.
#' @param beta Finite numeric vector of outcome coefficients (intercept
#'   first), length `p >= 1`.
#' @param gamma Finite numeric vector of selection coefficients (intercept
#'   first), length `q >= p`; the last `q - p` entries belong to the
#'   selection-only covariates `w`.
#' @param sigma Single positive finite number: standard deviation of the
#'   outcome heterogeneity.
#' @param rho Single number strictly between -1 and 1: correlation of the
#'   two error terms.
#'
#' @return A `data.frame` with columns `y` (counts, `NA` where `s = 0`),
#'   `s` (0/1), the shared covariates `x1, ...` (if `p > 1`) and the
#'   selection-only covariates `w1, ...` (if `q > p`).
#'
#' @examples
#' set.seed(1)
#' d <- simulate_poisselect(200,
#'   beta = c(1, 0.5), gamma = c(0.5, 1, -0.8),
#'   sigma = 0.5, rho = 0.6
#' )
#' head(d)
#' table(d$s)
#'
#' @export
simulate_poisselect <- function(n, beta, gamma, sigma, rho) {
  n <- checkmate::asInt(n, lower = 50)
  checkmate::assert_numeric(beta,
    finite = TRUE, any.missing = FALSE,
    min.len = 1
  )
  checkmate::assert_numeric(gamma,
    finite = TRUE, any.missing = FALSE,
    min.len = 1
  )
  if (length(gamma) < length(beta)) {
    stop("`gamma` must be at least as long as `beta`: the selection ",
      "equation contains all shared covariates plus the exclusion ",
      "variables.",
      call. = FALSE
    )
  }
  checkmate::assert_number(sigma, lower = 1e-8, finite = TRUE)
  checkmate::assert_number(rho, finite = TRUE)
  if (abs(rho) >= 1) {
    stop("`rho` must lie strictly between -1 and 1.", call. = FALSE)
  }
  p <- length(beta)
  q <- length(gamma)
  # Shared covariates enter both design matrices; the w block only z.
  x_shared <- matrix(stats::rnorm(n * (p - 1L)), nrow = n)
  w_exclusion <- matrix(stats::rnorm(n * (q - p)), nrow = n)
  x <- cbind(1, x_shared)
  z <- cbind(1, x_shared, w_exclusion)
  # Error construction: with independent v, e ~ N(0, 1), eps = sigma * v
  # and u = rho * v + sqrt(1 - rho^2) * e give Var(u) = 1 and
  # Cov(eps, u) = sigma * rho, i.e. Corr(eps, u) = rho exactly.
  v <- stats::rnorm(n)
  e <- stats::rnorm(n)
  eps <- sigma * v
  u <- rho * v + sqrt(1 - rho^2) * e
  y_full <- stats::rpois(n, lambda = exp(drop(x %*% beta) + eps))
  s <- as.integer(drop(z %*% gamma) + u > 0)
  # y is unobserved (NA) wherever the unit is not selected.
  result <- data.frame(y = ifelse(s == 1L, y_full, NA_integer_), s = s)
  if (p > 1L) {
    result[paste0("x", seq_len(p - 1L))] <- x_shared
  }
  if (q > p) {
    result[paste0("w", seq_len(q - p))] <- w_exclusion
  }
  result
}
