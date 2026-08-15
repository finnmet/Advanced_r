# Starting values for the maximum likelihood optimisation. Two cheap GLMs
# provide good neighbourhoods for beta and gamma; sigma and rho start at
# neutral values.

#' Fit one starting-value GLM with a safe fallback
#'
#' Wraps `stats::glm.fit()` in `tryCatch()` so a degenerate GLM can never
#' abort the estimation; `suppressWarnings()` silences glm.fit chatter
#' (e.g. fitted probabilities of 0 or 1) that is harmless for mere starting
#' values. Any failure - an error, NA or non-finite coefficients - falls
#' back to the supplied default.
#'
#' @noRd
compute_glm_start <- function(design, response, family, fallback) {
  estimate <- tryCatch(
    suppressWarnings(
      stats::glm.fit(design, response, family = family)$coefficients
    ),
    error = function(condition) NULL
  )
  if (is.null(estimate) || anyNA(estimate) || any(!is.finite(estimate))) {
    return(fallback)
  }
  unname(estimate)
}

#' Compute starting values for the poisselect optimisation
#'
#' `beta` starts at a Poisson GLM fit on the selected rows, `gamma` at a
#' probit GLM fit on all rows - `stats::glm()`/`glm.fit()` may be used for
#' starting values only (assignment rule).
#'
#' @param y,x,z,s Model pieces as built by `build_model_data()`.
#'
#' @return List with `beta`, `gamma`, `sigma` (= 1), `rho` (= 0).
#'
#' @noRd
compute_starting_values <- function(y, x, z, s) {
  selected <- s == 1L
  # Fallback for beta: log mean count for the intercept (offset by 0.5 so
  # an all-zero outcome cannot produce log(0)) and zero slopes.
  beta <- compute_glm_start(
    x[selected, , drop = FALSE], y[selected],
    family = stats::poisson(),
    fallback = c(log(mean(y[selected]) + 0.5), rep(0, ncol(x) - 1L))
  )
  gamma <- compute_glm_start(
    z, s,
    family = stats::binomial(link = "probit"),
    fallback = rep(0, ncol(z))
  )
  # Neutral starts for the error parameters: unit variance and rho = 0,
  # the point of ignorable selection.
  list(beta = beta, gamma = gamma, sigma = 1, rho = 0)
}
