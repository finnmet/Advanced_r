# summary method: full coefficient tables with z tests for both equations,
# following the standard R pattern of a summary object plus its own print
# method.

#' Build one coefficient table (Estimate, SE, z, p)
#'
#' Per coefficient the null hypothesis "coefficient = 0" is tested against
#' the two-sided alternative: z = Estimate/SE is asymptotically standard
#' normal under H0, hence p = 2 * Phi(-|z|).
#'
#' @noRd
build_coefficient_table <- function(estimates, standard_errors) {
  z_value <- estimates / standard_errors
  p_value <- 2 * stats::pnorm(-abs(z_value))
  cbind(
    "Estimate" = estimates, "Std. Error" = standard_errors,
    "z value" = z_value, "Pr(>|z|)" = p_value
  )
}

#' Summarise a fitted Poisson selection model
#'
#' Computes full coefficient tables (estimate, standard error, z value,
#' two-sided p value) for both the outcome and the selection equation and
#' collects rho, sigma, log-likelihood and AIC.
#'
#' @param object An object of class `"poisselect"`.
#' @param ... Ignored (present for compatibility with the generic).
#'
#' @return An object of class `"summary.poisselect"`: a list with the
#'   `call`, the matrices `coefficients_outcome` and
#'   `coefficients_selection` (columns `Estimate`, `Std. Error`, `z value`,
#'   `Pr(>|z|)`), `sigma`, `rho`, `loglik`, `aic`, `n`, `n_selected`,
#'   `convergence` and `k`.
#'
#' @examples
#' set.seed(1)
#' d <- simulate_poisselect(300,
#'   beta = c(1, 0.5), gamma = c(0.5, 1, -0.8),
#'   sigma = 0.5, rho = 0.5
#' )
#' fit <- poisselect(y ~ x1, s ~ x1 + w1, data = d, k = 8)
#' summary(fit)
#'
#' @export
summary.poisselect <- function(object, ...) {
  structure(
    list(
      call = object$call,
      coefficients_outcome = build_coefficient_table(
        object$coefficients$outcome, object$standard_errors$outcome
      ),
      coefficients_selection = build_coefficient_table(
        object$coefficients$selection, object$standard_errors$selection
      ),
      sigma = object$coefficients$sigma,
      rho = object$coefficients$rho,
      loglik = object$loglik,
      aic = object$aic,
      n = object$n,
      n_selected = object$n_selected,
      convergence = object$convergence,
      k = object$model$k
    ),
    class = "summary.poisselect"
  )
}

#' Print a poisselect summary
#'
#' @param x An object of class `"summary.poisselect"`.
#' @param ... Passed on to [stats::printCoefmat()].
#'
#' @return `x`, invisibly.
#'
#' @examples
#' set.seed(1)
#' d <- simulate_poisselect(300,
#'   beta = c(1, 0.5), gamma = c(0.5, 1, -0.8),
#'   sigma = 0.5, rho = 0.5
#' )
#' fit <- poisselect(y ~ x1, s ~ x1 + w1, data = d, k = 8)
#' print(summary(fit))
#'
#' @export
print.summary.poisselect <- function(x, ...) {
  cat("Poisson selection model (poisselect)\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\nOutcome equation (Poisson):\n")
  stats::printCoefmat(x$coefficients_outcome, signif.stars = TRUE, ...)
  cat("\nSelection equation (probit):\n")
  stats::printCoefmat(x$coefficients_selection, signif.stars = TRUE, ...)
  cat("\nsigma: ", format(x$sigma, digits = 4L),
    "   rho: ", format(x$rho, digits = 4L), "\n",
    sep = ""
  )
  cat("Log-likelihood: ", format(x$loglik, digits = 7L),
    "   AIC: ", format(x$aic, digits = 7L),
    "   (K = ", x$k, " quadrature nodes)\n",
    sep = ""
  )
  cat("Observations: ", x$n, " (selected: ", x$n_selected, ")\n", sep = "")
  status <- if (x$convergence == 0L) "successful" else "NOT converged"
  cat("Convergence: ", status, " (optim code ", x$convergence, ")\n",
    sep = ""
  )
  invisible(x)
}
