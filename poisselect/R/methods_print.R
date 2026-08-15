# print method: a compact overview of the fitted model.

#' Print a fitted Poisson selection model
#'
#' Shows the call, case counts, log-likelihood, convergence status and the
#' point estimates of both equations.
#'
#' @param x An object of class `"poisselect"`.
#' @param digits Number of significant digits for the coefficients.
#' @param ... Ignored (present for compatibility with the generic).
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
#' print(fit)
#'
#' @export
print.poisselect <- function(x, digits = max(3L, getOption("digits") - 3L),
                             ...) {
  cat("Poisson selection model (poisselect)\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\nObservations: ", x$n, " (selected: ", x$n_selected, ")\n",
    sep = ""
  )
  cat("Log-likelihood: ", format(x$loglik, digits = digits + 4L),
    " (K = ", x$model$k, " quadrature nodes)\n",
    sep = ""
  )
  status <- if (x$convergence == 0L) "successful" else "NOT converged"
  cat("Convergence: ", status, " (optim code ", x$convergence, ")\n\n",
    sep = ""
  )
  cat("Outcome coefficients (Poisson equation):\n")
  print(round(x$coefficients$outcome, digits))
  cat("\nSelection coefficients (probit equation):\n")
  print(round(x$coefficients$selection, digits))
  cat("\nsigma: ", format(x$coefficients$sigma, digits = digits),
    "   rho: ", format(x$coefficients$rho, digits = digits), "\n",
    sep = ""
  )
  invisible(x)
}
