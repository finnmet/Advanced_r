# Standard extractor methods: coef, vcov, logLik.

#' Extract all coefficients of a Poisson selection model
#'
#' @param object An object of class `"poisselect"`.
#' @param ... Ignored (present for compatibility with the generic).
#'
#' @return Named numeric vector: the outcome coefficients (prefixed
#'   `outcome:`), the selection coefficients (prefixed `selection:`), then
#'   `sigma` and `rho`.
#'
#' @examples
#' set.seed(1)
#' d <- simulate_poisselect(300,
#'   beta = c(1, 0.5), gamma = c(0.5, 1, -0.8),
#'   sigma = 0.5, rho = 0.5
#' )
#' fit <- poisselect(y ~ x1, s ~ x1 + w1, data = d, k = 8)
#' coef(fit)
#'
#' @export
coef.poisselect <- function(object, ...) {
  estimates <- object$coefficients
  c(
    stats::setNames(
      estimates$outcome,
      paste0("outcome:", names(estimates$outcome))
    ),
    stats::setNames(
      estimates$selection,
      paste0("selection:", names(estimates$selection))
    ),
    sigma = estimates$sigma, rho = estimates$rho
  )
}

#' Covariance matrix of the outcome and selection coefficients
#'
#' @param object An object of class `"poisselect"`.
#' @param ... Ignored (present for compatibility with the generic).
#'
#' @return The (p + q) x (p + q) covariance matrix of `(beta, gamma)` with
#'   `outcome:`/`selection:` prefixed dimnames (all `NA` if the Hessian
#'   could not be inverted).
#'
#' @examples
#' set.seed(1)
#' d <- simulate_poisselect(300,
#'   beta = c(1, 0.5), gamma = c(0.5, 1, -0.8),
#'   sigma = 0.5, rho = 0.5
#' )
#' fit <- poisselect(y ~ x1, s ~ x1 + w1, data = d, k = 8)
#' vcov(fit)
#'
#' @export
vcov.poisselect <- function(object, ...) {
  object$vcov
}

#' Log-likelihood of a fitted Poisson selection model
#'
#' @param object An object of class `"poisselect"`.
#' @param ... Ignored (present for compatibility with the generic).
#'
#' @return An object of class `"logLik"` with attributes `df` (number of
#'   estimated parameters) and `nobs` (number of units), so that e.g. the
#'   [stats::AIC()] generic works out of the box.
#'
#' @examples
#' set.seed(1)
#' d <- simulate_poisselect(300,
#'   beta = c(1, 0.5), gamma = c(0.5, 1, -0.8),
#'   sigma = 0.5, rho = 0.5
#' )
#' fit <- poisselect(y ~ x1, s ~ x1 + w1, data = d, k = 8)
#' logLik(fit)
#' AIC(fit)
#'
#' @export
logLik.poisselect <- function(object, ...) {
  structure(object$loglik,
    df = object$df, nobs = object$n,
    class = "logLik"
  )
}
