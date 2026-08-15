# Main fitting function and its private helpers: working-scale
# parameterisation, optim objective, covariance extraction, result assembly
# and fit-time warnings.

#' Fit a Poisson selection model by maximum likelihood
#'
#' Estimates the Poisson selection model
#' \deqn{\ln \mu_i = x_i'\beta + \varepsilon_i, \quad y_i \sim
#'   \mathrm{Poisson}(\mu_i)}
#' \deqn{s_i^* = z_i'\gamma + u_i, \quad s_i = 1\{s_i^* > 0\}}
#' where the outcome \eqn{y_i} is only observed for selected units
#' (\eqn{s_i = 1}) and the errors are jointly normal with
#' \eqn{\mathrm{Var}(\varepsilon_i) = \sigma^2},
#' \eqn{\mathrm{Var}(u_i) = 1} and
#' \eqn{\mathrm{Corr}(\varepsilon_i, u_i) = \rho}. A non-zero \eqn{\rho}
#' makes selection non-ignorable; joint maximum likelihood estimation of
#' both equations corrects the bias a plain Poisson GLM would have on the
#' selected subsample. The likelihood integral over the outcome
#' heterogeneity is approximated by Gauss-Hermite quadrature with `k` nodes
#' and evaluated on the log scale with the log-sum-exp trick. The range
#' restrictions \eqn{\sigma > 0} and \eqn{-1 < \rho < 1} are enforced by
#' optimising over \eqn{\log \sigma} and \eqn{\mathrm{atanh}\, \rho}.
#'
#' @param outcome Two-sided formula for the outcome equation, e.g.
#'   `y ~ x1 + x2`. The response are the observed counts; it may (and
#'   should) be `NA` wherever `s = 0`, those values are ignored.
#' @param selection Two-sided formula for the selection equation, e.g.
#'   `s ~ x1 + w1`. The response is the fully observed 0/1 (or logical)
#'   selection indicator. For credible identification the selection
#'   equation should contain at least one variable that is not in the
#'   outcome equation (exclusion restriction); a warning is issued
#'   otherwise.
#' @param data A `data.frame` containing all variables used in the two
#'   formulas.
#' @param k Single integerish value between 2 and 500: the number of
#'   Gauss-Hermite quadrature nodes. Default `20`.
#' @param start Optional numeric vector of starting values in the order
#'   `c(beta, gamma, sigma, rho)` (length: number of outcome coefficients
#'   + number of selection coefficients + 2), with `sigma > 0` and
#'   `abs(rho) < 1`. Default `NULL` computes starting values from a
#'   Poisson GLM (outcome) and a probit GLM (selection).
#' @param control List of control arguments passed on to
#'   [stats::optim()]; allowed entries are `maxit`, `reltol`, `abstol`,
#'   `trace`, `REPORT`, `fnscale`, `parscale` and `ndeps`. `maxit`
#'   defaults to 500.
#'
#' @return An object of class `"poisselect"`: a list with elements
#' \describe{
#'   \item{`coefficients`}{List with named vectors `outcome` (beta) and
#'     `selection` (gamma) plus scalars `sigma` and `rho`.}
#'   \item{`standard_errors`}{List with named vectors `outcome` and
#'     `selection`; `NA` if the Hessian could not be inverted.}
#'   \item{`vcov`}{Covariance matrix of the (beta, gamma) block with
#'     `outcome:`/`selection:` prefixed dimnames.}
#'   \item{`loglik`}{Approximated log-likelihood at the optimum.}
#'   \item{`df`}{Number of estimated parameters (p + q + 2).}
#'   \item{`aic`}{Akaike information criterion, `-2 loglik + 2 df`.}
#'   \item{`convergence`}{Convergence code from [stats::optim()] (0 means
#'     success).}
#'   \item{`optim_message`}{Message from [stats::optim()], or `NULL`.}
#'   \item{`n`}{Number of units.}
#'   \item{`n_selected`}{Number of selected units.}
#'   \item{`call`}{The matched call.}
#'   \item{`model`}{List with the model pieces (`y`, `s`, `x`, `z`, both
#'     `terms` objects, both `xlevels` lists, `k` and the quadrature
#'     rule `gh`) used by the methods.}
#' }
#'
#' @examples
#' set.seed(1)
#' d <- simulate_poisselect(400,
#'   beta = c(1, 0.5), gamma = c(0.5, 1, -0.8),
#'   sigma = 0.5, rho = 0.6
#' )
#' fit <- poisselect(y ~ x1, s ~ x1 + w1, data = d, k = 10)
#' print(fit)
#' summary(fit)
#' head(predict(fit))
#'
#' @export
poisselect <- function(outcome, selection, data, k = 20, start = NULL,
                       control = list()) {
  cl <- match.call()
  k <- check_poisselect_inputs(outcome, selection, data, k, start, control)
  model <- build_model_data(outcome, selection, data)
  p <- ncol(model$x)
  q <- ncol(model$z)
  if (!is.null(start)) {
    check_start_values(start, p, q)
  }
  gh <- compute_gh_rule(k)
  par0 <- build_start_parameters(start, model, p, q)
  objective <- build_objective(model, gh, p, q)
  fit <- stats::optim(
    par = par0, fn = objective, method = "BFGS",
    control = utils::modifyList(list(maxit = 500), control),
    hessian = TRUE
  )
  result <- build_poisselect_result(fit, model, gh, k, p, q, cl)
  warn_on_fit_issues(fit, result)
  result
}

#' Assemble the working-scale start vector
#'
#' The working scale replaces sigma by log(sigma) and rho by atanh(rho):
#' the back-transformations exp() and tanh() map the whole real line into
#' the admissible ranges sigma > 0 and -1 < rho < 1, so unconstrained BFGS
#' automatically respects the range restrictions of the assignment.
#'
#' @noRd
build_start_parameters <- function(start, model, p, q) {
  if (is.null(start)) {
    values <- compute_starting_values(model$y, model$x, model$z, model$s)
    start <- c(values$beta, values$gamma, values$sigma, values$rho)
  }
  c(start[seq_len(p + q)], log(start[p + q + 1L]), atanh(start[p + q + 2L]))
}

#' Build the optim objective (negative log-likelihood, working scale)
#'
#' @noRd
build_objective <- function(model, gh, p, q) {
  function(par) {
    beta <- par[seq_len(p)]
    gamma <- par[p + seq_len(q)]
    # Back-transform from the working scale (see build_start_parameters).
    sigma <- exp(par[p + q + 1L])
    rho <- tanh(par[p + q + 2L])
    loglik <- compute_poisselect_loglik(
      beta, gamma, sigma, rho, model$y,
      model$x, model$z, model$s, gh
    )
    # A large finite penalty instead of Inf/NaN lets optim's line search
    # back off from numerically degenerate regions without aborting.
    if (!is.finite(loglik)) {
      return(1e10)
    }
    -loglik
  }
}

#' Invert the Hessian and extract the (beta, gamma) covariance block
#'
#' optim() minimised the negative log-likelihood, so its Hessian at the
#' optimum is the observed information and the inverse is the covariance
#' matrix on the working scale. The (beta, gamma) block is invariant to the
#' exp/tanh reparameterisation of sigma and rho: the Jacobian of the
#' back-transformation is block-diagonal with an identity on the
#' coefficient block, so that block can be reported directly as the
#' natural-scale vcov. Any failure (singular or non-positive-definite
#' Hessian) degrades to NA results with a warning instead of an error.
#'
#' @noRd
compute_vcov_and_se <- function(hessian, coef_names, p, q) {
  block <- seq_len(p + q)
  vcov_working <- tryCatch(solve(hessian), error = function(condition) NULL)
  variances <- if (is.null(vcov_working)) {
    NULL
  } else {
    diag(vcov_working)[block]
  }
  hessian_usable <- !is.null(variances) && all(is.finite(variances)) &&
    all(variances > 0)
  if (hessian_usable) {
    vcov_block <- vcov_working[block, block, drop = FALSE]
  } else {
    warning("The Hessian at the optimum is not invertible or not positive ",
      "definite; standard errors and vcov are reported as NA.",
      call. = FALSE
    )
    vcov_block <- matrix(NA_real_, nrow = p + q, ncol = p + q)
  }
  dimnames(vcov_block) <- list(coef_names, coef_names)
  list(vcov = vcov_block, se = sqrt(diag(vcov_block)))
}

#' Assemble the S3 result object of class "poisselect"
#'
#' @noRd
build_poisselect_result <- function(fit, model, gh, k, p, q, cl) {
  beta <- stats::setNames(fit$par[seq_len(p)], colnames(model$x))
  gamma <- stats::setNames(fit$par[p + seq_len(q)], colnames(model$z))
  sigma <- exp(fit$par[p + q + 1L])
  rho <- tanh(fit$par[p + q + 2L])
  coef_names <- c(
    paste0("outcome:", colnames(model$x)),
    paste0("selection:", colnames(model$z))
  )
  inference <- compute_vcov_and_se(fit$hessian, coef_names, p, q)
  loglik <- -fit$value
  df <- p + q + 2L
  structure(
    list(
      coefficients = list(
        outcome = beta, selection = gamma,
        sigma = sigma, rho = rho
      ),
      standard_errors = list(
        outcome = stats::setNames(
          inference$se[seq_len(p)],
          colnames(model$x)
        ),
        selection = stats::setNames(
          inference$se[p + seq_len(q)],
          colnames(model$z)
        )
      ),
      vcov = inference$vcov,
      loglik = loglik,
      df = df,
      aic = -2 * loglik + 2 * df,
      convergence = fit$convergence,
      optim_message = fit$message,
      n = length(model$s),
      n_selected = sum(model$s == 1L),
      call = cl,
      model = c(model, list(k = k, gh = gh))
    ),
    class = "poisselect"
  )
}

#' Emit fit-time warnings (non-convergence, boundary estimates)
#'
#' These are warnings rather than errors: the returned object is still
#' useful for diagnosis, the user just must not trust it blindly.
#'
#' @noRd
warn_on_fit_issues <- function(fit, result) {
  if (fit$convergence != 0L) {
    suffix <- if (is.null(fit$message)) "" else paste0(": ", fit$message)
    warning("optim did not converge (code ", fit$convergence, ")", suffix,
      "; consider more iterations via control = list(maxit = ...).",
      call. = FALSE
    )
  }
  if (abs(result$coefficients$rho) > 0.99) {
    warning("The estimate of rho is at or near the boundary ",
      "(|rho| > 0.99); inference may be unreliable.",
      call. = FALSE
    )
  }
  if (result$coefficients$sigma < 1e-4) {
    warning("The estimate of sigma is close to zero (< 1e-4), i.e. at the ",
      "parameter boundary; inference may be unreliable.",
      call. = FALSE
    )
  }
  invisible(result)
}
