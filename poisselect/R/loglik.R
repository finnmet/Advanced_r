# Approximated log-likelihood of the Poisson selection model. Everything is
# vectorised over units x quadrature nodes (n x K matrices); there is no
# loop over observations anywhere.

#' Build the node-level matrices shared by likelihood and plot method
#'
#' Both the log-likelihood and the model-implied count distribution (plot 1)
#' need `log(mu_ik) = x_i'beta + sqrt(2) * sigma * t_k` and
#' `eta_ik = (z_i'gamma + sqrt(2) * rho * t_k) / sqrt(1 - rho^2)` at the same
#' parameter values, so the construction lives in one shared helper (DRY).
#'
#' @param beta,gamma,sigma,rho Parameters on the natural scale.
#' @param x,z Design matrices for the rows the matrices are needed for.
#' @param gh Gauss-Hermite rule from `compute_gh_rule()`.
#'
#' @return List with `log_mu` and `eta`, both `nrow(x) x K` matrices.
#'
#' @noRd
compute_node_matrices <- function(beta, gamma, sigma, rho, x, z, gh) {
  # outer() sweeps the scaled nodes across all units at once: row i, column
  # k holds the linear predictor of unit i evaluated at node t_k.
  log_mu <- outer(drop(x %*% beta), sqrt(2) * sigma * gh$nodes, `+`)
  eta <- outer(drop(z %*% gamma), sqrt(2) * rho * gh$nodes, `+`) /
    sqrt(1 - rho^2)
  list(log_mu = log_mu, eta = eta)
}

#' Approximated log-likelihood of the Poisson selection model
#'
#' Evaluates the Gauss-Hermite approximation of the log-likelihood entirely
#' on the log scale: the node sum for selected units is combined with the
#' log-sum-exp trick so that large counts cannot over- or underflow the
#' Poisson probabilities.
#'
#' @param beta,gamma,sigma,rho Parameters on the natural scale.
#' @param y Outcome counts (may be `NA` where `s = 0`; those entries are
#'   never touched).
#' @param x,z Full design matrices (all n rows).
#' @param s Selection indicator (integer 0/1, length n).
#' @param gh Gauss-Hermite rule from `compute_gh_rule()`.
#'
#' @return The scalar approximated log-likelihood.
#'
#' @noRd
compute_poisselect_loglik <- function(beta, gamma, sigma, rho, y, x, z, s,
                                      gh) {
  selected <- s == 1L
  lin_sel <- drop(z %*% gamma)
  # Unselected term ln(1 - Phi(z'gamma)): pnorm() returns the log of the
  # upper tail directly; log(1 - pnorm(.)) would underflow to -Inf for
  # large linear predictors.
  ll_unselected <- stats::pnorm(lin_sel[!selected],
    lower.tail = FALSE,
    log.p = TRUE
  )
  # Selected term: n_sel x K matrices of log(mu_ik) and eta_ik.
  node <- compute_node_matrices(
    beta, gamma, sigma, rho,
    x[selected, , drop = FALSE],
    z[selected, , drop = FALSE], gh
  )
  log_p <- stats::dpois(y[selected], lambda = exp(node$log_mu), log = TRUE)
  # y recycles down each column of the n_sel x K matrix (column-major), so
  # each row keeps its own count; the explicit dim<- guarantees the matrix
  # shape because the d/p/q functions do not reliably keep the dims of
  # their second argument.
  dim(log_p) <- dim(node$log_mu)
  log_phi <- stats::pnorm(node$eta, log.p = TRUE)
  # Log-scale summands a_ik = ln(w_k) - ln(pi)/2 + ln p(y|mu_ik)
  # + ln Phi(eta_ik), combined over nodes with the log-sum-exp trick:
  # subtracting the row maximum before exp() keeps every term in range.
  a <- sweep(log_p + log_phi,
    MARGIN = 2, STATS = gh$log_weights,
    FUN = `+`
  ) - 0.5 * log(pi)
  a_max <- apply(a, 1, max)
  ll_selected <- a_max + log(rowSums(exp(a - a_max)))
  sum(ll_selected) + sum(ll_unselected)
}
