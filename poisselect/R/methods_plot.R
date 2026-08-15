# plot method: (1) observed vs. model-implied count distribution of the
# selected units, (2) log-likelihood profile along rho.

#' Model-implied count distribution of the selected units
#'
#' For a single selected unit,
#' `Pr(y = m | x, z, s = 1)` is approximated by the ratio of the
#' quadrature sums `sum_k (w_k/sqrt(pi)) p(m | mu_ik) Phi(eta_ik)` and
#' `sum_k (w_k/sqrt(pi)) Phi(eta_ik)`; the reported frequency of m is the
#' mean over the selected units. Reuses the likelihood's node matrices
#' (DRY).
#'
#' @noRd
compute_implied_distribution <- function(object, count_values) {
  model <- object$model
  estimates <- object$coefficients
  selected <- model$s == 1L
  node <- compute_node_matrices(
    estimates$outcome, estimates$selection, estimates$sigma, estimates$rho,
    model$x[selected, , drop = FALSE], model$z[selected, , drop = FALSE],
    model$gh
  )
  mu <- exp(node$log_mu)
  # Fold the quadrature weights into the Phi factor once; numerator and
  # denominator then only differ by the Poisson probability factor.
  phi_weighted <- sweep(stats::pnorm(node$eta),
    MARGIN = 2,
    STATS = model$gh$weights / sqrt(pi), FUN = `*`
  )
  denominator <- rowSums(phi_weighted)
  # Small fixed-size loop over count values with a pre-allocated result.
  implied <- numeric(length(count_values))
  for (index in seq_along(count_values)) {
    poisson_prob <- stats::dpois(count_values[index], mu)
    numerator <- rowSums(poisson_prob * phi_weighted)
    implied[index] <- mean(numerator / denominator)
  }
  implied
}

#' Plot 1: observed vs. model-implied count distribution
#'
#' @noRd
plot_count_distribution <- function(object, ...) {
  model <- object$model
  y_selected <- model$y[model$s == 1L]
  count_values <- 0:max(y_selected)
  observed <- tabulate(y_selected + 1L, nbins = max(y_selected) + 1L) /
    length(y_selected)
  implied <- compute_implied_distribution(object, count_values)
  positions <- graphics::barplot(
    observed,
    names.arg = count_values, xlab = "Count",
    ylab = "Relative frequency", ylim = c(0, max(observed, implied) * 1.1),
    main = "Observed vs. model-implied\ncount distribution (selected)",
    ...
  )
  graphics::points(positions, implied, type = "b", pch = 19, col = "red3")
  graphics::legend("topright",
    legend = c("Observed", "Model-implied"),
    pch = c(22, 19), pt.bg = c("grey", NA),
    col = c("black", "red3"),
    lty = c(NA, 1), bty = "n"
  )
}

#' Plot 2: log-likelihood profile along rho
#'
#' @noRd
plot_loglik_profile <- function(object, rho_grid, ...) {
  model <- object$model
  estimates <- object$coefficients
  # All other parameters are held at their estimates; vapply pre-allocates
  # the result vector (no growing objects).
  profile <- vapply(
    rho_grid,
    function(rho) {
      compute_poisselect_loglik(
        estimates$outcome, estimates$selection, estimates$sigma, rho,
        model$y, model$x, model$z, model$s, model$gh
      )
    },
    numeric(1)
  )
  graphics::plot(rho_grid, profile,
    type = "l", xlab = expression(rho),
    ylab = "Log-likelihood",
    main = expression("Log-likelihood along" ~ rho), ...
  )
  graphics::abline(v = estimates$rho, lty = 2)
  graphics::points(rho_grid[which.max(profile)], max(profile), pch = 19)
}

#' Plot diagnostics for a fitted Poisson selection model
#'
#' Plot 1 compares the observed relative frequencies of the selected counts
#' with the model-implied count distribution (averaged over the selected
#' units). Plot 2 shows the approximated log-likelihood along a grid of
#' rho values with all other parameters held at their estimates; the dashed
#' line marks the estimate, the dot the maximum on the grid.
#'
#' @param x An object of class `"poisselect"`.
#' @param which Subset of `c(1, 2)` selecting the plots to draw.
#' @param rho_grid Numeric grid of rho values strictly inside (-1, 1) for
#'   plot 2.
#' @param ... Passed on to the underlying plotting calls.
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
#' plot(fit)
#'
#' @export
plot.poisselect <- function(x, which = c(1, 2),
                            rho_grid = seq(-0.95, 0.95, by = 0.05), ...) {
  checkmate::assert_subset(which, choices = c(1, 2), empty.ok = FALSE)
  checkmate::assert_numeric(rho_grid,
    finite = TRUE, any.missing = FALSE,
    min.len = 2
  )
  if (any(abs(rho_grid) >= 1)) {
    stop("`rho_grid` values must lie strictly between -1 and 1.",
      call. = FALSE
    )
  }
  if (length(unique(which)) == 2L) {
    # Lecture pattern: remember the full graphics state and restore it on
    # exit, no matter how this function is left.
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))
    graphics::par(mfrow = c(1, 2))
  }
  if (1 %in% which) {
    plot_count_distribution(x, ...)
  }
  if (2 %in% which) {
    plot_loglik_profile(x, rho_grid, ...)
  }
  invisible(x)
}
