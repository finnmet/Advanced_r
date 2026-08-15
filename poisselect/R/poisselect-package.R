# Package-level documentation. All functions reference other namespaces
# explicitly via :: (stats, graphics, utils, checkmate), so no importFrom
# directives are needed here.

#' poisselect: Poisson Selection Model for Count Data
#'
#' Maximum likelihood estimation of a Poisson selection model for count data
#' that are only observed for a selected subsample. The model couples a
#' Poisson outcome equation with log-normal heterogeneity,
#' `ln mu = x'beta + eps`, and a probit selection equation,
#' `s = 1{z'gamma + u > 0}`, through jointly normal errors `(eps, u)` with
#' `Corr(eps, u) = rho`. A non-zero `rho` makes selection non-ignorable; the
#' joint likelihood corrects the resulting bias of a plain Poisson GLM fit
#' to the selected subsample. The likelihood integral over the outcome
#' heterogeneity is approximated by Gauss-Hermite quadrature (Golub-Welsch
#' nodes and weights) and evaluated on the log scale with the log-sum-exp
#' trick for numerical stability.
#'
#' @keywords internal
"_PACKAGE"
