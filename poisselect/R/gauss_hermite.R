# Gauss-Hermite quadrature rule computed from first principles via the
# Golub-Welsch eigenvalue method, as required by the assignment (no
# quadrature packages allowed).

#' Compute a Gauss-Hermite quadrature rule via Golub-Welsch
#'
#' Nodes and weights for the weight function `exp(-t^2)` ("physicists'"
#' Hermite polynomials). The nodes are the eigenvalues of the symmetric
#' tridiagonal Jacobi matrix of the three-term Hermite recurrence and the
#' weights follow from the first components of its normalised eigenvectors.
#'
#' @param node_count Single integerish value `>= 1`, the number of nodes.
#'
#' @return A list with elements `nodes` (ascending), `weights`, and
#'   `log_weights` (`log(weights)`, precomputed once because the
#'   log-likelihood only ever needs the logs - DRY).
#'
#' @noRd
compute_gh_rule <- function(node_count) {
  node_count <- checkmate::asInt(node_count, lower = 1)
  # The Jacobi matrix of the Hermite recurrence has a zero diagonal and
  # off-diagonal entries sqrt(k / 2); fill both bands via index matrices so
  # no loop over entries is needed.
  jacobi <- matrix(0, nrow = node_count, ncol = node_count)
  band <- cbind(seq_len(node_count - 1L), seq_len(node_count - 1L) + 1L)
  jacobi[band] <- sqrt(seq_len(node_count - 1L) / 2)
  jacobi[band[, c(2L, 1L), drop = FALSE]] <- sqrt(seq_len(node_count - 1L) / 2)
  decomposition <- eigen(jacobi, symmetric = TRUE)
  # eigen() returns eigenvalues in decreasing order; the rule is reported
  # with ascending nodes, so reverse values and eigenvectors consistently.
  ascending <- rev(seq_len(node_count))
  nodes <- decomposition$values[ascending]
  weights <- sqrt(pi) * decomposition$vectors[1L, ascending]^2
  list(nodes = nodes, weights = weights, log_weights = log(weights))
}
