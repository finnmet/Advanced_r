# predict method: link, unconditional response expectation, and selection
# probability, on the stored data or on new data.

#' Build the design matrix a prediction type needs
#'
#' With `newdata = NULL` the stored model matrix is returned (all n units).
#' Otherwise the matrix is rebuilt from the stored terms and factor levels
#' so that factor codings match the fit exactly; only the variables of the
#' requested equation are required in `newdata`.
#'
#' @noRd
build_prediction_matrix <- function(object, newdata, equation) {
  if (is.null(newdata)) {
    return(object$model[[if (equation == "outcome") "x" else "z"]])
  }
  checkmate::assert_data_frame(newdata, min.rows = 1)
  terms_rhs <- stats::delete.response(
    object$model[[paste0("terms_", equation)]]
  )
  missing_vars <- setdiff(all.vars(terms_rhs), names(newdata))
  if (length(missing_vars) > 0L) {
    stop("`newdata` is missing variables required for this prediction ",
      "type: ", paste(missing_vars, collapse = ", "), ".",
      call. = FALSE
    )
  }
  # na.pass + explicit finiteness check below: na.omit would silently drop
  # rows, returning fewer predictions than newdata rows.
  frame <- stats::model.frame(
    terms_rhs,
    data = newdata, na.action = stats::na.pass,
    xlev = object$model[[paste0("xlevels_", equation)]]
  )
  design <- stats::model.matrix(terms_rhs, frame)
  bad <- colnames(design)[colSums(!is.finite(design)) > 0L]
  if (length(bad) > 0L) {
    stop("`newdata` contains NA/NaN/Inf in required covariates: ",
      paste(bad, collapse = ", "), ".",
      call. = FALSE
    )
  }
  design
}

#' Predict from a fitted Poisson selection model
#'
#' @param object An object of class `"poisselect"`.
#' @param newdata Optional `data.frame` with the covariates needed by the
#'   requested `type` (`"link"`/`"response"`: outcome covariates,
#'   `"pselect"`: selection covariates). Default `NULL` predicts for the
#'   data the model was fitted on (all n units).
#' @param type Type of prediction: `"response"` (default) returns the
#'   unconditional population expectation
#'   `E[Y | x] = exp(x'beta + sigma^2 / 2)` (averaged over the outcome
#'   heterogeneity eps - not the selected mean `E[Y | x, s = 1]`),
#'   `"link"` returns the linear predictor `x'beta`, and `"pselect"`
#'   returns the selection probability `Phi(z'gamma)`.
#' @param ... Ignored (present for compatibility with the generic).
#'
#' @return Named numeric vector of predictions (names are the row names of
#'   the data used).
#'
#' @examples
#' set.seed(1)
#' d <- simulate_poisselect(300,
#'   beta = c(1, 0.5), gamma = c(0.5, 1, -0.8),
#'   sigma = 0.5, rho = 0.5
#' )
#' fit <- poisselect(y ~ x1, s ~ x1 + w1, data = d, k = 8)
#' head(predict(fit))
#' predict(fit, newdata = data.frame(x1 = c(-1, 0, 1)), type = "link")
#'
#' @export
predict.poisselect <- function(object, newdata = NULL,
                               type = c("response", "link", "pselect"),
                               ...) {
  # match.arg() validates and resolves the default ("response").
  type <- match.arg(type)
  estimates <- object$coefficients
  # link/response live on the outcome equation, pselect on the selection
  # equation; only that equation's covariates are needed (and required).
  equation <- if (type == "pselect") "selection" else "outcome"
  design <- build_prediction_matrix(object, newdata, equation)
  linear_predictor <- drop(design %*% estimates[[equation]])
  # switch() for value dispatch instead of an if/else chain (style guide).
  switch(type,
    link = linear_predictor,
    response = exp(linear_predictor + estimates$sigma^2 / 2),
    pselect = stats::pnorm(linear_predictor)
  )
}
