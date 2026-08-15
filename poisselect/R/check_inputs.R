# All user-facing validation for poisselect(). Argument-level checks run
# before any data is touched ("fail fast"); data-dependent checks run on the
# built model pieces so their error messages can name the offending
# variables ("fail appropriately").

#' Assert that an argument is a two-sided formula
#'
#' Uses `assert_class()` rather than `assert_formula()` so that older
#' checkmate versions without `assert_formula()` also work; the class check
#' plus the length check is all that is needed here.
#'
#' @noRd
assert_two_sided_formula <- function(formula, arg_name) {
  checkmate::assert_class(formula, classes = "formula", .var.name = arg_name)
  if (length(formula) != 3L) {
    stop("`", arg_name, "` must be a two-sided formula of the form ",
      "response ~ covariates.",
      call. = FALSE
    )
  }
  invisible(formula)
}

#' Argument-level input checks for poisselect()
#'
#' Validates everything that can be validated without building the model
#' matrices. The length check for `start` needs the numbers of coefficients
#' and therefore lives in `check_start_values()`, which runs right after
#' the model matrices are built.
#'
#' @return `k` coerced to integer (for a clean pass-through to the
#'   quadrature rule).
#'
#' @noRd
check_poisselect_inputs <- function(outcome, selection, data, k, start,
                                    control) {
  checkmate::assert_data_frame(data, min.rows = 2, min.cols = 2)
  assert_two_sided_formula(outcome, "outcome")
  assert_two_sided_formula(selection, "selection")
  # Report all missing variables at once so the user can fix them in one go
  # instead of hitting them one by one.
  needed <- unique(c(all.vars(outcome), all.vars(selection)))
  missing_vars <- setdiff(needed, names(data))
  if (length(missing_vars) > 0L) {
    stop("The following variables are used in the formulas but missing ",
      "from `data`: ", paste(missing_vars, collapse = ", "), ".",
      call. = FALSE
    )
  }
  k <- checkmate::asInt(k, lower = 2, upper = 500)
  if (!is.null(start)) {
    checkmate::assert_numeric(start, finite = TRUE, any.missing = FALSE)
  }
  checkmate::assert_list(control)
  if (length(control) > 0L) {
    checkmate::assert_names(
      names(control),
      subset.of = c(
        "maxit", "reltol", "abstol", "trace", "REPORT",
        "fnscale", "parscale", "ndeps"
      )
    )
  }
  k
}

#' Check a user-supplied start vector once the dimensions are known
#'
#' @param start Numeric vector (already validated as finite and complete).
#' @param p,q Number of outcome and selection coefficients.
#'
#' @noRd
check_start_values <- function(start, p, q) {
  expected_length <- p + q + 2L
  if (length(start) != expected_length) {
    stop("`start` must have length ", expected_length, ": ", p,
      " outcome coefficients, ", q, " selection coefficients, then ",
      "sigma and rho (in this order).",
      call. = FALSE
    )
  }
  if (start[p + q + 1L] <= 0) {
    stop("`start`: the sigma entry (position ", p + q + 1L,
      ") must be strictly positive.",
      call. = FALSE
    )
  }
  if (abs(start[p + q + 2L]) >= 1) {
    stop("`start`: the rho entry (position ", p + q + 2L,
      ") must lie strictly between -1 and 1.",
      call. = FALSE
    )
  }
  invisible(start)
}

#' Validate the selection indicator and coerce it to integer 0/1
#'
#' @noRd
check_selection_indicator <- function(s) {
  if (anyNA(s)) {
    stop("The selection indicator must not contain missing values (NA); ",
      "s must be observed for every unit.",
      call. = FALSE
    )
  }
  # Logical indicators are meaningful (TRUE = selected), so coerce them
  # instead of rejecting.
  if (is.logical(s)) {
    s <- as.integer(s)
  }
  if (!is.numeric(s) || !all(s %in% c(0, 1))) {
    stop("The selection indicator must contain only 0 and 1 (or logical ",
      "values).",
      call. = FALSE
    )
  }
  s <- as.integer(s)
  if (all(s == 0L)) {
    stop("No selected observations (all s = 0): the outcome equation ",
      "cannot be estimated without any observed outcomes.",
      call. = FALSE
    )
  }
  if (all(s == 1L)) {
    stop("All observations are selected (all s = 1) - a plain Poisson GLM ",
      "would do; poisselect needs unselected units to identify the ",
      "selection equation.",
      call. = FALSE
    )
  }
  s
}

#' Validate the outcome counts on the selected rows
#'
#' `y` may be `NA` on unselected rows (that is the natural encoding for
#' this model); only the selected rows are checked.
#'
#' @noRd
check_outcome_response <- function(y_selected) {
  if (anyNA(y_selected)) {
    stop("The outcome variable contains NA on selected rows (s = 1); y ",
      "must be observed wherever s = 1 (NA is only allowed where ",
      "s = 0).",
      call. = FALSE
    )
  }
  if (!is.numeric(y_selected) || any(!is.finite(y_selected))) {
    stop("The outcome variable must be finite and numeric on selected ",
      "rows.",
      call. = FALSE
    )
  }
  if (any(y_selected < 0) || any(abs(y_selected - round(y_selected)) > 1e-8)) {
    stop("The outcome variable must contain non-negative integer counts ",
      "on selected rows.",
      call. = FALSE
    )
  }
  invisible(y_selected)
}

#' Data-dependent checks on the built model pieces
#'
#' @noRd
check_model_data <- function(y, s, x, z) {
  selected <- s == 1L
  check_outcome_response(y[selected])
  # The probit part of the likelihood uses all n rows, so the selection
  # covariates must be complete for everyone; name the offending columns.
  bad_z <- colnames(z)[colSums(!is.finite(z)) > 0L]
  if (length(bad_z) > 0L) {
    stop("Selection covariates contain NA/NaN/Inf values in: ",
      paste(bad_z, collapse = ", "),
      ". Selection covariates must be complete for all units.",
      call. = FALSE
    )
  }
  # The outcome equation only ever sees selected rows, so NA on unselected
  # rows is fine there - check the selected block only.
  x_selected <- x[selected, , drop = FALSE]
  bad_x <- colnames(x)[colSums(!is.finite(x_selected)) > 0L]
  if (length(bad_x) > 0L) {
    stop("Outcome covariates contain NA/NaN/Inf values on selected rows ",
      "in: ", paste(bad_x, collapse = ", "),
      ". Outcome covariates must be complete wherever s = 1.",
      call. = FALSE
    )
  }
  # Rank checks run after the finiteness checks above, otherwise qr() on a
  # matrix with NA would abort with an unhelpful low-level error.
  if (qr(x_selected)$rank < ncol(x)) {
    stop("The outcome model matrix has collinear columns on the selected ",
      "rows; drop redundant covariates from `outcome`.",
      call. = FALSE
    )
  }
  if (qr(z)$rank < ncol(z)) {
    stop("The selection model matrix has collinear columns; drop ",
      "redundant covariates from `selection`.",
      call. = FALSE
    )
  }
  n <- length(s)
  n_selected <- sum(selected)
  if (n_selected < ncol(x) + 2L) {
    stop("Too few selected observations: ", n_selected, " selected units ",
      "cannot support ", ncol(x), " outcome coefficients plus sigma ",
      "and rho.",
      call. = FALSE
    )
  }
  if (n < ncol(x) + ncol(z) + 2L) {
    stop("Too few observations: n = ", n, " units for ",
      ncol(x) + ncol(z) + 2L, " model parameters.",
      call. = FALSE
    )
  }
  # Identification: at least one selection covariate that is not part of
  # the outcome equation. Estimation still proceeds, hence a warning.
  selection_only <- setdiff(colnames(z), c("(Intercept)", colnames(x)))
  if (length(selection_only) == 0L) {
    warning("No exclusion restriction found: every selection covariate ",
      "also appears in the outcome equation, so the model is ",
      "identified only through functional form.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Build the model pieces from formulas and data, with all data checks
#'
#' @return List with `y`, `s`, `x`, `z`, the two `terms` objects and the
#'   two `xlevels` lists (both needed by `predict(newdata = ...)`).
#'
#' @noRd
build_model_data <- function(outcome, selection, data) {
  # na.pass is essential here: the default na.omit would silently drop NA
  # rows and the informative NA checks in check_model_data() could never
  # fire on the data the user actually passed.
  frame_outcome <- stats::model.frame(outcome,
    data = data,
    na.action = stats::na.pass
  )
  frame_selection <- stats::model.frame(selection,
    data = data,
    na.action = stats::na.pass
  )
  terms_outcome <- attr(frame_outcome, "terms")
  terms_selection <- attr(frame_selection, "terms")
  y <- stats::model.response(frame_outcome)
  s <- check_selection_indicator(stats::model.response(frame_selection))
  x <- stats::model.matrix(terms_outcome, frame_outcome)
  z <- stats::model.matrix(terms_selection, frame_selection)
  check_model_data(y, s, x, z)
  # The counts passed validation as integer-valued up to 1e-8; store them
  # rounded so dpois() never sees floating-point fuzz (which would trigger
  # its non-integer warning).
  list(
    y = round(as.numeric(y)), s = s, x = x, z = z,
    terms_outcome = terms_outcome, terms_selection = terms_selection,
    xlevels_outcome = stats::.getXlevels(terms_outcome, frame_outcome),
    xlevels_selection = stats::.getXlevels(terms_selection, frame_selection)
  )
}
