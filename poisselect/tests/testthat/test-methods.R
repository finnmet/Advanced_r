# Tests for the S3 methods on the shared helper fit (test_fit, fitted on
# test_data with p = 2 outcome and q = 3 selection coefficients).

p <- length(test_fit$coefficients$outcome)
q <- length(test_fit$coefficients$selection)

test_that("poisselect returns a classed object and print shows the info", {
  expect_s3_class(test_fit, "poisselect")
  expect_output(print(test_fit), "Observations")
  expect_output(print(test_fit), "selected")
  expect_output(print(test_fit), "Log-likelihood")
  expect_output(print(test_fit), "sigma")
  expect_output(print(test_fit), "rho")
})

test_that("summary builds correct coefficient tables", {
  fit_summary <- summary(test_fit)
  expect_s3_class(fit_summary, "summary.poisselect")
  expected_columns <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  expect_identical(dim(fit_summary$coefficients_outcome), c(p, 4L))
  expect_identical(
    colnames(fit_summary$coefficients_outcome),
    expected_columns
  )
  expect_identical(dim(fit_summary$coefficients_selection), c(q, 4L))
  expect_identical(
    colnames(fit_summary$coefficients_selection),
    expected_columns
  )
  tables <- rbind(
    fit_summary$coefficients_outcome,
    fit_summary$coefficients_selection
  )
  expect_true(all(tables[, "Pr(>|z|)"] >= 0 & tables[, "Pr(>|z|)"] <= 1))
  expect_equal(
    tables[, "z value"],
    tables[, "Estimate"] / tables[, "Std. Error"]
  )
  expect_output(print(fit_summary), "Outcome equation")
  expect_output(print(fit_summary), "Selection equation")
  expect_output(print(fit_summary), "AIC")
})

test_that("the AIC identity holds and AIC() works via logLik()", {
  expect_equal(test_fit$aic, -2 * test_fit$loglik + 2 * (p + q + 2))
  expect_equal(AIC(test_fit), test_fit$aic)
})

test_that("coef returns the full prefixed coefficient vector", {
  coefficients <- coef(test_fit)
  expect_length(coefficients, p + q + 2)
  expect_identical(
    names(coefficients),
    c(
      paste0("outcome:", names(test_fit$coefficients$outcome)),
      paste0("selection:", names(test_fit$coefficients$selection)),
      "sigma", "rho"
    )
  )
})

test_that("vcov is a symmetric matrix with the prefixed dimnames", {
  covariance <- vcov(test_fit)
  expect_identical(dim(covariance), c(p + q, p + q))
  expect_true(isSymmetric(covariance))
  expect_identical(
    rownames(covariance),
    names(coef(test_fit))[seq_len(p + q)]
  )
  expect_identical(colnames(covariance), rownames(covariance))
})

test_that("logLik carries the correct class and attributes", {
  loglik <- logLik(test_fit)
  expect_s3_class(loglik, "logLik")
  expect_equal(as.numeric(loglik), test_fit$loglik)
  expect_identical(attr(loglik, "df"), p + q + 2L)
  expect_identical(attr(loglik, "nobs"), test_fit$n)
})

test_that("predict computes all three types correctly", {
  beta <- test_fit$coefficients$outcome
  gamma <- test_fit$coefficients$selection
  sigma <- test_fit$coefficients$sigma
  link <- predict(test_fit, type = "link")
  expect_equal(link, drop(test_fit$model$x %*% beta))
  response <- predict(test_fit, type = "response")
  expect_equal(response, exp(link + sigma^2 / 2))
  # The default type must be "response" (assignment requirement).
  expect_equal(predict(test_fit), response)
  pselect <- predict(test_fit, type = "pselect")
  expect_equal(pselect, pnorm(drop(test_fit$model$z %*% gamma)))
  expect_true(all(pselect >= 0 & pselect <= 1))
})

test_that("predict works on newdata and validates it", {
  beta <- test_fit$coefficients$outcome
  newdata <- data.frame(x1 = c(-1, 0, 2), w1 = c(0.5, 0, -0.5))
  link <- predict(test_fit, newdata = newdata, type = "link")
  expect_equal(unname(link), beta[[1]] + beta[[2]] * newdata$x1)
  pselect <- predict(test_fit, newdata = newdata, type = "pselect")
  gamma <- test_fit$coefficients$selection
  expect_equal(
    unname(pselect),
    pnorm(gamma[[1]] + gamma[[2]] * newdata$x1 + gamma[[3]] * newdata$w1)
  )
  expect_error(
    predict(test_fit,
      newdata = data.frame(w1 = 1),
      type = "link"
    ),
    "x1"
  )
  expect_error(
    predict(test_fit,
      newdata = data.frame(x1 = NA_real_),
      type = "link"
    ),
    "NA"
  )
  expect_error(predict(test_fit, type = "nonsense"))
})

test_that("plot runs silently for valid `which` and errors otherwise", {
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(test_fit))
  expect_silent(plot(test_fit, which = 1))
  expect_silent(plot(test_fit, which = 2))
  expect_error(plot(test_fit, which = 3))
  expect_error(plot(test_fit, rho_grid = c(-1, 0, 1)), "between -1 and 1")
})
