# Tests for LOO-CV (v1.0.6)

# =============================================================================
# PSIS Helper Functions
# =============================================================================

test_that("estimate_pareto_k returns reasonable values", {
  set.seed(123)

  # Generate data from Pareto-like distribution
  x <- rexp(100) * 2

  k <- phaseR:::estimate_pareto_k(x)
  expect_true(is.finite(k))
  expect_true(k >= -0.5 && k <= 1.5)
})

test_that("estimate_pareto_k handles small samples", {
  x <- c(1, 2, 3)
  k <- phaseR:::estimate_pareto_k(x)
  expect_true(is.na(k))
})

test_that("smooth_pareto_tail returns correct length", {
  set.seed(456)
  x <- sort(rexp(20), decreasing = TRUE)
  k <- 0.3

  smoothed <- phaseR:::smooth_pareto_tail(x, k, 20)
  expect_length(smoothed, 20)
  expect_true(all(smoothed >= 0))
})

test_that("psis_smooth returns log_weights and k", {
  set.seed(789)
  log_ratios <- rnorm(100, 0, 1)

  result <- phaseR:::psis_smooth(log_ratios)

  expect_true("log_weights" %in% names(result))
  expect_true("k" %in% names(result))
  expect_length(result$log_weights, 100)
  expect_true(is.finite(result$k))
})

test_that("psis_loo computes estimates for each observation", {
  set.seed(101)

  # Create a mock log-likelihood matrix
  n_draws <- 100
  n_obs <- 20
  ll_matrix <- matrix(rnorm(n_draws * n_obs, -2, 0.5), nrow = n_draws)

  result <- phaseR:::psis_loo(ll_matrix)

  expect_length(result$elpd_loo, n_obs)
  expect_length(result$p_loo, n_obs)
  expect_length(result$pareto_k, n_obs)

  expect_true(all(is.finite(result$elpd_loo)))
  expect_true(all(is.finite(result$pareto_k)))
})


# =============================================================================
# loo() Function
# =============================================================================

test_that("loo() returns loo_phaseR object", {
  skip_on_cran()

  set.seed(111)
  n_units <- 15
  n_times <- 4

  data <- expand.grid(id = 1:n_units, time = 1:n_times)
  data$x <- rnorm(nrow(data))
  data$y <- 5 + 0.5 * data$x + rnorm(nrow(data), 0, 0.5)

  model <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ x),
    dynamics("pre", y ~ x),
    dynamics("post", y ~ x)
  )

  fit <- fit_phaseR(model, data, backend = "laplace", n_samples = 200, seed = 222)

  result <- loo(fit)

  expect_s3_class(result, "loo_phaseR")
  expect_true("elpd_loo" %in% names(result))
  expect_true("p_loo" %in% names(result))
  expect_true("looic" %in% names(result))
  expect_true("pointwise" %in% names(result))
  expect_true("diagnostics" %in% names(result))
})

test_that("loo() pointwise has correct dimensions", {
  skip_on_cran()

  set.seed(333)
  n_units <- 10
  n_times <- 4
  n_obs <- n_units * n_times

  data <- expand.grid(id = 1:n_units, time = 1:n_times)
  data$x <- rnorm(nrow(data))
  data$y <- 5 + 0.5 * data$x + rnorm(nrow(data), 0, 0.5)

  model <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ 1),
    dynamics("pre", y ~ x),
    dynamics("post", y ~ x)
  )

  fit <- fit_phaseR(model, data, backend = "laplace", n_samples = 200, seed = 444)

  result <- loo(fit)

  expect_equal(nrow(result$pointwise), n_obs)
  expect_true("elpd_loo" %in% names(result$pointwise))
  expect_true("p_loo" %in% names(result$pointwise))
  expect_true("pareto_k" %in% names(result$pointwise))
})

test_that("loo() diagnostics count k values correctly", {
  skip_on_cran()

  set.seed(555)
  n_units <- 12
  n_times <- 4

  data <- expand.grid(id = 1:n_units, time = 1:n_times)
  data$x <- rnorm(nrow(data))
  data$y <- 5 + 0.5 * data$x + rnorm(nrow(data), 0, 0.5)

  model <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ 1),
    dynamics("pre", y ~ x),
    dynamics("post", y ~ x)
  )

  fit <- fit_phaseR(model, data, backend = "laplace", n_samples = 200, seed = 666)

  result <- loo(fit)

  # Check that counts add up
  total <- result$diagnostics$n_good + result$diagnostics$n_okay + result$diagnostics$n_bad
  expect_equal(total, nrow(result$pointwise))

  # Check counts match pointwise values
  k_vals <- result$pointwise$pareto_k
  expect_equal(result$diagnostics$n_bad, sum(k_vals >= 0.7))
  expect_equal(result$diagnostics$n_okay, sum(k_vals >= 0.5 & k_vals < 0.7))
  expect_equal(result$diagnostics$n_good, sum(k_vals < 0.5))
})

test_that("loo() rejects non phase_fit objects", {
  expect_error(loo("not a fit"), "must be a phase_fit")
  expect_error(loo(list(a = 1)), "must be a phase_fit")
})

test_that("print.loo_phaseR works", {
  skip_on_cran()

  set.seed(777)
  data <- expand.grid(id = 1:10, time = 1:4)
  data$x <- rnorm(nrow(data))
  data$y <- 5 + 0.5 * data$x + rnorm(nrow(data), 0, 0.5)

  model <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ 1),
    dynamics("pre", y ~ 1),
    dynamics("post", y ~ 1)
  )

  fit <- fit_phaseR(model, data, backend = "laplace", n_samples = 200, seed = 888)

  result <- loo(fit)

  expect_output(print(result), "Leave-One-Out Cross-Validation")
  expect_output(print(result), "elpd_loo")
  expect_output(print(result), "Pareto k diagnostic")
})


# =============================================================================
# loo_compare()
# =============================================================================

test_that("loo_compare() compares two models", {
  skip_on_cran()

  set.seed(999)
  n_units <- 15
  n_times <- 4

  data <- expand.grid(id = 1:n_units, time = 1:n_times)
  data$x <- rnorm(nrow(data))
  data$y <- 5 + 0.5 * data$x + rnorm(nrow(data), 0, 0.5)

  model1 <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ 1),
    dynamics("pre", y ~ 1),
    dynamics("post", y ~ 1)
  )

  model2 <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ x),
    dynamics("pre", y ~ x),
    dynamics("post", y ~ x)
  )

  fit1 <- fit_phaseR(model1, data, backend = "laplace", n_samples = 200, seed = 1010)
  fit2 <- fit_phaseR(model2, data, backend = "laplace", n_samples = 200, seed = 1010)

  result <- loo_compare(fit1, fit2, model_names = c("Null", "Full"))

  expect_s3_class(result, "loo_compare_phaseR")
  expect_equal(nrow(result), 2)
  expect_true("model" %in% names(result))
  expect_true("elpd_loo" %in% names(result))
  expect_true("elpd_diff" %in% names(result))
  expect_true("se_diff" %in% names(result))

  # First model should have elpd_diff = 0
  expect_equal(result$elpd_diff[1], 0)
  expect_equal(result$se_diff[1], 0)
})

test_that("loo_compare() accepts pre-computed loo objects", {
  skip_on_cran()

  set.seed(1111)
  data <- expand.grid(id = 1:12, time = 1:4)
  data$x <- rnorm(nrow(data))
  data$y <- 5 + 0.5 * data$x + rnorm(nrow(data), 0, 0.5)

  model1 <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ 1),
    dynamics("pre", y ~ 1),
    dynamics("post", y ~ 1)
  )

  model2 <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ 1),
    dynamics("pre", y ~ x),
    dynamics("post", y ~ x)
  )

  fit1 <- fit_phaseR(model1, data, backend = "laplace", n_samples = 200, seed = 1212)
  fit2 <- fit_phaseR(model2, data, backend = "laplace", n_samples = 200, seed = 1212)

  loo1 <- loo(fit1)
  loo2 <- loo(fit2)

  result <- loo_compare(loo1, loo2)

  expect_s3_class(result, "loo_compare_phaseR")
  expect_equal(nrow(result), 2)
})

test_that("loo_compare() requires at least 2 models", {
  skip_on_cran()

  set.seed(1313)
  data <- expand.grid(id = 1:10, time = 1:4)
  data$x <- rnorm(nrow(data))
  data$y <- rnorm(nrow(data))

  model <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ 1),
    dynamics("pre", y ~ 1),
    dynamics("post", y ~ 1)
  )

  fit <- fit_phaseR(model, data, backend = "laplace", n_samples = 200, seed = 1414)

  expect_error(loo_compare(fit), "at least 2 models")
})

test_that("print.loo_compare_phaseR works", {
  skip_on_cran()

  set.seed(1515)
  data <- expand.grid(id = 1:12, time = 1:4)
  data$x <- rnorm(nrow(data))
  data$y <- 5 + 0.5 * data$x + rnorm(nrow(data), 0, 0.5)

  model1 <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ 1),
    dynamics("pre", y ~ 1),
    dynamics("post", y ~ 1)
  )

  model2 <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ 1),
    dynamics("pre", y ~ x),
    dynamics("post", y ~ x)
  )

  fit1 <- fit_phaseR(model1, data, backend = "laplace", n_samples = 200, seed = 1616)
  fit2 <- fit_phaseR(model2, data, backend = "laplace", n_samples = 200, seed = 1616)

  result <- loo_compare(fit1, fit2)

  expect_output(print(result), "Model Comparison")
  expect_output(print(result), "elpd_loo")
})


# =============================================================================
# LOO vs WAIC comparison
# =============================================================================

test_that("LOO and WAIC give similar rankings", {
  skip_on_cran()

  set.seed(1717)
  n_units <- 20
  n_times <- 5

  data <- expand.grid(id = 1:n_units, time = 1:n_times)
  data$x <- rnorm(nrow(data))
  # y depends on x, so model with x should be better
  data$y <- 5 + 0.8 * data$x + rnorm(nrow(data), 0, 0.3)

  model1 <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ 1),
    dynamics("pre", y ~ 1),
    dynamics("post", y ~ 1)
  )

  model2 <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ 1),
    dynamics("pre", y ~ x),
    dynamics("post", y ~ x)
  )

  fit1 <- fit_phaseR(model1, data, backend = "laplace", n_samples = 300, seed = 1818)
  fit2 <- fit_phaseR(model2, data, backend = "laplace", n_samples = 300, seed = 1818)

  loo1 <- loo(fit1)
  loo2 <- loo(fit2)

  waic1 <- waic(fit1)
  waic2 <- waic(fit2)

  # Both should prefer model2 (with x predictor)
  # LOO: higher elpd is better
  # WAIC: lower waic is better
  expect_gt(loo2$elpd_loo, loo1$elpd_loo)
  expect_lt(waic2$waic, waic1$waic)
})


# =============================================================================
# Edge cases
# =============================================================================

test_that("loo() handles small number of draws", {
  skip_on_cran()

  set.seed(1919)
  data <- expand.grid(id = 1:8, time = 1:4)
  data$x <- rnorm(nrow(data))
  data$y <- 5 + 0.5 * data$x + rnorm(nrow(data), 0, 0.5)

  model <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ 1),
    dynamics("pre", y ~ 1),
    dynamics("post", y ~ 1)
  )

  # Small number of samples
  fit <- fit_phaseR(model, data, backend = "laplace", n_samples = 50, seed = 2020)

  result <- loo(fit)

  expect_s3_class(result, "loo_phaseR")
  expect_true(is.finite(result$elpd_loo))
})
