test_that("posterior_predict returns correct dimensions", {

  skip_on_cran()

  dat <- sim_phaseR(n_units = 20, n_times = 3, seed = 123)

  model <- phase_model(
    phases("inactive", "active"),
    transition("inactive", "active", ~ 1),
    dynamics("inactive", y ~ 1),
    dynamics("active", y ~ 1)
  )

  fit <- fit_phaseR(model, data = dat, chains = 1, iter = 50, warmup = 25, seed = 456)

  preds <- posterior_predict(fit, n_sims = 10)

  expect_equal(nrow(preds), nrow(dat))
  expect_equal(ncol(preds), 10)
})


test_that("posterior_predict type='phase' returns probabilities", {

  skip_on_cran()

  dat <- sim_phaseR(n_units = 15, n_times = 3, seed = 111)

  model <- phase_model(
    phases("inactive", "active"),
    transition("inactive", "active", ~ 1),
    dynamics("inactive", y ~ 1),
    dynamics("active", y ~ 1)
  )

  fit <- fit_phaseR(model, data = dat, chains = 1, iter = 50, warmup = 25, seed = 222)

  probs <- posterior_predict(fit, type = "phase", n_sims = 10)

  expect_equal(nrow(probs), nrow(dat))
  expect_true(all(probs >= 0 & probs <= 1, na.rm = TRUE))

  # First time point should always be phase 0
  first_rows <- which(dat$time == 1)
  expect_true(all(probs[first_rows, ] == 0))
})


test_that("counterfactual returns expected structure", {

  skip_on_cran()

  dat <- sim_phaseR(n_units = 20, n_times = 3, seed = 333)
  dat$treatment <- rbinom(nrow(dat), 1, 0.5)

  model <- phase_model(
    phases("inactive", "active"),
    transition("inactive", "active", ~ treatment),
    dynamics("inactive", y ~ 1),
    dynamics("active", y ~ 1)
  )

  fit <- fit_phaseR(model, data = dat, chains = 1, iter = 50, warmup = 25, seed = 444)

  cf_data <- dat
  cf_data$treatment <- 1
  cf <- counterfactual(fit, newdata = cf_data, n_sims = 10)

  expect_true("mean" %in% names(cf))
  expect_true("sd" %in% names(cf))
  expect_true("samples" %in% names(cf))
  expect_equal(length(cf$mean), nrow(dat))
})


test_that("ate computes treatment effect", {

  skip_on_cran()

  dat <- sim_phaseR(n_units = 30, n_times = 3, seed = 555)
  dat$treatment <- rbinom(nrow(dat), 1, 0.5)

  model <- phase_model(
    phases("inactive", "active"),
    transition("inactive", "active", ~ treatment),
    dynamics("inactive", y ~ 1),
    dynamics("active", y ~ 1)
  )

  fit <- fit_phaseR(model, data = dat, chains = 1, iter = 50, warmup = 25, seed = 666)

  effect <- ate(fit, treatment_var = "treatment", n_sims = 10)

  expect_true("ate" %in% names(effect))
  expect_true("ate_sd" %in% names(effect))
  expect_true("ate_ci" %in% names(effect))
  expect_equal(length(effect$ate_ci), 2)
})


test_that("ite returns unit-level effects", {

  skip_on_cran()

  dat <- sim_phaseR(n_units = 15, n_times = 3, seed = 777)
  dat$treatment <- rbinom(nrow(dat), 1, 0.5)

  model <- phase_model(
    phases("inactive", "active"),
    transition("inactive", "active", ~ treatment),
    dynamics("inactive", y ~ 1),
    dynamics("active", y ~ 1)
  )

  fit <- fit_phaseR(model, data = dat, chains = 1, iter = 50, warmup = 25, seed = 888)

  effects <- ite(fit, treatment_var = "treatment", n_sims = 10)

  expect_s3_class(effects, "data.frame")
  expect_equal(nrow(effects), 15)
  expect_true(all(c("id", "ite_mean", "ite_sd") %in% names(effects)))
})


test_that("ate errors on missing treatment variable", {

  skip_on_cran()

  dat <- sim_phaseR(n_units = 10, n_times = 3, seed = 999)

  model <- phase_model(
    phases("inactive", "active"),
    transition("inactive", "active", ~ 1),
    dynamics("inactive", y ~ 1),
    dynamics("active", y ~ 1)
  )

  fit <- fit_phaseR(model, data = dat, chains = 1, iter = 50, warmup = 25, seed = 101)

  expect_error(ate(fit, treatment_var = "nonexistent"), "not found")
})


test_that("log_sum_exp handles edge cases", {

  expect_equal(log_sum_exp(-Inf, 0), 0)
  expect_equal(log_sum_exp(0, -Inf), 0)
  expect_equal(log_sum_exp(-Inf, -Inf), -Inf)

  # Standard case
  result <- log_sum_exp(log(0.3), log(0.7))
  expect_equal(exp(result), 1.0, tolerance = 1e-10)
})
