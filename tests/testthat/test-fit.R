test_that("sim_phaseR generates valid data", {

  dat <- sim_phaseR(n_units = 10, n_times = 4, seed = 123)

  expect_true(is.data.frame(dat))
  expect_true(all(c("id", "time", "x", "y") %in% names(dat)))
  expect_equal(nrow(dat), 10 * 4)
  expect_equal(length(unique(dat$id)), 10)
})

test_that("fit_phaseR runs on simulated data", {

  skip_on_cran()

  # Simulate data
  dat <- sim_phaseR(n_units = 30, n_times = 4, seed = 123)

  # Define model
  model <- phase_model(
    phases("inactive", "active"),
    transition("inactive", "active", ~ x),
    dynamics("inactive", y ~ x),
    dynamics("active", y ~ x)
  )

  # Fit with minimal iterations
  fit <- fit_phaseR(model, data = dat, chains = 1, iter = 50, warmup = 25, seed = 456)

  expect_s3_class(fit, "phase_fit")
  expect_equal(fit$n_chains, 1)
  expect_true(nrow(fit$draws) > 0)
})

test_that("summary.phase_fit produces valid output", {

  skip_on_cran()

  dat <- sim_phaseR(n_units = 20, n_times = 3, seed = 456)
  model <- phase_model(
    phases("inactive", "active"),
    transition("inactive", "active", ~ 1),
    dynamics("inactive", y ~ 1),
    dynamics("active", y ~ 1)
  )
  fit <- fit_phaseR(model, data = dat, chains = 1, iter = 50, warmup = 25, seed = 789)

  s <- summary(fit)
  expect_s3_class(s, "summary.phase_fit")
  expect_true("mean" %in% colnames(s$summaries))
})

test_that("predict.phase_fit returns predictions", {

  skip_on_cran()

  dat <- sim_phaseR(n_units = 20, n_times = 3, seed = 789)
  model <- phase_model(
    phases("inactive", "active"),
    transition("inactive", "active", ~ 1),
    dynamics("inactive", y ~ 1),
    dynamics("active", y ~ 1)
  )
  fit <- fit_phaseR(model, data = dat, chains = 1, iter = 50, warmup = 25, seed = 111)

  preds <- predict(fit, type = "phase")
  expect_true(is.data.frame(preds))
  expect_equal(nrow(preds), nrow(dat))
})

test_that("coef.phase_fit returns coefficients", {

  skip_on_cran()

  dat <- sim_phaseR(n_units = 20, n_times = 3, seed = 222)
  model <- phase_model(
    phases("inactive", "active"),
    transition("inactive", "active", ~ 1),
    dynamics("inactive", y ~ 1),
    dynamics("active", y ~ 1)
  )
  fit <- fit_phaseR(model, data = dat, chains = 1, iter = 50, warmup = 25, seed = 333)

  coefs <- coef(fit)
  expect_true(is.numeric(coefs))
  expect_true(length(coefs) > 0)
})
