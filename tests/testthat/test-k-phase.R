test_that("sim_phaseR_k generates valid data for 3 phases", {

  dat <- sim_phaseR_k(n_units = 20, n_times = 6, k = 3, seed = 123)

  expect_s3_class(dat, "data.frame")
  expect_true(all(c("id", "time", "y", "true_phase") %in% names(dat)))
  expect_equal(nrow(dat), 20 * 6)
  expect_true(all(dat$true_phase %in% 0:2))
})


test_that("sim_phaseR_k generates valid data for 4 phases", {

  dat <- sim_phaseR_k(n_units = 30, n_times = 10, k = 4, seed = 456)

  expect_equal(nrow(dat), 30 * 10)
  expect_true(all(dat$true_phase %in% 0:3))
})


test_that("phase_model accepts 3 phases with correct transitions", {

  model <- phase_model(
    phases("A", "B", "C"),
    transition("A", "B", ~ 1),
    transition("B", "C", ~ 1),
    dynamics("A", y ~ 1),
    dynamics("B", y ~ 1),
    dynamics("C", y ~ 1)
  )

  expect_s3_class(model, "phase_model")
  expect_equal(model$phases$n_phases, 3)
  expect_equal(length(model$transitions), 2)
  expect_equal(length(model$dynamics), 3)
})


test_that("phase_model rejects invalid transition order", {

  expect_error(
    phase_model(
      phases("A", "B", "C"),
      transition("B", "A", ~ 1),  # backwards
      transition("B", "C", ~ 1),
      dynamics("A", y ~ 1),
      dynamics("B", y ~ 1),
      dynamics("C", y ~ 1)
    ),
    "violates phase ordering"
  )
})


test_that("prepare_stan_data_k builds correct structure", {

  dat <- sim_phaseR_k(n_units = 10, n_times = 5, k = 3, seed = 789)

  model <- phase_model(
    phases("A", "B", "C"),
    transition("A", "B", ~ 1),
    transition("B", "C", ~ 1),
    dynamics("A", y ~ 1),
    dynamics("B", y ~ 1),
    dynamics("C", y ~ 1)
  )

  stan_data <- prepare_stan_data(model, dat)

  expect_equal(stan_data$k_phases, 3)
  expect_equal(length(stan_data$X_trans_list), 2)
  expect_equal(length(stan_data$X_dyn_list), 3)
  expect_equal(stan_data$n_units, 10)
})


test_that("fit_phaseR works with 3-phase model", {

  skip_on_cran()

  dat <- sim_phaseR_k(n_units = 20, n_times = 5, k = 3, seed = 111)

  model <- phase_model(
    phases("A", "B", "C"),
    transition("A", "B", ~ 1),
    transition("B", "C", ~ 1),
    dynamics("A", y ~ 1),
    dynamics("B", y ~ 1),
    dynamics("C", y ~ 1)
  )

  fit <- fit_phaseR(model, data = dat, chains = 1, iter = 50, warmup = 25, seed = 222)

  expect_s3_class(fit, "phase_fit")
  expect_equal(fit$k_phases, 3)
  expect_true(nrow(fit$draws) > 0)
})


test_that("k-phase param names are correctly structured", {

  dat <- sim_phaseR_k(n_units = 10, n_times = 4, k = 3, seed = 333)

  model <- phase_model(
    phases("early", "mid", "late"),
    transition("early", "mid", ~ 1),
    transition("mid", "late", ~ 1),
    dynamics("early", y ~ 1),
    dynamics("mid", y ~ 1),
    dynamics("late", y ~ 1)
  )

  stan_data <- prepare_stan_data(model, dat)

  # Check transition param names
  expect_true(any(grepl("beta_trans_early_to_mid", stan_data$param_names)))
  expect_true(any(grepl("beta_trans_mid_to_late", stan_data$param_names)))

  # Check dynamics param names
  expect_true(any(grepl("beta_early_", stan_data$param_names)))
  expect_true(any(grepl("beta_mid_", stan_data$param_names)))
  expect_true(any(grepl("beta_late_", stan_data$param_names)))

  # Check sigma names
  expect_true("log_sigma_early" %in% stan_data$param_names)
  expect_true("log_sigma_mid" %in% stan_data$param_names)
  expect_true("log_sigma_late" %in% stan_data$param_names)
})
