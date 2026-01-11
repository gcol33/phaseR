test_that("encode_family correctly encodes families", {

  expect_equal(encode_family(gaussian()), 0L)
  expect_equal(encode_family(poisson()), 1L)
  expect_equal(encode_family(binomial()), 2L)
  expect_equal(encode_family(NULL), 0L)
})


test_that("sim_phaseR_poisson generates count data", {

  dat <- sim_phaseR_poisson(n_units = 20, n_times = 4, seed = 123)

  expect_s3_class(dat, "data.frame")
  expect_true(all(dat$y >= 0))
  expect_true(all(dat$y == floor(dat$y)))  # integers
})


test_that("prepare_stan_data detects Poisson family", {

  dat <- sim_phaseR_poisson(n_units = 10, n_times = 3, seed = 456)

  model <- phase_model(
    phases("inactive", "active"),
    transition("inactive", "active", ~ 1),
    dynamics("inactive", y ~ 1, family = poisson()),
    dynamics("active", y ~ 1, family = poisson())
  )

  stan_data <- prepare_stan_data(model, dat)

  expect_true(stan_data$is_glm)
  expect_equal(stan_data$family_0, 1L)
  expect_equal(stan_data$family_1, 1L)
})


test_that("prepare_stan_data handles mixed families", {

  dat <- sim_phaseR(n_units = 10, n_times = 3, seed = 789)

  model <- phase_model(
    phases("inactive", "active"),
    transition("inactive", "active", ~ 1),
    dynamics("inactive", y ~ 1, family = gaussian()),
    dynamics("active", y ~ 1, family = gaussian())
  )

  stan_data <- prepare_stan_data(model, dat)

  expect_false(stan_data$is_glm)
  expect_equal(stan_data$family_0, 0L)
  expect_equal(stan_data$family_1, 0L)
})


test_that("fit_phaseR works with Poisson family", {

  skip_on_cran()

  dat <- sim_phaseR_poisson(n_units = 25, n_times = 4, seed = 111)

  model <- phase_model(
    phases("inactive", "active"),
    transition("inactive", "active", ~ 1),
    dynamics("inactive", y ~ 1, family = poisson()),
    dynamics("active", y ~ 1, family = poisson())
  )

  fit <- fit_phaseR(model, data = dat, chains = 1, iter = 50, warmup = 25, seed = 222)

  expect_s3_class(fit, "phase_fit")
  expect_true(nrow(fit$draws) > 0)
})


test_that("Poisson model recovers reasonable parameters", {

  skip_on_cran()

  # Simulate with known parameters
  dat <- sim_phaseR_poisson(
    n_units = 40,
    n_times = 5,
    beta_trans = c(-0.5),
    beta_0 = c(2),    # log(lambda) ~ 7.4
    beta_1 = c(1.5),  # log(lambda) ~ 4.5
    seed = 333
  )

  model <- phase_model(
    phases("inactive", "active"),
    transition("inactive", "active", ~ 1),
    dynamics("inactive", y ~ 1, family = poisson()),
    dynamics("active", y ~ 1, family = poisson())
  )

  fit <- fit_phaseR(model, data = dat, chains = 1, iter = 100, warmup = 50, seed = 444)

  # Get posterior means
  beta_inactive <- mean(fit$draws[, "beta_inactive_(Intercept)"])
  beta_active <- mean(fit$draws[, "beta_active_(Intercept)"])

  # Should be in reasonable range (not testing exact recovery with short chain)
  expect_true(beta_inactive > 1 && beta_inactive < 3)
  expect_true(beta_active > 0.5 && beta_active < 2.5)
})
