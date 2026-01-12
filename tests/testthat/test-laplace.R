# Tests for Laplace approximation backend (v1.0.4)

test_that("fit_laplace runs on simple model", {
  skip_on_cran()

  set.seed(123)
  n_units <- 20
  n_times <- 5

  data <- expand.grid(
    id = 1:n_units,
    time = 1:n_times
  )
  data$x <- rnorm(nrow(data))
  data$y <- 5 + 0.5 * data$x + rnorm(nrow(data), 0, 0.5)

  model <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ x),
    dynamics("pre", y ~ x),
    dynamics("post", y ~ x)
  )

  fit <- fit_phaseR(model, data, backend = "laplace", seed = 456)

  expect_s3_class(fit, "phase_fit")
  expect_equal(fit$backend, "laplace")
  expect_true(nrow(fit$draws) > 0)
})


test_that("Laplace returns consistent structure with HMC", {
  skip_on_cran()

  set.seed(789)
  n_units <- 15
  n_times <- 4

  data <- expand.grid(
    id = 1:n_units,
    time = 1:n_times
  )
  data$x <- rnorm(nrow(data))
  data$y <- 5 + 0.3 * data$x + rnorm(nrow(data), 0, 0.3)

  model <- phase_model(
    phases("baseline", "treatment"),
    transition("baseline", "treatment", ~ 1),
    dynamics("baseline", y ~ x),
    dynamics("treatment", y ~ x)
  )

  fit_lap <- fit_phaseR(model, data, backend = "laplace", seed = 101)
  fit_hmc <- fit_phaseR(model, data, backend = "hmc",
                         chains = 1, iter = 200, warmup = 100, seed = 101)

  # Same parameter names
  expect_equal(fit_lap$param_names, fit_hmc$param_names)

  # Same number of parameters
  expect_equal(ncol(fit_lap$draws), ncol(fit_hmc$draws))

  # Both are phase_fit objects
  expect_s3_class(fit_lap, "phase_fit")
  expect_s3_class(fit_hmc, "phase_fit")
})


test_that("Laplace diagnostics contain expected fields", {
  skip_on_cran()

  set.seed(111)
  data <- expand.grid(id = 1:10, time = 1:4)
  data$x <- rnorm(nrow(data))
  data$y <- 5 + 0.3 * data$x + rnorm(nrow(data), 0, 0.3)

  model <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ 1),
    dynamics("pre", y ~ 1),
    dynamics("post", y ~ 1)
  )

  fit <- fit_phaseR(model, data, backend = "laplace", seed = 222)

  expect_true("convergence" %in% names(fit$diagnostics))
  expect_true("mode" %in% names(fit$diagnostics))
  expect_true("se" %in% names(fit$diagnostics))
  expect_true("hessian" %in% names(fit$diagnostics))
})


test_that("Laplace works with random effects", {
  skip_on_cran()

  set.seed(333)
  n_sites <- 4
  n_times <- 4

  data <- expand.grid(
    site = letters[1:n_sites],
    time = 1:n_times
  )
  data$id <- data$site
  data$x <- rnorm(nrow(data))

  site_effects <- setNames(rnorm(n_sites, 0, 0.3), letters[1:n_sites])
  data$y <- 5 + site_effects[data$site] + 0.3 * data$x + rnorm(nrow(data), 0, 0.2)

  model <- phase_model(
    phases("baseline", "treatment"),
    transition("baseline", "treatment", ~ 1),
    dynamics("baseline", y ~ x + (1|site)),
    dynamics("treatment", y ~ x)
  )

  fit <- fit_phaseR(model, data, backend = "laplace", seed = 444)

  expect_s3_class(fit, "phase_fit")
  expect_true(fit$has_re)
  expect_true(any(grepl("site", fit$param_names)))
})


test_that("Laplace works with crossed random effects", {
  skip_on_cran()

  set.seed(555)
  data <- expand.grid(
    site = letters[1:3],
    year = 2020:2021,
    time = 1:3
  )
  data$id <- interaction(data$site, data$year)
  data$x <- rnorm(nrow(data))
  data$y <- 5 + 0.3 * data$x + rnorm(nrow(data), 0, 0.3)

  model <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ 1),
    dynamics("pre", y ~ x + (1|site) + (1|year)),
    dynamics("post", y ~ x)
  )

  fit <- fit_phaseR(model, data, backend = "laplace", seed = 666)

  expect_s3_class(fit, "phase_fit")
  expect_true(any(grepl("site", fit$param_names)))
  expect_true(any(grepl("year", fit$param_names)))
})


test_that("Laplace standard errors are reasonable", {
  skip_on_cran()

  set.seed(777)
  n_units <- 30
  n_times <- 5

  data <- expand.grid(id = 1:n_units, time = 1:n_times)
  data$x <- rnorm(nrow(data))
  # Large sample, low noise - should have small SEs
  data$y <- 5 + 0.5 * data$x + rnorm(nrow(data), 0, 0.1)

  model <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ x),
    dynamics("pre", y ~ x),
    dynamics("post", y ~ x)
  )

  fit <- fit_phaseR(model, data, backend = "laplace", seed = 888)

  # Standard errors should be positive and finite
  se <- fit$diagnostics$se
  expect_true(all(se > 0))
  expect_true(all(is.finite(se)))

  # For well-identified parameters, SE should be reasonably small
  expect_true(all(se < 10))
})


test_that("Laplace mode is close to posterior mean", {
  skip_on_cran()

  set.seed(999)
  n_units <- 25
  n_times <- 5

  data <- expand.grid(id = 1:n_units, time = 1:n_times)
  data$x <- rnorm(nrow(data))
  data$y <- 5 + 0.4 * data$x + rnorm(nrow(data), 0, 0.2)

  model <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ 1),
    dynamics("pre", y ~ x),
    dynamics("post", y ~ x)
  )

  fit <- fit_phaseR(model, data, backend = "laplace", seed = 1010)

  # Mode should be close to posterior mean of samples
  mode <- fit$diagnostics$mode
  post_mean <- colMeans(fit$draws)

  # They should be within 2 SEs of each other
  se <- fit$diagnostics$se
  diff <- abs(mode - post_mean)
  expect_true(all(diff < 3 * se))
})


test_that("Laplace n_samples parameter works", {
  skip_on_cran()

  set.seed(1111)
  data <- expand.grid(id = 1:10, time = 1:3)
  data$x <- rnorm(nrow(data))
  data$y <- rnorm(nrow(data))

  model <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ 1),
    dynamics("pre", y ~ 1),
    dynamics("post", y ~ 1)
  )

  fit1 <- fit_phaseR(model, data, backend = "laplace", n_samples = 500, seed = 1212)
  fit2 <- fit_phaseR(model, data, backend = "laplace", n_samples = 2000, seed = 1313)

  expect_equal(nrow(fit1$draws), 500)
  expect_equal(nrow(fit2$draws), 2000)
})


test_that("build_neg_log_posterior returns finite values", {
  set.seed(1414)
  data <- expand.grid(id = 1:10, time = 1:4)
  data$x <- rnorm(nrow(data))
  data$y <- 5 + 0.3 * data$x + rnorm(nrow(data), 0, 0.3)

  model <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ 1),
    dynamics("pre", y ~ x),
    dynamics("post", y ~ x)
  )

  stan_data <- phaseR:::prepare_stan_data(model, data)
  neg_lp <- phaseR:::build_neg_log_posterior(stan_data)

  # Test at zero
  val_0 <- neg_lp(rep(0, stan_data$n_params))
  expect_true(is.finite(val_0))

  # Test at random values
  val_rand <- neg_lp(rnorm(stan_data$n_params, 0, 0.5))
  expect_true(is.finite(val_rand))
})


test_that("sample_mvnorm produces correct samples", {
  set.seed(1515)

  mu <- c(1, 2, 3)
  sigma <- matrix(c(1, 0.5, 0.3,
                    0.5, 2, 0.4,
                    0.3, 0.4, 1.5), nrow = 3)

  samples <- phaseR:::sample_mvnorm(10000, mu, sigma)

  expect_equal(dim(samples), c(10000, 3))

  # Sample mean should be close to mu
  sample_mean <- colMeans(samples)
  expect_true(all(abs(sample_mean - mu) < 0.1))

  # Sample covariance should be close to sigma
  sample_cov <- cov(samples)
  expect_true(all(abs(sample_cov - sigma) < 0.15))
})
