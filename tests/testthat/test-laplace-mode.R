readme_model <- function() {
  phase_model(
    phases("inactive", "active"),
    transition("inactive", "active", ~ x),
    dynamics("inactive", y ~ x),
    dynamics("active", y ~ x)
  )
}

test_that("Laplace finds the mode NUTS samples around when a zero start does not", {
  skip_on_cran()
  dat <- sim_phaseR(n_units = 60, n_times = 4, seed = 8)
  model <- readme_model()

  # From all-zero parameters the optimiser stops in a mode in which no unit
  # enters the active phase.
  objective <- phaseR:::build_neg_log_posterior(
    phaseR:::prepare_stan_data(model, dat, prior()))
  from_zero <- optim(rep(0, 8), objective$fn, objective$gr, method = "L-BFGS-B")

  lap <- fit_phaseR(model, dat, backend = "laplace", seed = 3)
  nuts <- fit_phaseR(model, dat, chains = 4, iter = 1000, seed = 3)
  expect_lt(-from_zero$value, max(lap$diagnostics$start_values) - 20)

  nuts_mean <- colMeans(nuts$draws)
  nuts_sd <- apply(nuts$draws, 2, sd)
  lap_sd <- lap$diagnostics$se
  expect_true(all(abs(lap$diagnostics$mode - nuts_mean) < 0.5 * nuts_sd),
              info = paste(round((lap$diagnostics$mode - nuts_mean) / nuts_sd, 2),
                           collapse = " "))
  expect_true(all(lap_sd / nuts_sd > 0.7 & lap_sd / nuts_sd < 1.3),
              info = paste(round(lap_sd / nuts_sd, 2), collapse = " "))
  expect_gt(lap$diagnostics$mode[["beta_active_(Intercept)"]], 7.5)
})

test_that("every NUTS chain starts in the basin of the main mode", {
  skip_on_cran()
  dat <- sim_phaseR(n_units = 60, n_times = 4, seed = 8)
  nuts <- fit_phaseR(readme_model(), dat, chains = 4, iter = 600, seed = 3)
  active <- split(nuts$draws[, "beta_active_(Intercept)"], nuts$chain_id)
  expect_true(all(vapply(active, mean, numeric(1)) > 7.5))
})

test_that("Laplace fits are reproducible from the seed", {
  dat <- sim_phaseR(n_units = 30, n_times = 4, seed = 2)
  fit_a <- fit_phaseR(readme_model(), dat, backend = "laplace", seed = 7)
  fit_b <- fit_phaseR(readme_model(), dat, backend = "laplace", seed = 7)
  expect_identical(fit_a$draws, fit_b$draws)
  expect_identical(fit_a$diagnostics$start_values,
                   fit_b$diagnostics$start_values)
})

test_that("starting paths begin in the first phase and move one phase at a time", {
  dat <- sim_phaseR_k(n_units = 20, n_times = 6, k = 4, seed = 1)
  model <- phase_model(
    phases("A", "B", "C", "D"),
    transition("A", "B", ~ 1), transition("B", "C", ~ 1),
    transition("C", "D", ~ 1),
    dynamics("A", y ~ 1), dynamics("B", y ~ 1), dynamics("C", y ~ 1),
    dynamics("D", y ~ 1)
  )
  target <- phaseR:::posterior_target(phaseR:::prepare_stan_data(model, dat))
  set.seed(4)
  for (random in c(FALSE, TRUE, TRUE, TRUE)) {
    phase <- phaseR:::phase_paths(target, random)
    for (i in seq_along(target$unit_start)) {
      path <- phase[(target$unit_start[i]:target$unit_end[i]) + 1L]
      expect_equal(path[1], 0L)
      expect_true(all(diff(path) %in% c(0L, 1L)))
    }
  }
  split <- phaseR:::phase_paths(target, random = FALSE)
  expect_equal(split[1:6], c(0L, 0L, 1L, 2L, 2L, 3L))
})

test_that("ridge_glm recovers the coefficients of each family", {
  set.seed(11)
  X <- cbind(1, rnorm(2000))
  beta <- c(0.5, -0.8)
  eta <- drop(X %*% beta)
  w <- rep(1, 2000)
  expect_equal(phaseR:::ridge_glm(X, eta + rnorm(2000, 0, 0.1), w, "gaussian"),
               beta, tolerance = 0.02)
  expect_equal(phaseR:::ridge_glm(X, rpois(2000, exp(eta)), w, "poisson"),
               beta, tolerance = 0.1)
  expect_equal(phaseR:::ridge_glm(X, rbinom(2000, 1, plogis(eta)), w, "binomial"),
               beta, tolerance = 0.15)
})
