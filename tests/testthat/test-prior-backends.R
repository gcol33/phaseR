prior_test_data <- function() {
  dat <- sim_phaseR(n_units = 40, n_times = 4, seed = 21)
  dat$site <- factor((dat$id - 1) %% 5)
  dat$year <- factor(dat$time)
  dat$count <- rpois(nrow(dat), 3)
  dat
}

two_phase_model <- function(dyn_0, dyn_1, family = gaussian()) {
  phase_model(
    phases("inactive", "active"),
    transition("inactive", "active", ~ x),
    dynamics("inactive", dyn_0, family = family),
    dynamics("active", dyn_1, family = family)
  )
}

prior_target_cases <- function() {
  list(
    two_phase = two_phase_model(y ~ x, y ~ x),
    glm = two_phase_model(count ~ x, count ~ 1, family = poisson()),
    re = two_phase_model(y ~ x + (1 | site), y ~ x),
    multi_re = two_phase_model(y ~ x + (1 | site) + (1 | year), y ~ x),
    slopes = two_phase_model(y ~ x + (1 + x | site), y ~ x),
    k_phase = phase_model(
      phases("A", "B", "C"),
      transition("A", "B", ~ 1),
      transition("B", "C", ~ 1),
      dynamics("A", y ~ 1),
      dynamics("B", y ~ x),
      dynamics("C", y ~ 1)
    )
  )
}

# Coefficients and log residual SDs lead every parameter vector.
fixed_block_sizes <- function(stan_data) {
  if (stan_data$k_phases > 2) {
    return(c(beta = sum(stan_data$n_trans_coef) + sum(stan_data$n_dyn_coef),
             sigma = stan_data$k_phases))
  }
  c(beta = stan_data$n_trans_coef + stan_data$n_dyn_coef_0 +
      stan_data$n_dyn_coef_1, sigma = 2)
}

test_that("every model family evaluates the prior it is given", {
  dat <- prior_test_data()
  prior_a <- prior()
  prior_b <- prior(beta ~ normal(1, 0.3), sigma ~ half_cauchy(2))

  for (kind in names(prior_target_cases())) {
    model <- prior_target_cases()[[kind]]
    sd_a <- phaseR:::prepare_stan_data(model, dat, prior_a)
    sd_b <- phaseR:::prepare_stan_data(model, dat, prior_b)
    target_a <- phaseR:::posterior_target(sd_a)
    target_b <- phaseR:::posterior_target(sd_b)

    set.seed(3)
    theta <- rnorm(sd_a$n_params, 0, 0.3)
    sizes <- fixed_block_sizes(sd_a)
    beta <- theta[seq_len(sizes[["beta"]])]
    sigma <- exp(theta[sizes[["beta"]] + seq_len(sizes[["sigma"]])])
    expected <- sum(dnorm(beta, 0, 2.5, log = TRUE) -
                      dnorm(beta, 1, 0.3, log = TRUE)) +
      sum(dnorm(sigma, 0, 1, log = TRUE) - dcauchy(sigma, 0, 2, log = TRUE))

    diff <- phaseR:::log_posterior(target_a, theta) -
      phaseR:::log_posterior(target_b, theta)
    expect_equal(diff, expected, tolerance = 1e-8, info = kind)
  }
})

test_that("the Laplace objective is the shared log posterior", {
  dat <- prior_test_data()
  model <- prior_target_cases()$re
  stan_data <- phaseR:::prepare_stan_data(
    model, dat, prior(beta ~ student_t(3, 0, 1), re ~ half_cauchy(1)))
  theta <- seq(-0.5, 0.5, length.out = stan_data$n_params)
  expect_equal(
    phaseR:::build_neg_log_posterior(stan_data)$fn(theta),
    -phaseR:::log_posterior(phaseR:::posterior_target(stan_data), theta)
  )
})

test_that("a strong prior moves NUTS and Laplace posteriors alike", {
  skip_on_cran()
  dat <- sim_phaseR(n_units = 80, n_times = 5, seed = 42)
  model <- two_phase_model(y ~ x, y ~ x)
  strong <- prior(beta ~ normal(3, 0.05))
  betas <- c("beta_trans_(Intercept)", "beta_trans_x",
             "beta_inactive_(Intercept)", "beta_inactive_x",
             "beta_active_(Intercept)", "beta_active_x")

  nuts_default <- fit_phaseR(model, data = dat, chains = 2, iter = 600,
                             warmup = 300, seed = 4)
  nuts_strong <- fit_phaseR(model, data = dat, chains = 2, iter = 600,
                            warmup = 300, seed = 4, prior = strong)
  lap_default <- fit_phaseR(model, data = dat, backend = "laplace", seed = 4)
  lap_strong <- fit_phaseR(model, data = dat, backend = "laplace", seed = 4,
                           prior = strong)

  mean_of <- function(fit) colMeans(fit$draws[, betas])
  sd_of <- function(fit) apply(fit$draws[, betas], 2, sd)

  # The default-prior posterior is far from the strong prior's centre ...
  expect_gt(max(abs(mean_of(nuts_default) - 3)), 2)
  # ... and the strong prior pulls both backends onto it.
  expect_true(all(abs(mean_of(nuts_strong) - 3) < 0.25))
  expect_true(all(abs(mean_of(lap_strong) - 3) < 0.25))
  expect_true(all(abs(mean_of(nuts_strong) - mean_of(lap_strong)) <
                    3 * sd_of(nuts_strong) + 0.01))
  shift_nuts <- mean_of(nuts_strong) - mean_of(nuts_default)
  shift_lap <- mean_of(lap_strong) - mean_of(lap_default)
  expect_true(all(sign(shift_nuts) == sign(shift_lap)))
  expect_lt(max(abs(shift_nuts - shift_lap)), 0.25)
  expect_equal(nuts_strong$diagnostics$n_divergent, 0L)
})
