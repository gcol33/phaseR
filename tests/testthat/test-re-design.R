# A random-effect term and the fixed-effect design with one coefficient per
# group level give the same likelihood when the coefficients equal the
# effects. Changing the effects changes both log posteriors by the same
# likelihood difference plus their own (closed-form) prior differences.

likelihood_change <- function(model, dat, theta, theta_new, prior_change) {
  stan_data <- phaseR:::prepare_stan_data(model, dat, prior())
  target <- phaseR:::posterior_target(stan_data)
  phaseR:::log_posterior(target, theta_new) -
    phaseR:::log_posterior(target, theta) - prior_change
}

re_vs_fixed_case <- function(re_formula, fixed_formula, fixed_coef, dat,
                             n_levels) {
  model_re <- reference_two_phase(~ x, re_formula, y ~ x)
  model_fixed <- reference_two_phase(~ x, fixed_formula, y ~ x)
  sd_re <- phaseR:::prepare_stan_data(model_re, dat, prior())
  sd_fixed <- phaseR:::prepare_stan_data(model_fixed, dat, prior())

  # Shared leading parameters: transition (2), then phase-0 coefficients.
  set.seed(9)
  trans <- c(-1, 0.4)
  dyn_1 <- c(8, -0.5)
  log_sigma <- c(0, 0.4)
  effects <- rnorm(n_levels, 0, 0.6)
  effects_new <- effects + rnorm(n_levels, 0, 0.3)
  log_sd_re <- log(0.6)

  # Random-effect model: phase-0 intercept 10, then the term's effects;
  # fixed_coef() gives the fixed-effect coefficients with the same
  # predictor.
  theta_re <- function(u) c(trans, 10, dyn_1, log_sigma, u, log_sd_re)
  theta_fixed <- function(u) c(trans, fixed_coef(u), dyn_1, log_sigma)
  expect_equal(sd_re$n_params, length(theta_re(effects)))
  expect_equal(sd_fixed$n_params, length(theta_fixed(effects)))

  beta_sd <- prior()$beta$sd
  change_re <- likelihood_change(
    model_re, dat, theta_re(effects), theta_re(effects_new),
    sum(dnorm(effects_new, 0, 0.6, log = TRUE) -
          dnorm(effects, 0, 0.6, log = TRUE)))
  change_fixed <- likelihood_change(
    model_fixed, dat, theta_fixed(effects), theta_fixed(effects_new),
    sum(dnorm(fixed_coef(effects_new), 0, beta_sd, log = TRUE)) -
      sum(dnorm(fixed_coef(effects), 0, beta_sd, log = TRUE)))
  c(re = change_re, fixed = change_fixed)
}

test_that("a grouping factor that varies within units indexes each observation", {
  dat <- reference_data()
  dat$year <- factor(dat$time)
  change <- re_vs_fixed_case(y ~ 1 + (1 | year), y ~ 0 + year,
                             function(u) 10 + u, dat, nlevels(dat$year))
  expect_equal(change[["re"]], change[["fixed"]], tolerance = 1e-10)
  expect_gt(abs(change[["re"]]), 1)
})

test_that("a random slope without intercept multiplies its covariate", {
  dat <- reference_data()
  change <- re_vs_fixed_case(y ~ 1 + (0 + x | site), y ~ 1 + site:x,
                             function(u) c(10, u), dat, nlevels(dat$site))
  expect_equal(change[["re"]], change[["fixed"]], tolerance = 1e-10)
  expect_gt(abs(change[["re"]]), 1)
})

test_that("Cholesky factor names follow the order the likelihood reads", {
  dat <- reference_data()
  model <- reference_two_phase(~ x, y ~ x + z + (1 + x + z | region), y ~ x)
  stan_data <- phaseR:::prepare_stan_data(model, dat, prior(re ~ half_normal(0.8)))
  L_names <- grep("^L_", stan_data$param_names, value = TRUE)
  index <- substring(L_names, nchar(L_names) - 2)
  expect_equal(index, c("1_1", "2_1", "3_1", "2_2", "3_2", "3_3"))
  # initial_values() sets the diagonal of the factor, located by the C++
  # packing, to the scale of the random-effect prior.
  start <- phaseR:::initial_values(phaseR:::posterior_target(stan_data), 1)[[1]]
  names(start) <- stan_data$param_names
  diagonal <- substring(index, 1, 1) == substring(index, 3, 3)
  expect_equal(unname(start[L_names[diagonal]]), rep(0.8, 3))
  expect_equal(unname(start[L_names[!diagonal]]), rep(0, 3))
})
