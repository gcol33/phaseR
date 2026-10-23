fake_fit <- function(model, data, draws) {
  stan_data <- phaseR:::prepare_stan_data(model, data, prior())
  draws <- matrix(draws, ncol = stan_data$n_params,
                  dimnames = list(NULL, stan_data$param_names))
  structure(list(model = model, data = data, draws = draws,
                 param_names = stan_data$param_names,
                 k_phases = stan_data$k_phases),
            class = "phase_fit")
}

two_phase_draws <- function(n) {
  set.seed(2)
  cbind(rnorm(n, -1, 0.2), rnorm(n, 0.5, 0.2), rnorm(n, 10, 0.2),
        rnorm(n, 0, 0.2), rnorm(n, 8, 0.2), rnorm(n, -0.5, 0.2),
        rnorm(n, 0, 0.1), rnorm(n, 0.4, 0.1))
}

test_that("two-phase Gaussian predictions agree with the direct recursions", {
  dat <- sim_phaseR(n_units = 10, n_times = 5, seed = 4)
  dat <- dat[sample(nrow(dat)), ]
  draws <- two_phase_draws(4)
  fit <- fake_fit(reference_two_phase(~ x, y ~ x, y ~ x), dat, draws)

  filtered <- posterior_predict(fit, type = "phase")
  moves <- predict(fit, type = "transition", summary = FALSE)
  for (d in 1:4) {
    th <- draws[d, ]
    p <- plogis(th[1] + th[2] * dat$x)
    expect_equal(moves[d, dat$time > 1], p[dat$time > 1], tolerance = 1e-12)
    for (unit in unique(dat$id)) {
      rows <- which(dat$id == unit)
      rows <- rows[order(dat$time[rows])]
      e0 <- dnorm(dat$y[rows], th[3] + th[4] * dat$x[rows], exp(th[7]))
      e1 <- dnorm(dat$y[rows], th[5] + th[6] * dat$x[rows], exp(th[8]))
      a <- c(e0[1], 0)
      prob_1 <- 0
      for (t in seq_along(rows)[-1]) {
        a <- c(a[1] * (1 - p[rows[t]]) * e0[t],
               (a[1] * p[rows[t]] + a[2]) * e1[t])
        a <- a / sum(a)
        prob_1 <- c(prob_1, a[2])
      }
      expect_equal(filtered[rows, d], prob_1, tolerance = 1e-10)
    }
  }

  # With one draw, phase probabilities given the covariates are the
  # transition-only recursion and the expected outcome mixes the phase means.
  fit_1 <- fake_fit(fit$model, dat, draws[1, ])
  th <- draws[1, ]
  phase <- predict(fit_1, type = "phase")
  response <- predict(fit_1, type = "response", summary = FALSE)
  p <- plogis(th[1] + th[2] * dat$x)
  for (unit in unique(dat$id)) {
    rows <- which(dat$id == unit)
    rows <- rows[order(dat$time[rows])]
    stay <- cumprod(c(1, 1 - p[rows[-1]]))
    expect_equal(phase$inactive[rows], stay, tolerance = 1e-12)
    expect_equal(phase$active[rows], 1 - stay, tolerance = 1e-12)
    expect_equal(response[1, rows],
                 stay * (th[3] + th[4] * dat$x[rows]) +
                   (1 - stay) * (th[5] + th[6] * dat$x[rows]),
                 tolerance = 1e-12)
  }
})

# Direct simulation of the generative model of a two-phase model with a
# Poisson first phase, a Gaussian second phase and a site random intercept
# in the first phase.
simulate_truth <- function(dat, th, u_site, reps) {
  out <- matrix(0, reps, nrow(dat))
  for (unit in unique(dat$id)) {
    rows <- which(dat$id == unit)
    rows <- rows[order(dat$time[rows])]
    phase <- rep(0, reps)
    for (t in seq_along(rows)) {
      r <- rows[t]
      if (t > 1) {
        phase <- pmax(phase, stats::runif(reps) < plogis(th[1] + th[2] * dat$x[r]))
      }
      mu_0 <- exp(th[3] + th[4] * dat$x[r] + u_site[as.integer(dat$site[r])])
      mu_1 <- th[5] + th[6] * dat$x[r]
      out[, r] <- ifelse(phase == 0, rpois(reps, mu_0),
                         rnorm(reps, mu_1, exp(th[8])))
    }
  }
  out
}

test_that("predictions of a Poisson random-effect model match direct simulation", {
  dat <- reference_data()
  dat <- dat[dat$id <= 12, ]
  model <- reference_two_phase(~ x, count ~ x + (1 | site), count ~ x,
                               family_0 = poisson(), family_1 = gaussian())
  set.seed(8)
  u_site <- rnorm(8, 0, 0.4)
  th <- c(-0.8, 0.6, 1.2, 0.3, 4, -0.5, 0, -0.3)
  fit <- fake_fit(model, dat, c(th, u_site, log(0.4)))

  set.seed(9)
  truth <- simulate_truth(dat, th, u_site, 40000)
  expected <- predict(fit, type = "response", summary = FALSE)[1, ]
  mc_se <- apply(truth, 2, sd) / sqrt(nrow(truth))
  expect_true(all(abs(colMeans(truth) - expected) < 4 * mc_se))

  set.seed(10)
  sims <- posterior_predict(fit, n_sims = 1)
  expect_equal(dim(sims), c(nrow(dat), 1L))
  big <- fake_fit(model, dat, matrix(c(th, u_site, log(0.4)), 20000, 17,
                                     byrow = TRUE))
  sims <- posterior_predict(big)
  expect_true(all(abs(rowMeans(sims) - expected) < 4 * mc_se * sqrt(2)))

  # Treatment effect of x = 1 against x = 0 on the expected outcome.
  effect <- ate(fit, treatment_var = "x")
  set.seed(11)
  treat <- simulate_truth(transform(dat, x = 1), th, u_site, 40000)
  control <- simulate_truth(transform(dat, x = 0), th, u_site, 40000)
  true_ate <- mean(treat) - mean(control)
  ate_se <- sqrt(var(rowMeans(treat)) + var(rowMeans(control))) /
    sqrt(40000)
  expect_lt(abs(effect$ate - true_ate), 4 * ate_se)
})

test_that("predictions cover k-phase models and unseen random-effect levels", {
  dat <- sim_phaseR_k(n_units = 8, n_times = 6, k = 3, seed = 3)
  model <- phase_model(
    phases("A", "B", "C"),
    transition("A", "B", ~ x), transition("B", "C", ~ 1),
    dynamics("A", y ~ 1), dynamics("B", y ~ x), dynamics("C", y ~ 1)
  )
  th <- c(-0.5, 0.4, -1, 12, 10, -0.5, 8, 0, 0, 0)
  fit <- fake_fit(model, dat, th)
  phase <- predict(fit, type = "phase")
  expect_equal(names(phase), c("A", "B", "C"))
  expect_equal(rowSums(phase), rep(1, nrow(dat)), tolerance = 1e-12)
  moves <- predict(fit, type = "transition", summary = FALSE)
  expect_named(moves, c("A_to_B", "B_to_C"))
  expect_equal(moves$B_to_C[1, dat$time > 1],
               rep(plogis(-1), sum(dat$time > 1)), tolerance = 1e-12)
  filtered_c <- posterior_predict(fit, type = "phase", phase = "C")
  expect_true(all(filtered_c[dat$time <= 2, ] == 0))

  re_dat <- reference_data()
  re_model <- reference_two_phase(~ x, y ~ x + (1 | site), y ~ x)
  stan_data <- phaseR:::prepare_stan_data(re_model, re_dat, prior())
  theta <- seq(-0.3, 0.3, length.out = stan_data$n_params)
  re_fit <- fake_fit(re_model, re_dat, theta)
  new <- re_dat[re_dat$id == 1, ]
  new$site <- factor("unseen")
  no_effect <- predict(re_fit, newdata = new, summary = FALSE)
  site_free <- fake_fit(reference_two_phase(~ x, y ~ x, y ~ x), re_dat,
                        theta[1:8])
  expect_equal(no_effect, predict(site_free, newdata = new, summary = FALSE))
})

test_that("fitted values are one-step predictive means and logLik is the marginal likelihood", {
  dat <- sim_phaseR(n_units = 10, n_times = 5, seed = 4)
  draws <- two_phase_draws(1)
  fit <- fake_fit(reference_two_phase(~ x, y ~ x, y ~ x), dat, draws)
  target <- phaseR:::posterior_target(
    phaseR:::prepare_stan_data(fit$model, dat, prior()))
  expect_equal(as.numeric(logLik(fit)),
               sum(phaseR:::compute_pointwise_ll(fit)), tolerance = 1e-12)

  th <- draws[1, ]
  p <- plogis(th[1] + th[2] * dat$x)
  mu_0 <- th[3] + th[4] * dat$x
  mu_1 <- th[5] + th[6] * dat$x
  rows <- which(dat$id == 1)
  e0 <- dnorm(dat$y[rows], mu_0[rows], exp(th[7]))
  e1 <- dnorm(dat$y[rows], mu_1[rows], exp(th[8]))
  a <- c(1, 0)
  one_step <- mu_0[rows[1]]
  for (t in 2:5) {
    a <- a * c(e0[t - 1], e1[t - 1])
    a <- a / sum(a)
    pred <- c(a[1] * (1 - p[rows[t]]), a[1] * p[rows[t]] + a[2])
    one_step <- c(one_step, sum(pred * c(mu_0[rows[t]], mu_1[rows[t]])))
    a <- pred
  }
  expect_equal(fitted(fit)[rows], one_step, tolerance = 1e-10)
  expect_equal(residuals(fit)[rows], dat$y[rows] - one_step, tolerance = 1e-10)
})
