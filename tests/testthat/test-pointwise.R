# Normalised forward filter of a two-phase Gaussian model: the log
# predictive density of each observation given the earlier ones of its unit.
filter_log_lik <- function(y, x, theta) {
  n <- length(y)
  mu_0 <- theta[3] + theta[4] * x
  mu_1 <- theta[5] + theta[6] * x
  p <- plogis(theta[1] + theta[2] * x)
  e_0 <- dnorm(y, mu_0, exp(theta[7]))
  e_1 <- dnorm(y, mu_1, exp(theta[8]))
  ll <- numeric(n)
  a <- c(1, 0)
  ll[1] <- log(e_0[1])
  for (t in seq_len(n)[-1]) {
    new <- c(a[1] * (1 - p[t]) * e_0[t], (a[1] * p[t] + a[2]) * e_1[t])
    ll[t] <- log(sum(new))
    a <- new / sum(new)
  }
  ll
}

test_that("pointwise log likelihoods follow the forward filter in data order", {
  dat <- sim_phaseR(n_units = 12, n_times = 5, seed = 4)
  dat <- dat[sample(nrow(dat)), ]
  model <- reference_two_phase(~ x, y ~ x, y ~ x)
  set.seed(2)
  draws <- cbind(rnorm(3, -1, 0.2), rnorm(3, 0.5, 0.2), rnorm(3, 10, 0.2),
                 rnorm(3, 0, 0.2), rnorm(3, 8, 0.2), rnorm(3, -0.5, 0.2),
                 rnorm(3, 0, 0.1), rnorm(3, 0.4, 0.1))
  fit <- list(model = model, data = dat, draws = draws)

  ll <- phaseR:::compute_pointwise_ll(fit)
  expect_equal(dim(ll), c(3L, nrow(dat)))
  for (d in 1:3) {
    for (unit in unique(dat$id)) {
      rows <- which(dat$id == unit)
      rows <- rows[order(dat$time[rows])]
      expect_equal(ll[d, rows],
                   filter_log_lik(dat$y[rows], dat$x[rows], draws[d, ]),
                   tolerance = 1e-10)
    }
  }
})

test_that("pointwise log likelihoods use the drawn random effects", {
  dat <- reference_data()
  model <- reference_two_phase(~ x, y ~ x + (1 | site), y ~ x)
  stan_data <- phaseR:::prepare_stan_data(model, dat, prior())
  set.seed(3)
  theta <- rnorm(stan_data$n_params, 0, 0.3)
  shifted <- theta
  site_3 <- which(stan_data$param_names == "v_inactive_site_3")
  shifted[site_3] <- theta[site_3] + 1
  fit <- list(model = model, data = dat, draws = rbind(theta, shifted))

  ll <- phaseR:::compute_pointwise_ll(fit)
  changed <- abs(ll[2, ] - ll[1, ]) > 1e-12
  # Only units of site 3 change; their first observation, certainly in the
  # inactive phase, always does.
  expect_false(any(changed[dat$site != "3"]))
  expect_true(all(changed[dat$site == "3" & dat$time == 1]))
})
