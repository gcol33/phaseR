sim_phase_designs <- function(f0, f1, beta_0, beta_1, beta_trans = -0.5,
                              sigma = c(0.5, 0.5), n_units = 150,
                              n_times = 5, seed = 1) {
  set.seed(seed)
  dat <- expand.grid(time = seq_len(n_times), id = seq_len(n_units))
  dat <- dat[, c("id", "time")]
  dat$x <- rnorm(nrow(dat))
  dat$z <- rnorm(nrow(dat))
  mu_0 <- drop(model.matrix(f0, dat) %*% beta_0)
  mu_1 <- drop(model.matrix(f1, dat) %*% beta_1)
  p_trans <- plogis(beta_trans)

  dat$phase <- 0L
  for (i in seq_len(n_units)) {
    rows <- which(dat$id == i)
    for (k in seq_along(rows)[-1]) {
      prev <- dat$phase[rows[k - 1]]
      dat$phase[rows[k]] <- if (prev == 1L) 1L else rbinom(1, 1, p_trans)
    }
  }
  dat$y <- ifelse(dat$phase == 0L,
                  rnorm(nrow(dat), mu_0, sigma[1]),
                  rnorm(nrow(dat), mu_1, sigma[2]))
  dat
}

expect_dynamics_recovered <- function(fit, truth) {
  draws <- fit$draws[, names(truth), drop = FALSE]
  z <- (colMeans(draws) - truth) / apply(draws, 2, sd)
  expect_true(all(abs(z) < 4), info = paste(names(truth), round(z, 2),
                                            collapse = ", "))
}

phase_design_cases <- list(
  different_column_counts = list(
    f0 = ~ 1, f1 = ~ x + z, beta_0 = 2, beta_1 = c(5, 1.5, -1),
    truth = c("beta_inactive_(Intercept)" = 2,
              "beta_active_(Intercept)" = 5,
              "beta_active_x" = 1.5,
              "beta_active_z" = -1)
  ),
  same_count_different_variables = list(
    f0 = ~ x, f1 = ~ z, beta_0 = c(2, 1), beta_1 = c(5, -1.5),
    truth = c("beta_inactive_(Intercept)" = 2,
              "beta_inactive_x" = 1,
              "beta_active_(Intercept)" = 5,
              "beta_active_z" = -1.5)
  )
)

for (case_name in names(phase_design_cases)) {
  case <- phase_design_cases[[case_name]]

  test_that(paste("each phase uses its own design matrix:", case_name), {
    skip_on_cran()
    dat <- sim_phase_designs(case$f0, case$f1, case$beta_0, case$beta_1)
    model <- phase_model(
      phases("inactive", "active"),
      transition("inactive", "active", ~ 1),
      dynamics("inactive", update(case$f0, y ~ .)),
      dynamics("active", update(case$f1, y ~ .))
    )

    expect_no_warning(
      fit_laplace <- fit_phaseR(model, data = dat, backend = "laplace",
                                seed = 1)
    )
    expect_dynamics_recovered(fit_laplace, case$truth)

    expect_no_warning(
      fit_nuts <- fit_phaseR(model, data = dat, chains = 2, iter = 500,
                             warmup = 250, seed = 1)
    )
    expect_equal(fit_nuts$diagnostics$n_divergent, 0L)
    expect_dynamics_recovered(fit_nuts, case$truth)
  })
}
