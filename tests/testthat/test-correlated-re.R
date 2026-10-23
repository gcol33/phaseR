correlated_case <- function(prior_spec = prior(re ~ half_normal(0.8),
                                                cor ~ lkj(2.5))) {
  dat <- reference_data()
  model <- reference_two_phase(~ x, y ~ x + z + (1 + x + z | region), y ~ x)
  stan_data <- phaseR:::prepare_stan_data(model, dat, prior_spec)
  term <- phaseR:::posterior_target(stan_data)$dynamics[[1]]$re[[1]]
  list(stan_data = stan_data, target = phaseR:::posterior_target(stan_data),
       effects = term$offset + seq_len(term$n_groups * term$n_coefs),
       log_sd = term$offset + term$n_groups * term$n_coefs + 1:3,
       cpc = term$offset + term$n_groups * term$n_coefs + 4:6)
}

# Correlation matrix from the unconstrained values, written out for three
# coefficients.
correlation_3 <- function(y) {
  z <- tanh(y)
  L <- diag(3)
  L[2, 1] <- z[1]
  L[2, 2] <- sqrt(1 - z[1]^2)
  L[3, 1] <- z[2]
  L[3, 2] <- z[3] * sqrt(1 - z[2]^2)
  L[3, 3] <- sqrt(1 - L[3, 1]^2 - L[3, 2]^2)
  L %*% t(L)
}

test_that("correlated terms are named log SDs and atanh partial correlations", {
  cc <- correlated_case()
  names <- cc$stan_data$param_names
  expect_equal(names[cc$log_sd],
               paste0("log_sigma_inactive_re_region_", c("intercept", "x", "z")))
  expect_equal(names[cc$cpc],
               paste0("atanh_cpc_inactive_re_region_", c("2_1", "3_1", "3_2")))
  start <- phaseR:::initial_values(cc$target, 1)[[1]]
  expect_equal(start[cc$log_sd], rep(log(0.8), 3))
  expect_equal(start[cc$cpc], rep(0, 3))
})

test_that("the correlation prior is LKJ(eta) on the correlation matrix", {
  cc <- correlated_case()
  eta <- 2.5
  # With the standardised effects at zero the likelihood does not depend on
  # the correlations, so log posterior differences are prior differences:
  # (eta - 1) log det(Omega) plus the log Jacobian of y -> Omega.
  log_density <- function(y) {
    offdiag <- function(v) {
      O <- correlation_3(v)
      c(O[2, 1], O[3, 1], O[3, 2])
    }
    J <- vapply(1:3, function(m) {
      e <- replace(numeric(3), m, 1e-6)
      (offdiag(y + e) - offdiag(y - e)) / 2e-6
    }, numeric(3))
    (eta - 1) * log(det(correlation_3(y))) + log(abs(det(J)))
  }
  set.seed(6)
  theta <- rnorm(cc$target$n_params, 0, 0.3)
  theta[cc$effects] <- 0
  lp <- function(y) phaseR:::log_posterior(cc$target, replace(theta, cc$cpc, y))
  for (i in 1:4) {
    y_a <- rnorm(3, 0, 0.6)
    y_b <- rnorm(3, 0, 0.6)
    expect_equal(lp(y_a) - lp(y_b), log_density(y_a) - log_density(y_b),
                 tolerance = 1e-7)
  }
})

test_that("correlated coefficients are diag(sigma) L z and VarCorr reads them", {
  cc <- correlated_case()
  set.seed(7)
  theta <- rnorm(cc$target$n_params, 0, 0.4)
  fit <- structure(list(model = reference_two_phase(~ x, y ~ x + z +
                                                      (1 + x + z | region),
                                                    y ~ x),
                        data = reference_data(),
                        draws = matrix(theta, 1,
                                       dimnames = list(NULL, cc$stan_data$param_names)),
                        param_names = cc$stan_data$param_names, k_phases = 2L),
                   class = "phase_fit")
  sd <- exp(theta[cc$log_sd])
  Omega <- correlation_3(theta[cc$cpc])
  L <- t(chol(Omega))
  z <- matrix(theta[cc$effects], ncol = 3, byrow = TRUE)
  b <- z %*% t(L) %*% diag(sd)

  re <- ranef(fit, summary = FALSE)[["dynamics_inactive:region"]]
  expect_equal(re[1, , ], b, ignore_attr = TRUE, tolerance = 1e-12)
  vc <- VarCorr(fit)[["dynamics_inactive:region"]]
  expect_equal(vc$sd$mean, sd, tolerance = 1e-12)
  expect_equal(vc$cor, Omega, ignore_attr = TRUE, tolerance = 1e-12)
  summary_rows <- ranef(fit)
  expect_equal(nrow(summary_rows), 3 * 3)
})
