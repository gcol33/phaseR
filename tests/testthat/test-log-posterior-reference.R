test_that("log posteriors match the reference values of every model kind", {
  ref <- readRDS(test_path("fixtures", "log_posterior_reference.rds"))
  dat <- reference_data()
  cases <- reference_cases()
  for (name in names(ref)) {
    target <- reference_target(cases[[name]], dat)
    value <- apply(ref[[name]]$params, 1, function(p) {
      phaseR:::phaseR_log_posterior(target, p)
    })
    expect_equal(value, ref[[name]]$value, tolerance = 1e-10, info = name)
  }
})

test_that("uncorrelated random slopes get one Normal(0, sigma_j) per coefficient", {
  dat <- reference_data()
  case <- reference_cases()$slopes_uncorrelated
  stan_data <- phaseR:::prepare_stan_data(case$model, dat, case$prior)
  target <- phaseR:::posterior_target(stan_data)
  term <- target$dynamics[[1]]$re[[1]]
  n_effects <- term$n_groups * term$n_coefs
  effects <- term$offset + seq_len(n_effects)
  log_sd <- term$offset + n_effects + seq_len(term$n_coefs)

  set.seed(5)
  theta <- rnorm(target$n_params, 0, 0.4)
  theta_2 <- theta
  theta_2[log_sd] <- theta[log_sd] + c(0.3, -0.2)

  u <- matrix(theta[effects], ncol = term$n_coefs, byrow = TRUE)
  re_prior <- function(s) {
    sum(dnorm(u, 0, rep(exp(s), each = nrow(u)), log = TRUE)) +
      sum(dnorm(exp(s), 0, 0.5, log = TRUE) + log(2) + s)
  }
  expect_equal(
    phaseR:::phaseR_log_posterior(target, theta_2) -
      phaseR:::phaseR_log_posterior(target, theta),
    re_prior(theta_2[log_sd]) - re_prior(theta[log_sd]),
    tolerance = 1e-10
  )
})
