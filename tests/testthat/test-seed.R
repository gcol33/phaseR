seed_test_setup <- function() {
  dat <- sim_phaseR(n_units = 20, n_times = 3, seed = 11)
  model <- phase_model(
    phases("inactive", "active"),
    transition("inactive", "active", ~ x),
    dynamics("inactive", y ~ x),
    dynamics("active", y ~ x)
  )
  list(dat = dat, model = model)
}

test_that("HMC fits with the same seed return identical draws", {
  skip_on_cran()
  s <- seed_test_setup()

  fit_a <- fit_phaseR(s$model, data = s$dat, chains = 2, iter = 60,
                      warmup = 30, seed = 2024)
  fit_b <- fit_phaseR(s$model, data = s$dat, chains = 2, iter = 60,
                      warmup = 30, seed = 2024)
  fit_c <- fit_phaseR(s$model, data = s$dat, chains = 2, iter = 60,
                      warmup = 30, seed = 2025)

  expect_identical(fit_a$draws, fit_b$draws)
  expect_identical(fit_a$seed, 2024L)
  expect_false(isTRUE(all.equal(fit_a$draws, fit_c$draws)))

  chain_1 <- fit_a$draws[fit_a$chain_id == 1, ]
  chain_2 <- fit_a$draws[fit_a$chain_id == 2, ]
  expect_false(isTRUE(all.equal(chain_1, chain_2)))
})

test_that("seed = NULL draws the seed from R's RNG", {
  skip_on_cran()
  s <- seed_test_setup()

  set.seed(99)
  fit_a <- fit_phaseR(s$model, data = s$dat, chains = 1, iter = 40,
                      warmup = 20)
  set.seed(99)
  fit_b <- fit_phaseR(s$model, data = s$dat, chains = 1, iter = 40,
                      warmup = 20)

  expect_identical(fit_a$draws, fit_b$draws)
  expect_identical(fit_a$seed, fit_b$seed)

  fit_c <- fit_phaseR(s$model, data = s$dat, chains = 1, iter = 40,
                      warmup = 20, seed = fit_a$seed)
  expect_identical(fit_a$draws, fit_c$draws)
})

test_that("Laplace fits with the same seed return identical draws", {
  s <- seed_test_setup()

  fit_a <- fit_phaseR(s$model, data = s$dat, backend = "laplace", seed = 7)
  fit_b <- fit_phaseR(s$model, data = s$dat, backend = "laplace", seed = 7)
  fit_c <- fit_phaseR(s$model, data = s$dat, backend = "laplace", seed = 8)

  expect_identical(fit_a$draws, fit_b$draws)
  expect_false(isTRUE(all.equal(fit_a$draws, fit_c$draws)))
})

test_that("seed must be a single integer", {
  s <- seed_test_setup()
  expect_error(
    fit_phaseR(s$model, data = s$dat, backend = "laplace", seed = 1.5),
    "single integer"
  )
  expect_error(
    fit_phaseR(s$model, data = s$dat, backend = "laplace", seed = c(1, 2)),
    "single integer"
  )
})
