run_mvnormal <- function(mu, Sigma, chains = 4, iter = 2000, warmup = 1000,
                         seed = 1L, adapt_delta = 0.8, max_treedepth = 10L) {
  target <- list(
    kind = "mvnormal",
    n_params = length(mu),
    data = list(mu = mu, precision = chol2inv(chol(Sigma))),
    args = list()
  )
  runs <- lapply(seq_len(chains), function(chain) {
    phaseR:::phaseR_nuts_sampler(target, as.integer(iter), as.integer(warmup),
                                 as.integer(seed), as.integer(chain),
                                 adapt_delta, as.integer(max_treedepth))
  })
  draws <- do.call(rbind, lapply(runs, `[[`, "draws"))
  chain_id <- rep(seq_len(chains), each = iter - warmup)
  list(draws = draws, chain_id = chain_id, runs = runs,
       ess = phaseR:::compute_ess(draws, chain_id))
}

expect_moments_recovered <- function(res, mu, Sigma) {
  sds <- sqrt(diag(Sigma))
  est_mean <- colMeans(res$draws)
  est_sd <- apply(res$draws, 2, sd)
  mcse_mean <- sds / sqrt(res$ess)
  centered_sq <- sweep(res$draws, 2, colMeans(res$draws))^2
  ess_sq <- phaseR:::compute_ess(centered_sq, res$chain_id)
  mcse_sd <- sds / sqrt(2 * ess_sq)
  expect_true(all(abs(est_mean - mu) < 4 * mcse_mean),
              info = paste(round((est_mean - mu) / mcse_mean, 2), collapse = " "))
  expect_true(all(abs(est_sd - sds) < 4 * mcse_sd + 0.02 * sds),
              info = paste(round(est_sd / sds, 3), collapse = " "))
}

test_that("NUTS recovers a correlated bivariate normal", {
  skip_on_cran()
  mu <- c(1, -2)
  sds <- c(1, 3)
  R <- matrix(c(1, 0.95, 0.95, 1), 2)
  Sigma <- diag(sds) %*% R %*% diag(sds)

  res <- run_mvnormal(mu, Sigma, seed = 11L)

  expect_equal(sum(vapply(res$runs, `[[`, integer(1), "n_divergent")), 0L)
  expect_true(all(res$ess > 400))
  expect_moments_recovered(res, mu, Sigma)
  expect_lt(abs(cor(res$draws)[1, 2] - 0.95), 0.02)
})

test_that("NUTS adapts the diagonal metric to badly scaled targets", {
  skip_on_cran()
  mu <- c(0, 5, -1)
  sds <- c(10, 1, 0.1)
  R <- matrix(c(1, 0.5, 0, 0.5, 1, -0.3, 0, -0.3, 1), 3)
  Sigma <- diag(sds) %*% R %*% diag(sds)

  res <- run_mvnormal(mu, Sigma, seed = 12L)

  expect_equal(sum(vapply(res$runs, `[[`, integer(1), "n_divergent")), 0L)
  expect_true(all(res$ess > 400))
  expect_moments_recovered(res, mu, Sigma)

  for (run in res$runs) {
    ratio <- run$inv_metric / sds^2
    expect_true(all(ratio > 0.6 & ratio < 1.6),
                info = paste(round(ratio, 2), collapse = " "))
    expect_gt(run$accept_prob, 0.65)
  }
})

test_that("NUTS recovers a 20-dimensional standard normal", {
  skip_on_cran()
  d <- 20
  res <- run_mvnormal(rep(0, d), diag(d), seed = 13L)

  expect_equal(sum(vapply(res$runs, `[[`, integer(1), "n_divergent")), 0L)
  expect_true(all(res$ess > 1000))
  expect_moments_recovered(res, rep(0, d), diag(d))
  treedepth <- unlist(lapply(res$runs, `[[`, "treedepth"))
  expect_true(all(treedepth >= 2 & treedepth <= 5))
})

test_that("NUTS respects max_treedepth and reports its diagnostics", {
  skip_on_cran()
  res <- run_mvnormal(c(0, 0), diag(c(1, 100)), chains = 1, iter = 300,
                      warmup = 150, seed = 14L, max_treedepth = 2L)
  run <- res$runs[[1]]
  expect_true(all(run$treedepth <= 2))
  expect_true(all(run$n_leapfrog <= 2^2 - 1))
  expect_length(run$accept_stat, 150)
  expect_length(run$energy, 150)
  expect_true(all(run$accept_stat >= 0 & run$accept_stat <= 1))
})

test_that("NUTS streams are reproducible per seed and chain", {
  skip_on_cran()
  a <- run_mvnormal(c(0, 0), diag(2), chains = 2, iter = 200, warmup = 100, seed = 5L)
  b <- run_mvnormal(c(0, 0), diag(2), chains = 2, iter = 200, warmup = 100, seed = 5L)
  c <- run_mvnormal(c(0, 0), diag(2), chains = 2, iter = 200, warmup = 100, seed = 6L)
  expect_identical(a$draws, b$draws)
  expect_false(isTRUE(all.equal(a$draws, c$draws)))
  expect_false(isTRUE(all.equal(a$runs[[1]]$draws, a$runs[[2]]$draws)))
})

test_that("NUTS rejects invalid settings", {
  target <- list(kind = "mvnormal", n_params = 1L,
                 data = list(mu = 0, precision = matrix(1)), args = list())
  expect_error(phaseR:::phaseR_nuts_sampler(target, 10L, 10L, 1L, 1L, 0.8, 10L),
               "warmup")
  expect_error(phaseR:::phaseR_nuts_sampler(target, 10L, 5L, 1L, 1L, 1, 10L),
               "adapt_delta")
  expect_error(phaseR:::phaseR_nuts_sampler(target, 10L, 5L, 1L, 1L, 0.8, 0L),
               "max_treedepth")
})

test_that("compute_ess matches the variance of independent draws", {
  set.seed(1)
  x <- matrix(rnorm(4000), ncol = 1)
  ess <- phaseR:::compute_ess(x, rep(1:4, each = 1000))
  expect_gt(ess, 3400)
  expect_lt(ess, 4600)

  ar <- as.numeric(stats::filter(rnorm(4000), 0.9, method = "recursive"))
  ess_ar <- phaseR:::compute_ess(matrix(ar, ncol = 1), rep(1:4, each = 1000))
  expect_lt(ess_ar, 4000 * (1 - 0.9) / (1 + 0.9) * 2)
})
