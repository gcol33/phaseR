laplace_at <- function(target) {
  laplace <- phaseR:::phaseR_laplace(target)
  layout <- phaseR:::phaseR_laplace_layout(laplace)
  list(laplace = laplace, layout = layout,
       value = function(theta) {
         phaseR:::phaseR_laplace_objective(laplace, theta)$value
       })
}

test_that("the Laplace marginal gradient matches finite differences", {
  gc <- gradient_cases()
  for (name in names(gc$cases)) {
    target <- reference_target(gc$cases[[name]], gc$dat)
    lap <- laplace_at(target)
    # Near the data-informed start, where the effects given theta have one
    # mode.
    set.seed(1)
    theta <- phaseR:::initial_values(target, 1)[[1]][lap$layout$theta] +
      rnorm(length(lap$layout$theta), 0, 0.05)
    exact <- phaseR:::phaseR_laplace_objective(lap$laplace, theta)
    expect_true(exact$converged, info = name)
    numeric <- richardson_gradient(lap$value, theta, h = 1e-2, levels = 2)
    rel_error <- max(abs(exact$gradient - numeric)) / max(1, abs(numeric))
    expect_lt(rel_error, 1e-6, label = name)
  }
})

test_that("without random effects the Laplace objective is the log posterior", {
  case <- reference_cases()$two_phase
  target <- reference_target(case)
  lap <- laplace_at(target)
  expect_length(lap$layout$u, 0)
  theta <- seq(-0.4, 0.4, length.out = target$n_params)
  expect_equal(phaseR:::phaseR_laplace_objective(lap$laplace, theta)[1:2],
               phaseR:::log_posterior(target, theta, gradient = TRUE),
               tolerance = 1e-12)
})

test_that("the joint posterior is unbounded as a random-effect SD shrinks; the marginal is not", {
  case <- reference_cases()$re
  target <- reference_target(case)
  lap <- laplace_at(target)
  theta <- phaseR:::initial_values(target, 1)[[1]]
  term <- target$dynamics[[1]]$re[[1]]
  effects <- term$offset + seq_len(term$n_groups)
  log_sd <- term$offset + term$n_groups + 1L

  log_sds <- c(-2, -6, -10, -14)
  joint <- vapply(log_sds, function(s) {
    phaseR:::log_posterior(target, replace(replace(theta, effects, 0),
                                           log_sd, s))
  }, numeric(1))
  marginal <- vapply(log_sds, function(s) {
    lap$value(replace(theta, log_sd, s)[lap$layout$theta])
  }, numeric(1))
  # Joint: each Normal(0, sd) density of a zero effect adds -log(sd), the
  # prior on log(sd) only +log(sd) - sd^2 / (2 s^2).
  expect_true(all(diff(joint) > 0))
  expect_gt(joint[4] - joint[1], 50)
  expect_true(all(diff(marginal) < 0))
})

test_that("Laplace matches NUTS on a random-intercept model", {
  skip_on_cran()
  set.seed(10)
  dat <- sim_phaseR(n_units = 60, n_times = 5, seed = 31)
  dat$site <- factor((dat$id - 1) %% 10)
  dat$y <- dat$y + rnorm(10, 0, 0.7)[as.integer(dat$site)]
  model <- phase_model(
    phases("inactive", "active"),
    transition("inactive", "active", ~ x),
    dynamics("inactive", y ~ x + (1 | site)),
    dynamics("active", y ~ x + (1 | site))
  )

  expect_no_warning(lap <- fit_phaseR(model, dat, backend = "laplace",
                                      seed = 2))
  nuts <- fit_phaseR(model, dat, chains = 2, iter = 1500, seed = 2)

  z <- (colMeans(lap$draws) - colMeans(nuts$draws)) / apply(nuts$draws, 2, sd)
  ratio <- apply(lap$draws, 2, sd) / apply(nuts$draws, 2, sd)
  sd_params <- grepl("_re_site$", names(z))
  # The posterior of a log random-effect SD is skewed, which a Gaussian
  # matches less closely.
  expect_true(all(abs(z[!sd_params]) < 0.5),
              info = paste(round(z, 2), collapse = " "))
  expect_true(all(abs(z[sd_params]) < 0.75),
              info = paste(round(z[sd_params], 2), collapse = " "))
  expect_true(all(ratio[!sd_params] > 0.8 & ratio[!sd_params] < 1.25),
              info = paste(round(ratio, 2), collapse = " "))
  expect_true(all(ratio[sd_params] > 0.6 & ratio[sd_params] < 1.5),
              info = paste(round(ratio[sd_params], 2), collapse = " "))
  expect_setequal(lap$diagnostics$marginalised,
                  grep("^v_", lap$param_names, value = TRUE))
})
