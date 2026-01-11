test_that("parse_formula_re detects random effects", {

  parsed <- parse_formula_re(~ x + (1 | site))
  expect_true(parsed$has_re)
  expect_equal(parsed$re_groups, "site")
  expect_equal(deparse(parsed$fixed_formula), "~x")
})

test_that("parse_formula_re handles no random effects", {

  parsed <- parse_formula_re(~ x + z)
  expect_false(parsed$has_re)
  expect_equal(parsed$re_groups, character(0))
})

test_that("parse_formula_re handles intercept-only after removing RE", {

  parsed <- parse_formula_re(~ (1 | site))
  expect_true(parsed$has_re)
  expect_equal(parsed$re_groups, "site")
  expect_equal(deparse(parsed$fixed_formula), "~1")
})

test_that("build_re_index creates correct structure", {

  data <- data.frame(
    id = 1:10,
    site = rep(c("A", "B", "C"), length.out = 10)
  )

  idx <- build_re_index(data, "site")
  expect_equal(idx$n_groups, 3)
  expect_equal(length(idx$idx), 10)
  expect_equal(idx$levels, c("A", "B", "C"))
})

test_that("build_re_index errors on missing variable", {

  data <- data.frame(id = 1:5)
  expect_error(build_re_index(data, "site"), "not found")
})

test_that("fit_phaseR works with random intercepts in transition", {

  skip_on_cran()

  # Simulate data with site
  dat <- sim_phaseR(n_units = 30, n_times = 3, seed = 123)
  dat$site <- sample(c("A", "B", "C"), nrow(dat), replace = TRUE)

  model <- phase_model(
    phases("inactive", "active"),
    transition("inactive", "active", ~ (1 | site)),
    dynamics("inactive", y ~ 1),
    dynamics("active", y ~ 1)
  )

  fit <- fit_phaseR(model, data = dat, chains = 1, iter = 50, warmup = 25, seed = 456)

  expect_s3_class(fit, "phase_fit")
  expect_true(fit$has_re)
  expect_true(fit$has_trans_re)
})

test_that("ranef returns NULL for models without RE", {

  skip_on_cran()

  dat <- sim_phaseR(n_units = 20, n_times = 3, seed = 789)

  model <- phase_model(
    phases("inactive", "active"),
    transition("inactive", "active", ~ 1),
    dynamics("inactive", y ~ 1),
    dynamics("active", y ~ 1)
  )

  fit <- fit_phaseR(model, data = dat, chains = 1, iter = 50, warmup = 25, seed = 111)

  expect_message(ranef(fit), "no random effects")
})
