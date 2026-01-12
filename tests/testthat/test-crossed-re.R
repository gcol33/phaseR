# Tests for crossed random effects (v1.0.1)

test_that("formula parser extracts multiple RE terms", {
  f <- y ~ x + (1|site) + (1|year)
  parsed <- parse_formula_re(f)

  expect_true(parsed$has_re)
  expect_equal(parsed$n_re_terms, 2)
  expect_equal(names(parsed$re_terms), c("site", "year"))
  expect_equal(parsed$re_groups, c("site", "year"))

  # Fixed formula should not contain RE terms
  expect_equal(deparse(parsed$fixed_formula), "y ~ x")
})


test_that("formula parser handles single RE term",
  {
  f <- y ~ x + (1|site)
  parsed <- parse_formula_re(f)

  expect_true(parsed$has_re)
  expect_equal(parsed$n_re_terms, 1)
  expect_equal(names(parsed$re_terms), "site")
})


test_that("formula parser handles no RE terms", {
  f <- y ~ x + z
  parsed <- parse_formula_re(f)

  expect_false(parsed$has_re)
  expect_equal(parsed$n_re_terms, 0)
  expect_length(parsed$re_terms, 0)
})


test_that("build_multi_re creates correct structure", {
  data <- data.frame(
    y = rnorm(20),
    site = rep(letters[1:4], 5),
    year = rep(2020:2024, 4)
  )

  re_terms <- list(
    site = list(group = "site", type = "intercept", term_idx = 1),
    year = list(group = "year", type = "intercept", term_idx = 2)
  )

  result <- build_multi_re(data, re_terms)

  expect_equal(result$n_re_terms, 2)
  expect_equal(names(result$re_list), c("site", "year"))
  expect_equal(result$re_list$site$n_groups, 4)
  expect_equal(result$re_list$year$n_groups, 5)
  # Total params: 4 + 1 (site) + 5 + 1 (year) = 11
  expect_equal(result$total_re_params, 11)
})


test_that("model with crossed RE in dynamics compiles", {
  skip_on_cran()

  # Create data with crossed structure
  set.seed(123)
  n_sites <- 3
  n_years <- 2
  n_times <- 4

  data <- expand.grid(
    site = letters[1:n_sites],
    year = 2020:(2020 + n_years - 1),
    time = 1:n_times
  )
  data$id <- interaction(data$site, data$year)
  data$x <- rnorm(nrow(data))
  data$y <- 5 + 0.5 * data$x + rnorm(nrow(data), 0, 0.5)

  # Model with crossed RE in dynamics
  model <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ x),
    dynamics("pre", y ~ x + (1|site) + (1|year)),
    dynamics("post", y ~ x)
  )

  # Should be able to prepare data
  expect_no_error({
    stan_data <- phaseR:::prepare_stan_data(model, data)
  })

  # Check multi-RE structure is created
  stan_data <- phaseR:::prepare_stan_data(model, data)
  expect_equal(stan_data$dyn_re_multi_0$n_re_terms, 2)
})


test_that("crossed RE model fits and produces valid output", {
  skip_on_cran()

  # Create simple crossed data
  set.seed(456)
  n_sites <- 3
  n_years <- 2
  n_times <- 3

  data <- expand.grid(
    site = letters[1:n_sites],
    year = 2020:(2020 + n_years - 1),
    time = 1:n_times
  )
  data$id <- interaction(data$site, data$year)
  data$x <- rnorm(nrow(data))

  # Add site and year effects
  site_effects <- setNames(rnorm(n_sites, 0, 0.5), letters[1:n_sites])
  year_effects <- setNames(rnorm(n_years, 0, 0.3), as.character(2020:(2020 + n_years - 1)))

  data$y <- 5 + site_effects[data$site] + year_effects[as.character(data$year)] +
            0.3 * data$x + rnorm(nrow(data), 0, 0.2)

  model <- phase_model(
    phases("baseline", "treatment"),
    transition("baseline", "treatment", ~ 1),
    dynamics("baseline", y ~ x + (1|site) + (1|year)),
    dynamics("treatment", y ~ x)
  )

  # Fit with minimal iterations
  fit <- fit_phaseR(model, data, chains = 1, iter = 200, warmup = 100, seed = 789)

  expect_s3_class(fit, "phase_fit")
  expect_true(fit$has_re)

  # Check parameter names include both RE terms
  expect_true(any(grepl("site", fit$param_names)))
  expect_true(any(grepl("year", fit$param_names)))

  # Check draws have right dimensions
  # 1 trans beta + 2 baseline betas + 2 treatment betas + 2 sigmas = 7
  # + 3 site RE + 1 sigma + 2 year RE + 1 sigma = 7
  n_expected <- 14
  expect_equal(ncol(fit$draws), n_expected)
})
