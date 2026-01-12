# Tests for random slopes (v1.0.2)

test_that("formula parser extracts correlated random slopes", {
  f <- y ~ x + (1 + x | site)
  parsed <- parse_formula_re(f)

  expect_true(parsed$has_re)
  expect_true(parsed$has_slopes)
  expect_equal(parsed$n_re_terms, 1)

  term <- parsed$re_terms$site
  expect_equal(term$type, "intercept_slope")
  expect_true(term$has_intercept)
  expect_equal(term$slope_vars, "x")
  expect_true(term$correlated)
  expect_equal(term$n_coefs, 2)

  expect_equal(deparse(parsed$fixed_formula), "y ~ x")
})


test_that("formula parser extracts uncorrelated random slopes", {
  f <- y ~ x + (1 + x || site)
  parsed <- parse_formula_re(f)

  expect_true(parsed$has_re)
  expect_true(parsed$has_slopes)

  term <- parsed$re_terms$site
  expect_equal(term$type, "intercept_slope")
  expect_false(term$correlated)
  expect_equal(term$n_coefs, 2)
})


test_that("formula parser extracts slope-only RE", {
  f <- y ~ x + (0 + x | site)
  parsed <- parse_formula_re(f)

  expect_true(parsed$has_re)
  expect_true(parsed$has_slopes)

  term <- parsed$re_terms$site
  expect_equal(term$type, "slope")
  expect_false(term$has_intercept)
  expect_equal(term$slope_vars, "x")
  expect_equal(term$n_coefs, 1)
})


test_that("formula parser handles multiple slope variables", {
  f <- y ~ x + z + (1 + x + z | site)
  parsed <- parse_formula_re(f)

  term <- parsed$re_terms$site
  expect_equal(term$type, "intercept_slope")
  expect_equal(term$slope_vars, c("x", "z"))
  expect_equal(term$n_coefs, 3)
})


test_that("build_re_z_matrix creates correct Z matrix", {
  data <- data.frame(
    y = rnorm(10),
    x = 1:10,
    site = rep(letters[1:2], 5)
  )

  term <- list(
    group = "site",
    type = "intercept_slope",
    has_intercept = TRUE,
    slope_vars = "x",
    correlated = TRUE,
    n_coefs = 2
  )

  Z <- build_re_z_matrix(data, term)

  expect_equal(dim(Z), c(10, 2))
  expect_equal(Z[, 1], rep(1, 10))  # Intercept column
  expect_equal(Z[, 2], 1:10)        # Slope column
})


test_that("build_multi_re calculates correct parameter counts for correlated", {
  data <- data.frame(
    y = rnorm(20),
    x = rnorm(20),
    site = rep(letters[1:4], 5)
  )

  re_terms <- list(
    site = list(
      group = "site",
      type = "intercept_slope",
      has_intercept = TRUE,
      slope_vars = "x",
      correlated = TRUE,
      n_coefs = 2
    )
  )

  result <- build_multi_re(data, re_terms)

  # 4 groups * 2 coefs = 8 RE values
  # Cholesky: 2*(2+1)/2 = 3 parameters
  # Total: 11
  expect_equal(result$total_re_params, 11)
  expect_true(result$has_slopes)
  expect_equal(result$re_list$site$n_coefs, 2)
})


test_that("build_multi_re calculates correct parameter counts for uncorrelated", {
  data <- data.frame(
    y = rnorm(20),
    x = rnorm(20),
    site = rep(letters[1:4], 5)
  )

  re_terms <- list(
    site = list(
      group = "site",
      type = "intercept_slope",
      has_intercept = TRUE,
      slope_vars = "x",
      correlated = FALSE,
      n_coefs = 2
    )
  )

  result <- build_multi_re(data, re_terms)

  # 4 groups * 2 coefs = 8 RE values
  # Diagonal: 2 log_sigma parameters
  # Total: 10
  expect_equal(result$total_re_params, 10)
})


test_that("build_re_param_names generates correct names for correlated slopes", {
  re_info <- list(
    n_groups = 3,
    n_coefs = 2,
    levels = c("a", "b", "c"),
    has_intercept = TRUE,
    slope_vars = "x",
    correlated = TRUE
  )

  names <- build_re_param_names("baseline", "site", re_info)

  # Should have: 3 groups * 2 coefs = 6 RE values + 3 Cholesky params
  expect_equal(length(names), 9)

  # Check RE value names
  expect_true("v_baseline_site_a_intercept" %in% names)
  expect_true("v_baseline_site_a_x" %in% names)
  expect_true("v_baseline_site_c_x" %in% names)

  # Check Cholesky names
  expect_true("L_baseline_re_site_1_1" %in% names)
  expect_true("L_baseline_re_site_2_1" %in% names)
  expect_true("L_baseline_re_site_2_2" %in% names)
})


test_that("build_re_param_names generates correct names for uncorrelated slopes", {
  re_info <- list(
    n_groups = 3,
    n_coefs = 2,
    levels = c("a", "b", "c"),
    has_intercept = TRUE,
    slope_vars = "x",
    correlated = FALSE
  )

  names <- build_re_param_names("baseline", "site", re_info)

  # Should have: 3 groups * 2 coefs = 6 RE values + 2 log_sigma params
  expect_equal(length(names), 8)

  # Check variance names
  expect_true("log_sigma_baseline_re_site_intercept" %in% names)
  expect_true("log_sigma_baseline_re_site_x" %in% names)
})


test_that("model with random slopes in dynamics compiles", {
  skip_on_cran()

  set.seed(123)
  n_sites <- 3
  n_times <- 4

  data <- expand.grid(
    site = letters[1:n_sites],
    time = 1:n_times
  )
  data$id <- data$site
  data$x <- rnorm(nrow(data))
  data$y <- 5 + 0.5 * data$x + rnorm(nrow(data), 0, 0.5)

  # Model with correlated random slope in dynamics
  model <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ x),
    dynamics("pre", y ~ x + (1 + x | site)),
    dynamics("post", y ~ x)
  )

  # Should be able to prepare data
  expect_no_error({
    stan_data <- phaseR:::prepare_stan_data(model, data)
  })

  # Check slope structure is created
  stan_data <- phaseR:::prepare_stan_data(model, data)
  expect_true(stan_data$has_slopes)
  expect_true(stan_data$dyn_re_multi_0$has_slopes)
  expect_equal(stan_data$dyn_re_multi_0$re_list$site$n_coefs, 2)
})


test_that("model with uncorrelated random slopes compiles", {
  skip_on_cran()

  set.seed(456)
  n_sites <- 3
  n_times <- 4

  data <- expand.grid(
    site = letters[1:n_sites],
    time = 1:n_times
  )
  data$id <- data$site
  data$x <- rnorm(nrow(data))
  data$y <- 5 + 0.5 * data$x + rnorm(nrow(data), 0, 0.5)

  # Model with uncorrelated random slope
  model <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ x),
    dynamics("pre", y ~ x + (1 + x || site)),
    dynamics("post", y ~ x)
  )

  stan_data <- phaseR:::prepare_stan_data(model, data)
  expect_true(stan_data$has_slopes)
  expect_false(stan_data$dyn_re_multi_0$re_list$site$correlated)
})


test_that("random slopes model fits and produces valid output", {
  skip_on_cran()

  # Create data with random slopes
  set.seed(789)
  n_sites <- 4
  n_times <- 5

  data <- expand.grid(
    site = letters[1:n_sites],
    time = 1:n_times
  )
  data$id <- data$site
  data$x <- rnorm(nrow(data))

  # True random effects: intercept and slope vary by site
  site_int <- setNames(rnorm(n_sites, 0, 0.5), letters[1:n_sites])
  site_slope <- setNames(rnorm(n_sites, 0, 0.3), letters[1:n_sites])

  data$y <- 5 + site_int[data$site] +
            (0.5 + site_slope[data$site]) * data$x +
            rnorm(nrow(data), 0, 0.2)

  model <- phase_model(
    phases("baseline", "treatment"),
    transition("baseline", "treatment", ~ 1),
    dynamics("baseline", y ~ x + (1 + x | site)),
    dynamics("treatment", y ~ x)
  )

  # Fit with minimal iterations
  fit <- fit_phaseR(model, data, chains = 1, iter = 200, warmup = 100, seed = 101)

  expect_s3_class(fit, "phase_fit")
  expect_true(fit$has_re)

  # Check parameter names include slope components
  expect_true(any(grepl("_intercept", fit$param_names)))
  expect_true(any(grepl("L_baseline_re_site", fit$param_names)))
})
