# Tests for nested random effects (v1.0.3)

test_that("expand_nested_re expands simple nesting", {
  formula_str <- "y ~ x + (1|site/plot)"
  result <- expand_nested_re(formula_str)

  expect_true(result$has_nested)
  expect_equal(result$interactions, "site:plot")
  expect_true(grepl("(1|site)", result$formula_str, fixed = TRUE))
  expect_true(grepl("(1|site:plot)", result$formula_str, fixed = TRUE))
})


test_that("expand_nested_re expands three-level nesting", {
  formula_str <- "y ~ x + (1|site/plot/subplot)"
  result <- expand_nested_re(formula_str)

  expect_true(result$has_nested)
  expect_equal(sort(result$interactions), c("site:plot", "site:plot:subplot"))
  expect_true(grepl("(1|site)", result$formula_str, fixed = TRUE))
  expect_true(grepl("(1|site:plot)", result$formula_str, fixed = TRUE))
  expect_true(grepl("(1|site:plot:subplot)", result$formula_str, fixed = TRUE))
})


test_that("expand_nested_re handles non-nested terms unchanged", {
  formula_str <- "y ~ x + (1|site)"
  result <- expand_nested_re(formula_str)

  expect_false(result$has_nested)
  expect_length(result$interactions, 0)
  expect_equal(result$formula_str, formula_str)
})


test_that("expand_nested_re handles mixed nested and non-nested", {
  formula_str <- "y ~ x + (1|year) + (1|site/plot)"
  result <- expand_nested_re(formula_str)

  expect_true(result$has_nested)
  expect_equal(result$interactions, "site:plot")
  expect_true(grepl("(1|year)", result$formula_str, fixed = TRUE))
  expect_true(grepl("(1|site)", result$formula_str, fixed = TRUE))
  expect_true(grepl("(1|site:plot)", result$formula_str, fixed = TRUE))
})


test_that("parse_formula_re handles nested syntax", {
  f <- y ~ x + (1|site/plot)
  parsed <- parse_formula_re(f)

  expect_true(parsed$has_re)
  expect_true(parsed$has_nested)
  expect_equal(parsed$nested_interactions, "site:plot")
  expect_equal(parsed$n_re_terms, 2)
  expect_true("site" %in% names(parsed$re_terms))
  expect_true("site:plot" %in% names(parsed$re_terms))
})


test_that("create_nested_interactions creates correct columns", {
  data <- data.frame(
    y = rnorm(12),
    site = rep(c("A", "B"), each = 6),
    plot = rep(c("p1", "p2", "p3"), 4)
  )

  result <- create_nested_interactions(data, "site:plot")

  expect_true("site:plot" %in% names(result))
  expect_equal(nlevels(factor(result[["site:plot"]])), 6)  # 2 sites * 3 plots
})


test_that("create_nested_interactions handles multiple interactions", {
  data <- data.frame(
    y = rnorm(24),
    site = rep(c("A", "B"), each = 12),
    plot = rep(c("p1", "p2", "p3"), 8),
    subplot = rep(c("s1", "s2"), 12)
  )

  result <- create_nested_interactions(data, c("site:plot", "site:plot:subplot"))

  expect_true("site:plot" %in% names(result))
  expect_true("site:plot:subplot" %in% names(result))
})


test_that("create_nested_interactions errors on missing variables", {
  data <- data.frame(y = rnorm(10), site = rep(c("A", "B"), 5))

  expect_error(
    create_nested_interactions(data, "site:plot"),
    "Grouping variable.*plot.*not found"
  )
})


test_that("validate_nesting accepts valid nesting", {
  # Plots are truly nested within sites
  data <- data.frame(
    site = c("A", "A", "B", "B"),
    plot = c("A_p1", "A_p2", "B_p1", "B_p2")
  )

  expect_true(validate_nesting(data, "site", "plot"))
})


test_that("validate_nesting rejects invalid nesting", {
  # Plot "shared" appears in both sites - invalid nesting
  data <- data.frame(
    site = c("A", "A", "B", "B"),
    plot = c("p1", "shared", "shared", "p2")
  )

  expect_error(
    validate_nesting(data, "site", "plot"),
    "Invalid nesting"
  )
})


test_that("model with nested RE in dynamics compiles", {
  skip_on_cran()

  set.seed(123)
  n_sites <- 2
  n_plots <- 3
  n_times <- 4

  data <- expand.grid(
    site = letters[1:n_sites],
    plot = 1:n_plots,
    time = 1:n_times
  )
  # Make plot names unique within site
  data$plot <- paste(data$site, data$plot, sep = "_p")
  data$id <- data$plot
  data$x <- rnorm(nrow(data))
  data$y <- 5 + 0.5 * data$x + rnorm(nrow(data), 0, 0.5)

  # Model with nested RE
  model <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ x),
    dynamics("pre", y ~ x + (1|site/plot)),
    dynamics("post", y ~ x)
  )

  # Should be able to prepare data
  expect_no_error({
    stan_data <- phaseR:::prepare_stan_data(model, data)
  })

  # Check nested structure is created
  stan_data <- phaseR:::prepare_stan_data(model, data)
  expect_equal(stan_data$dyn_re_multi_0$n_re_terms, 2)
  expect_true("site" %in% names(stan_data$dyn_re_multi_0$re_list))
  expect_true("site:plot" %in% names(stan_data$dyn_re_multi_0$re_list))
})


test_that("nested RE model fits and produces valid output", {
  skip_on_cran()

  set.seed(456)
  n_sites <- 3
  n_plots <- 2
  n_times <- 3

  data <- expand.grid(
    site = letters[1:n_sites],
    plot = 1:n_plots,
    time = 1:n_times
  )
  data$plot <- paste(data$site, data$plot, sep = "_p")
  data$id <- data$plot
  data$x <- rnorm(nrow(data))

  # Add site and plot effects
  site_effects <- setNames(rnorm(n_sites, 0, 0.5), letters[1:n_sites])
  plot_effects <- setNames(rnorm(n_sites * n_plots, 0, 0.3), unique(data$plot))

  data$y <- 5 + site_effects[data$site] + plot_effects[data$plot] +
            0.3 * data$x + rnorm(nrow(data), 0, 0.2)

  model <- phase_model(
    phases("baseline", "treatment"),
    transition("baseline", "treatment", ~ 1),
    dynamics("baseline", y ~ x + (1|site/plot)),
    dynamics("treatment", y ~ x)
  )

  # Fit with minimal iterations
  fit <- fit_phaseR(model, data, chains = 1, iter = 200, warmup = 100, seed = 789)

  expect_s3_class(fit, "phase_fit")
  expect_true(fit$has_re)

  # Check parameter names include both RE terms
  expect_true(any(grepl("site", fit$param_names)))
  expect_true(any(grepl("site:plot", fit$param_names) | grepl("site_p", fit$param_names)))
})


test_that("three-level nested RE works", {
  skip_on_cran()

  set.seed(789)
  data <- expand.grid(
    region = c("N", "S"),
    site = 1:2,
    plot = 1:2,
    time = 1:3
  )
  data$site <- paste(data$region, data$site, sep = "_s")
  data$plot <- paste(data$site, data$plot, sep = "_p")
  data$id <- data$plot
  data$x <- rnorm(nrow(data))
  data$y <- 5 + 0.5 * data$x + rnorm(nrow(data), 0, 0.5)

  model <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ 1),
    dynamics("pre", y ~ x + (1|region/site/plot)),
    dynamics("post", y ~ x)
  )

  # Should parse with 3 RE terms
  dyn_parsed <- parse_formula_re(model$dynamics[[1]]$formula)
  expect_equal(dyn_parsed$n_re_terms, 3)
  expect_true(dyn_parsed$has_nested)
  expect_equal(sort(dyn_parsed$nested_interactions), c("region:site", "region:site:plot"))
})
