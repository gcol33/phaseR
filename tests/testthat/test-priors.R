# Tests for custom prior specification (v1.0.5)

# =============================================================================
# Distribution Constructors
# =============================================================================

test_that("normal() creates prior_spec object", {
  p <- normal(0, 2.5)
  expect_s3_class(p, "prior_spec")
  expect_equal(p$type, "normal")
  expect_equal(p$mean, 0)
  expect_equal(p$sd, 2.5)
})

test_that("normal() validates sd > 0", {
  expect_error(normal(0, 0), "sd must be positive")
  expect_error(normal(0, -1), "sd must be positive")
})

test_that("student_t() creates prior_spec object", {
  p <- student_t(3, 0, 2.5)
  expect_s3_class(p, "prior_spec")
  expect_equal(p$type, "student_t")
  expect_equal(p$df, 3)
  expect_equal(p$mean, 0)
  expect_equal(p$scale, 2.5)
})

test_that("student_t() validates parameters", {
  expect_error(student_t(0, 0, 2.5), "df must be positive")
  expect_error(student_t(3, 0, 0), "scale must be positive")
})

test_that("cauchy() creates prior_spec object", {
  p <- cauchy(0, 2.5)
  expect_s3_class(p, "prior_spec")
  expect_equal(p$type, "cauchy")
  expect_equal(p$location, 0)
  expect_equal(p$scale, 2.5)
})

test_that("cauchy() validates scale > 0", {
  expect_error(cauchy(0, 0), "scale must be positive")
})

test_that("half_normal() creates prior_spec object", {
  p <- half_normal(1)
  expect_s3_class(p, "prior_spec")
  expect_equal(p$type, "half_normal")
  expect_equal(p$sd, 1)
})

test_that("half_normal() validates sd > 0", {
  expect_error(half_normal(0), "sd must be positive")
})

test_that("half_cauchy() creates prior_spec object", {
  p <- half_cauchy(1)
  expect_s3_class(p, "prior_spec")
  expect_equal(p$type, "half_cauchy")
  expect_equal(p$scale, 1)
})

test_that("half_cauchy() validates scale > 0", {
  expect_error(half_cauchy(0), "scale must be positive")
})

test_that("lkj() creates prior_spec object", {
  p <- lkj(2)
  expect_s3_class(p, "prior_spec")
  expect_equal(p$type, "lkj")
  expect_equal(p$eta, 2)
})

test_that("lkj() validates eta > 0", {
  expect_error(lkj(0), "eta must be positive")
})


# =============================================================================
# prior() Function
# =============================================================================

test_that("prior() with defaults creates phase_prior", {
  p <- prior()
  expect_s3_class(p, "phase_prior")
  expect_equal(p$beta$type, "normal")
  expect_equal(p$sigma$type, "half_normal")
  expect_equal(p$re$type, "half_normal")
  expect_equal(p$cor$type, "lkj")
})

test_that("prior() with SD syntax (backward compatible)", {
  p <- prior(beta_sd = 1, sigma_sd = 0.5, re_sd = 0.3)
  expect_s3_class(p, "phase_prior")
  expect_equal(p$beta_sd, 1)
  expect_equal(p$sigma_sd, 0.5)
  expect_equal(p$re_sd, 0.3)
  expect_equal(p$beta$sd, 1)
})

test_that("prior() with formula syntax", {
  p <- prior(
    beta ~ normal(0, 1),
    sigma ~ half_cauchy(2)
  )
  expect_s3_class(p, "phase_prior")
  expect_equal(p$beta$type, "normal")
  expect_equal(p$beta$sd, 1)
  expect_equal(p$sigma$type, "half_cauchy")
  expect_equal(p$sigma$scale, 2)
})

test_that("prior() validates parameter classes", {
  expect_error(prior(foo ~ normal(0, 1)), "Unknown parameter class 'foo'")
})

test_that("prior() validates distribution for parameter class", {
  # sigma requires positive distribution
  expect_error(prior(sigma ~ normal(0, 1)), "requires a positive distribution")

  # cor requires lkj
  expect_error(prior(cor ~ half_normal(1)), "requires lkj")

  # beta requires unbounded distribution
  expect_error(prior(beta ~ half_normal(1)), "requires an unbounded distribution")
})

test_that("prior() allows student_t for beta", {
  p <- prior(beta ~ student_t(3, 0, 2.5))
  expect_equal(p$beta$type, "student_t")
  expect_equal(p$beta$df, 3)
})

test_that("prior() allows cauchy for beta", {
  p <- prior(beta ~ cauchy(0, 2.5))
  expect_equal(p$beta$type, "cauchy")
})

test_that("prior() allows half_cauchy for sigma", {
  p <- prior(sigma ~ half_cauchy(1))
  expect_equal(p$sigma$type, "half_cauchy")
  expect_equal(p$sigma$scale, 1)
})


# =============================================================================
# print.phase_prior
# =============================================================================

test_that("print.phase_prior works", {
  p <- prior()
  expect_output(print(p), "phaseR prior specification")
  expect_output(print(p), "Coefficients.*Normal")
  expect_output(print(p), "Residual SD.*Half-Normal")
})


# =============================================================================
# is_default_prior
# =============================================================================

test_that("is_default_prior identifies defaults", {
  expect_true(phaseR:::is_default_prior(NULL))
  expect_true(phaseR:::is_default_prior(prior()))
  expect_false(phaseR:::is_default_prior(prior(beta_sd = 1)))
  expect_false(phaseR:::is_default_prior(prior(beta ~ normal(0, 1))))
})


# =============================================================================
# prior_to_cpp
# =============================================================================

test_that("prior_to_cpp converts default priors", {
  cpp <- phaseR:::prior_to_cpp(prior())

  expect_equal(cpp$beta_type, 0L)  # normal
  expect_equal(cpp$beta_params[1], 0)   # mean
  expect_equal(cpp$beta_params[2], 2.5) # sd

  expect_equal(cpp$sigma_type, 3L)  # half_normal
  expect_equal(cpp$sigma_params[1], 1)  # sd

  expect_equal(cpp$re_type, 3L)  # half_normal
  expect_equal(cpp$re_params[1], 0.5)  # sd

  expect_equal(cpp$cor_type, 5L)  # lkj
  expect_equal(cpp$cor_params[1], 2)  # eta
})

test_that("prior_to_cpp converts student_t", {
  p <- prior(beta ~ student_t(5, 1, 3))
  cpp <- phaseR:::prior_to_cpp(p)

  expect_equal(cpp$beta_type, 1L)  # student_t
  expect_equal(cpp$beta_params[1], 5)  # df
  expect_equal(cpp$beta_params[2], 1)  # mean
  expect_equal(cpp$beta_params[3], 3)  # scale
})

test_that("prior_to_cpp converts half_cauchy", {
  p <- prior(sigma ~ half_cauchy(2))
  cpp <- phaseR:::prior_to_cpp(p)

  expect_equal(cpp$sigma_type, 4L)  # half_cauchy
  expect_equal(cpp$sigma_params[1], 2)  # scale
})


# =============================================================================
# Integration Tests with fit_phaseR
# =============================================================================

test_that("fit_phaseR accepts custom priors (Laplace backend)", {
  skip_on_cran()

  set.seed(123)
  n_units <- 15
  n_times <- 4

  data <- expand.grid(
    id = 1:n_units,
    time = 1:n_times
  )
  data$x <- rnorm(nrow(data))
  data$y <- 5 + 0.5 * data$x + rnorm(nrow(data), 0, 0.5)

  model <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ x),
    dynamics("pre", y ~ x),
    dynamics("post", y ~ x)
  )

  # Fit with custom priors
  custom_prior <- prior(
    beta ~ normal(0, 1),
    sigma ~ half_normal(0.5)
  )

  fit <- fit_phaseR(model, data, backend = "laplace",
                    prior = custom_prior, seed = 456)

  expect_s3_class(fit, "phase_fit")
  expect_true(nrow(fit$draws) > 0)
})


test_that("fit_phaseR with student_t prior (Laplace backend)", {
  skip_on_cran()

  set.seed(789)
  n_units <- 15
  n_times <- 4

  data <- expand.grid(
    id = 1:n_units,
    time = 1:n_times
  )
  data$x <- rnorm(nrow(data))
  data$y <- 5 + 0.5 * data$x + rnorm(nrow(data), 0, 0.5)

  model <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ x),
    dynamics("pre", y ~ x),
    dynamics("post", y ~ x)
  )

  # Fit with robust student_t prior
  robust_prior <- prior(beta ~ student_t(3, 0, 2.5))

  fit <- fit_phaseR(model, data, backend = "laplace",
                    prior = robust_prior, seed = 101)

  expect_s3_class(fit, "phase_fit")
  expect_true(nrow(fit$draws) > 0)
})


test_that("fit_phaseR with half_cauchy sigma prior (Laplace backend)", {
  skip_on_cran()

  set.seed(111)
  n_units <- 15
  n_times <- 4

  data <- expand.grid(
    id = 1:n_units,
    time = 1:n_times
  )
  data$x <- rnorm(nrow(data))
  data$y <- 5 + 0.5 * data$x + rnorm(nrow(data), 0, 0.5)

  model <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ x),
    dynamics("pre", y ~ x),
    dynamics("post", y ~ x)
  )

  # Fit with half-Cauchy prior on sigma
  hc_prior <- prior(sigma ~ half_cauchy(1))

  fit <- fit_phaseR(model, data, backend = "laplace",
                    prior = hc_prior, seed = 222)

  expect_s3_class(fit, "phase_fit")
  expect_true(nrow(fit$draws) > 0)
})


test_that("Custom priors with random effects (Laplace backend)", {
  skip_on_cran()

  set.seed(333)
  n_sites <- 4
  n_times <- 4

  data <- expand.grid(
    site = letters[1:n_sites],
    time = 1:n_times
  )
  data$id <- data$site
  data$x <- rnorm(nrow(data))
  data$y <- 5 + 0.3 * data$x + rnorm(nrow(data), 0, 0.2)

  model <- phase_model(
    phases("baseline", "treatment"),
    transition("baseline", "treatment", ~ 1),
    dynamics("baseline", y ~ x + (1|site)),
    dynamics("treatment", y ~ x)
  )

  # Custom priors including RE prior
  custom_prior <- prior(
    beta ~ normal(0, 1),
    sigma ~ half_normal(0.5),
    re ~ half_cauchy(0.5)
  )

  fit <- fit_phaseR(model, data, backend = "laplace",
                    prior = custom_prior, seed = 444)

  expect_s3_class(fit, "phase_fit")
  expect_true(fit$has_re)
})


test_that("Tighter priors affect estimates", {
  skip_on_cran()

  set.seed(555)
  n_units <- 20
  n_times <- 5

  data <- expand.grid(
    id = 1:n_units,
    time = 1:n_times
  )
  data$x <- rnorm(nrow(data))
  # True intercept is 10
  data$y <- 10 + 0.5 * data$x + rnorm(nrow(data), 0, 0.5)

  model <- phase_model(
    phases("pre", "post"),
    transition("pre", "post", ~ 1),
    dynamics("pre", y ~ x),
    dynamics("post", y ~ x)
  )

  # Weak prior
  fit_weak <- fit_phaseR(model, data, backend = "laplace",
                          prior = prior(beta ~ normal(0, 10)),
                          seed = 666, n_samples = 500)

  # Strong prior centered at 0
  fit_strong <- fit_phaseR(model, data, backend = "laplace",
                            prior = prior(beta ~ normal(0, 0.5)),
                            seed = 666, n_samples = 500)

  # With weak prior, intercept should be closer to true value (10)
  # With strong prior centered at 0, intercept should be pulled toward 0
  weak_intercept <- mean(fit_weak$draws[, "beta_pre_(Intercept)"])
  strong_intercept <- mean(fit_strong$draws[, "beta_pre_(Intercept)"])

  expect_gt(abs(weak_intercept), abs(strong_intercept))
})
