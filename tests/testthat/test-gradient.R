# Central differences with one Richardson extrapolation step: the error is
# O(h^4) in the step and O(eps |f| / h) in rounding.
richardson_gradient <- function(f, x, h = 1e-4) {
  vapply(seq_along(x), function(i) {
    step <- h * max(1, abs(x[i]))
    central <- function(s) {
      e <- replace(numeric(length(x)), i, s)
      (f(x + e) - f(x - e)) / (2 * s)
    }
    (4 * central(step / 2) - central(step)) / 3
  }, numeric(1))
}

gradient_cases <- function() {
  dat <- reference_data()
  dat$year <- factor(dat$time)
  cases <- reference_cases()
  cases$crossed_year <- list(
    model = reference_two_phase(~ x + (1 | year), y ~ x + (1 | site),
                                y ~ x + (1 | year)),
    prior = prior(re ~ half_cauchy(1)))
  cases$slope_only <- list(
    model = reference_two_phase(~ x, y ~ 1 + (0 + x | site), y ~ x),
    prior = prior())
  cases$correlated_three <- list(
    model = reference_two_phase(~ x, y ~ x + z + (1 + x + z | region), y ~ x),
    prior = prior(cor ~ lkj(3), re ~ half_normal(1)))
  cases$binomial_re <- list(
    model = reference_two_phase(~ x, bin ~ x + (1 | site), bin ~ 1,
                                family_0 = binomial()),
    prior = prior())
  list(dat = dat, cases = cases)
}

test_that("exact gradients match finite differences for every model kind", {
  gc <- gradient_cases()
  for (name in names(gc$cases)) {
    target <- reference_target(gc$cases[[name]], gc$dat)
    f <- function(p) phaseR:::log_posterior(target, p)
    set.seed(nchar(name))
    for (point in 1:3) {
      theta <- rnorm(target$n_params, 0, 0.4)
      exact <- phaseR:::log_posterior(target, theta, gradient = TRUE)
      expect_equal(exact$value, f(theta), tolerance = 1e-12, info = name)
      numeric <- richardson_gradient(f, theta)
      rel_error <- max(abs(exact$gradient - numeric) / pmax(1, abs(numeric)))
      expect_lt(rel_error, 1e-6, label = paste(name, "point", point))
    }
  }
})

test_that("the multivariate normal target has gradient -P (q - mu)", {
  P <- matrix(c(2, 0.5, 0.5, 1), 2)
  target <- list(kind = "mvnormal", n_params = 2L,
                 data = list(mu = c(1, -1), precision = P))
  q <- c(0.3, 0.2)
  out <- phaseR:::log_posterior(target, q, gradient = TRUE)
  expect_equal(out$gradient, -drop(P %*% (q - c(1, -1))))
  expect_equal(out$value, -0.5 * drop(t(q - c(1, -1)) %*% P %*% (q - c(1, -1))))
})
