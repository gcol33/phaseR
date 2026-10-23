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
