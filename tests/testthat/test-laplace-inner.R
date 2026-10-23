# The Laplace objective written out with a dense numerical Hessian in the
# effects at the mode the C++ search returns.
dense_laplace_value <- function(target, theta) {
  lap <- phaseR:::phaseR_laplace(target)
  layout <- phaseR:::phaseR_laplace_layout(lap)
  cond <- phaseR:::phaseR_laplace_conditional(lap, theta)
  q <- numeric(target$n_params)
  q[layout$theta] <- theta
  q[layout$u] <- cond$mode
  density <- phaseR:::phaseR_density(target)
  grad <- function(p) phaseR:::phaseR_log_density_gradient(density, p)$gradient
  H <- vapply(layout$u, function(i) {
    e <- replace(numeric(length(q)), i, 1e-5)
    (grad(q + e) - grad(q - e)) / 2e-5
  }, numeric(length(q)))[layout$u, , drop = FALSE]
  list(
    value = phaseR:::phaseR_log_density(density, q) +
      length(layout$u) / 2 * log(2 * pi) -
      0.5 * determinant(-(H + t(H)) / 2)$modulus[1],
    max_grad_u = max(abs(grad(q)[layout$u])),
    objective = phaseR:::phaseR_laplace_objective(lap, theta)$value
  )
}

test_that("the block-diagonal Laplace objective equals the dense formula", {
  gc <- gradient_cases()
  unit_re <- list(
    model = reference_two_phase(~ x + (1 | id), y ~ x + (1 | id), y ~ x),
    prior = prior())
  cases <- c(gc$cases[c("re", "multi_re", "crossed_year", "slopes_correlated",
                        "k_phase_glm_re")], list(unit_re = unit_re))
  for (name in names(cases)) {
    target <- reference_target(cases[[name]], gc$dat)
    lap <- phaseR:::phaseR_laplace(target)
    layout <- phaseR:::phaseR_laplace_layout(lap)
    theta <- phaseR:::initial_values(target, 1)[[1]][layout$theta]
    out <- dense_laplace_value(target, theta)
    expect_lt(out$max_grad_u, 1e-6, label = name)
    expect_equal(out$objective, out$value, tolerance = 1e-7, info = name)
  }
})

test_that("the Laplace objective is -Inf where the effects' mode is not finite", {
  case <- reference_cases()$re
  target <- reference_target(case)
  lap <- phaseR:::phaseR_laplace(target)
  layout <- phaseR:::phaseR_laplace_layout(lap)
  far <- rep(c(800, -800), length.out = length(layout$theta))
  out <- phaseR:::phaseR_laplace_objective(lap, far)
  expect_equal(out$value, -Inf)
  expect_true(all(is.na(out$gradient)))
})

test_that("starts that reach the same point share one search", {
  starts <- list(c(1, 2), c(1 + 1e-6, 2), c(3, 2), c(1, 2 - 1e-7), c(3, 2))
  expect_equal(phaseR:::distinct_starts(starts), c(1L, 1L, 3L, 1L, 3L))
})
