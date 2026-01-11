test_that("dynamics() creates valid phase_dynamics", {

  dyn <- dynamics("active", y ~ x)
  expect_s3_class(dyn, "phase_dynamics")
  expect_equal(dyn$phase, "active")
  expect_equal(dyn$response, "y")
})

test_that("dynamics() requires two-sided formula", {
  expect_error(dynamics("A", ~ x), "two-sided")
})

test_that("dynamics() rejects difference terms", {
  expect_error(dynamics("A", diff(y) ~ x), "difference|delta")
  expect_error(dynamics("A", delta_y ~ x), "difference|delta")
})

test_that("print.phase_dynamics works", {
  dyn <- dynamics("active", y ~ x, family = poisson())
  expect_output(print(dyn), "active")
  expect_output(print(dyn), "poisson")
})
