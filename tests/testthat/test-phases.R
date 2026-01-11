test_that("phases() creates valid phase_spec", {

  p <- phases("A", "B")
  expect_s3_class(p, "phase_spec")
  expect_equal(p$n_phases, 2)
  expect_equal(p$names, c("A", "B"))
  expect_equal(p$absorbing, "B")
})

test_that("phases() requires at least 2 phases", {
  expect_error(phases("A"), "At least 2 phases")
})

test_that("phases() rejects duplicate names", {
  expect_error(phases("A", "A", "B"), "unique")
})

test_that("phases() validates absorbing phases", {
  expect_error(phases("A", "B", absorbing = "C"), "declared phases")
})

test_that("print.phase_spec works", {
  p <- phases("inactive", "active")
  expect_output(print(p), "Phase specification")
  expect_output(print(p), "inactive")
  expect_output(print(p), "active.*absorbing")
})
