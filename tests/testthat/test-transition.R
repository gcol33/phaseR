test_that("transition() creates valid phase_transition", {

  tr <- transition("A", "B", ~ x)
  expect_s3_class(tr, "phase_transition")
  expect_equal(tr$from, "A")
  expect_equal(tr$to, "B")
})

test_that("transition() defaults to intercept-only", {
  tr <- transition("A", "B")
  expect_equal(deparse(tr$formula), "~1")
})

test_that("transition() rejects two-sided formula", {
  expect_error(transition("A", "B", y ~ x), "one-sided")
})

test_that("print.phase_transition works", {
  tr <- transition("inactive", "active", ~ treatment)
  expect_output(print(tr), "inactive -> active")
  expect_output(print(tr), "treatment")
})
