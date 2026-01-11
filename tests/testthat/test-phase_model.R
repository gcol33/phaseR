test_that("phase_model() composes valid model", {

  model <- phase_model(
    phases("inactive", "active"),
    transition("inactive", "active", ~ x),
    dynamics("inactive", y ~ 1),
    dynamics("active", y ~ x)
  )

  expect_s3_class(model, "phase_model")
  expect_length(model$transitions, 1)
  expect_length(model$dynamics, 2)
})

test_that("phase_model() validates transition phases", {
  expect_error(
    phase_model(
      phases("A", "B"),
      transition("A", "C", ~ 1),
      dynamics("A", y ~ 1),
      dynamics("B", y ~ 1)
    ),
    "not in declared phases"
  )
})

test_that("phase_model() rejects backward transitions", {
  expect_error(
    phase_model(
      phases("A", "B"),
      transition("B", "A", ~ 1),
      dynamics("A", y ~ 1)
    ),
    "violates phase ordering"
  )
})

test_that("phase_model() requires dynamics for non-absorbing phases", {
  expect_error(
    phase_model(
      phases("A", "B", "C"),
      transition("A", "B", ~ 1),
      transition("B", "C", ~ 1),
      dynamics("A", y ~ 1)
      # Missing dynamics for "B"
    ),
    "Missing dynamics"
  )
})

test_that("print.phase_model works", {
  model <- phase_model(
    phases("inactive", "active"),
    transition("inactive", "active", ~ x),
    dynamics("inactive", y ~ 1),
    dynamics("active", y ~ x)
  )
  expect_output(print(model), "Phase Model")
  expect_output(print(model), "inactive -> active")
})
