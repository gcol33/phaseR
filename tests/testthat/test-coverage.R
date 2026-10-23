# Coverage of 90% central posterior intervals over data sets simulated from
# known parameters. The light version runs in the suite (not on CRAN); the
# full study over every model kind and both backends runs with
# PHASER_COVERAGE=true and prints its table.

test_that("90% intervals cover the truth at about the nominal rate", {
  skip_on_cran()
  kind <- coverage_kinds()$gaussian
  laplace <- run_coverage(kind, "laplace", n_rep = 40, seed = 1)
  nuts <- run_coverage(kind, "hmc", n_rep = 15, seed = 2, chains = 2,
                       iter = 600, warmup = 300)
  # Pooled over parameters, 40 (15) data sets x 8 parameters: the binomial
  # SD of the coverage is about 0.015 (0.025) under independence.
  expect_gt(mean(laplace$covered), 0.82)
  expect_gt(mean(nuts$covered), 0.80)
  expect_lt(mean(laplace$covered), 0.97)
  expect_lt(mean(nuts$covered), 0.98)
})

test_that("coverage study over model kinds and backends", {
  skip_if_not(identical(Sys.getenv("PHASER_COVERAGE"), "true"),
              "set PHASER_COVERAGE=true for the full coverage study")
  n_rep <- as.integer(Sys.getenv("PHASER_COVERAGE_REPS", "100"))
  table <- list()
  for (name in names(coverage_kinds())) {
    kind <- coverage_kinds()[[name]]
    for (backend in c("laplace", "hmc")) {
      res <- if (backend == "hmc") {
        run_coverage(kind, backend, n_rep, seed = 3, chains = 2, iter = 800,
                     warmup = 400)
      } else {
        run_coverage(kind, backend, n_rep, seed = 3)
      }
      table[[length(table) + 1]] <- cbind(kind = name, backend = backend,
                                          coverage_by_class(res))
    }
  }
  table <- do.call(rbind, table)
  print(table, row.names = FALSE, digits = 3)
  out <- Sys.getenv("PHASER_COVERAGE_OUT")
  if (nzchar(out)) utils::write.csv(table, out, row.names = FALSE)
  # A 90% interval covers in [0.80, 0.97] for every class with n >= 100.
  expect_true(all(table$coverage[table$n >= 100] > 0.80))
})
