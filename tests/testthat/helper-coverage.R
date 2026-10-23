# Coverage of posterior intervals over data sets simulated from known
# parameters. Each kind gives a model, the data design (covariates) and the
# true parameter vector; outcomes are simulated from the truth with
# posterior_predict() on a fit whose only draw is the truth (random effects
# drawn afresh from their distribution for every data set).

coverage_kinds <- function() {
  design <- function(n_units, n_times) {
    d <- expand.grid(time = seq_len(n_times), id = seq_len(n_units))[, c("id", "time")]
    d$x <- stats::rnorm(nrow(d))
    d$site <- factor((d$id - 1) %% 12)
    d$y <- 0
    d$count <- 0
    d
  }
  two_phase <- function(dyn_0, dyn_1, family = gaussian()) {
    phase_model(
      phases("inactive", "active"),
      transition("inactive", "active", ~ x),
      dynamics("inactive", dyn_0, family = family),
      dynamics("active", dyn_1, family = family)
    )
  }
  list(
    gaussian = list(
      model = two_phase(y ~ x, y ~ x),
      design = function() design(60, 5),
      truth = c(-1, 0.5, 10, 0, 8, -0.5, 0, log(1.5))),
    poisson = list(
      model = two_phase(count ~ x, count ~ x, family = poisson()),
      design = function() design(60, 5),
      truth = c(-1, 0.5, 1.5, 0.3, 0.5, -0.4, 0, 0)),
    random_intercept = list(
      model = two_phase(y ~ x + (1 | site), y ~ x),
      design = function() design(60, 5),
      truth = c(-1, 0.5, 10, 0, 8, -0.5, 0, log(1.5),
                rep(0, 12), log(0.5))),
    three_phase = list(
      model = phase_model(
        phases("A", "B", "C"),
        transition("A", "B", ~ x),
        transition("B", "C", ~ 1),
        dynamics("A", y ~ 1), dynamics("B", y ~ x), dynamics("C", y ~ 1)
      ),
      design = function() design(60, 6),
      truth = c(-1, 0.5, -1, 12, 10, -0.5, 8, 0, 0, 0))
  )
}

# Class of each parameter: transition or dynamics coefficient, log residual
# SD, log random-effect SD; random effects themselves are "effect".
parameter_class <- function(names) {
  ifelse(grepl("^beta_trans", names), "transition coefficient",
  ifelse(grepl("^beta_", names), "dynamics coefficient",
  ifelse(grepl("^log_sigma_.*_re_", names), "log random-effect SD",
  ifelse(grepl("^log_sigma_", names), "log residual SD", "effect"))))
}

simulate_coverage_data <- function(kind) {
  dat <- kind$design()
  stan_data <- phaseR:::prepare_stan_data(kind$model, dat, prior())
  truth <- kind$truth
  names(truth) <- stan_data$param_names
  # Fresh random effects from their distribution.
  effects <- grepl("^v_|^u_", names(truth))
  if (any(effects)) {
    truth[effects] <- stats::rnorm(sum(effects), 0,
                                   exp(truth[grepl("_re_", names(truth))][1]))
  }
  fit <- structure(list(model = kind$model, data = dat,
                        draws = matrix(truth, 1, dimnames = list(NULL, names(truth))),
                        param_names = names(truth),
                        k_phases = stan_data$k_phases),
                   class = "phase_fit")
  response <- kind$model$dynamics[[1]]$response
  dat[[response]] <- drop(posterior_predict(fit))
  list(data = dat, truth = truth)
}

# One row per data set and parameter: whether the central `level` interval
# of the posterior draws contains the truth. Parameters that enter only
# their prior (residual SDs of non-Gaussian phases) and the random effects
# are left out.
run_coverage <- function(kind, backend, n_rep, level = 0.9, seed = 1, ...) {
  gaussian_phase <- vapply(kind$model$dynamics, function(d) {
    d$family$family == "gaussian"
  }, logical(1))
  out <- vector("list", n_rep)
  for (r in seq_len(n_rep)) {
    set.seed(seed * 1000 + r)
    sim <- simulate_coverage_data(kind)
    fit <- suppressMessages(fit_phaseR(kind$model, sim$data, backend = backend,
                                       seed = seed * 1000 + r, ...))
    lower <- apply(fit$draws, 2, stats::quantile, (1 - level) / 2)
    upper <- apply(fit$draws, 2, stats::quantile, 1 - (1 - level) / 2)
    keep <- parameter_class(names(sim$truth)) != "effect"
    sigma_names <- paste0("log_sigma_", kind$model$phases$names[!gaussian_phase])
    keep <- keep & !(names(sim$truth) %in% sigma_names)
    out[[r]] <- data.frame(
      rep = r,
      parameter = names(sim$truth)[keep],
      class = parameter_class(names(sim$truth))[keep],
      covered = sim$truth[keep] >= lower[keep] & sim$truth[keep] <= upper[keep],
      row.names = NULL
    )
  }
  do.call(rbind, out)
}

coverage_by_class <- function(results) {
  agg <- stats::aggregate(covered ~ class, data = results,
                          FUN = function(x) c(coverage = mean(x), n = length(x)))
  data.frame(class = agg$class, coverage = agg$covered[, "coverage"],
             n = agg$covered[, "n"])
}
