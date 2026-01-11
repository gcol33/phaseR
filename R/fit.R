#' Fit a phase model
#'
#' @param model A `phase_model` object
#' @param data A data frame with columns: id, time, outcome(s), covariates
#' @param chains Number of MCMC chains (default 4)
#' @param iter Total iterations per chain (default 2000)
#' @param warmup Warmup iterations per chain (default iter/2)
#' @param backend Estimation backend: "hmc" (default) or "laplace"
#' @param cores Number of cores for parallel chains (default 1)
#' @param seed Random seed for reproducibility
#' @param ... Additional arguments passed to backend
#'
#' @return A `phase_fit` object
#' @export
#'
#' @examples
#' \dontrun{
#' model <- phase_model(
#'   phases("inactive", "active"),
#'   transition("inactive", "active", ~ x),
#'   dynamics("inactive", y ~ 1),
#'   dynamics("active", y ~ x)
#' )
#' dat <- sim_phaseR(n_units = 50, n_times = 5)
#' fit <- fit_phaseR(model, data = dat)
#' }
#'
fit_phaseR <- function(model,
                       data,
                       chains = 4,
                       iter = 2000,
                       warmup = floor(iter / 2),
                       backend = c("hmc", "laplace"),
                       cores = 1,
                       seed = NULL,
                       ...) {

  backend <- match.arg(backend)

 # Validate
  validate_model_data(model, data)

  # Prepare data for C++
  stan_data <- prepare_stan_data(model, data)

  # Set seed
  if (!is.null(seed)) set.seed(seed)

  # Dispatch to backend
  fit_raw <- switch(backend,
    hmc = fit_hmc(stan_data, chains, iter, warmup, cores, ...),
    laplace = fit_laplace(stan_data, ...),
    stop("Unknown backend")
  )

  # Convert to phase_fit
  convert_to_phase_fit(fit_raw, model, data, backend, chains, iter, warmup)
}


#' Prepare data for C++ likelihood
#' @keywords internal
prepare_stan_data <- function(model, data) {

  # Sort data by id, then time
  data <- data[order(data$id, data$time), ]

  # Get unit indices
  unit_ids <- unique(data$id)
  n_units <- length(unit_ids)

  unit_start <- integer(n_units)
  unit_end <- integer(n_units)

  for (i in seq_along(unit_ids)) {
    idx <- which(data$id == unit_ids[i])
    unit_start[i] <- min(idx) - 1L  # 0-indexed for C++
    unit_end[i] <- max(idx) - 1L
  }

  # Build design matrices
  # Transition: use first transition's formula
  trans_formula <- model$transitions[[1]]$formula
  X_trans <- model.matrix(trans_formula, data)

  # Dynamics: get formula for each phase
  response_var <- model$dynamics[[1]]$response

  # For simplicity in v0.2.0: assume same predictors for both phases
  # (can be relaxed in later versions)
  dyn_formula_0 <- model$dynamics[[1]]$formula
  dyn_formula_1 <- if (length(model$dynamics) > 1) model$dynamics[[2]]$formula else dyn_formula_0

  # Build dynamics design matrix (use RHS only)
  X_dyn_0 <- model.matrix(update(dyn_formula_0, NULL ~ .), data)
  X_dyn_1 <- model.matrix(update(dyn_formula_1, NULL ~ .), data)

  # Parameter counts
  n_trans_coef <- ncol(X_trans)
  n_dyn_coef_0 <- ncol(X_dyn_0)
  n_dyn_coef_1 <- ncol(X_dyn_1)

  # Total parameters: trans + dyn_0 + dyn_1 + 2 sigmas
 n_params <- n_trans_coef + n_dyn_coef_0 + n_dyn_coef_1 + 2

  # Parameter names
  param_names <- c(
    paste0("beta_trans_", colnames(X_trans)),
    paste0("beta_", model$phases$names[1], "_", colnames(X_dyn_0)),
    paste0("beta_", model$phases$names[2], "_", colnames(X_dyn_1)),
    paste0("log_sigma_", model$phases$names[1]),
    paste0("log_sigma_", model$phases$names[2])
  )

  list(
    id = as.integer(factor(data$id)),
    time = data$time,
    y = data[[response_var]],
    X_trans = X_trans,
    X_dyn_0 = X_dyn_0,
    X_dyn_1 = X_dyn_1,
    n_units = n_units,
    unit_start = unit_start,
    unit_end = unit_end,
    n_trans_coef = n_trans_coef,
    n_dyn_coef_0 = n_dyn_coef_0,
    n_dyn_coef_1 = n_dyn_coef_1,
    n_params = n_params,
    param_names = param_names
  )
}


#' Convert raw fit to phase_fit object
#' @keywords internal
convert_to_phase_fit <- function(fit_raw, model, data, backend, chains, iter, warmup) {

  structure(
    list(
      draws = fit_raw$draws,
      chain_id = fit_raw$chain_id,
      param_names = fit_raw$param_names,
      model = model,
      data = data,
      backend = backend,
      n_chains = chains,
      n_iter = iter,
      n_warmup = warmup,
      n_units = length(unique(data$id)),
      n_times = max(data$time),
      diagnostics = fit_raw$diagnostics
    ),
    class = "phase_fit"
  )
}
