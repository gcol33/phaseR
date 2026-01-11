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
  convert_to_phase_fit(fit_raw, model, data, backend, chains, iter, warmup, stan_data)
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

  # Parse formulas for random effects
  trans_formula <- model$transitions[[1]]$formula
  trans_parsed <- parse_formula_re(trans_formula)

  dyn_formula_0 <- model$dynamics[[1]]$formula
  dyn_formula_1 <- if (length(model$dynamics) > 1) model$dynamics[[2]]$formula else dyn_formula_0
  response_var <- model$dynamics[[1]]$response

  dyn_parsed_0 <- parse_formula_re(dyn_formula_0)
  dyn_parsed_1 <- parse_formula_re(dyn_formula_1)

  # Build design matrices using fixed-effects formulas
  X_trans <- model.matrix(trans_parsed$fixed_formula, data)
  X_dyn_0 <- model.matrix(update(dyn_parsed_0$fixed_formula, NULL ~ .), data)
  X_dyn_1 <- model.matrix(update(dyn_parsed_1$fixed_formula, NULL ~ .), data)

  # Parameter counts
  n_trans_coef <- ncol(X_trans)
  n_dyn_coef_0 <- ncol(X_dyn_0)
  n_dyn_coef_1 <- ncol(X_dyn_1)

  # Check for random effects
  has_trans_re <- trans_parsed$has_re
  has_dyn_re_0 <- dyn_parsed_0$has_re
  has_dyn_re_1 <- dyn_parsed_1$has_re
  has_re <- has_trans_re || has_dyn_re_0 || has_dyn_re_1

  # Initialize RE structures
  n_trans_re <- 0L
  n_dyn_re_0 <- 0L
  n_dyn_re_1 <- 0L
  trans_re_idx <- integer(0)
  dyn_re_idx_0 <- integer(0)
  dyn_re_idx_1 <- integer(0)
  trans_re_levels <- character(0)
  dyn_re_levels_0 <- character(0)
  dyn_re_levels_1 <- character(0)

  if (has_trans_re) {
    group_var <- trans_parsed$re_groups[1]
    re_info <- build_re_index(data, group_var)
    trans_re_idx <- re_info$idx
    n_trans_re <- re_info$n_groups
    trans_re_levels <- re_info$levels
  }

  if (has_dyn_re_0) {
    group_var <- dyn_parsed_0$re_groups[1]
    re_info <- build_re_index(data, group_var)
    dyn_re_idx_0 <- re_info$idx
    n_dyn_re_0 <- re_info$n_groups
    dyn_re_levels_0 <- re_info$levels
  }

  if (has_dyn_re_1) {
    group_var <- dyn_parsed_1$re_groups[1]
    re_info <- build_re_index(data, group_var)
    dyn_re_idx_1 <- re_info$idx
    n_dyn_re_1 <- re_info$n_groups
    dyn_re_levels_1 <- re_info$levels
  }

  # Total parameters
  # Fixed: trans + dyn_0 + dyn_1 + 2 sigmas
  # RE: u_trans(n_trans_re) + log_sigma_trans_re + v_0(n_dyn_re_0) + log_sigma_dyn_re_0 + ...
  n_params <- n_trans_coef + n_dyn_coef_0 + n_dyn_coef_1 + 2
  if (has_trans_re) n_params <- n_params + n_trans_re + 1
  if (has_dyn_re_0) n_params <- n_params + n_dyn_re_0 + 1
  if (has_dyn_re_1) n_params <- n_params + n_dyn_re_1 + 1

  # Parameter names
  param_names <- c(
    paste0("beta_trans_", colnames(X_trans)),
    paste0("beta_", model$phases$names[1], "_", colnames(X_dyn_0)),
    paste0("beta_", model$phases$names[2], "_", colnames(X_dyn_1)),
    paste0("log_sigma_", model$phases$names[1]),
    paste0("log_sigma_", model$phases$names[2])
  )

  if (has_trans_re) {
    param_names <- c(param_names,
                     paste0("u_trans_", trans_re_levels),
                     "log_sigma_trans_re")
  }
  if (has_dyn_re_0) {
    param_names <- c(param_names,
                     paste0("v_dyn_", model$phases$names[1], "_", dyn_re_levels_0),
                     paste0("log_sigma_dyn_re_", model$phases$names[1]))
  }
  if (has_dyn_re_1) {
    param_names <- c(param_names,
                     paste0("v_dyn_", model$phases$names[2], "_", dyn_re_levels_1),
                     paste0("log_sigma_dyn_re_", model$phases$names[2]))
  }

  list(
    id = as.integer(factor(data$id)),
    time = data$time,
    y = data[[response_var]],
    X_trans = X_trans,
    X_dyn = X_dyn_0,  # For backward compatibility
    X_dyn_0 = X_dyn_0,
    X_dyn_1 = X_dyn_1,
    n_units = n_units,
    unit_start = unit_start,
    unit_end = unit_end,
    n_trans_coef = n_trans_coef,
    n_dyn_coef_0 = n_dyn_coef_0,
    n_dyn_coef_1 = n_dyn_coef_1,
    n_params = n_params,
    param_names = param_names,
    # RE info
    has_re = has_re,
    has_trans_re = has_trans_re,
    has_dyn_re_0 = has_dyn_re_0,
    has_dyn_re_1 = has_dyn_re_1,
    n_trans_re = n_trans_re,
    n_dyn_re_0 = n_dyn_re_0,
    n_dyn_re_1 = n_dyn_re_1,
    trans_re_idx = trans_re_idx,
    dyn_re_idx_0 = dyn_re_idx_0,
    dyn_re_idx_1 = dyn_re_idx_1,
    trans_re_levels = trans_re_levels,
    dyn_re_levels_0 = dyn_re_levels_0,
    dyn_re_levels_1 = dyn_re_levels_1
  )
}


#' Convert raw fit to phase_fit object
#' @keywords internal
convert_to_phase_fit <- function(fit_raw, model, data, backend, chains, iter, warmup, stan_data) {

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
      diagnostics = fit_raw$diagnostics,
      # RE info
      has_re = stan_data$has_re,
      has_trans_re = stan_data$has_trans_re,
      has_dyn_re = stan_data$has_dyn_re_0 || stan_data$has_dyn_re_1,
      trans_re_levels = stan_data$trans_re_levels,
      dyn_re_levels = list(
        phase_0 = stan_data$dyn_re_levels_0,
        phase_1 = stan_data$dyn_re_levels_1
      )
    ),
    class = "phase_fit"
  )
}
