#' HMC backend for phase models
#'
#' @param stan_data Prepared data list
#' @param chains Number of chains
#' @param iter Iterations per chain
#' @param warmup Warmup iterations
#' @param cores Number of cores
#' @param ... Additional HMC parameters
#'
#' @return Raw fit object
#' @keywords internal
#'
fit_hmc <- function(stan_data, chains, iter, warmup, cores, ...) {

  n_params <- stan_data$n_params

  # Initialize chains
  inits <- lapply(seq_len(chains), function(i) {
    rnorm(n_params, 0, 0.1)
  })

  # Run chains (parallel if cores > 1)
  if (cores > 1 && chains > 1) {
    draws_list <- parallel::mclapply(seq_len(chains), function(i) {
      run_nuts_chain(stan_data, inits[[i]], iter, warmup, ...)
    }, mc.cores = min(cores, chains))
  } else {
    draws_list <- lapply(seq_len(chains), function(i) {
      run_nuts_chain(stan_data, inits[[i]], iter, warmup, ...)
    })
  }

  # Combine chains
  all_draws <- do.call(rbind, lapply(draws_list, function(x) x$draws))
  colnames(all_draws) <- stan_data$param_names

  list(
    draws = all_draws,
    chain_id = rep(seq_len(chains), each = iter - warmup),
    diagnostics = combine_diagnostics(draws_list),
    n_params = n_params,
    param_names = stan_data$param_names
  )
}


#' Run a single NUTS chain
#' @keywords internal
run_nuts_chain <- function(stan_data, init, iter, warmup, ...) {

  # Call C++ NUTS sampler
  result <- phaseR_nuts_sampler(
    data = list(
      id = stan_data$id,
      time = stan_data$time,
      y = stan_data$y,
      X_trans = stan_data$X_trans,
      X_dyn = stan_data$X_dyn_0,  # Use phase 0 design matrix
      n_units = stan_data$n_units,
      unit_start = stan_data$unit_start,
      unit_end = stan_data$unit_end
    ),
    init = init,
    n_iter = as.integer(iter),
    n_warmup = as.integer(warmup),
    n_trans_coef = stan_data$n_trans_coef,
    n_dyn_coef_0 = stan_data$n_dyn_coef_0,
    n_dyn_coef_1 = stan_data$n_dyn_coef_1
  )

  list(
    draws = result$draws,
    diagnostics = list(
      n_divergent = result$n_divergent,
      accept_prob = result$accept_prob,
      step_size = result$step_size
    )
  )
}


#' Combine diagnostics from multiple chains
#' @keywords internal
combine_diagnostics <- function(draws_list) {
  list(
    n_divergent = sum(vapply(draws_list, function(x) x$diagnostics$n_divergent, integer(1))),
    mean_accept_prob = mean(vapply(draws_list, function(x) x$diagnostics$accept_prob, numeric(1))),
    step_sizes = vapply(draws_list, function(x) x$diagnostics$step_size, numeric(1))
  )
}


#' Laplace approximation backend (placeholder)
#' @keywords internal
fit_laplace <- function(stan_data, ...) {
  stop("Laplace backend not yet implemented", call. = FALSE)
}
