#' HMC backend for phase models
#'
#' Runs the No-U-Turn sampler on the log posterior of the model, one chain
#' after another or in parallel.
#'
#' @param stan_data Prepared data list
#' @param chains Number of chains
#' @param iter Iterations per chain
#' @param warmup Warmup iterations
#' @param cores Number of cores
#' @param seed Integer seed; chain `c` samples from the stream seeded by
#'   `(seed, c)`, which also draws its initial values
#' @param adapt_delta Target average acceptance statistic of the step-size
#'   adaptation, in (0, 1)
#' @param max_treedepth Maximum tree depth; a transition takes at most
#'   `2^max_treedepth - 1` leapfrog steps
#' @param ... Unused
#'
#' @return Raw fit object
#' @keywords internal
#'
fit_hmc <- function(stan_data, chains, iter, warmup, cores, seed,
                    adapt_delta = 0.8, max_treedepth = 10L, ...) {

  target <- nuts_target(stan_data)
  run_chain <- function(chain) {
    phaseR_nuts_sampler(
      target = target,
      n_iter = as.integer(iter),
      n_warmup = as.integer(warmup),
      seed = as.integer(seed),
      chain = as.integer(chain),
      adapt_delta = as.numeric(adapt_delta),
      max_treedepth = as.integer(max_treedepth)
    )
  }

  if (cores > 1 && chains > 1) {
    chain_list <- parallel::mclapply(seq_len(chains), run_chain,
                                     mc.cores = min(cores, chains))
  } else {
    chain_list <- lapply(seq_len(chains), run_chain)
  }

  all_draws <- do.call(rbind, lapply(chain_list, function(x) x$draws))
  colnames(all_draws) <- stan_data$param_names

  list(
    draws = all_draws,
    chain_id = rep(seq_len(chains), each = iter - warmup),
    diagnostics = combine_diagnostics(chain_list, max_treedepth),
    n_params = stan_data$n_params,
    param_names = stan_data$param_names,
    k_phases = stan_data$k_phases
  )
}


#' Describe the log posterior of a phase model for the NUTS sampler
#'
#' @param stan_data Prepared data list
#' @return List with `kind` (which C++ log posterior to use), `n_params`,
#'   `data` and `args` (the dimensions and random-effect structures that log
#'   posterior takes)
#' @keywords internal
nuts_target <- function(stan_data) {

  base_data <- list(
    id = stan_data$id,
    time = stan_data$time,
    y = stan_data$y,
    n_units = stan_data$n_units,
    unit_start = stan_data$unit_start,
    unit_end = stan_data$unit_end
  )
  two_phase_data <- c(base_data, list(X_trans = stan_data$X_trans,
                                      X_dyn_0 = stan_data$X_dyn_0,
                                      X_dyn_1 = stan_data$X_dyn_1))
  dims <- list(
    n_trans_coef = stan_data$n_trans_coef,
    n_dyn_coef_0 = stan_data$n_dyn_coef_0,
    n_dyn_coef_1 = stan_data$n_dyn_coef_1
  )
  re_info <- function() {
    list(
      trans_re_info = build_re_info_for_cpp(stan_data$trans_re_multi),
      dyn_re_info_0 = build_re_info_for_cpp(stan_data$dyn_re_multi_0),
      dyn_re_info_1 = build_re_info_for_cpp(stan_data$dyn_re_multi_1)
    )
  }
  target <- function(kind, data, args) {
    list(kind = kind, n_params = as.integer(stan_data$n_params),
         data = data, args = args)
  }

  if (stan_data$k_phases > 2) {
    return(target(
      "k_phase",
      c(base_data, list(X_trans_list = stan_data$X_trans_list,
                        X_dyn_list = stan_data$X_dyn_list)),
      list(n_dyn_coef = as.integer(stan_data$n_dyn_coef),
           n_trans_coef = as.integer(stan_data$n_trans_coef),
           k_phases = as.integer(stan_data$k_phases))
    ))
  }

  if (isTRUE(stan_data$is_glm)) {
    return(target(
      "glm", two_phase_data,
      c(dims, list(family_0 = stan_data$family_0,
                   family_1 = stan_data$family_1,
                   n_trials_0 = stan_data$n_trials_0,
                   n_trials_1 = stan_data$n_trials_1))
    ))
  }

  if (isTRUE(stan_data$has_re)) {
    if (isTRUE(stan_data$has_slopes)) {
      return(target("slopes", two_phase_data, c(dims, re_info())))
    }
    has_multi_re <- (stan_data$trans_re_multi$n_re_terms > 1) ||
                    (stan_data$dyn_re_multi_0$n_re_terms > 1) ||
                    (stan_data$dyn_re_multi_1$n_re_terms > 1)
    if (has_multi_re) {
      return(target("multi_re", two_phase_data, c(dims, re_info())))
    }
    return(target(
      "re",
      c(two_phase_data, list(trans_re_idx = stan_data$trans_re_idx,
                             dyn_re_idx_0 = stan_data$dyn_re_idx_0,
                             dyn_re_idx_1 = stan_data$dyn_re_idx_1)),
      c(dims, list(n_trans_re = stan_data$n_trans_re,
                   n_dyn_re_0 = stan_data$n_dyn_re_0,
                   n_dyn_re_1 = stan_data$n_dyn_re_1,
                   has_trans_re = stan_data$has_trans_re,
                   has_dyn_re_0 = stan_data$has_dyn_re_0,
                   has_dyn_re_1 = stan_data$has_dyn_re_1))
    ))
  }

  target("two_phase", two_phase_data, dims)
}


#' Build RE info structure for C++
#'
#' Converts R multi-RE structure to format expected by C++ sampler
#'
#' @param re_multi Multi-RE structure from build_multi_re()
#' @return List with n_terms, n_groups_vec, idx_list, and slope info
#' @keywords internal
build_re_info_for_cpp <- function(re_multi) {
  n_terms <- re_multi$n_re_terms

  if (n_terms == 0) {
    return(list(
      n_terms = 0L,
      n_groups_vec = integer(0),
      n_coefs_vec = integer(0),
      idx_list = list(),
      Z_list = list(),
      correlated_vec = logical(0),
      has_slopes = FALSE
    ))
  }

  n_groups_vec <- vapply(re_multi$re_list, function(x) x$n_groups, integer(1))
  n_coefs_vec <- vapply(re_multi$re_list, function(x) x$n_coefs, integer(1))
  idx_list <- lapply(re_multi$re_list, function(x) x$idx)
  Z_list <- lapply(re_multi$re_list, function(x) {
    if (is.null(x$Z)) {
      # For intercept-only, create a matrix of 1s (will be ignored but needed)
      NULL
    } else {
      x$Z
    }
  })
  correlated_vec <- vapply(re_multi$re_list, function(x) x$correlated, logical(1))

  list(
    n_terms = as.integer(n_terms),
    n_groups_vec = as.integer(n_groups_vec),
    n_coefs_vec = as.integer(n_coefs_vec),
    idx_list = idx_list,
    Z_list = Z_list,
    correlated_vec = correlated_vec,
    has_slopes = re_multi$has_slopes
  )
}


#' Combine diagnostics from multiple chains
#'
#' @param chain_list List of per-chain results of `phaseR_nuts_sampler()`
#' @param max_treedepth Maximum tree depth the chains ran with
#' @return List with the post-warmup number of divergent transitions, the
#'   number of transitions that hit `max_treedepth`, the mean acceptance
#'   statistic, and per chain the step size, mean number of leapfrog steps,
#'   E-BFMI and adapted inverse metric
#' @keywords internal
combine_diagnostics <- function(chain_list, max_treedepth) {
  ebfmi <- function(energy) {
    sum(diff(energy)^2) / sum((energy - mean(energy))^2)
  }
  list(
    n_divergent = sum(vapply(chain_list, function(x) x$n_divergent, integer(1))),
    n_max_treedepth = sum(vapply(chain_list, function(x) {
      sum(x$treedepth >= max_treedepth)
    }, integer(1))),
    mean_accept_prob = mean(vapply(chain_list, function(x) x$accept_prob, numeric(1))),
    step_sizes = vapply(chain_list, function(x) x$step_size, numeric(1)),
    mean_leapfrog = vapply(chain_list, function(x) mean(x$n_leapfrog), numeric(1)),
    ebfmi = vapply(chain_list, function(x) ebfmi(x$energy), numeric(1)),
    inv_metric = lapply(chain_list, function(x) x$inv_metric)
  )
}
