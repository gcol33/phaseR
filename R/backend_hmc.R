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

  target <- posterior_target(stan_data)
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
