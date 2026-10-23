#' Describe the log posterior of a phase model
#'
#' Both backends evaluate the log posterior through this description: the
#' NUTS sampler in `phaseR_nuts_sampler()` and the Laplace approximation.
#' A model with `K` phases has `K - 1` transition linear predictors (the logit
#' of moving from phase `k` to `k + 1`) and `K` dynamics linear predictors
#' (the mean of each phase on the scale of its link). Each is a fixed-effect
#' design matrix plus random-effect terms, whose parameters sit at the
#' offsets of the parameter vector that `prepare_stan_data()` names.
#'
#' @param stan_data Prepared data list
#' @return List with `kind = "phase"`, `n_params`, `n_phases`, the outcome
#'   `y`, the 0-based observation range of each unit (`unit_start`,
#'   `unit_end`), the outcome family code and number of trials of each phase,
#'   the 0-based offset of the log residual SDs (`sigma_offset`), the linear
#'   predictors (`transitions`, `dynamics`; each a list with design matrix
#'   `X`, 0-based coefficient `offset` and random-effect terms `re`) and the
#'   prior of each parameter class (`prior`, from [prior_to_cpp()])
#' @keywords internal
posterior_target <- function(stan_data) {

  k <- stan_data$k_phases
  if (k > 2) {
    X_trans <- stan_data$X_trans_list
    X_dyn <- stan_data$X_dyn_list
    family <- rep(0L, k)
    re <- rep(list(NULL), 2L * k - 1L)
  } else {
    X_trans <- list(stan_data$X_trans)
    X_dyn <- list(stan_data$X_dyn_0, stan_data$X_dyn_1)
    family <- c(stan_data$family_0, stan_data$family_1)
    re <- list(stan_data$trans_re_multi, stan_data$dyn_re_multi_0,
               stan_data$dyn_re_multi_1)
  }

  unit_size <- stan_data$unit_end - stan_data$unit_start + 1L
  unit_first_obs <- rep(stan_data$unit_start + 1L, unit_size)

  offset <- 0L
  predictor <- function(X) {
    out <- list(X = X, offset = offset, re = list())
    offset <<- offset + ncol(X)
    out
  }
  predictors <- lapply(c(X_trans, X_dyn), predictor)
  sigma_offset <- offset
  offset <- offset + k

  re_term <- function(term) {
    kind <- if (term$n_coefs == 1L) 0L else if (term$correlated) 2L else 1L
    n_hyper <- switch(kind + 1L, 1L, term$n_coefs,
                      (term$n_coefs * (term$n_coefs + 1L)) %/% 2L)
    out <- list(
      kind = kind,
      n_groups = as.integer(term$n_groups),
      n_coefs = as.integer(term$n_coefs),
      offset = offset,
      group = as.integer(term$idx[unit_first_obs]) - 1L,
      Z = if (term$n_coefs > 1L) term$Z else NULL
    )
    offset <<- offset + term$n_groups * term$n_coefs + n_hyper
    out
  }
  for (c in seq_along(predictors)) {
    if (!is.null(re[[c]]) && re[[c]]$n_re_terms > 0) {
      predictors[[c]]$re <- lapply(unname(re[[c]]$re_list), re_term)
    }
  }
  if (offset != stan_data$n_params) {
    stop("The model reads ", offset, " parameters but names ",
         stan_data$n_params, ".", call. = FALSE)
  }

  list(
    kind = "phase",
    n_params = as.integer(stan_data$n_params),
    n_phases = as.integer(k),
    y = as.numeric(stan_data$y),
    unit_start = as.integer(stan_data$unit_start),
    unit_end = as.integer(stan_data$unit_end),
    family = as.integer(family),
    n_trials = rep(1L, k),
    sigma_offset = sigma_offset,
    transitions = predictors[seq_len(k - 1L)],
    dynamics = predictors[k - 1L + seq_len(k)],
    prior = stan_data$prior_cpp
  )
}
