#' Describe the log posterior of a phase model
#'
#' Both backends evaluate the log posterior through this description: the
#' NUTS sampler in `phaseR_nuts_sampler()` and the Laplace approximation
#' through `phaseR_log_posterior()`.
#'
#' @param stan_data Prepared data list
#' @return List with `kind` (which C++ log posterior to use), `n_params`,
#'   `data`, `args` (the dimensions and random-effect structures that log
#'   posterior takes) and `prior` (the prior of each parameter class, from
#'   [prior_to_cpp()])
#' @keywords internal
posterior_target <- function(stan_data) {

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
         data = data, args = args, prior = stan_data$prior_cpp)
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
