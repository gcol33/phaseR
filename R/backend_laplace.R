#' Laplace approximation backend for phase models
#'
#' Fast approximate inference using the Laplace approximation.
#' Finds the posterior mode and approximates the posterior as
#' Gaussian with covariance equal to the inverse Hessian.
#'
#' @param stan_data Prepared data list
#' @param n_samples Number of samples to draw from Gaussian approximation
#' @param ... Additional arguments (ignored)
#'
#' @return Raw fit object with same structure as HMC backend
#' @keywords internal
#'
fit_laplace <- function(stan_data, n_samples = 1000, ...) {

  n_params <- stan_data$n_params

  # Build objective function (negative log posterior)
  neg_log_post <- build_neg_log_posterior(stan_data)

  # Initial values
  init <- rep(0, n_params)

  # Find the mode using L-BFGS-B
  opt_result <- stats::optim(
    par = init,
    fn = neg_log_post,
    method = "L-BFGS-B",
    control = list(maxit = 1000, factr = 1e7),
    hessian = TRUE
  )

  if (opt_result$convergence != 0) {
    warning("Optimization did not converge: ", opt_result$message)
  }

  # Extract mode and Hessian
  mode <- opt_result$par
  hessian <- opt_result$hessian

  # Compute covariance matrix (inverse Hessian)
  # Add small ridge for numerical stability
  ridge <- 1e-6 * diag(n_params)
  cov_matrix <- tryCatch({
    solve(hessian + ridge)
  }, error = function(e) {
    warning("Hessian not invertible, using pseudo-inverse")
    MASS_ginv(hessian + ridge)
  })

  # Ensure covariance is symmetric and positive definite
  cov_matrix <- (cov_matrix + t(cov_matrix)) / 2
  eigenvals <- eigen(cov_matrix, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigenvals <= 0)) {
    # Make positive definite
    min_eig <- min(eigenvals)
    cov_matrix <- cov_matrix + (-min_eig + 1e-6) * diag(n_params)
  }

  # Sample from Gaussian approximation
  draws <- sample_mvnorm(n_samples, mode, cov_matrix)
  colnames(draws) <- stan_data$param_names

  # Compute standard errors
  se <- sqrt(diag(cov_matrix))

  list(
    draws = draws,
    chain_id = rep(1L, n_samples),
    diagnostics = list(
      convergence = opt_result$convergence,
      value = opt_result$value,
      counts = opt_result$counts,
      mode = mode,
      se = se,
      hessian = hessian,
      cov = cov_matrix
    ),
    n_params = n_params,
    param_names = stan_data$param_names,
    k_phases = stan_data$k_phases
  )
}


#' Build negative log posterior function for optimization
#'
#' Creates a function that computes -log p(theta | y) for use with optim().
#'
#' @param stan_data Prepared data list
#' @return Function of parameter vector returning negative log posterior
#' @keywords internal
#'
build_neg_log_posterior <- function(stan_data) {

  k_phases <- stan_data$k_phases

  if (k_phases > 2) {
    # k-phase model
    function(params) {
      -phase_log_posterior_k(
        params = params,
        data = list(
          id = stan_data$id,
          time = stan_data$time,
          y = stan_data$y,
          X_trans_list = stan_data$X_trans_list,
          X_dyn_list = stan_data$X_dyn_list,
          n_units = stan_data$n_units,
          unit_start = stan_data$unit_start,
          unit_end = stan_data$unit_end
        ),
        n_dyn_coef = as.integer(stan_data$n_dyn_coef),
        n_trans_coef = as.integer(stan_data$n_trans_coef),
        k_phases = as.integer(k_phases)
      )
    }
  } else if (isTRUE(stan_data$is_glm)) {
    # 2-phase GLM
    function(params) {
      -phase_log_posterior_glm(
        params = params,
        data = list(
          id = stan_data$id,
          time = stan_data$time,
          y = stan_data$y,
          X_trans = stan_data$X_trans,
          X_dyn = stan_data$X_dyn,
          n_units = stan_data$n_units,
          unit_start = stan_data$unit_start,
          unit_end = stan_data$unit_end
        ),
        n_trans_coef = stan_data$n_trans_coef,
        n_dyn_coef_0 = stan_data$n_dyn_coef_0,
        n_dyn_coef_1 = stan_data$n_dyn_coef_1,
        family_0 = stan_data$family_0,
        family_1 = stan_data$family_1,
        n_trials_0 = stan_data$n_trials_0,
        n_trials_1 = stan_data$n_trials_1
      )
    }
  } else if (isTRUE(stan_data$has_slopes)) {
    # Random slopes model
    function(params) {
      -phase_log_posterior_slopes(
        params = params,
        data = list(
          id = stan_data$id,
          time = stan_data$time,
          y = stan_data$y,
          X_trans = stan_data$X_trans,
          X_dyn = stan_data$X_dyn,
          n_units = stan_data$n_units,
          unit_start = stan_data$unit_start,
          unit_end = stan_data$unit_end
        ),
        n_trans_coef = stan_data$n_trans_coef,
        n_dyn_coef_0 = stan_data$n_dyn_coef_0,
        n_dyn_coef_1 = stan_data$n_dyn_coef_1,
        trans_re_info = build_re_info_for_cpp(stan_data$trans_re_multi),
        dyn_re_info_0 = build_re_info_for_cpp(stan_data$dyn_re_multi_0),
        dyn_re_info_1 = build_re_info_for_cpp(stan_data$dyn_re_multi_1)
      )
    }
  } else if (isTRUE(stan_data$has_re)) {
    # Check for multi-RE
    has_multi_re <- (stan_data$trans_re_multi$n_re_terms > 1) ||
                    (stan_data$dyn_re_multi_0$n_re_terms > 1) ||
                    (stan_data$dyn_re_multi_1$n_re_terms > 1)

    if (has_multi_re) {
      # Multi-RE model
      function(params) {
        -phase_log_posterior_multi_re(
          params = params,
          data = list(
            id = stan_data$id,
            time = stan_data$time,
            y = stan_data$y,
            X_trans = stan_data$X_trans,
            X_dyn = stan_data$X_dyn,
            n_units = stan_data$n_units,
            unit_start = stan_data$unit_start,
            unit_end = stan_data$unit_end
          ),
          n_trans_coef = stan_data$n_trans_coef,
          n_dyn_coef_0 = stan_data$n_dyn_coef_0,
          n_dyn_coef_1 = stan_data$n_dyn_coef_1,
          trans_re_info = build_re_info_for_cpp(stan_data$trans_re_multi),
          dyn_re_info_0 = build_re_info_for_cpp(stan_data$dyn_re_multi_0),
          dyn_re_info_1 = build_re_info_for_cpp(stan_data$dyn_re_multi_1)
        )
      }
    } else {
      # Single-RE model with custom priors
      prior_cpp <- stan_data$prior_cpp
      function(params) {
        -phase_log_posterior_re_custom(
          params = params,
          data = list(
            id = stan_data$id,
            time = stan_data$time,
            y = stan_data$y,
            X_trans = stan_data$X_trans,
            X_dyn = stan_data$X_dyn,
            n_units = stan_data$n_units,
            unit_start = stan_data$unit_start,
            unit_end = stan_data$unit_end,
            trans_re_idx = stan_data$trans_re_idx,
            dyn_re_idx_0 = stan_data$dyn_re_idx_0,
            dyn_re_idx_1 = stan_data$dyn_re_idx_1
          ),
          n_trans_coef = stan_data$n_trans_coef,
          n_dyn_coef_0 = stan_data$n_dyn_coef_0,
          n_dyn_coef_1 = stan_data$n_dyn_coef_1,
          n_trans_re = stan_data$n_trans_re,
          n_dyn_re_0 = stan_data$n_dyn_re_0,
          n_dyn_re_1 = stan_data$n_dyn_re_1,
          has_trans_re = stan_data$has_trans_re,
          has_dyn_re_0 = stan_data$has_dyn_re_0,
          has_dyn_re_1 = stan_data$has_dyn_re_1,
          beta_type = prior_cpp$beta_type,
          beta_params = prior_cpp$beta_params,
          sigma_type = prior_cpp$sigma_type,
          sigma_params = prior_cpp$sigma_params,
          re_type = prior_cpp$re_type,
          re_params = prior_cpp$re_params
        )
      }
    }
  } else {
    # 2-phase Gaussian without RE - use custom priors
    prior_cpp <- stan_data$prior_cpp
    function(params) {
      -phase_log_posterior_custom(
        params = params,
        data = list(
          id = stan_data$id,
          time = stan_data$time,
          y = stan_data$y,
          X_trans = stan_data$X_trans,
          X_dyn = stan_data$X_dyn,
          n_units = stan_data$n_units,
          unit_start = stan_data$unit_start,
          unit_end = stan_data$unit_end
        ),
        n_trans_coef = stan_data$n_trans_coef,
        n_dyn_coef_0 = stan_data$n_dyn_coef_0,
        n_dyn_coef_1 = stan_data$n_dyn_coef_1,
        beta_type = prior_cpp$beta_type,
        beta_params = prior_cpp$beta_params,
        sigma_type = prior_cpp$sigma_type,
        sigma_params = prior_cpp$sigma_params
      )
    }
  }
}


#' Sample from multivariate normal
#'
#' @param n Number of samples
#' @param mu Mean vector
#' @param sigma Covariance matrix
#' @return Matrix of samples (n x length(mu))
#' @keywords internal
#'
sample_mvnorm <- function(n, mu, sigma) {
  p <- length(mu)

  # Cholesky decomposition
  L <- tryCatch({
    chol(sigma)
  }, error = function(e) {
    # If Cholesky fails, use eigendecomposition
    eig <- eigen(sigma, symmetric = TRUE)
    eig$vectors %*% diag(sqrt(pmax(eig$values, 0))) %*% t(eig$vectors)
  })

  # Generate standard normal samples
  z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)

  # Transform to target distribution
  # For Cholesky L (upper triangular): X = Z %*% L + mu
  # For eigendecomposition: same formula works
  t(t(z %*% L) + mu)
}


#' Generalized inverse (Moore-Penrose pseudo-inverse)
#'
#' Simple implementation without MASS dependency.
#'
#' @param X Matrix to invert
#' @param tol Tolerance for singular values
#' @return Pseudo-inverse of X
#' @keywords internal
#'
MASS_ginv <- function(X, tol = sqrt(.Machine$double.eps)) {
  svd_X <- svd(X)
  pos <- svd_X$d > max(tol * svd_X$d[1], 0)
  if (all(pos)) {
    svd_X$v %*% (1 / svd_X$d * t(svd_X$u))
  } else if (!any(pos)) {
    array(0, dim(X)[2:1])
  } else {
    svd_X$v[, pos, drop = FALSE] %*% ((1 / svd_X$d[pos]) *
      t(svd_X$u[, pos, drop = FALSE]))
  }
}
