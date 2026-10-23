#' Laplace approximation backend for phase models
#'
#' Fast approximate inference using the Laplace approximation.
#' Finds the posterior mode and approximates the posterior as
#' Gaussian with covariance equal to the inverse Hessian.
#'
#' @param stan_data Prepared data list
#' @param n_samples Number of samples to draw from Gaussian approximation
#' @param seed Integer seed for the draws from the Gaussian approximation
#' @param ... Additional arguments (ignored)
#'
#' @return Raw fit object with same structure as HMC backend
#' @keywords internal
#'
fit_laplace <- function(stan_data, n_samples = 1000, seed, ...) {

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
  set.seed(seed)
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
#' Creates a function that computes -log p(theta | y) for use with optim(),
#' from the same model description and prior the NUTS backend samples.
#'
#' @param stan_data Prepared data list
#' @return Function of parameter vector returning negative log posterior
#' @keywords internal
#'
build_neg_log_posterior <- function(stan_data) {
  target <- posterior_target(stan_data)
  function(params) -phaseR_log_posterior(target, params)
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
