#' Laplace approximation backend for phase models
#'
#' Approximates the posterior by a Gaussian at its mode. In a model with
#' random effects the effects u are integrated out first: for the other
#' parameters theta (coefficients, log SDs, Cholesky factors) the marginal
#' posterior is approximated by the Laplace formula
#' `p(theta | y) ~ p(theta, u*(theta) | y) (2 pi)^(n_u / 2) det(A)^(-1/2)`,
#' with `u*(theta)` the mode of the effects given theta and `A` the negative
#' Hessian of the log posterior in u there (Tierney and Kadane 1986; the
#' approach of INLA and TMB). The joint mode of (theta, u) does not exist:
#' as a random-effect SD shrinks with its effects at zero, the Normal
#' density of the effects grows without bound, faster than any proper SD
#' prior falls. The marginal posterior of theta is bounded, so its mode is
#' well defined. theta is drawn from the Gaussian at that mode, with
#' covariance the inverse negative Hessian of the log marginal (central
#' differences of its exact gradient), and u given theta from
#' `u*(theta_hat) + du*/dtheta (theta - theta_hat) + N(0, A^-1)`. In a model
#' without random effects this is the Gaussian approximation at the joint
#' mode.
#'
#' The log posterior of a phase model is multimodal: a mode in which a phase
#' is never entered leaves that phase's parameters at their prior. The mode
#' is therefore searched from `n_starts` starting points built from
#' provisional phase assignments ([initial_values()]), and the start reaching
#' the highest log marginal posterior is kept.
#'
#' @param stan_data Prepared data list
#' @param n_samples Number of samples to draw from Gaussian approximation
#' @param seed Integer seed for the random starting points and the draws
#'   from the Gaussian approximation
#' @param n_starts Number of starting points of the mode search
#' @param ... Additional arguments (ignored)
#'
#' @return Raw fit object with same structure as HMC backend
#' @keywords internal
#'
fit_laplace <- function(stan_data, n_samples = 1000, seed, n_starts = 5L,
                        ...) {

  n_params <- stan_data$n_params
  target <- posterior_target(stan_data)
  laplace <- phaseR_laplace(target)
  layout <- phaseR_laplace_layout(laplace)
  objective <- laplace_objective(laplace)

  starts <- with_seed(seed, initial_values(target, n_starts))
  runs <- lapply(starts, function(init) {
    stats::optim(
      par = init[layout$theta],
      fn = objective$fn,
      gr = objective$gr,
      method = "L-BFGS-B",
      control = list(maxit = 1000, factr = 1e7)
    )
  })
  values <- vapply(runs, function(r) r$value, numeric(1))
  opt_result <- runs[[which.min(values)]]
  theta_hat <- opt_result$par

  if (opt_result$convergence != 0) {
    warning("Optimization did not converge: ", opt_result$message,
            call. = FALSE)
  }
  if (!phaseR_laplace_objective(laplace, theta_hat)$converged) {
    warning("The mode of the random effects was not reached to full ",
            "precision at the posterior mode.", call. = FALSE)
  }

  # Gaussian approximation of theta
  hessian <- -phaseR_laplace_hessian(laplace, theta_hat)
  R_theta <- tryCatch(chol(hessian), error = function(e) {
    stop("The negative Hessian of the log posterior is not positive ",
         "definite at the mode found, so the Laplace approximation is not ",
         "defined there. Use backend = \"hmc\".", call. = FALSE)
  })
  cov_theta <- chol2inv(R_theta)

  # Effects given theta, linearised in theta
  cond <- phaseR_laplace_conditional(laplace, theta_hat)
  n_u <- length(layout$u)

  mode <- numeric(n_params)
  mode[layout$theta] <- theta_hat
  mode[layout$u] <- cond$mode
  names(mode) <- stan_data$param_names

  cov_matrix <- matrix(0, n_params, n_params,
                       dimnames = list(stan_data$param_names,
                                       stan_data$param_names))
  cov_matrix[layout$theta, layout$theta] <- cov_theta
  draws <- with_seed(seed, {
    z_theta <- matrix(stats::rnorm(n_samples * length(theta_hat)),
                      nrow = n_samples)
    theta_dev <- t(backsolve(R_theta, t(z_theta)))
    out <- matrix(0, n_samples, n_params)
    out[, layout$theta] <- sweep(theta_dev, 2, theta_hat, "+")
    if (n_u > 0) {
      R_u <- chol(cond$precision)
      z_u <- matrix(stats::rnorm(n_samples * n_u), nrow = n_samples)
      out[, layout$u] <- sweep(theta_dev %*% t(cond$slope) +
                                 t(backsolve(R_u, t(z_u))),
                               2, cond$mode, "+")
    }
    out
  })
  colnames(draws) <- stan_data$param_names
  if (n_u > 0) {
    cov_u_theta <- cond$slope %*% cov_theta
    cov_matrix[layout$u, layout$theta] <- cov_u_theta
    cov_matrix[layout$theta, layout$u] <- t(cov_u_theta)
    cov_matrix[layout$u, layout$u] <- chol2inv(chol(cond$precision)) +
      cov_u_theta %*% t(cond$slope)
  }
  dimnames(hessian) <- list(stan_data$param_names[layout$theta],
                            stan_data$param_names[layout$theta])

  list(
    draws = draws,
    chain_id = rep(1L, n_samples),
    diagnostics = list(
      convergence = opt_result$convergence,
      value = opt_result$value,
      counts = opt_result$counts,
      start_values = -values,
      mode = mode,
      se = sqrt(diag(cov_matrix)),
      hessian = hessian,
      cov = cov_matrix,
      marginalised = stan_data$param_names[layout$u]
    ),
    n_params = n_params,
    param_names = stan_data$param_names,
    k_phases = stan_data$k_phases
  )
}


#' Negative Laplace log marginal posterior for optimization
#'
#' @param laplace External pointer from `phaseR_laplace()`
#' @return List of two functions of the parameters outside the random
#'   effects: `fn`, the negative log marginal posterior with the random
#'   effects integrated out by the Laplace approximation, and `gr`, its exact
#'   gradient. The C++ object caches the last evaluation, so `fn` and `gr`
#'   at the same point cost one evaluation.
#' @keywords internal
laplace_objective <- function(laplace) {
  list(
    fn = function(theta) -phaseR_laplace_objective(laplace, theta)$value,
    gr = function(theta) -phaseR_laplace_objective(laplace, theta)$gradient
  )
}


#' Build negative log posterior function for optimization
#'
#' Creates -log p(theta | y) and its exact gradient for use with optim(),
#' from the same model description and prior the NUTS backend samples.
#'
#' @param stan_data Prepared data list
#' @return List of two functions of the parameter vector: `fn`, the negative
#'   log posterior, and `gr`, its gradient
#' @keywords internal
#'
build_neg_log_posterior <- function(stan_data) {
  density <- phaseR_density(posterior_target(stan_data))
  list(
    fn = function(params) -phaseR_log_density(density, params),
    gr = function(params) -phaseR_log_density_gradient(density, params)$gradient
  )
}
