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
#' Given theta the log posterior splits over blocks of effects that share no
#' unit (one block per unit for unit-level effects, one block per level of a
#' grouping factor nested in the units), so `A` is block diagonal and the
#' mode of the effects is found by Newton's method from `u = 0` in each block
#' separately. When the effects of a block have several modes given theta,
#' for example unit-level effects in phases whose timing is uncertain, the
#' approximation can become degenerate where two modes merge; the fit then
#' stops or, if several modes remain at the optimum, warns, and the HMC
#' backend is the reliable choice.
#'
#' The log posterior of a phase model is multimodal: a mode in which a phase
#' is never entered leaves that phase's parameters at their prior. The mode
#' is therefore searched from `n_starts` starting points built from
#' provisional phase assignments ([initial_values()]), and the start reaching
#' the highest log marginal posterior is kept. With random effects, each
#' start is first moved to the mode of the fixed effects with the effects at
#' zero ([fixed_effects_mode()]), and starts that arrive at the same point
#' share one search.
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
  if (length(layout$u) > 0) {
    starts <- lapply(starts, fixed_effects_mode, target = target)
  }
  # Starts that reached the same point give the same run.
  run_of <- distinct_starts(lapply(starts, function(s) s[layout$theta]))
  # A search that steps where the approximation is not defined (the mode of
  # the effects is not a strict maximum or not finite) ends without a value.
  runs <- lapply(unique(run_of), function(i) {
    tryCatch(
      stats::optim(
        par = starts[[i]][layout$theta],
        fn = objective$fn,
        gr = objective$gr,
        method = "L-BFGS-B",
        control = list(maxit = 1000, factr = 1e7)
      ),
      error = function(e) list(value = Inf)
    )
  })
  run_values <- vapply(runs, function(r) r$value, numeric(1))
  if (!any(is.finite(run_values))) {
    stop("The Laplace search failed from every starting point: along it the ",
         "mode of the random effects given the other parameters became ",
         "degenerate or not finite, where the approximation is not defined. ",
         "This happens when the effects of a group have several modes, for ",
         "example unit-level effects in phases whose timing is uncertain. ",
         "Use backend = \"hmc\".", call. = FALSE)
  }
  values <- run_values[match(run_of, unique(run_of))]
  opt_result <- runs[[which.min(run_values)]]
  theta_hat <- opt_result$par

  if (opt_result$convergence != 0) {
    warning("Optimization did not converge: ", opt_result$message,
            call. = FALSE)
  }
  if (!phaseR_laplace_objective(laplace, theta_hat)$converged) {
    warning("The mode of the random effects was not reached to full ",
            "precision at the posterior mode.", call. = FALSE)
  }
  n_multimodal <- phaseR_laplace_multimodal(laplace, theta_hat)
  if (n_multimodal > 0) {
    warning("At the posterior mode the random effects of ", n_multimodal,
            " groups have several modes given the other parameters; the ",
            "Laplace approximation keeps one of them. Compare with ",
            "backend = \"hmc\".", call. = FALSE)
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

  draws <- with_seed(seed, {
    z_theta <- matrix(stats::rnorm(n_samples * length(theta_hat)),
                      nrow = n_samples)
    theta_dev <- t(backsolve(R_theta, t(z_theta)))
    out <- matrix(0, n_samples, n_params)
    out[, layout$theta] <- sweep(theta_dev, 2, theta_hat, "+")
    if (n_u > 0) {
      z_u <- matrix(stats::rnorm(n_samples * n_u), nrow = n_samples)
      noise <- phaseR_laplace_effect_noise(laplace, theta_hat, z_u)
      out[, layout$u] <- sweep(theta_dev %*% t(cond$slope) + noise,
                               2, cond$mode, "+")
    }
    out
  })
  colnames(draws) <- stan_data$param_names
  se <- numeric(n_params)
  se[layout$theta] <- sqrt(diag(cov_theta))
  if (n_u > 0) {
    se[layout$u] <- sqrt(cond$variance +
                           rowSums((cond$slope %*% cov_theta) * cond$slope))
  }
  names(se) <- stan_data$param_names
  dimnames(cov_theta) <- list(stan_data$param_names[layout$theta],
                              stan_data$param_names[layout$theta])
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
      se = se,
      hessian = hessian,
      cov = cov_theta,
      marginalised = stan_data$param_names[layout$u],
      multimodal_groups = n_multimodal
    ),
    n_params = n_params,
    param_names = stan_data$param_names,
    k_phases = stan_data$k_phases
  )
}


#' Mode of the fixed effects with the random effects at zero
#'
#' Maximises the log posterior over the coefficients and residual SDs with
#' every random effect at zero and its covariance parameters at their
#' starting values. Each step costs one evaluation of the log posterior and
#' its gradient instead of an inner mode search, and it brings a starting
#' point of the Laplace search close to the mode of the phase structure it
#' leads to.
#'
#' @param start Parameter vector from [initial_values()]
#' @param target Model description from [posterior_target()]
#' @return `start` with the coefficients and log residual SDs replaced by
#'   their mode
#' @keywords internal
fixed_effects_mode <- function(start, target) {
  free <- seq_len(target$sigma_offset + target$n_phases)
  density <- phaseR_density(target)
  full <- function(x) replace(start, free, x)
  fit <- stats::optim(
    par = start[free],
    fn = function(x) -phaseR_log_density(density, full(x)),
    gr = function(x) -phaseR_log_density_gradient(density, full(x))$gradient[free],
    method = "L-BFGS-B",
    control = list(maxit = 1000, factr = 1e7)
  )
  full(fit$par)
}


#' Index of the first equal start for each start
#'
#' @param starts List of parameter vectors
#' @param tol Relative tolerance
#' @return Integer vector: for each start, the index of the first start
#'   within `tol` of it
#' @keywords internal
distinct_starts <- function(starts, tol = 1e-3) {
  out <- seq_along(starts)
  for (i in seq_along(starts)[-1]) {
    for (j in seq_len(i - 1L)) {
      if (out[j] == j &&
          max(abs(starts[[i]] - starts[[j]])) <= tol * (1 + max(abs(starts[[j]])))) {
        out[i] <- j
        break
      }
    }
  }
  out
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
