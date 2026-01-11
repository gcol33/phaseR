#' @export
print.phase_fit <- function(x, ...) {

  cat("phaseR model fit\n")
  cat(sprintf("  Phases: %s\n", paste(x$model$phases$names, collapse = " -> ")))
  cat(sprintf("  Units: %d\n", x$n_units))
  cat(sprintf("  Time points: %d\n", x$n_times))
  cat(sprintf("  Backend: %s\n", x$backend))
  cat(sprintf("  Chains: %d, Iterations: %d (warmup: %d)\n",
              x$n_chains, x$n_iter, x$n_warmup))

  if (x$backend == "hmc" && !is.null(x$diagnostics)) {
    cat(sprintf("  Divergences: %d\n", x$diagnostics$n_divergent))
  }

  cat("\nUse summary() for parameter estimates, plot() for diagnostics.\n")
  invisible(x)
}


#' @export
summary.phase_fit <- function(object, prob = 0.95, ...) {

  draws <- object$draws
  param_names <- object$param_names

  alpha <- (1 - prob) / 2
  probs <- c(alpha, 0.5, 1 - alpha)

  # Compute summaries
  summaries <- t(apply(draws, 2, function(x) {
    c(
      mean = mean(x),
      sd = sd(x),
      quantile(x, probs)
    )
  }))

  colnames(summaries) <- c("mean", "sd",
                           sprintf("%.1f%%", probs[1] * 100),
                           "median",
                           sprintf("%.1f%%", probs[3] * 100))
  rownames(summaries) <- param_names

  # Add Rhat and ESS if multiple chains
  if (object$n_chains > 1) {
    summaries <- cbind(summaries,
                       Rhat = compute_rhat(draws, object$chain_id),
                       ESS = compute_ess(draws, object$chain_id))
  }

  structure(
    list(
      summaries = summaries,
      model = object$model,
      prob = prob
    ),
    class = "summary.phase_fit"
  )
}


#' @export
print.summary.phase_fit <- function(x, digits = 2, ...) {

  cat("Phase Model Summary\n")
  cat(sprintf("Phases: %s\n\n", paste(x$model$phases$names, collapse = " -> ")))

  # Print transition coefficients
  cat("Transition coefficients:\n")
  trans_idx <- grep("^beta_trans", rownames(x$summaries))
  if (length(trans_idx) > 0) {
    print(round(x$summaries[trans_idx, , drop = FALSE], digits))
  }

  # Print dynamics coefficients by phase
  for (phase in x$model$phases$names) {
    cat(sprintf("\nDynamics coefficients (%s):\n", phase))
    phase_pattern <- sprintf("^beta_%s_", phase)
    phase_idx <- grep(phase_pattern, rownames(x$summaries))
    if (length(phase_idx) > 0) {
      print(round(x$summaries[phase_idx, , drop = FALSE], digits))
    }
  }

  # Print variance parameters
  cat("\nVariance parameters:\n")
  var_idx <- grep("^log_sigma", rownames(x$summaries))
  if (length(var_idx) > 0) {
    print(round(x$summaries[var_idx, , drop = FALSE], digits))
  }

  invisible(x)
}


#' @export
plot.phase_fit <- function(x, type = c("trace", "density", "both"),
                           pars = NULL, ...) {

  type <- match.arg(type)

  # Select parameters
  if (is.null(pars)) {
    # Default: main parameters, not auxiliary
    pars <- grep("^(beta|log_sigma)", x$param_names, value = TRUE)
  }

  par_idx <- which(x$param_names %in% pars)
  draws <- x$draws[, par_idx, drop = FALSE]

  # Base R plotting
  plot_base_diagnostics(draws, x$chain_id, type, pars)
}


#' Base R diagnostic plots
#' @keywords internal
plot_base_diagnostics <- function(draws, chain_id, type, pars) {

  n_pars <- ncol(draws)
  n_chains <- length(unique(chain_id))

  if (type %in% c("trace", "both")) {
    # Trace plots
    old_par <- par(mfrow = c(min(n_pars, 4), 1), mar = c(2, 4, 2, 1))
    on.exit(par(old_par))

    colors <- rainbow(n_chains)

    for (i in seq_len(min(n_pars, 4))) {
      plot(draws[, i], type = "n", main = colnames(draws)[i],
           xlab = "Iteration", ylab = "Value")
      for (chain in seq_len(n_chains)) {
        chain_idx <- which(chain_id == chain)
        lines(chain_idx, draws[chain_idx, i], col = colors[chain])
      }
    }
  }

  if (type == "density") {
    old_par <- par(mfrow = c(min(n_pars, 4), 1), mar = c(2, 4, 2, 1))
    on.exit(par(old_par))

    for (i in seq_len(min(n_pars, 4))) {
      d <- density(draws[, i])
      plot(d, main = colnames(draws)[i], xlab = "Value")
    }
  }
}


#' Compute R-hat convergence diagnostic
#' @keywords internal
compute_rhat <- function(draws, chain_id) {

  n_params <- ncol(draws)
  chains <- unique(chain_id)
  n_chains <- length(chains)

  rhats <- numeric(n_params)

  for (j in seq_len(n_params)) {
    chain_means <- numeric(n_chains)
    chain_vars <- numeric(n_chains)
    chain_lengths <- numeric(n_chains)

    for (i in seq_along(chains)) {
      chain_draws <- draws[chain_id == chains[i], j]
      chain_means[i] <- mean(chain_draws)
      chain_vars[i] <- var(chain_draws)
      chain_lengths[i] <- length(chain_draws)
    }

    n <- mean(chain_lengths)
    grand_mean <- mean(chain_means)

    # Between-chain variance
    B <- n * var(chain_means)

    # Within-chain variance
    W <- mean(chain_vars)

    # Pooled variance estimate
    var_hat <- ((n - 1) / n) * W + (1 / n) * B

    rhats[j] <- sqrt(var_hat / W)
  }

  rhats
}


#' Compute effective sample size
#' @keywords internal
compute_ess <- function(draws, chain_id) {

  n_params <- ncol(draws)
  n_total <- nrow(draws)

  ess <- numeric(n_params)

  for (j in seq_len(n_params)) {
    # Simple ESS based on autocorrelation
    x <- draws[, j]
    ac <- acf(x, plot = FALSE, lag.max = min(100, n_total / 2))$acf[-1]

    # Find cutoff where autocorrelation goes negative
    cutoff <- which(ac < 0)[1]
    if (is.na(cutoff)) cutoff <- length(ac)

    rho_sum <- sum(ac[seq_len(cutoff - 1)])
    ess[j] <- n_total / (1 + 2 * rho_sum)
  }

  pmax(1, round(ess))
}


#' Extract coefficients from phase_fit
#'
#' Returns posterior means of all parameters.
#'
#' @param object A `phase_fit` object
#' @param ... Additional arguments (unused)
#'
#' @return Named numeric vector of posterior means
#' @export
coef.phase_fit <- function(object, ...) {
  colMeans(object$draws)
}


#' Extract fitted values from phase_fit
#'
#' Returns predicted response values at each observation.
#'
#' @param object A `phase_fit` object
#' @param ... Additional arguments (unused)
#'
#' @return Numeric vector of fitted values
#' @export
fitted.phase_fit <- function(object, ...) {
  predict(object, type = "response")
}


#' Extract residuals from phase_fit
#'
#' Computes residuals (observed - fitted) marginalized over phase uncertainty.
#'
#' @param object A `phase_fit` object
#' @param type Type of residuals: "response" (default) or "pearson"
#' @param ... Additional arguments (unused)
#'
#' @return Numeric vector of residuals
#' @export
residuals.phase_fit <- function(object, type = c("response", "pearson"), ...) {

  type <- match.arg(type)

  y <- object$data[[object$model$dynamics[[1]]$response]]
  fitted_vals <- fitted(object)

  if (is.data.frame(fitted_vals)) {
    fitted_vals <- fitted_vals$mean
  }

  resid <- y - fitted_vals

  if (type == "pearson") {
    # Get residual SD (weighted average of phase-specific SDs)
    sigma_idx <- grep("^log_sigma_", object$param_names)
    sigmas <- exp(colMeans(object$draws[, sigma_idx, drop = FALSE]))

    # Use average sigma for simplicity
    sigma_avg <- mean(sigmas)
    resid <- resid / sigma_avg
  }

  resid
}


#' Number of observations in phase_fit
#'
#' @param object A `phase_fit` object
#' @param ... Additional arguments (unused)
#'
#' @return Integer count of observations
#' @export
nobs.phase_fit <- function(object, ...) {
  nrow(object$data)
}


#' Log-likelihood of phase_fit
#'
#' Returns the log-likelihood at the posterior mean.
#'
#' @param object A `phase_fit` object
#' @param ... Additional arguments (unused)
#'
#' @return A logLik object
#' @export
logLik.phase_fit <- function(object, ...) {

  # Get posterior mean parameters
  params <- colMeans(object$draws)

  # Compute log-likelihood at posterior mean
  y <- object$data[[object$model$dynamics[[1]]$response]]
  fitted_vals <- fitted(object)
  if (is.data.frame(fitted_vals)) {
    fitted_vals <- fitted_vals$mean
  }

  # Get sigma (average of phase-specific)
  sigma_idx <- grep("^log_sigma_", object$param_names)
  sigmas <- exp(colMeans(object$draws[, sigma_idx, drop = FALSE]))
  sigma_avg <- mean(sigmas)

  # Compute log-likelihood assuming Gaussian
  ll <- sum(dnorm(y, fitted_vals, sigma_avg, log = TRUE))

  # Number of parameters (excluding random effects for df)
  n_fixed <- sum(!grepl("^u_|^v_", object$param_names))

  structure(ll,
            df = n_fixed,
            nobs = nobs(object),
            class = "logLik")
}


#' Confidence intervals for phase_fit parameters
#'
#' Returns credible intervals from the posterior distribution.
#'
#' @param object A `phase_fit` object
#' @param parm Character vector of parameter names (NULL for all)
#' @param level Credible interval level (default 0.95)
#' @param ... Additional arguments (unused)
#'
#' @return Matrix with lower and upper bounds
#' @export
confint.phase_fit <- function(object, parm = NULL, level = 0.95, ...) {

  alpha <- (1 - level) / 2
  probs <- c(alpha, 1 - alpha)

  draws <- object$draws

  if (!is.null(parm)) {
    idx <- which(object$param_names %in% parm)
    if (length(idx) == 0) {
      stop("No matching parameters found", call. = FALSE)
    }
    draws <- draws[, idx, drop = FALSE]
  }

  ci <- t(apply(draws, 2, quantile, probs = probs))
  colnames(ci) <- paste0(probs * 100, "%")
  ci
}
