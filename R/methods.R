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
    cat(sprintf("  Transitions at max tree depth: %d\n",
                x$diagnostics$n_max_treedepth))
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
#'
#' Multi-chain effective sample size of Geyer's initial monotone sequence
#' estimator, as in Stan: autocovariances are averaged over chains, combined
#' with the between-chain variance, and summed in pairs of lags until a pair
#' turns negative.
#'
#' @param draws Matrix of draws, one column per parameter
#' @param chain_id Chain of each row of `draws`
#' @return Numeric vector of effective sample sizes, one per column
#' @keywords internal
compute_ess <- function(draws, chain_id) {
  draws <- as.matrix(draws)
  apply(draws, 2, ess_geyer, chain_id = chain_id)
}


#' Effective sample size of one parameter
#' @keywords internal
ess_geyer <- function(x, chain_id) {
  chains <- split(x, chain_id)
  m <- length(chains)
  n <- min(lengths(chains))
  if (n < 4) return(NA_real_)
  X <- vapply(chains, function(v) v[seq_len(n)], numeric(n))
  acov <- apply(X, 2, function(v) {
    stats::acf(v, lag.max = n - 1, type = "covariance", demean = TRUE,
               plot = FALSE)$acf[, 1, 1]
  })
  mean_var <- mean(acov[1, ]) * n / (n - 1)
  var_plus <- mean_var * (n - 1) / n
  if (m > 1) var_plus <- var_plus + stats::var(colMeans(X))
  if (!is.finite(var_plus) || var_plus <= 0) return(NA_real_)

  rho <- 1 - (mean_var - rowMeans(acov)) / var_plus
  rho[1] <- 1
  pair_sums <- numeric(0)
  t <- 1
  while (t + 1 <= n) {
    p <- rho[t] + rho[t + 1]
    if (p <= 0) break
    pair_sums <- c(pair_sums, p)
    t <- t + 2
  }
  tau <- -1 + 2 * sum(cummin(pair_sums))
  m * n / max(tau, 1 / log10(m * n))
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
#' The posterior mean of the one-step predictive mean of each observation,
#' E(y_t | y_0, ..., y_{t-1}): the outcome mean of each phase weighted by
#' the probability of the phase given the unit's earlier outcomes.
#'
#' @param object A `phase_fit` object
#' @param ... Additional arguments (unused)
#'
#' @return Numeric vector of fitted values, one per row of the data
#' @export
fitted.phase_fit <- function(object, ...) {
  one_step_moments(object)$mean
}


#' Extract residuals from phase_fit
#'
#' One-step prediction residuals: the observed outcome minus [fitted()],
#' divided for `type = "pearson"` by the SD of the one-step predictive
#' distribution (phase mixture of the family variances, plus the posterior
#' variance of the predictive mean).
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
  m <- one_step_moments(object)
  resid <- y - m$mean
  if (type == "pearson") resid <- resid / sqrt(m$var)
  resid
}


#' Moments of the one-step predictive distribution
#'
#' @param object A `phase_fit` object
#' @return List with `mean` and `var` per row of the data: the predictive
#'   mean and variance of y_t given y_0, ..., y_{t-1}, over the posterior
#' @keywords internal
one_step_moments <- function(object) {
  s <- phase_summaries(object, object$data, object$draws)
  family_var <- s$mean
  for (k in seq_along(s$family)) {
    family_var[, , k] <- switch(s$family[k] + 1L,
      s$sigma[, k]^2,
      s$mean[, , k],
      s$mean[, , k] * (1 - s$mean[, , k])
    )
  }
  mean_d <- apply(s$predictive * s$mean, c(1, 2), sum)
  second_d <- apply(s$predictive * (family_var + s$mean^2), c(1, 2), sum)
  mean <- colMeans(mean_d)
  list(mean = mean, var = colMeans(second_d) - mean^2)
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
#' The log-likelihood at the posterior mean of the parameters, marginalised
#' over the latent phases with the forward recursion the backends fit with.
#'
#' @param object A `phase_fit` object
#' @param ... Additional arguments (unused)
#'
#' @return A logLik object
#' @export
logLik.phase_fit <- function(object, ...) {
  at_mean <- object
  at_mean$draws <- matrix(colMeans(object$draws), nrow = 1,
                          dimnames = list(NULL, colnames(object$draws)))
  ll <- sum(compute_pointwise_ll(at_mean))

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
