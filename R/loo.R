#' Leave-One-Out Cross-Validation via PSIS
#'
#' Computes approximate leave-one-out cross-validation using Pareto-smoothed
#' importance sampling (PSIS-LOO). This provides an estimate of out-of-sample
#' predictive accuracy that is more robust than WAIC.
#'
#' @param object A `phase_fit` object
#' @param ... Additional arguments (unused)
#'
#' @return A `loo_phaseR` object with components:
#' \itemize{
#'   \item `elpd_loo`: Expected log pointwise predictive density for LOO
#'   \item `p_loo`: Effective number of parameters
#'   \item `looic`: LOO information criterion (-2 * elpd_loo)
#'   \item `se_elpd_loo`: Standard error of elpd_loo
#'   \item `pointwise`: Data frame with per-observation values
#'   \item `diagnostics`: List with k-hat values and convergence info
#' }
#'
#' @export
#'
#' @details
#' PSIS-LOO approximates leave-one-out cross-validation without refitting
#' the model. It uses importance sampling with Pareto-smoothed weights to
#' estimate the predictive density for each held-out observation.
#'
#' The k-hat diagnostic indicates the reliability of the estimate:
#' \itemize{
#'   \item k < 0.5: Good, estimates are reliable
#'   \item 0.5 <= k < 0.7: Okay, estimates may have some error
#'   \item k >= 0.7: Bad, estimates are unreliable for these observations
#' }
#'
#' If many observations have k > 0.7, consider using K-fold cross-validation
#' instead.
#'
#' @references
#' Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC.
#' Statistics and Computing 27, 1413-1432.
#'
#' @seealso [waic()], [loo_compare()]
#'
#' @examples
#' \dontrun{
#' fit <- fit_phaseR(model, data)
#' loo_result <- loo(fit)
#' print(loo_result)
#'
#' # Check diagnostics
#' plot(loo_result$diagnostics$pareto_k)
#' }
#'
loo <- function(object, ...) {

 if (!inherits(object, "phase_fit")) {
    stop("Object must be a phase_fit", call. = FALSE)
  }

  # Get pointwise log-likelihood matrix (draws x observations)
  ll_matrix <- compute_pointwise_ll(object)
  n_draws <- nrow(ll_matrix)
  n_obs <- ncol(ll_matrix)

  # Compute PSIS for each observation
  psis_result <- psis_loo(ll_matrix)

  # Extract results
  elpd_loo_i <- psis_result$elpd_loo
  p_loo_i <- psis_result$p_loo
  pareto_k <- psis_result$pareto_k

  # Aggregate
elpd_loo <- sum(elpd_loo_i)
  p_loo <- sum(p_loo_i)
  looic <- -2 * elpd_loo

  # Standard errors
  se_elpd_loo <- sqrt(n_obs * var(elpd_loo_i))
  se_looic <- 2 * se_elpd_loo

  # Count problematic observations
  n_bad <- sum(pareto_k > 0.7)
  n_okay <- sum(pareto_k > 0.5 & pareto_k <= 0.7)

  if (n_bad > 0) {
    warning(sprintf(
      "Found %d observation(s) with Pareto k > 0.7. These estimates may be unreliable.",
      n_bad
    ), call. = FALSE)
  }

  structure(
    list(
      elpd_loo = elpd_loo,
      p_loo = p_loo,
      looic = looic,
      se_elpd_loo = se_elpd_loo,
      se_looic = se_looic,
      pointwise = data.frame(
        elpd_loo = elpd_loo_i,
        p_loo = p_loo_i,
        looic = -2 * elpd_loo_i,
        pareto_k = pareto_k
      ),
      diagnostics = list(
        pareto_k = pareto_k,
        n_bad = n_bad,
        n_okay = n_okay,
        n_good = n_obs - n_bad - n_okay
      )
    ),
    class = "loo_phaseR"
  )
}


#' @export
print.loo_phaseR <- function(x, digits = 1, ...) {
  cat("Leave-One-Out Cross-Validation (PSIS-LOO)\n\n")
  cat(sprintf("  elpd_loo:  %.1f (SE: %.1f)\n", x$elpd_loo, x$se_elpd_loo))
  cat(sprintf("  p_loo:     %.1f\n", x$p_loo))
  cat(sprintf("  looic:     %.1f (SE: %.1f)\n", x$looic, x$se_looic))
  cat("\n")
  cat("Pareto k diagnostic:\n")
  cat(sprintf("  Good (k < 0.5):       %d\n", x$diagnostics$n_good))
  cat(sprintf("  Okay (0.5 <= k < 0.7): %d\n", x$diagnostics$n_okay))
  cat(sprintf("  Bad (k >= 0.7):        %d\n", x$diagnostics$n_bad))
  invisible(x)
}


#' Compare models using LOO-CV
#'
#' Compare multiple fitted phase models using leave-one-out cross-validation.
#'
#' @param ... `phase_fit` objects to compare, or `loo_phaseR` objects
#' @param model_names Optional names for the models
#'
#' @return A `loo_compare` object (data frame with comparison)
#' @export
#'
#' @details
#' Models are ranked by expected log pointwise predictive density (elpd_loo),
#' with higher values indicating better predictive performance. The difference
#' in elpd and its standard error are computed relative to the best model.
#'
#' @examples
#' \dontrun{
#' fit1 <- fit_phaseR(model1, data)
#' fit2 <- fit_phaseR(model2, data)
#' loo_compare(fit1, fit2, model_names = c("Simple", "Complex"))
#' }
#'
loo_compare <- function(..., model_names = NULL) {

  objects <- list(...)

  if (length(objects) < 2) {
    stop("Need at least 2 models to compare", call. = FALSE)
  }

  if (is.null(model_names)) {
    model_names <- paste0("Model_", seq_along(objects))
  }

  # Compute LOO for each object (or use if already computed)
  loos <- lapply(objects, function(obj) {
    if (inherits(obj, "loo_phaseR")) {
      obj
    } else if (inherits(obj, "phase_fit")) {
      loo(obj)
    } else {
      stop("Objects must be phase_fit or loo_phaseR", call. = FALSE)
    }
  })

  # Extract elpd_loo values
  elpd_values <- vapply(loos, function(l) l$elpd_loo, numeric(1))
  se_values <- vapply(loos, function(l) l$se_elpd_loo, numeric(1))

  # Sort by elpd_loo (higher is better)
  ord <- order(elpd_values, decreasing = TRUE)

  result <- data.frame(
    model = model_names[ord],
    elpd_loo = elpd_values[ord],
    se_elpd_loo = se_values[ord],
    p_loo = vapply(loos, function(l) l$p_loo, numeric(1))[ord],
    looic = vapply(loos, function(l) l$looic, numeric(1))[ord]
  )

  # Compute differences from best model
  result$elpd_diff <- result$elpd_loo - result$elpd_loo[1]

  # SE of difference (requires pointwise values)
  # For simplicity, use sqrt(sum of variances) as approximation
  # More accurate would be to compute from pointwise differences
  best_loo <- loos[[ord[1]]]
  result$se_diff <- NA_real_
  result$se_diff[1] <- 0

  for (i in 2:nrow(result)) {
    this_loo <- loos[[ord[i]]]
    # Pointwise difference
    diff_i <- best_loo$pointwise$elpd_loo - this_loo$pointwise$elpd_loo
    n_obs <- length(diff_i)
    result$se_diff[i] <- sqrt(n_obs * var(diff_i))
  }

  rownames(result) <- NULL

  structure(result, class = c("loo_compare_phaseR", "data.frame"))
}


#' @export
print.loo_compare_phaseR <- function(x, digits = 1, ...) {
  cat("Model Comparison (LOO-CV)\n\n")
  cat("Ranked by elpd_loo (higher is better):\n\n")
  print.data.frame(x, digits = digits, row.names = FALSE)
  invisible(x)
}


# =============================================================================
# PSIS Implementation
# =============================================================================

#' Pareto-Smoothed Importance Sampling for LOO
#'
#' Core PSIS-LOO algorithm. Computes smoothed importance weights and
#' LOO estimates from a log-likelihood matrix.
#'
#' @param ll_matrix Matrix of log-likelihoods (draws x observations)
#' @return List with elpd_loo, p_loo, and pareto_k vectors
#' @keywords internal
#'
psis_loo <- function(ll_matrix) {

  n_draws <- nrow(ll_matrix)
  n_obs <- ncol(ll_matrix)

  elpd_loo <- numeric(n_obs)
  p_loo <- numeric(n_obs)
  pareto_k <- numeric(n_obs)

  for (i in seq_len(n_obs)) {
    ll_i <- ll_matrix[, i]

    # Raw importance ratios: 1/p(y_i|theta) = exp(-ll_i)
    # We work with log ratios for stability
    log_ratios <- -ll_i

    # Apply PSIS smoothing
    psis_out <- psis_smooth(log_ratios)

    # Normalized weights (on log scale, then exponentiate)
    log_weights <- psis_out$log_weights
    weights <- exp(log_weights - max(log_weights))
    weights <- weights / sum(weights)

    # LOO log predictive density
    # elpd_loo_i = log(sum(w * p(y_i|theta))) = log(sum(w * exp(ll_i)))
    # Using log-sum-exp trick
    max_ll <- max(ll_i)
    elpd_loo[i] <- log(sum(weights * exp(ll_i - max_ll))) + max_ll

    # p_loo: difference between full posterior and LOO
    # Approximated as: log(mean(exp(ll_i))) - elpd_loo_i
    lppd_i <- log(mean(exp(ll_i - max_ll))) + max_ll
    p_loo[i] <- lppd_i - elpd_loo[i]

    pareto_k[i] <- psis_out$k
  }

  list(
    elpd_loo = elpd_loo,
    p_loo = p_loo,
    pareto_k = pareto_k
  )
}


#' Apply Pareto smoothing to importance ratios
#'
#' Fits a generalized Pareto distribution to the tail of the importance
#' ratios and replaces extreme values with expected order statistics.
#'
#' @param log_ratios Log importance ratios
#' @return List with smoothed log_weights and Pareto k estimate
#' @keywords internal
#'
psis_smooth <- function(log_ratios) {

  n <- length(log_ratios)

  # Stabilize by subtracting max
  log_ratios <- log_ratios - max(log_ratios)

  # Number of tail samples to use for Pareto fit
  # Following Vehtari et al. (2017): min(n/5, 3*sqrt(n))
  m <- min(floor(n / 5), ceiling(3 * sqrt(n)))
  m <- max(m, 5)  # At least 5 samples

  # Sort ratios (descending on original scale = ascending on log scale when negated)
  ord <- order(log_ratios, decreasing = TRUE)
  sorted_log <- log_ratios[ord]

  # Get the m largest values
  tail_log <- sorted_log[1:m]

  # Threshold (smallest of the tail values)
  threshold <- tail_log[m]

  # Fit generalized Pareto to exceedances
  exceedances <- exp(tail_log) - exp(threshold)
  exceedances <- pmax(exceedances, .Machine$double.eps)

  # Estimate Pareto shape parameter k using probability-weighted moments
  k <- estimate_pareto_k(exceedances)

  # If k is reasonable, smooth the tail
  if (k < 0.7 && k > -0.5) {
    # Replace tail values with expected order statistics from fitted Pareto
    smoothed_tail <- smooth_pareto_tail(exceedances, k, m)

    # Put back on log scale
    sorted_log[1:m] <- log(smoothed_tail + exp(threshold))
  }

  # Restore original order
  log_weights <- numeric(n)
  log_weights[ord] <- sorted_log

  list(
    log_weights = log_weights,
    k = k
  )
}


#' Estimate Pareto shape parameter
#'
#' Uses probability-weighted moments estimator for the generalized Pareto
#' distribution shape parameter.
#'
#' @param x Exceedances (positive values)
#' @return Estimated shape parameter k
#' @keywords internal
#'
estimate_pareto_k <- function(x) {

  n <- length(x)

  if (n < 5) {
    return(NA_real_)
  }

  # Sort in ascending order
  x_sorted <- sort(x)

  # Probability-weighted moments estimator
  # Following Zhang & Stephens (2009) approach
  m1 <- mean(x)

  # Weighted mean with probability weights
  probs <- (seq_len(n) - 0.5) / n
  m2 <- sum(x_sorted * probs) / n

  # Shape estimate
  if (m1 > 0 && m2 > 0) {
    k <- m1 / (m1 - 2 * m2)
    # Constrain to reasonable range
    k <- min(max(k, -0.5), 1.5)
  } else {
    k <- 0.7  # Default to "bad" if estimation fails
  }

  k
}


#' Smooth Pareto tail values
#'
#' Replaces extreme values with expected order statistics from fitted
#' generalized Pareto distribution.
#'
#' @param x Exceedances
#' @param k Shape parameter
#' @param m Number of tail samples
#' @return Smoothed exceedances
#' @keywords internal
#'
smooth_pareto_tail <- function(x, k, m) {

  n <- length(x)

  # Scale parameter estimate
  sigma <- mean(x) * (1 + k)

  if (sigma <= 0 || !is.finite(sigma)) {
    return(x)  # Return original if estimation fails
  }

  # Expected order statistics for Pareto(k, sigma)
  # E[X_(j)] where j is the order statistic index
  # For generalized Pareto: E[X_(j:m)] = sigma * (B((m-j+1)/m, 1-k) - 1) / k
  # where B is the beta function

  # Simplified: use quantile function
  # F^{-1}(p) = sigma * ((1-p)^(-k) - 1) / k  for k != 0

  probs <- (seq_len(m) - 0.5) / m
  probs <- 1 - probs  # For upper tail

  if (abs(k) < 1e-10) {
    # k ≈ 0: exponential case
    smoothed <- sigma * (-log(probs))
  } else {
    # General case
    smoothed <- sigma * (probs^(-k) - 1) / k
  }

  # Sort descending (largest first)
  smoothed <- sort(smoothed, decreasing = TRUE)

  # Truncate any negative values
  smoothed <- pmax(smoothed, 0)

  smoothed
}
