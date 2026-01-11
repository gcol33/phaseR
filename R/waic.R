#' Compute WAIC for model comparison
#'
#' Computes the Widely Applicable Information Criterion (WAIC) for a
#' fitted phase model. This provides an estimate of out-of-sample
#' predictive accuracy.
#'
#' @param object A `phase_fit` object
#' @param ... Additional arguments (unused)
#'
#' @return A list with components:
#' \itemize{
#'   \item `waic`: The WAIC value (lower is better)
#'   \item `lppd`: Log pointwise predictive density
#'   \item `p_waic`: Effective number of parameters
#'   \item `se`: Standard error of WAIC
#' }
#'
#' @export
#'
#' @details
#' WAIC is computed as:
#' \deqn{WAIC = -2 * (lppd - p_{WAIC})}
#'
#' where lppd is the log pointwise predictive density and p_WAIC is the
#' effective number of parameters (variance of the log-likelihood).
#'
#' For phase models, the log-likelihood is computed by marginalizing over
#' phase uncertainty at the posterior mean parameters.
#'
#' @references
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and
#' widely applicable information criterion in singular learning theory.
#' Journal of Machine Learning Research 11, 3571-3594.
#'
#' @examples
#' \dontrun{
#' fit1 <- fit_phaseR(model1, data)
#' fit2 <- fit_phaseR(model2, data)
#' waic(fit1)$waic
#' waic(fit2)$waic
#' }
#'
waic <- function(object, ...) {

  if (!inherits(object, "phase_fit")) {
    stop("Object must be a phase_fit", call. = FALSE)
  }

  # Get data
  y <- object$data[[object$model$dynamics[[1]]$response]]
  n_obs <- length(y)
  n_draws <- nrow(object$draws)

  # Compute pointwise log-likelihood for each draw
  ll_matrix <- compute_pointwise_ll(object)

  # lppd: log pointwise predictive density
  # For each observation, compute log(mean(exp(ll)))
  lppd_i <- apply(ll_matrix, 2, function(ll_obs) {
    max_ll <- max(ll_obs)
    log(mean(exp(ll_obs - max_ll))) + max_ll
  })
  lppd <- sum(lppd_i)

  # p_waic: effective number of parameters
  # Variance of log-likelihood for each observation
  p_waic_i <- apply(ll_matrix, 2, var)
  p_waic <- sum(p_waic_i)

  # WAIC
  waic_val <- -2 * (lppd - p_waic)

  # Standard error (based on observation-level contributions)
  waic_i <- -2 * (lppd_i - p_waic_i)
  se <- sqrt(n_obs * var(waic_i))

  structure(
    list(
      waic = waic_val,
      lppd = lppd,
      p_waic = p_waic,
      se = se,
      pointwise = data.frame(
        lppd = lppd_i,
        p_waic = p_waic_i,
        waic = waic_i
      )
    ),
    class = "waic_phaseR"
  )
}


#' @export
print.waic_phaseR <- function(x, digits = 1, ...) {
  cat("Widely Applicable Information Criterion (WAIC)\n\n")
  cat(sprintf("  WAIC:    %.1f (SE: %.1f)\n", x$waic, x$se))
  cat(sprintf("  lppd:    %.1f\n", x$lppd))
  cat(sprintf("  p_WAIC:  %.1f\n", x$p_waic))
  invisible(x)
}


#' Compute pointwise log-likelihood
#' @keywords internal
compute_pointwise_ll <- function(fit) {

  draws <- fit$draws
  n_draws <- nrow(draws)
  data <- fit$data
  model <- fit$model
  n_obs <- nrow(data)

  # Build design matrices
  trans_formula <- model$transitions[[1]]$formula
  dyn_formula_0 <- model$dynamics[[1]]$formula
  dyn_formula_1 <- if (length(model$dynamics) > 1) model$dynamics[[2]]$formula else dyn_formula_0

  X_trans <- model.matrix(trans_formula, data)
  X_dyn_0 <- model.matrix(update(dyn_formula_0, NULL ~ .), data)
  X_dyn_1 <- model.matrix(update(dyn_formula_1, NULL ~ .), data)

  # Get parameter indices
  param_names <- fit$param_names
  trans_idx <- grep("^beta_trans", param_names)
  beta_0_idx <- grep(paste0("^beta_", model$phases$names[1], "_"), param_names)
  beta_1_idx <- grep(paste0("^beta_", model$phases$names[2], "_"), param_names)
  sigma_0_idx <- grep(paste0("^log_sigma_", model$phases$names[1], "$"), param_names)
  sigma_1_idx <- grep(paste0("^log_sigma_", model$phases$names[2], "$"), param_names)

  y <- data[[model$dynamics[[1]]$response]]
  unit_ids <- unique(data$id)

  # Matrix: draws x observations

  ll_matrix <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)

  for (d in seq_len(n_draws)) {
    beta_trans <- draws[d, trans_idx]
    beta_0 <- draws[d, beta_0_idx]
    beta_1 <- draws[d, beta_1_idx]
    sigma_0 <- exp(draws[d, sigma_0_idx])
    sigma_1 <- exp(draws[d, sigma_1_idx])

    # Compute phase probabilities via forward algorithm
    for (uid in unit_ids) {
      idx <- which(data$id == uid)
      n_t <- length(idx)

      alpha_0 <- 1
      alpha_1 <- 0

      for (t in seq_len(n_t)) {
        obs_idx <- idx[t]

        # Emission probabilities
        mu_0 <- sum(X_dyn_0[obs_idx, ] * beta_0)
        mu_1 <- sum(X_dyn_1[obs_idx, ] * beta_1)
        emit_0 <- dnorm(y[obs_idx], mu_0, sigma_0)
        emit_1 <- dnorm(y[obs_idx], mu_1, sigma_1)

        if (t == 1) {
          # Initial: always phase 0
          ll_matrix[d, obs_idx] <- log(emit_0 + 1e-300)
        } else {
          # Transition probability
          p_trans <- plogis(sum(X_trans[obs_idx, ] * beta_trans))

          # Forward step
          new_alpha_0 <- alpha_0 * (1 - p_trans) * emit_0
          new_alpha_1 <- (alpha_0 * p_trans + alpha_1) * emit_1

          # Log-likelihood for this observation (marginalized)
          ll_matrix[d, obs_idx] <- log(new_alpha_0 + new_alpha_1 + 1e-300)

          # Normalize for next step
          total <- new_alpha_0 + new_alpha_1
          if (total > 0) {
            alpha_0 <- new_alpha_0 / total
            alpha_1 <- new_alpha_1 / total
          }
        }
      }
    }
  }

  ll_matrix
}


#' Compare models using WAIC
#'
#' Compare multiple fitted phase models using WAIC.
#'
#' @param ... `phase_fit` objects to compare
#' @param model_names Optional names for the models
#'
#' @return A data frame with WAIC comparisons
#' @export
#'
#' @examples
#' \dontrun{
#' compare_waic(fit1, fit2, model_names = c("Simple", "Complex"))
#' }
#'
compare_waic <- function(..., model_names = NULL) {

  fits <- list(...)

  if (length(fits) < 2) {
    stop("Need at least 2 models to compare", call. = FALSE)
  }

  if (is.null(model_names)) {
    model_names <- paste0("Model_", seq_along(fits))
  }

  waics <- lapply(fits, waic)

  result <- data.frame(
    model = model_names,
    waic = vapply(waics, function(w) w$waic, numeric(1)),
    se = vapply(waics, function(w) w$se, numeric(1)),
    lppd = vapply(waics, function(w) w$lppd, numeric(1)),
    p_waic = vapply(waics, function(w) w$p_waic, numeric(1))
  )

  # Sort by WAIC (best first)
  result <- result[order(result$waic), ]
  result$delta_waic <- result$waic - result$waic[1]
  rownames(result) <- NULL

  structure(result, class = c("compare_waic", "data.frame"))
}


#' @export
print.compare_waic <- function(x, digits = 1, ...) {
  cat("Model Comparison (WAIC)\n\n")
  print.data.frame(x, digits = digits, row.names = FALSE)
  invisible(x)
}
