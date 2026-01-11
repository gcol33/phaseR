#' Predict from a phase model
#'
#' @param object A `phase_fit` object
#' @param newdata Optional new data for prediction
#' @param type Type of prediction: "response" (default), "phase", or "transition"
#' @param summary Return summary statistics (default TRUE) or full draws
#' @param ... Additional arguments (unused)
#'
#' @return A data frame with predictions and credible intervals
#' @export
#'
predict.phase_fit <- function(object,
                              newdata = NULL,
                              type = c("response", "phase", "transition"),
                              summary = TRUE,
                              ...) {

  type <- match.arg(type)

  if (is.null(newdata)) {
    newdata <- object$data
  }

  result <- switch(type,
    response = predict_response(object, newdata),
    phase = predict_phase(object, newdata),
    transition = predict_transition(object, newdata)
  )

  if (summary && is.matrix(result)) {
    summarize_predictions(result)
  } else {
    result
  }
}


#' Predict response marginalizing over phases
#' @keywords internal
predict_response <- function(fit, newdata) {

  # Get phase probabilities
  phase_probs <- predict_phase_internal(fit, newdata)

  # Get draws
  draws <- fit$draws
  n_draws <- nrow(draws)
  n_obs <- nrow(newdata)

  # Build design matrices
  model <- fit$model
  dyn_formula_0 <- model$dynamics[[1]]$formula
  dyn_formula_1 <- if (length(model$dynamics) > 1) model$dynamics[[2]]$formula else dyn_formula_0

  X_dyn_0 <- model.matrix(update(dyn_formula_0, NULL ~ .), newdata)
  X_dyn_1 <- model.matrix(update(dyn_formula_1, NULL ~ .), newdata)

  # Extract parameter indices
  param_names <- fit$param_names
  beta_0_idx <- grep(paste0("^beta_", model$phases$names[1], "_"), param_names)
  beta_1_idx <- grep(paste0("^beta_", model$phases$names[2], "_"), param_names)

  # Compute predictions for each draw
  pred_draws <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)

  for (d in seq_len(n_draws)) {
    beta_0 <- draws[d, beta_0_idx]
    beta_1 <- draws[d, beta_1_idx]

    mu_0 <- X_dyn_0 %*% beta_0
    mu_1 <- X_dyn_1 %*% beta_1

    # Marginalize over phase uncertainty
    pred_draws[d, ] <- phase_probs[, 1] * mu_0 + phase_probs[, 2] * mu_1
  }

  pred_draws
}


#' Predict phase probabilities
#' @keywords internal
predict_phase <- function(fit, newdata) {
  phase_probs <- predict_phase_internal(fit, newdata)
  colnames(phase_probs) <- fit$model$phases$names
  as.data.frame(phase_probs)
}


#' Internal function for phase probability computation
#' @keywords internal
predict_phase_internal <- function(fit, newdata) {

  # This uses a simplified approach for v0.2.0
  # Full forward-backward algorithm would be used for exact inference

  draws <- fit$draws
  model <- fit$model
  n_obs <- nrow(newdata)

  # Get transition parameters (posterior mean for simplicity)
  param_names <- fit$param_names
  trans_idx <- grep("^beta_trans", param_names)
  beta_trans <- colMeans(draws[, trans_idx, drop = FALSE])

  # Build transition design matrix
  trans_formula <- model$transitions[[1]]$formula
  X_trans <- model.matrix(trans_formula, newdata)

  # Compute transition probabilities
  eta <- X_trans %*% beta_trans
  p_trans <- 1 / (1 + exp(-eta))

  # Simple forward pass for phase probabilities
  # Group by unit
  unit_ids <- unique(newdata$id)
  phase_probs <- matrix(0, nrow = n_obs, ncol = 2)

  for (uid in unit_ids) {
    idx <- which(newdata$id == uid)
    n_t <- length(idx)

    # Initial: always phase 0
    alpha_0 <- 1
    alpha_1 <- 0

    for (t in seq_len(n_t)) {
      obs_idx <- idx[t]

      if (t == 1) {
        phase_probs[obs_idx, ] <- c(1, 0)
      } else {
        p_t <- p_trans[obs_idx]
        new_alpha_0 <- alpha_0 * (1 - p_t)
        new_alpha_1 <- alpha_0 * p_t + alpha_1

        # Normalize
        total <- new_alpha_0 + new_alpha_1
        phase_probs[obs_idx, ] <- c(new_alpha_0, new_alpha_1) / total

        alpha_0 <- new_alpha_0 / total
        alpha_1 <- new_alpha_1 / total
      }
    }
  }

  phase_probs
}


#' Predict transition probabilities
#' @keywords internal
predict_transition <- function(fit, newdata) {

  draws <- fit$draws
  model <- fit$model

  # Get transition parameters
  param_names <- fit$param_names
  trans_idx <- grep("^beta_trans", param_names)

  # Build design matrix
  trans_formula <- model$transitions[[1]]$formula
  X_trans <- model.matrix(trans_formula, newdata)

  n_draws <- nrow(draws)
  n_obs <- nrow(newdata)

  # Compute for each draw
  trans_draws <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)

  for (d in seq_len(n_draws)) {
    beta_trans <- draws[d, trans_idx]
    eta <- X_trans %*% beta_trans
    trans_draws[d, ] <- 1 / (1 + exp(-eta))
  }

  trans_draws
}


#' Summarize prediction draws
#' @keywords internal
summarize_predictions <- function(draws) {
  data.frame(
    mean = colMeans(draws),
    sd = apply(draws, 2, sd),
    q2.5 = apply(draws, 2, quantile, 0.025),
    q50 = apply(draws, 2, median),
    q97.5 = apply(draws, 2, quantile, 0.975)
  )
}
