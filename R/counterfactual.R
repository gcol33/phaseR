#' Posterior predictive simulation
#'
#' Simulate outcomes from the posterior predictive distribution.
#'
#' @param object A `phase_fit` object
#' @param newdata Optional new data frame. If NULL, uses original data.
#' @param n_sims Number of posterior draws to use (default: all)
#' @param type "response" for outcome predictions, "phase" for phase probabilities
#'
#' @return A matrix of predictions (rows = observations, cols = posterior samples)
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- fit_phaseR(model, data)
#' preds <- posterior_predict(fit)
#' }
#'
posterior_predict <- function(object, newdata = NULL, n_sims = NULL,
                               type = c("response", "phase")) {

  type <- match.arg(type)

  if (!inherits(object, "phase_fit")) {
    stop("object must be a phase_fit", call. = FALSE)
  }

  if (object$k_phases > 2) {
    stop("posterior_predict not yet implemented for k > 2 phases", call. = FALSE)
  }

  # Get data
  data <- if (is.null(newdata)) object$data else newdata
  data <- data[order(data$id, data$time), ]

  # Get posterior draws
  draws <- object$draws
  if (!is.null(n_sims) && n_sims < nrow(draws)) {
    idx <- sample(nrow(draws), n_sims)
    draws <- draws[idx, , drop = FALSE]
  }
  n_draws <- nrow(draws)

  # Get model info
  model <- object$model
  response_var <- model$dynamics[[1]]$response

  # Build design matrices
  trans_formula <- model$transitions[[1]]$formula
  trans_parsed <- parse_formula_re(trans_formula)
  dyn_formula_0 <- model$dynamics[[1]]$formula
  dyn_formula_1 <- if (length(model$dynamics) > 1) model$dynamics[[2]]$formula else dyn_formula_0
  dyn_parsed_0 <- parse_formula_re(dyn_formula_0)
  dyn_parsed_1 <- parse_formula_re(dyn_formula_1)

  X_trans <- model.matrix(trans_parsed$fixed_formula, data)
  X_dyn_0 <- model.matrix(update(dyn_parsed_0$fixed_formula, NULL ~ .), data)
  X_dyn_1 <- model.matrix(update(dyn_parsed_1$fixed_formula, NULL ~ .), data)

  n_trans <- ncol(X_trans)
  n_dyn_0 <- ncol(X_dyn_0)
  n_dyn_1 <- ncol(X_dyn_1)

  # Get unit structure
  unit_ids <- unique(data$id)
  n_units <- length(unit_ids)

  N <- nrow(data)

  if (type == "response") {
    # Simulate outcomes
    pred_matrix <- matrix(NA_real_, nrow = N, ncol = n_draws)

    for (s in seq_len(n_draws)) {
      params <- draws[s, ]

      # Extract parameters
      idx <- 1
      beta_trans <- params[idx:(idx + n_trans - 1)]
      idx <- idx + n_trans
      beta_0 <- params[idx:(idx + n_dyn_0 - 1)]
      idx <- idx + n_dyn_0
      beta_1 <- params[idx:(idx + n_dyn_1 - 1)]
      idx <- idx + n_dyn_1
      sigma_0 <- exp(params[idx])
      sigma_1 <- exp(params[idx + 1])

      # Simulate for each unit
      for (u in seq_along(unit_ids)) {
        unit_rows <- which(data$id == unit_ids[u])

        # First time point: always phase 0
        t1 <- unit_rows[1]
        mu_0 <- sum(X_dyn_0[t1, ] * beta_0)
        pred_matrix[t1, s] <- rnorm(1, mu_0, sigma_0)

        current_phase <- 0

        # Subsequent time points
        for (i in 2:length(unit_rows)) {
          t_i <- unit_rows[i]

          if (current_phase == 0) {
            # Check for transition
            p_trans <- plogis(sum(X_trans[t_i, ] * beta_trans))
            if (runif(1) < p_trans) {
              current_phase <- 1
            }
          }

          # Generate outcome
          if (current_phase == 0) {
            mu <- sum(X_dyn_0[t_i, ] * beta_0)
            pred_matrix[t_i, s] <- rnorm(1, mu, sigma_0)
          } else {
            mu <- sum(X_dyn_1[t_i, ] * beta_1)
            pred_matrix[t_i, s] <- rnorm(1, mu, sigma_1)
          }
        }
      }
    }

    return(pred_matrix)

  } else {
    # Return phase probabilities using forward algorithm
    phase_probs <- matrix(NA_real_, nrow = N, ncol = n_draws)

    for (s in seq_len(n_draws)) {
      params <- draws[s, ]

      idx <- 1
      beta_trans <- params[idx:(idx + n_trans - 1)]
      idx <- idx + n_trans
      beta_0 <- params[idx:(idx + n_dyn_0 - 1)]
      idx <- idx + n_dyn_0
      beta_1 <- params[idx:(idx + n_dyn_1 - 1)]
      idx <- idx + n_dyn_1
      sigma_0 <- exp(params[idx])
      sigma_1 <- exp(params[idx + 1])

      y <- data[[response_var]]

      for (u in seq_along(unit_ids)) {
        unit_rows <- which(data$id == unit_ids[u])

        # Forward algorithm
        t1 <- unit_rows[1]
        mu_0 <- sum(X_dyn_0[t1, ] * beta_0)
        log_emit_0 <- dnorm(y[t1], mu_0, sigma_0, log = TRUE)

        alpha_0 <- log_emit_0
        alpha_1 <- -Inf

        phase_probs[t1, s] <- 0  # First point always phase 0

        for (i in 2:length(unit_rows)) {
          t_i <- unit_rows[i]

          p_trans <- plogis(sum(X_trans[t_i, ] * beta_trans))
          mu_0 <- sum(X_dyn_0[t_i, ] * beta_0)
          mu_1 <- sum(X_dyn_1[t_i, ] * beta_1)

          log_emit_0 <- dnorm(y[t_i], mu_0, sigma_0, log = TRUE)
          log_emit_1 <- dnorm(y[t_i], mu_1, sigma_1, log = TRUE)

          new_alpha_0 <- alpha_0 + log(1 - p_trans) + log_emit_0
          new_alpha_1 <- log_sum_exp(
            alpha_0 + log(p_trans) + log_emit_1,
            alpha_1 + log_emit_1
          )

          alpha_0 <- new_alpha_0
          alpha_1 <- new_alpha_1

          # P(phase = 1 | data)
          total <- log_sum_exp(alpha_0, alpha_1)
          phase_probs[t_i, s] <- exp(alpha_1 - total)
        }
      }
    }

    return(phase_probs)
  }
}


#' Log-sum-exp helper
#' @keywords internal
log_sum_exp <- function(a, b) {
  if (is.infinite(a) && a < 0) return(b)
  if (is.infinite(b) && b < 0) return(a)
  if (a > b) {
    a + log1p(exp(b - a))
  } else {
    b + log1p(exp(a - b))
  }
}


#' Counterfactual prediction
#'
#' Predict outcomes under hypothetical covariate values, integrating over
#' phase uncertainty.
#'
#' @param object A `phase_fit` object
#' @param newdata Data frame with counterfactual covariate values
#' @param n_sims Number of posterior draws to use
#'
#' @return A list with:
#'   - `mean`: Mean predicted outcome per observation
#'   - `sd`: SD of predicted outcome per observation
#'   - `samples`: Matrix of samples (rows = obs, cols = draws)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Predict under treatment = 1 for all
#' cf_data <- data
#' cf_data$treatment <- 1
#' cf <- counterfactual(fit, newdata = cf_data)
#' }
#'
counterfactual <- function(object, newdata, n_sims = NULL) {

  if (!inherits(object, "phase_fit")) {
    stop("object must be a phase_fit", call. = FALSE)
  }

  samples <- posterior_predict(object, newdata = newdata, n_sims = n_sims,
                                type = "response")

  list(
    mean = rowMeans(samples),
    sd = apply(samples, 1, sd),
    lower = apply(samples, 1, quantile, 0.025),
    upper = apply(samples, 1, quantile, 0.975),
    samples = samples
  )
}


#' Average treatment effect
#'
#' Compute the average treatment effect by comparing counterfactual outcomes
#' under treatment vs control conditions.
#'
#' @param object A `phase_fit` object
#' @param treatment_var Name of the treatment variable
#' @param treat_value Value representing treatment (default 1)
#' @param control_value Value representing control (default 0)
#' @param n_sims Number of posterior draws
#'
#' @return A list with:
#'   - `ate`: Posterior mean of ATE
#'   - `ate_sd`: Posterior SD of ATE
#'   - `ate_ci`: 95% credible interval
#'   - `samples`: Vector of ATE samples
#'
#' @export
#'
#' @examples
#' \dontrun{
#' effect <- ate(fit, treatment_var = "treatment")
#' }
#'
ate <- function(object, treatment_var, treat_value = 1, control_value = 0,
                n_sims = NULL) {

  if (!inherits(object, "phase_fit")) {
    stop("object must be a phase_fit", call. = FALSE)
  }

  if (!treatment_var %in% names(object$data)) {
    stop(sprintf("Variable '%s' not found in data", treatment_var), call. = FALSE)
  }

  # Create counterfactual datasets
  data_treat <- object$data
  data_control <- object$data
  data_treat[[treatment_var]] <- treat_value
  data_control[[treatment_var]] <- control_value

  # Get predictions
  cf_treat <- counterfactual(object, newdata = data_treat, n_sims = n_sims)
  cf_control <- counterfactual(object, newdata = data_control, n_sims = n_sims)

  # Compute ATE for each draw
  ate_samples <- colMeans(cf_treat$samples - cf_control$samples)

  list(
    ate = mean(ate_samples),
    ate_sd = sd(ate_samples),
    ate_ci = quantile(ate_samples, c(0.025, 0.975)),
    samples = ate_samples
  )
}


#' Individual treatment effects
#'
#' Compute treatment effects for each unit, integrating over phase uncertainty.
#'
#' @param object A `phase_fit` object
#' @param treatment_var Name of the treatment variable
#' @param treat_value Value representing treatment (default 1)
#' @param control_value Value representing control (default 0)
#' @param n_sims Number of posterior draws
#'
#' @return A data frame with unit-level treatment effects
#'
#' @export
#'
ite <- function(object, treatment_var, treat_value = 1, control_value = 0,
                n_sims = NULL) {

  if (!inherits(object, "phase_fit")) {
    stop("object must be a phase_fit", call. = FALSE)
  }

  if (!treatment_var %in% names(object$data)) {
    stop(sprintf("Variable '%s' not found in data", treatment_var), call. = FALSE)
  }

  # Create counterfactual datasets
  data_treat <- object$data
  data_control <- object$data
  data_treat[[treatment_var]] <- treat_value
  data_control[[treatment_var]] <- control_value

  # Get predictions
  cf_treat <- counterfactual(object, newdata = data_treat, n_sims = n_sims)
  cf_control <- counterfactual(object, newdata = data_control, n_sims = n_sims)

  # Compute ITE per observation
  ite_samples <- cf_treat$samples - cf_control$samples

  # Aggregate by unit
  unit_ids <- unique(object$data$id)
  results <- lapply(unit_ids, function(uid) {
    rows <- which(object$data$id == uid)
    unit_ite <- colMeans(ite_samples[rows, , drop = FALSE])
    data.frame(
      id = uid,
      ite_mean = mean(unit_ite),
      ite_sd = sd(unit_ite),
      ite_lower = quantile(unit_ite, 0.025),
      ite_upper = quantile(unit_ite, 0.975)
    )
  })

  do.call(rbind, results)
}
