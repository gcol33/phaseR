#' Posterior predictive simulation
#'
#' Simulates outcomes from the posterior predictive distribution given the
#' covariates of `newdata`, or returns phase probabilities given its
#' outcomes.
#'
#' For `type = "response"`, each posterior draw simulates every unit
#' forward: the unit starts in the first phase, moves to the next phase
#' with the transition probability its covariates imply at each later
#' observation, and draws each outcome from the dynamics of its current
#' phase (Gaussian, Poisson or binomial, random effects included; levels
#' the fit has not seen get no effect). The observed outcomes are not used,
#' so changing a covariate in `newdata` changes both who transitions and
#' the outcome in each phase.
#'
#' For `type = "phase"`, the probability that each observation is in
#' `phase` given the unit's outcomes up to and including that observation
#' (the forward filter), which needs the outcome in `newdata`.
#'
#' @param object A `phase_fit` object
#' @param newdata Optional new data frame. If NULL, uses original data.
#' @param n_sims Number of posterior draws to use (default: all)
#' @param type "response" for outcome predictions, "phase" for phase
#'   probabilities
#' @param phase For `type = "phase"`, the name of the phase whose
#'   probability is returned (default: the last phase)
#'
#' @return A matrix with one row per row of `newdata` (in its order) and one
#'   column per posterior draw
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- fit_phaseR(model, data)
#' preds <- posterior_predict(fit)
#' }
#'
posterior_predict <- function(object, newdata = NULL, n_sims = NULL,
                              type = c("response", "phase"), phase = NULL) {

  type <- match.arg(type)

  if (!inherits(object, "phase_fit")) {
    stop("object must be a phase_fit", call. = FALSE)
  }

  data <- if (is.null(newdata)) object$data else newdata
  draws <- object$draws
  if (!is.null(n_sims) && n_sims < nrow(draws)) {
    draws <- draws[sample(nrow(draws), n_sims), , drop = FALSE]
  }
  s <- phase_summaries(object, data, draws)
  phase_names <- object$model$phases$names

  if (type == "phase") {
    if (!s$has_response) {
      stop("Phase probabilities given the outcomes need the outcome in ",
           "`newdata`.", call. = FALSE)
    }
    if (is.null(phase)) phase <- phase_names[length(phase_names)]
    k <- match(phase, phase_names)
    if (is.na(k)) stop("Unknown phase '", phase, "'.", call. = FALSE)
    return(t(matrix(s$filtered[, , k], nrow = nrow(draws))))
  }

  t(simulate_outcomes(s, data))
}


#' Simulate outcomes forward through the phases
#'
#' @param s Phase summaries from [phase_summaries()]
#' @param data The data the summaries were computed on
#' @return Matrix of simulated outcomes, one row per draw and one column per
#'   row of `data`
#' @keywords internal
simulate_outcomes <- function(s, data) {
  n_draws <- dim(s$mean)[1]
  draw <- seq_len(n_draws)
  K <- dim(s$mean)[3]
  out <- matrix(NA_real_, n_draws, nrow(data))
  for (rows in split(seq_len(nrow(data)), data$id)) {
    rows <- rows[order(data$time[rows])]
    phase <- rep(1L, n_draws)
    for (i in seq_along(rows)) {
      t <- rows[i]
      if (i > 1L) {
        can_move <- phase < K
        p_move <- numeric(n_draws)
        p_move[can_move] <- s$move[cbind(draw, t, phase)[can_move, , drop = FALSE]]
        phase <- phase + (stats::runif(n_draws) < p_move)
      }
      mu <- s$mean[cbind(draw, t, phase)]
      family <- s$family[phase]
      y <- numeric(n_draws)
      g <- family == 0L
      y[g] <- stats::rnorm(sum(g), mu[g], s$sigma[cbind(draw, phase)][g])
      g <- family == 1L
      y[g] <- stats::rpois(sum(g), mu[g])
      g <- family == 2L
      y[g] <- stats::rbinom(sum(g), 1L, mu[g])
      out[, t] <- y
    }
  }
  out
}


#' Counterfactual prediction
#'
#' Simulates outcomes under hypothetical covariate values with
#' [posterior_predict()]: every unit moves through the phases with the
#' transition probabilities of `newdata` and draws its outcomes from the
#' dynamics of its current phase.
#'
#' @param object A `phase_fit` object
#' @param newdata Data frame with counterfactual covariate values
#' @param n_sims Number of posterior draws to use
#'
#' @return A list with:
#'   - `mean`: Mean predicted outcome per observation
#'   - `sd`: SD of predicted outcome per observation
#'   - `lower`, `upper`: 2.5% and 97.5% quantiles per observation
#'   - `samples`: Matrix of samples (rows = rows of `newdata`, cols = draws)
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
#' The effect of setting `treatment_var` to `treat_value` against
#' `control_value` for every observation, on the expected outcome given the
#' covariates: for each posterior draw, the expected outcome integrates over
#' the phase each observation is in, so a treatment that changes the
#' transition probabilities acts through the phase as well as through the
#' dynamics. The effect is averaged over the observations of the data.
#'
#' @param object A `phase_fit` object
#' @param treatment_var Name of the treatment variable
#' @param treat_value Value representing treatment (default 1)
#' @param control_value Value representing control (default 0)
#' @param n_sims Number of posterior draws (default: all)
#'
#' @return A list with:
#'   - `ate`: Posterior mean of ATE
#'   - `ate_sd`: Posterior SD of ATE
#'   - `ate_ci`: 95% credible interval
#'   - `samples`: Vector of ATE samples, one per draw
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
  effects <- treatment_effects(object, treatment_var, treat_value,
                               control_value, n_sims)
  ate_samples <- rowMeans(effects)
  list(
    ate = mean(ate_samples),
    ate_sd = sd(ate_samples),
    ate_ci = quantile(ate_samples, c(0.025, 0.975)),
    samples = ate_samples
  )
}


#' Individual treatment effects
#'
#' The effect of the treatment on the expected outcome, as in [ate()],
#' averaged over the observations of each unit.
#'
#' @inheritParams ate
#'
#' @return A data frame with one row per unit: `id` and the posterior mean,
#'   SD and 95% interval of its effect
#'
#' @export
#'
ite <- function(object, treatment_var, treat_value = 1, control_value = 0,
                n_sims = NULL) {
  effects <- treatment_effects(object, treatment_var, treat_value,
                               control_value, n_sims)
  unit_ids <- unique(object$data$id)
  results <- lapply(unit_ids, function(uid) {
    unit_ite <- rowMeans(effects[, object$data$id == uid, drop = FALSE])
    data.frame(
      id = uid,
      ite_mean = mean(unit_ite),
      ite_sd = sd(unit_ite),
      ite_lower = quantile(unit_ite, 0.025),
      ite_upper = quantile(unit_ite, 0.975)
    )
  })
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}


#' Effects of a treatment on the expected outcome
#'
#' @inheritParams ate
#' @return Matrix with one row per posterior draw and one column per row of
#'   the data: expected outcome under treatment minus under control
#' @keywords internal
treatment_effects <- function(object, treatment_var, treat_value,
                              control_value, n_sims) {
  if (!inherits(object, "phase_fit")) {
    stop("object must be a phase_fit", call. = FALSE)
  }
  if (!treatment_var %in% names(object$data)) {
    stop(sprintf("Variable '%s' not found in data", treatment_var), call. = FALSE)
  }
  draws <- object$draws
  if (!is.null(n_sims) && n_sims < nrow(draws)) {
    draws <- draws[sample(nrow(draws), n_sims), , drop = FALSE]
  }
  expected <- function(value) {
    newdata <- object$data
    newdata[[treatment_var]] <- value
    s <- phase_summaries(object, newdata, draws)
    apply(s$marginal * s$mean, c(1, 2), sum)
  }
  expected(treat_value) - expected(control_value)
}
