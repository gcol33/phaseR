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
#' For phase models, the log-likelihood of each observation is that given
#' the earlier observations of its unit, marginalised over the latent
#' phases, at every posterior draw.
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
#'
#' The log-likelihood of each observation given the earlier observations of
#' its unit, marginalised over the latent phases, for every posterior draw
#' (random effects at their drawn values). The forward recursion is the one
#' both backends fit with, so every model family is covered.
#'
#' @param fit A `phase_fit` object
#' @return Matrix with one row per draw and one column per row of
#'   `fit$data`
#' @keywords internal
compute_pointwise_ll <- function(fit) {
  stan_data <- prepare_stan_data(fit$model, fit$data, prior())
  ll <- phaseR_pointwise_log_likelihood(posterior_target(stan_data),
                                        fit$draws)
  # prepare_stan_data() sorts the rows by unit and time.
  sorted <- order(fit$data$id, fit$data$time)
  ll[, match(seq_len(nrow(fit$data)), sorted), drop = FALSE]
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
