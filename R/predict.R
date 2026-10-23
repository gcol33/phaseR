#' Predict from a phase model
#'
#' Predictions given the covariates of `newdata` alone, averaged over the
#' posterior draws: each unit starts in the first phase and moves through
#' the phases with the transition probabilities its covariates imply, and
#' each phase's outcome follows its dynamics, random effects included.
#' Random-effect levels that the fit has not seen contribute no effect.
#'
#' @param object A `phase_fit` object
#' @param newdata Optional new data for prediction, with columns `id`,
#'   `time` and the covariates of the model; defaults to the data of the fit
#' @param type `"response"` (default): the expected outcome, marginalised
#'   over the phase; `"phase"`: the probability of each phase;
#'   `"transition"`: the probability of each move between consecutive
#'   phases at each observation (zero at a unit's first observation)
#' @param summary Return summary statistics (default `TRUE`) or the draws
#'   (one row per draw, one column per row of `newdata`)
#' @param ... Additional arguments (unused)
#'
#' @return For `"response"`, a data frame of summaries (or the matrix of
#'   draws). For `"phase"`, a data frame with one column per phase holding
#'   the posterior mean probability. For `"transition"`, the summaries (or
#'   draws) of the single move of a two-phase model, or a named list of them
#'   with one element per move of a model with more phases.
#' @export
#'
predict.phase_fit <- function(object,
                              newdata = NULL,
                              type = c("response", "phase", "transition"),
                              summary = TRUE,
                              ...) {

  type <- match.arg(type)
  if (is.null(newdata)) newdata <- object$data
  s <- phase_summaries(object, newdata, object$draws)
  phase_names <- object$model$phases$names
  K <- length(phase_names)
  summarise <- function(draws) if (summary) summarize_predictions(draws) else draws

  switch(type,
    response = summarise(apply(s$marginal * s$mean, c(1, 2), sum)),
    phase = {
      probs <- apply(s$marginal, c(2, 3), mean)
      colnames(probs) <- phase_names
      as.data.frame(probs)
    },
    transition = {
      moves <- lapply(seq_len(K - 1L), function(j) {
        summarise(matrix(s$move[, , j], nrow = dim(s$move)[1]))
      })
      if (K == 2L) {
        moves[[1]]
      } else {
        stats::setNames(moves, paste0(phase_names[-K], "_to_", phase_names[-1]))
      }
    }
  )
}


#' Phase summaries of a fit on new data
#'
#' Evaluates the shared C++ model at every draw on the rows of `newdata`.
#'
#' @param fit A `phase_fit` object
#' @param newdata Data frame with `id`, `time` and the model's covariates;
#'   the outcome is needed only for the outcome-conditional summaries
#' @param draws Matrix of parameter draws (one row per draw)
#' @return List of arrays indexed `[draw, row of newdata, phase]`: `mean`
#'   (the outcome mean of each phase, on the response scale), `marginal`
#'   (the phase probability given the covariates), `predictive` (given the
#'   unit's earlier outcomes), `filtered` (given the unit's outcomes up to
#'   and including the row) and `move` (the probability of moving from each
#'   phase to the next at the row), plus `sigma` (`[draw, phase]`, the
#'   residual SDs), `family` (the family code of each phase) and
#'   `has_response`
#' @keywords internal
phase_summaries <- function(fit, newdata, draws) {
  pt <- prediction_target(fit, newdata)
  out <- phaseR_phase_summaries(pt$target, draws)
  rows <- match(seq_len(nrow(newdata)), pt$order)
  out <- lapply(out, function(a) a[, rows, , drop = FALSE])
  K <- pt$target$n_phases
  out$sigma <- exp(draws[, pt$target$sigma_offset + seq_len(K), drop = FALSE])
  out$family <- pt$target$family
  out$has_response <- pt$has_response
  out
}


#' Model description of a fit evaluated on new data
#'
#' The parameter layout is the fit's; the outcome, the units and the design
#' matrices come from `newdata`. Random-effect levels are matched to the
#' fit's by name; a level the fit has not seen gets no effect.
#'
#' @param fit A `phase_fit` object
#' @param newdata Data frame with `id`, `time` and the model's covariates
#' @return List with `target` (as from [posterior_target()]), `order` (the
#'   rows of `newdata` in the order of the target) and `has_response`
#' @keywords internal
prediction_target <- function(fit, newdata) {
  model <- fit$model
  if (!all(c("id", "time") %in% names(newdata))) {
    stop("`newdata` needs columns `id` and `time`.", call. = FALSE)
  }
  response <- model$dynamics[[1]]$response
  has_response <- response %in% names(newdata)
  if (!has_response) newdata[[response]] <- 0

  target <- posterior_target(prepare_stan_data(model, fit$data, prior()))
  new <- posterior_target(prepare_stan_data(model, newdata, prior()))

  target$y <- new$y
  target$unit_start <- new$unit_start
  target$unit_end <- new$unit_end
  for (part in c("transitions", "dynamics")) {
    for (c in seq_along(target[[part]])) {
      pred <- target[[part]][[c]]
      new_pred <- new[[part]][[c]]
      if (!identical(colnames(pred$X), colnames(new_pred$X))) {
        stop("`newdata` gives the design matrix columns ",
             paste(colnames(new_pred$X), collapse = ", "),
             " where the fit has ", paste(colnames(pred$X), collapse = ", "),
             ".", call. = FALSE)
      }
      pred$X <- new_pred$X
      for (r in seq_along(pred$re)) {
        new_term <- new_pred$re[[r]]
        level <- new_term$levels[new_term$group + 1L]
        group <- match(level, pred$re[[r]]$levels) - 1L
        pred$re[[r]]$group <- ifelse(is.na(group), -1L, group)
        pred$re[[r]]["Z"] <- list(new_term$Z)
      }
      target[[part]][[c]] <- pred
    }
  }
  list(target = target, order = order(newdata$id, newdata$time),
       has_response = has_response)
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
