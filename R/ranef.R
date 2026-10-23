#' Extract fixed effects
#'
#' Returns posterior summaries of fixed effect coefficients.
#'
#' @param object A `phase_fit` object
#' @param component Which component: "all" (default), "transition", or "dynamics"
#' @param phase For dynamics, which phase (NULL for all)
#'
#' @return A data frame with fixed effect estimates
#' @export
#'
fixef <- function(object, component = c("all", "transition", "dynamics"),
                  phase = NULL) {

  component <- match.arg(component)

  if (!inherits(object, "phase_fit")) {
    stop("Object must be a phase_fit", call. = FALSE)
  }

  draws <- object$draws
  param_names <- object$param_names
  phase_names <- object$model$phases$names

  results <- list()

  # Transition coefficients
  if (component %in% c("all", "transition")) {
    trans_idx <- grep("^beta_trans", param_names)
    if (length(trans_idx) > 0) {
      results$transition <- summarize_fixef(draws[, trans_idx, drop = FALSE])
    }
  }

  # Dynamics coefficients
  if (component %in% c("all", "dynamics")) {
    phases_to_use <- if (is.null(phase)) phase_names else phase
    for (p in phases_to_use) {
      pattern <- sprintf("^beta_%s_", p)
      dyn_idx <- grep(pattern, param_names)
      if (length(dyn_idx) > 0) {
        results[[paste0("dynamics_", p)]] <- summarize_fixef(draws[, dyn_idx, drop = FALSE])
      }
    }
  }

  if (length(results) == 1) {
    return(results[[1]])
  }
  results
}


#' Summarize fixed effects draws
#' @keywords internal
summarize_fixef <- function(draws) {
  summaries <- t(apply(draws, 2, function(x) {
    c(mean = mean(x),
      sd = sd(x),
      `2.5%` = unname(quantile(x, 0.025)),
      `50%` = unname(median(x)),
      `97.5%` = unname(quantile(x, 0.975)))
  }))

  data.frame(
    parameter = rownames(summaries),
    summaries,
    row.names = NULL
  )
}


#' Extract random effects
#'
#' Posterior summaries of the random-effect coefficients of every level:
#' the effects themselves for intercept and uncorrelated terms, and
#' `diag(sigma) L z` for correlated terms.
#'
#' @param object A `phase_fit` object
#' @param component Which component: "transition", "dynamics", or "all" (default)
#' @param phase For dynamics RE, which phase (ignored for transition)
#' @param summary If TRUE (default), return posterior summary; if FALSE, return draws
#'
#' @return With `summary = TRUE`, a data frame with one row per predictor,
#'   grouping term, level and coefficient and the posterior mean, SD and
#'   quantiles. With `summary = FALSE`, a named list of arrays of draws
#'   indexed `[draw, level, coefficient]`, one per predictor and term.
#' @export
#'
#' @examples
#' \dontrun{
#' ranef(fit)
#' ranef(fit, component = "transition")
#' ranef(fit, component = "dynamics", phase = "active")
#' }
#'
ranef <- function(object, component = c("all", "transition", "dynamics"),
                  phase = NULL, summary = TRUE) {

  component <- match.arg(component)

  if (!inherits(object, "phase_fit")) {
    stop("Object must be a phase_fit", call. = FALSE)
  }

  terms <- random_effect_terms(object)
  if (length(terms) == 0) {
    message("Model has no random effects")
    return(NULL)
  }
  keep <- vapply(terms, function(t) {
    (component == "all" || t$component == component) &&
      (is.null(phase) || t$component != "dynamics" || t$phase %in% phase)
  }, logical(1))
  terms <- terms[keep]
  if (length(terms) == 0) return(NULL)

  if (!summary) {
    return(stats::setNames(lapply(terms, function(t) t$coefficients),
                           vapply(terms, function(t) t$label, character(1))))
  }
  out <- do.call(rbind, lapply(terms, function(t) {
    b <- t$coefficients
    cells <- expand.grid(level = dimnames(b)[[2]], coef = dimnames(b)[[3]],
                         stringsAsFactors = FALSE)
    values <- matrix(b, nrow = dim(b)[1])
    data.frame(
      predictor = t$predictor,
      term = t$name,
      cells,
      summarize_draws(values),
      row.names = NULL,
      check.names = FALSE
    )
  }))
  out
}


#' Random-effect terms of a fit with their posterior draws
#'
#' @param object A `phase_fit` object
#' @return List with one element per term, in the order of the parameter
#'   vector: `predictor` (e.g. "transition", "dynamics_active"),
#'   `component`, `phase`, `name` (grouping factor), `label`,
#'   `coefficients` (draws `[draw, level, coefficient]`), `sd`
#'   (`[draw, coefficient]`) and `cor` (`[draw, coefficient, coefficient]`)
#' @keywords internal
random_effect_terms <- function(object) {
  stan_data <- prepare_stan_data(object$model, object$data, prior())
  target <- posterior_target(stan_data)
  values <- phaseR_random_effects(target, object$draws)
  phase_names <- stan_data$phase_names
  K <- length(phase_names)
  predictors <- c(
    lapply(seq_len(K - 1L), function(j) {
      list(pred = target$transitions[[j]], component = "transition",
           phase = phase_names[j],
           predictor = if (K == 2L) "transition" else
             paste0("transition_", phase_names[j], "_to_", phase_names[j + 1L]))
    }),
    lapply(seq_len(K), function(k) {
      list(pred = target$dynamics[[k]], component = "dynamics",
           phase = phase_names[k], predictor = paste0("dynamics_", phase_names[k]))
    })
  )
  out <- list()
  for (p in predictors) {
    for (term in p$pred$re) {
      v <- values[[length(out) + 1L]]
      dimnames(v$coefficients) <- list(NULL, term$levels, term$coef_names)
      dimnames(v$sd) <- list(NULL, term$coef_names)
      dimnames(v$cor) <- list(NULL, term$coef_names, term$coef_names)
      out[[length(out) + 1L]] <- c(
        list(predictor = p$predictor, component = p$component,
             phase = p$phase, name = term$name,
             label = paste0(p$predictor, ":", term$name)),
        v)
    }
  }
  out
}


#' Posterior summaries of the columns of a draws matrix
#' @keywords internal
summarize_draws <- function(draws) {
  data.frame(
    mean = colMeans(draws),
    sd = apply(draws, 2, stats::sd),
    `2.5%` = apply(draws, 2, stats::quantile, 0.025),
    `50%` = apply(draws, 2, stats::median),
    `97.5%` = apply(draws, 2, stats::quantile, 0.975),
    check.names = FALSE,
    row.names = NULL
  )
}


#' Extract variance components
#'
#' Posterior summaries of the SD of every random-effect coefficient and the
#' posterior mean of the correlation matrix of each term (the identity for
#' intercept and uncorrelated terms).
#'
#' @param object A `phase_fit` object
#'
#' @return A list with one element per predictor and term, each with `sd`
#'   (a data frame with one row per coefficient) and `cor` (a matrix)
#' @export
VarCorr <- function(object) {

  if (!inherits(object, "phase_fit")) {
    stop("Object must be a phase_fit", call. = FALSE)
  }
  terms <- random_effect_terms(object)
  if (length(terms) == 0) {
    message("Model has no random effects")
    return(NULL)
  }
  out <- lapply(terms, function(t) {
    list(
      sd = data.frame(coef = colnames(t$sd), summarize_draws(t$sd),
                      check.names = FALSE),
      cor = apply(t$cor, c(2, 3), mean)
    )
  })
  names(out) <- vapply(terms, function(t) t$label, character(1))
  structure(out, class = "VarCorr_phaseR")
}


#' Print method for VarCorr
#'
#' @param x A `VarCorr_phaseR` object
#' @param ... Additional arguments (unused)
#'
#' @return Invisibly returns the input object
#' @export
print.VarCorr_phaseR <- function(x, ...) {
  for (name in names(x)) {
    cat("Random effect:", name, "\n")
    print(x[[name]]$sd, digits = 3, row.names = FALSE)
    if (nrow(x[[name]]$sd) > 1) {
      cat("  Correlations (posterior mean):\n")
      print(round(x[[name]]$cor, 3))
    }
    cat("\n")
  }
  invisible(x)
}
