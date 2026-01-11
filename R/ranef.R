#' Extract random effects
#'
#' @param object A `phase_fit` object
#' @param component Which component: "transition", "dynamics", or "all" (default)
#' @param phase For dynamics RE, which phase (ignored for transition)
#' @param summary If TRUE (default), return posterior summary; if FALSE, return draws
#'
#' @return A data frame with random effect estimates
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

  if (!isTRUE(object$has_re)) {
    message("Model has no random effects")
    return(NULL)
  }

  results <- list()

  if (component %in% c("all", "transition") && isTRUE(object$has_trans_re)) {
    re_draws <- extract_re_draws(object, "trans")
    if (summary) {
      results$transition <- summarize_re(re_draws, object$trans_re_levels)
    } else {
      results$transition <- re_draws
    }
  }

  if (component %in% c("all", "dynamics") && isTRUE(object$has_dyn_re)) {
    phases_to_use <- if (is.null(phase)) object$model$phases$names else phase
    for (p in phases_to_use) {
      re_draws <- extract_re_draws(object, paste0("dyn_", p))
      name <- sprintf("dynamics_%s", p)
      if (!is.null(re_draws)) {
        if (summary) {
          results[[name]] <- summarize_re(re_draws, object$dyn_re_levels[[p]])
        } else {
          results[[name]] <- re_draws
        }
      }
    }
  }

  if (length(results) == 0) {
    return(NULL)
  }
  if (length(results) == 1) {
    return(results[[1]])
  }
  results
}


#' Extract RE draws from fitted model
#' @keywords internal
extract_re_draws <- function(fit, type) {
  pattern <- sprintf("^u_%s_|^v_%s_", type, type)
  idx <- grep(pattern, fit$param_names)
  if (length(idx) == 0) return(NULL)
  fit$draws[, idx, drop = FALSE]
}


#' Summarize random effects draws
#' @keywords internal
summarize_re <- function(draws, levels) {
  if (is.null(draws) || ncol(draws) == 0) return(NULL)

  summaries <- t(apply(draws, 2, function(x) {
    c(mean = mean(x),
      sd = sd(x),
      `2.5%` = unname(quantile(x, 0.025)),
      `50%` = unname(median(x)),
      `97.5%` = unname(quantile(x, 0.975)))
  }))

  if (is.null(levels) || length(levels) != nrow(summaries)) {
    levels <- paste0("group_", seq_len(nrow(summaries)))
  }

  data.frame(
    group = levels,
    summaries,
    row.names = NULL
  )
}


#' Extract variance components
#'
#' @param object A `phase_fit` object
#'
#' @return A data frame with variance component estimates
#' @export
#'
VarCorr <- function(object) {

  if (!inherits(object, "phase_fit")) {
    stop("Object must be a phase_fit", call. = FALSE)
  }

  if (!isTRUE(object$has_re)) {
    message("Model has no random effects")
    return(NULL)
  }

  # Find variance parameters
  sigma_idx <- grep("^log_sigma_.*_re", object$param_names)
  if (length(sigma_idx) == 0) return(NULL)

  draws <- object$draws[, sigma_idx, drop = FALSE]

  # Convert from log to actual SD
  sd_draws <- exp(draws)

  summaries <- t(apply(sd_draws, 2, function(x) {
    c(sd = mean(x),
      sd.sd = sd(x),
      `2.5%` = unname(quantile(x, 0.025)),
      `97.5%` = unname(quantile(x, 0.975)))
  }))

  # Clean up names
  comp_names <- gsub("log_sigma_(.*)_re", "\\1", colnames(draws))

  data.frame(
    component = comp_names,
    summaries,
    row.names = NULL
  )
}
