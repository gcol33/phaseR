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
  # Handle both old naming (u_trans_level) and new naming (u_trans_groupvar_level)
  # Exclude log_sigma parameters
  if (grepl("^trans", type)) {
    pattern <- "^u_trans_[^l]"  # Match u_trans_ but not log_sigma
  } else if (grepl("^dyn_", type)) {
    # Extract phase name
    phase_name <- sub("^dyn_", "", type)
    pattern <- sprintf("^v_%s_[^l]", phase_name)
  } else {
    pattern <- sprintf("^u_%s_|^v_%s_", type, type)
  }

  idx <- grep(pattern, fit$param_names)
  # Also exclude log_sigma parameters
  sigma_idx <- grep("^log_sigma", fit$param_names)
  idx <- setdiff(idx, sigma_idx)

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
#' Returns variance-covariance information for random effects.
#' For random intercepts, returns the SD. For random slopes,
#' returns the full covariance matrix (correlated) or diagonal
#' variances (uncorrelated).
#'
#' @param object A `phase_fit` object
#'
#' @return A list with variance component estimates by term
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

  results <- list()

  # Find variance parameters - both log_sigma and L (Cholesky)
  sigma_idx <- grep("^log_sigma_.*_re", object$param_names)
  L_idx <- grep("^L_", object$param_names)

  # Process log_sigma parameters (intercept-only or uncorrelated)
  if (length(sigma_idx) > 0) {
    draws <- object$draws[, sigma_idx, drop = FALSE]
    sd_draws <- exp(draws)

    for (i in seq_along(sigma_idx)) {
      param_name <- object$param_names[sigma_idx[i]]
      comp_name <- gsub("log_sigma_(.*)_re.*", "\\1", param_name)
      group_name <- gsub(".*_re_(.*)$", "\\1", param_name)
      # Check if this is a slope-specific sigma
      if (grepl("_re_.*_", param_name)) {
        coef_name <- gsub(".*_re_[^_]+_(.*)$", "\\1", param_name)
        full_name <- paste(comp_name, group_name, coef_name, sep = "_")
      } else {
        full_name <- paste(comp_name, group_name, sep = "_")
      }

      sd_vec <- sd_draws[, i]
      results[[full_name]] <- list(
        type = "sd",
        sd = mean(sd_vec),
        sd.sd = sd(sd_vec),
        `2.5%` = unname(quantile(sd_vec, 0.025)),
        `97.5%` = unname(quantile(sd_vec, 0.975))
      )
    }
  }

  # Process Cholesky parameters (correlated slopes)
  if (length(L_idx) > 0) {
    L_names <- object$param_names[L_idx]

    # Group by term (e.g., L_baseline_re_site_1_1, L_baseline_re_site_2_1, ...)
    # Extract unique term identifiers
    term_patterns <- unique(gsub("^(L_[^_]+_re_[^_]+)_.*", "\\1", L_names))

    for (term_pat in term_patterns) {
      term_L_idx <- grep(paste0("^", term_pat, "_"), object$param_names)
      n_L <- length(term_L_idx)

      # Determine dimension from number of Cholesky elements
      # n_L = dim * (dim + 1) / 2, so dim = (-1 + sqrt(1 + 8*n_L)) / 2
      dim <- as.integer((-1 + sqrt(1 + 8 * n_L)) / 2)

      # Extract L draws and compute covariance: Sigma = L * L'
      L_draws <- object$draws[, term_L_idx, drop = FALSE]
      n_samples <- nrow(L_draws)

      # Compute covariance matrix for each sample
      cov_samples <- array(0, dim = c(dim, dim, n_samples))

      for (s in seq_len(n_samples)) {
        # Reconstruct L matrix
        L <- matrix(0, dim, dim)
        k <- 1
        for (j in seq_len(dim)) {
          for (i in j:dim) {
            L[i, j] <- L_draws[s, k]
            k <- k + 1
          }
        }
        cov_samples[, , s] <- L %*% t(L)
      }

      # Compute posterior mean covariance
      cov_mean <- apply(cov_samples, c(1, 2), mean)
      cov_sd <- apply(cov_samples, c(1, 2), sd)

      # Compute correlation matrix
      sds <- sqrt(diag(cov_mean))
      cor_mean <- cov_mean / outer(sds, sds)

      comp_name <- gsub("^L_([^_]+)_re_.*", "\\1", term_pat)
      group_name <- gsub("^L_[^_]+_re_(.*)$", "\\1", term_pat)
      full_name <- paste(comp_name, group_name, sep = "_")

      results[[full_name]] <- list(
        type = "cov",
        cov = cov_mean,
        cov.sd = cov_sd,
        cor = cor_mean,
        sd = sds
      )
    }
  }

  if (length(results) == 0) return(NULL)

  structure(results, class = "VarCorr_phaseR")
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
    item <- x[[name]]

    if (item$type == "sd") {
      cat(sprintf("  StdDev: %.4f (SE: %.4f)\n", item$sd, item$sd.sd))
      cat(sprintf("  95%% CI: [%.4f, %.4f]\n", item$`2.5%`, item$`97.5%`))
    } else if (item$type == "cov") {
      cat("  StdDevs:\n")
      for (i in seq_along(item$sd)) {
        cat(sprintf("    [%d]: %.4f\n", i, item$sd[i]))
      }
      cat("  Correlations:\n")
      print(round(item$cor, 3))
    }
    cat("\n")
  }

  invisible(x)
}
