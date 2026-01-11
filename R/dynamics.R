#' Specify phase-specific dynamics
#'
#' @param phase Character name of the phase
#' @param formula A two-sided formula specifying the outcome model
#'   (e.g., `y ~ x1 + x2`).
#' @param family A family object (e.g., `gaussian()`, `poisson()`).
#'   Default is `gaussian()`.
#'
#' @return A `phase_dynamics` object
#' @export
#'
#' @examples
#' # Gaussian outcome
#' dynamics("active", y ~ treatment + time)
#'
#' # Poisson outcome
#' dynamics("infected", cases ~ age + region, family = poisson())
#'
dynamics <- function(phase, formula, family = gaussian()) {

  if (!is.character(phase) || length(phase) != 1) {
    stop("`phase` must be a single character string", call. = FALSE)
  }
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula", call. = FALSE)
  }
  if (length(formula) != 3) {
    stop("Dynamics formula must be two-sided (response ~ predictors)", call. = FALSE)
  }

  # Check for delta/difference terms (forbidden)
  formula_text <- deparse(formula)
  if (grepl("diff\\(|Diff\\(|delta|Delta|\\s-\\s", formula_text)) {
    stop("Dynamics formula cannot contain difference/delta terms. ",
         "phaseR models outcomes, not changes.", call. = FALSE)
  }

  structure(
    list(
      phase = phase,
      formula = formula,
      family = family,
      response = all.vars(formula[[2]])
    ),
    class = "phase_dynamics"
  )
}

#' @export
print.phase_dynamics <- function(x, ...) {
  cat(sprintf("Dynamics for phase '%s':\n", x$phase))
  cat(sprintf("  Formula: %s\n", deparse(x$formula)))
  cat(sprintf("  Family: %s\n", x$family$family))
  invisible(x)
}
