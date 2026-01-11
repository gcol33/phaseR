#' Specify a phase transition
#'
#' @param from Character name of the origin phase
#' @param to Character name of the destination phase
#' @param formula A one-sided formula specifying covariates for the transition
#'   probability (e.g., `~ x1 + x2`). Default is `~ 1` (intercept only).
#'
#' @return A `phase_transition` object
#' @export
#'
#' @examples
#' # Transition with covariates
#' transition("inactive", "active", ~ treatment + age)
#'
#' # Intercept-only transition
#' transition("stage1", "stage2")
#'
transition <- function(from, to, formula = ~ 1) {

  if (!is.character(from) || length(from) != 1) {
    stop("`from` must be a single character string", call. = FALSE)
  }
  if (!is.character(to) || length(to) != 1) {
    stop("`to` must be a single character string", call. = FALSE)
  }
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula", call. = FALSE)
  }

  if (length(formula) == 3) {
    stop("Transition formula must be one-sided (no response)", call. = FALSE)
  }

  structure(
    list(
      from = from,
      to = to,
      formula = formula
    ),
    class = "phase_transition"
  )
}

#' @export
print.phase_transition <- function(x, ...) {
  cat(sprintf("Transition: %s -> %s\n", x$from, x$to))
  cat(sprintf("  Formula: %s\n", deparse(x$formula)))
  invisible(x)
}
