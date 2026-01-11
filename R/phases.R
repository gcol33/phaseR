#' Declare ordered phases
#'
#' @param ... Character names of phases in order (e.g., "inactive", "active")
#' @param absorbing Character vector of phase names that are absorbing states.
#'   Defaults to the last declared phase.
#'
#' @return A `phase_spec` object
#' @export
#'
#' @examples
#' # Two-phase model
#' phases("inactive", "active")
#'
#' # Multi-phase with explicit absorbing state
#' phases("susceptible", "infected", "recovered", absorbing = "recovered")
#'
phases <- function(..., absorbing = NULL) {

  phase_names <- c(...)

  # Validation
  if (length(phase_names) < 2) {
    stop("At least 2 phases required", call. = FALSE)
  }
  if (any(duplicated(phase_names))) {
    stop("Phase names must be unique", call. = FALSE)
  }
  if (!all(vapply(phase_names, is.character, logical(1)))) {
    stop("Phase names must be character strings", call. = FALSE)
  }

  # Default absorbing: last phase
  if (is.null(absorbing)) {
    absorbing <- phase_names[length(phase_names)]
  }
  if (!all(absorbing %in% phase_names)) {
    stop("Absorbing phases must be declared phases", call. = FALSE)
  }

  structure(
    list(
      names = phase_names,
      n_phases = length(phase_names),
      absorbing = absorbing,
      order = setNames(seq_along(phase_names) - 1L, phase_names)
    ),
    class = "phase_spec"
  )
}

#' @export
print.phase_spec <- function(x, ...) {
  cat("Phase specification:\n")
  for (i in seq_along(x$names)) {
    marker <- if (x$names[i] %in% x$absorbing) " [absorbing]" else ""
    cat(sprintf("  %d. %s%s\n", i, x$names[i], marker))
  }
  invisible(x)
}
