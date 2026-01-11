#' Compose a phase model
#'
#' @param phases A `phase_spec` object from `phases()`
#' @param ... `phase_transition` and `phase_dynamics` objects
#'
#' @return A `phase_model` object
#' @export
#'
#' @examples
#' model <- phase_model(
#'   phases("inactive", "active"),
#'   transition("inactive", "active", ~ treatment),
#'   dynamics("inactive", y ~ 1),
#'   dynamics("active", y ~ treatment)
#' )
#'
phase_model <- function(phases, ...) {

  if (!inherits(phases, "phase_spec")) {
    stop("`phases` must be a phase_spec object from phases()", call. = FALSE)
  }

  components <- list(...)

  transitions <- Filter(function(x) inherits(x, "phase_transition"), components)
  dynamics_list <- Filter(function(x) inherits(x, "phase_dynamics"), components)

  # Validate transitions reference declared phases
  for (tr in transitions) {
    if (!tr$from %in% phases$names) {
      stop(sprintf("Transition 'from' phase '%s' not in declared phases", tr$from),
           call. = FALSE)
    }
    if (!tr$to %in% phases$names) {
      stop(sprintf("Transition 'to' phase '%s' not in declared phases", tr$to),
           call. = FALSE)
    }
    # Check ordering (no backward transitions)
    if (phases$order[tr$from] >= phases$order[tr$to]) {
      stop(sprintf("Transition %s -> %s violates phase ordering", tr$from, tr$to),
           call. = FALSE)
    }
  }

  # Validate dynamics reference declared phases
  for (dyn in dynamics_list) {
    if (!dyn$phase %in% phases$names) {
      stop(sprintf("Dynamics phase '%s' not in declared phases", dyn$phase),
           call. = FALSE)
    }
  }

  # Check all non-absorbing phases have dynamics
  non_absorbing <- setdiff(phases$names, phases$absorbing)
  phases_with_dynamics <- vapply(dynamics_list, function(x) x$phase, character(1))
  missing_dynamics <- setdiff(non_absorbing, phases_with_dynamics)
  if (length(missing_dynamics) > 0) {
    stop(sprintf("Missing dynamics for non-absorbing phases: %s",
                 paste(missing_dynamics, collapse = ", ")), call. = FALSE)
  }

  structure(
    list(
      phases = phases,
      transitions = transitions,
      dynamics = dynamics_list
    ),
    class = "phase_model"
  )
}

#' @export
print.phase_model <- function(x, ...) {
  cat("Phase Model\n")
  cat(sprintf("  Phases: %s\n", paste(x$phases$names, collapse = " -> ")))
  cat(sprintf("  Transitions: %d\n", length(x$transitions)))
  cat(sprintf("  Dynamics: %d\n", length(x$dynamics)))
  invisible(x)
}
