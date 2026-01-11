#' Validate phase model against data
#'
#' @param model A `phase_model` object
#' @param data A data frame
#'
#' @return Invisibly returns TRUE if valid, otherwise stops with error
#' @keywords internal
#'
validate_model_data <- function(model, data) {

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame", call. = FALSE)
  }

  # Check required columns
  required <- c("id", "time")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(sprintf("Data missing required columns: %s",
                 paste(missing, collapse = ", ")), call. = FALSE)
  }

  # Check id and time are not NA
  if (any(is.na(data$id))) {
    stop("Column 'id' contains NA values", call. = FALSE)
  }
  if (any(is.na(data$time))) {
    stop("Column 'time' contains NA values", call. = FALSE)
  }

  # Check minimum observations per unit
  obs_per_unit <- table(data$id)
  if (any(obs_per_unit < 2)) {
    bad_units <- names(obs_per_unit)[obs_per_unit < 2]
    stop(sprintf("Units with fewer than 2 observations: %s. ",
                 "Phase models require at least 2 time points per unit.",
                 paste(head(bad_units, 5), collapse = ", ")), call. = FALSE)
  }

  # Check response variables exist
  for (dyn in model$dynamics) {
    if (!dyn$response %in% names(data)) {
      stop(sprintf("Response variable '%s' not found in data", dyn$response),
           call. = FALSE)
    }
  }

  # Check for NA in response
  response_var <- model$dynamics[[1]]$response
  if (any(is.na(data[[response_var]]))) {
    n_na <- sum(is.na(data[[response_var]]))
    stop(sprintf("Response variable '%s' contains %d NA values. ",
                 "Remove or impute before fitting.",
                 response_var, n_na), call. = FALSE)
  }

  # Check transition covariates exist
  for (tr in model$transitions) {
    covars <- all.vars(tr$formula)
    missing <- setdiff(covars, names(data))
    if (length(missing) > 0) {
      stop(sprintf("Transition covariates not found in data: %s",
                   paste(missing, collapse = ", ")), call. = FALSE)
    }

    # Check for NA in covariates
    for (cov in covars) {
      if (any(is.na(data[[cov]]))) {
        n_na <- sum(is.na(data[[cov]]))
        stop(sprintf("Transition covariate '%s' contains %d NA values",
                     cov, n_na), call. = FALSE)
      }
    }
  }

  # Check dynamics covariates exist
  for (dyn in model$dynamics) {
    covars <- all.vars(dyn$formula)
    covars <- setdiff(covars, dyn$response)  # Exclude response
    missing <- setdiff(covars, names(data))
    if (length(missing) > 0) {
      stop(sprintf("Dynamics covariates not found in data: %s",
                   paste(missing, collapse = ", ")), call. = FALSE)
    }

    # Check for NA in covariates
    for (cov in covars) {
      if (any(is.na(data[[cov]]))) {
        n_na <- sum(is.na(data[[cov]]))
        stop(sprintf("Dynamics covariate '%s' contains %d NA values",
                     cov, n_na), call. = FALSE)
      }
    }
  }

  # Check time is sorted within units
  data_sorted <- data[order(data$id, data$time), ]
  if (!identical(data$id, data_sorted$id) || !identical(data$time, data_sorted$time)) {
    message("Note: Data will be sorted by id and time for fitting.")
  }

  # Warn if response appears in transition formula
  response_vars <- unique(vapply(model$dynamics, function(x) x$response, character(1)))
  for (tr in model$transitions) {
    tr_vars <- all.vars(tr$formula)
    overlap <- intersect(response_vars, tr_vars)
    if (length(overlap) > 0) {
      warning(sprintf("Response variable '%s' appears in transition formula. ",
                      "This may indicate outcome-dependent phase definition.",
                      overlap[1]), call. = FALSE)
    }
  }

  invisible(TRUE)
}


#' Validate phase model specification
#'
#' @param model A `phase_model` object
#'
#' @return Invisibly returns TRUE if valid, otherwise stops with error
#' @keywords internal
#'
validate_model <- function(model) {

  if (!inherits(model, "phase_model")) {
    stop("Expected a phase_model object. Use phase_model() to create one.",
         call. = FALSE)
  }

  # Check phases
  if (is.null(model$phases) || model$phases$n_phases < 2) {
    stop("Model must have at least 2 phases", call. = FALSE)
  }

  # Check transitions
  if (length(model$transitions) == 0) {
    stop("Model must have at least one transition", call. = FALSE)
  }

  # Check dynamics
  if (length(model$dynamics) == 0) {
    stop("Model must have at least one dynamics specification", call. = FALSE)
  }

  # Check all phases have dynamics
  phase_names <- model$phases$names
  dynamics_phases <- vapply(model$dynamics, function(d) d$phase, character(1))
  missing_dyn <- setdiff(phase_names, dynamics_phases)
  if (length(missing_dyn) > 0) {
    stop(sprintf("Missing dynamics specification for phase(s): %s",
                 paste(missing_dyn, collapse = ", ")), call. = FALSE)
  }

  invisible(TRUE)
}
