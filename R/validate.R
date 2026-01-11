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

  # Check response variables exist
  for (dyn in model$dynamics) {
    if (!dyn$response %in% names(data)) {
      stop(sprintf("Response variable '%s' not found in data", dyn$response),
           call. = FALSE)
    }
  }

  # Check transition covariates exist
  for (tr in model$transitions) {
    covars <- all.vars(tr$formula)
    missing <- setdiff(covars, names(data))
    if (length(missing) > 0) {
      stop(sprintf("Transition covariates not found in data: %s",
                   paste(missing, collapse = ", ")), call. = FALSE)
    }
  }

  # Check dynamics covariates exist
  for (dyn in model$dynamics) {
    covars <- all.vars(dyn$formula)
    missing <- setdiff(covars, names(data))
    if (length(missing) > 0) {
      stop(sprintf("Dynamics covariates not found in data: %s",
                   paste(missing, collapse = ", ")), call. = FALSE)
    }
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
