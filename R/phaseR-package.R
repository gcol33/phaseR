#' @keywords internal
"_PACKAGE"

#' phaseR: Phase-Structured Inference for Longitudinal Change
#'
#' A framework for analysing change in longitudinal data that treats
#' "before" and "after" as regimes, not time points.
#'
#' @section Core functions:
#' \itemize{
#'   \item \code{\link{phases}}: Declare ordered phases
#'   \item \code{\link{transition}}: Specify transition model
#'   \item \code{\link{dynamics}}: Specify phase-specific dynamics
#'   \item \code{\link{phase_model}}: Compose a complete model
#'   \item \code{\link{fit_phaseR}}: Fit the model via HMC
#'   \item \code{\link{sim_phaseR}}: Simulate data from a phase model
#' }
#'
#' @section Design principles:
#' \itemize{
#'   \item Change (delta) is never modelled directly
#'   \item Phases are regimes, not time points
#'   \item Phase membership is latent or probabilistic
#'   \item Transitions are modelled, not observed
#' }
#'
#' @useDynLib phaseR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @docType package
#' @name phaseR-package
NULL
