#' phaseR: Phase-Structured Inference for Longitudinal Change
#'
#' A framework for analysing change in longitudinal data that treats
#' "before" and "after" as regimes (phases), not time points. The key insight
#' is that **inference conditional on realized change is not a valid estimand**.
#'
#' @section The Problem:
#' Conventional analyses often:
#' 1. Observe when a unit "changes" (e.g., disease onset, treatment adoption)
#' 2. Compare outcomes before vs. after the change
#' 3. Attribute the difference to the "effect" of the change
#'
#' This approach conflates the true phase effect with selection effects from
#' the transition mechanism.
#'
#' @section The Solution:
#' phaseR models the complete generative process:
#' \itemize{
#'   \item **Transition model**: How units move between phases
#'   \item **Dynamics model**: How outcomes behave within each phase
#' }
#'
#' By jointly modeling both components and integrating over phase uncertainty,
#' phaseR produces valid causal estimates.
#'
#' @section Core Functions:
#' **Model specification:**
#' \itemize{
#'   \item \code{\link{phases}}: Declare ordered phases
#'   \item \code{\link{transition}}: Specify transition model
#'   \item \code{\link{dynamics}}: Specify phase-specific dynamics
#'   \item \code{\link{phase_model}}: Compose a complete model
#' }
#'
#' **Fitting:**
#' \itemize{
#'   \item \code{\link{fit_phaseR}}: Fit the model via HMC
#' }
#'
#' **Inference:**
#' \itemize{
#'   \item \code{\link{posterior_predict}}: Posterior predictive simulation
#'   \item \code{\link{counterfactual}}: Counterfactual prediction
#'   \item \code{\link{ate}}: Average treatment effects
#'   \item \code{\link{ite}}: Individual treatment effects
#' }
#'
#' **Simulation:**
#' \itemize{
#'   \item \code{\link{sim_phaseR}}: Simulate 2-phase Gaussian data
#'   \item \code{\link{sim_phaseR_k}}: Simulate k-phase data
#'   \item \code{\link{sim_phaseR_poisson}}: Simulate Poisson data
#' }
#'
#' @section Design Principles:
#' \itemize{
#'   \item Change (delta) is never modelled directly
#'   \item Phases are regimes, not time points
#'   \item Phase membership is latent or probabilistic
#'   \item Transitions are modelled, not conditioned on
#' }
#'
#' @useDynLib phaseR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats acf as.formula density dnorm gaussian median model.matrix
#'   plogis predict quantile rbinom rnorm rpois runif sd setNames update var
#' @importFrom graphics lines par
#' @importFrom grDevices rainbow
#' @keywords internal
"_PACKAGE"
