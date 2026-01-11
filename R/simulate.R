#' Simulate data from a phase model
#'
#' @param n_units Number of units
#' @param n_times Number of time points per unit
#' @param beta_trans Transition coefficients
#' @param beta_0 Phase 0 dynamics coefficients
#' @param beta_1 Phase 1 dynamics coefficients
#' @param sigma_0 Phase 0 residual SD
#' @param sigma_1 Phase 1 residual SD
#' @param seed Random seed
#'
#' @return A data frame suitable for `fit_phaseR()`
#' @export
#'
#' @examples
#' dat <- sim_phaseR(n_units = 100, n_times = 5)
#' head(dat)
#'
sim_phaseR <- function(n_units = 100,
                       n_times = 5,
                       beta_trans = c(-1, 0.5),  # intercept, x effect
                       beta_0 = c(10, 0),         # phase 0: intercept, x effect
                       beta_1 = c(8, -0.5),       # phase 1: intercept, x effect
                       sigma_0 = 1,
                       sigma_1 = 1.5,
                       seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Generate covariates
  data_list <- lapply(seq_len(n_units), function(i) {

    x <- rnorm(n_times)
    phase <- integer(n_times)
    y <- numeric(n_times)

    # Initial phase always 0
    phase[1] <- 0
    y[1] <- beta_0[1] + beta_0[2] * x[1] + rnorm(1, 0, sigma_0)

    for (t in 2:n_times) {
      if (phase[t - 1] == 0) {
        # Transition probability
        p_trans <- plogis(beta_trans[1] + beta_trans[2] * x[t])
        phase[t] <- rbinom(1, 1, p_trans)
      } else {
        # Absorbing
        phase[t] <- 1
      }

      # Generate outcome
      if (phase[t] == 0) {
        y[t] <- beta_0[1] + beta_0[2] * x[t] + rnorm(1, 0, sigma_0)
      } else {
        y[t] <- beta_1[1] + beta_1[2] * x[t] + rnorm(1, 0, sigma_1)
      }
    }

    data.frame(
      id = i,
      time = seq_len(n_times),
      x = x,
      y = y,
      true_phase = phase
    )
  })

  do.call(rbind, data_list)
}
