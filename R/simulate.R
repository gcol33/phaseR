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


#' Simulate data from a k-phase model
#'
#' @param n_units Number of units
#' @param n_times Number of time points per unit
#' @param k Number of phases (must be >= 2)
#' @param beta_trans List of transition coefficient vectors (length k-1)
#' @param beta_dyn List of dynamics coefficient vectors (length k)
#' @param sigma Vector of residual SDs (length k)
#' @param seed Random seed
#'
#' @return A data frame suitable for `fit_phaseR()`
#' @export
#'
#' @examples
#' # 3-phase model
#' dat <- sim_phaseR_k(n_units = 50, n_times = 8, k = 3)
#' table(dat$true_phase)
#'
sim_phaseR_k <- function(n_units = 100,
                          n_times = 8,
                          k = 3,
                          beta_trans = NULL,
                          beta_dyn = NULL,
                          sigma = NULL,
                          seed = NULL) {

  if (k < 2) stop("k must be at least 2", call. = FALSE)

  if (!is.null(seed)) set.seed(seed)

  # Default parameters if not provided
  if (is.null(beta_trans)) {
    # k-1 transitions, each with intercept only that decreases
    beta_trans <- lapply(seq_len(k - 1), function(i) c(-1.5 + 0.3 * i))
  }
  if (is.null(beta_dyn)) {
    # k phases, intercept decreases by 2 per phase
    beta_dyn <- lapply(seq_len(k), function(p) c(12 - 2 * (p - 1)))
  }
  if (is.null(sigma)) {
    sigma <- rep(1, k)
  }

  # Generate data
  data_list <- lapply(seq_len(n_units), function(i) {

    x <- rnorm(n_times)
    phase <- integer(n_times)
    y <- numeric(n_times)

    # Initial phase always 0
    phase[1] <- 0
    mu_1 <- sum(beta_dyn[[1]] * c(1, x[1])[seq_along(beta_dyn[[1]])])
    y[1] <- mu_1 + rnorm(1, 0, sigma[1])

    for (t in 2:n_times) {
      current_phase <- phase[t - 1]

      # Check for transition (if not in absorbing state)
      if (current_phase < k - 1) {
        # Get transition coefficients for current phase
        b_trans <- beta_trans[[current_phase + 1]]
        X_trans <- c(1, x[t])[seq_along(b_trans)]
        p_trans <- plogis(sum(b_trans * X_trans))

        if (runif(1) < p_trans) {
          phase[t] <- current_phase + 1
        } else {
          phase[t] <- current_phase
        }
      } else {
        # Absorbing state
        phase[t] <- current_phase
      }

      # Generate outcome for current phase
      p <- phase[t] + 1  # 1-indexed
      b_dyn <- beta_dyn[[p]]
      X_dyn <- c(1, x[t])[seq_along(b_dyn)]
      y[t] <- sum(b_dyn * X_dyn) + rnorm(1, 0, sigma[p])
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


#' Simulate Poisson data from a phase model
#'
#' @param n_units Number of units
#' @param n_times Number of time points per unit
#' @param beta_trans Transition coefficients
#' @param beta_0 Phase 0 dynamics coefficients (log scale)
#' @param beta_1 Phase 1 dynamics coefficients (log scale)
#' @param seed Random seed
#'
#' @return A data frame suitable for `fit_phaseR()` with `family = poisson()`
#' @export
#'
#' @examples
#' dat <- sim_phaseR_poisson(n_units = 50, n_times = 5)
#' head(dat)
#'
sim_phaseR_poisson <- function(n_units = 100,
                                n_times = 5,
                                beta_trans = c(-1, 0.5),
                                beta_0 = c(2, 0),      # log(lambda) ~ 7.4 at x=0
                                beta_1 = c(1.5, -0.3), # log(lambda) ~ 4.5 at x=0
                                seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Ensure coefficients have x effect (default 0 if not provided)
  if (length(beta_trans) == 1) beta_trans <- c(beta_trans, 0)
  if (length(beta_0) == 1) beta_0 <- c(beta_0, 0)
  if (length(beta_1) == 1) beta_1 <- c(beta_1, 0)

  data_list <- lapply(seq_len(n_units), function(i) {

    x <- rnorm(n_times)
    phase <- integer(n_times)
    y <- integer(n_times)

    # Initial phase always 0
    phase[1] <- 0
    lambda_0 <- exp(beta_0[1] + beta_0[2] * x[1])
    y[1] <- rpois(1, lambda_0)

    for (t in 2:n_times) {
      if (phase[t - 1] == 0) {
        p_trans <- plogis(beta_trans[1] + beta_trans[2] * x[t])
        phase[t] <- rbinom(1, 1, p_trans)
      } else {
        phase[t] <- 1
      }

      if (phase[t] == 0) {
        lambda <- exp(beta_0[1] + beta_0[2] * x[t])
        y[t] <- rpois(1, lambda)
      } else {
        lambda <- exp(beta_1[1] + beta_1[2] * x[t])
        y[t] <- rpois(1, lambda)
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
