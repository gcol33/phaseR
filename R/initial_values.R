#' Data-informed starting points
#'
#' Starting points of the Laplace mode search and initial values of the NUTS
#' chains. A phase model's posterior can have a mode in which a phase is
#' never entered and its parameters stay at their prior; starting all
#' parameters near zero, far from the data, can end there. Each starting
#' point is built from a provisional phase assignment: every
#' unit follows a path that starts in the first phase and moves at most one
#' phase forward per observation. Given the paths, each phase's dynamics are
#' fitted by a weighted ridge regression (Gaussian, Poisson or logistic) of
#' the outcome on that phase's design matrix, using the observations assigned
#' to the phase (weight 1) and all others (weight `1e-3`, which keeps the fit
#' defined for a phase no unit enters), and each transition by a ridge
#' logistic regression of "moved" on the transition design matrix over the
#' observations at risk of that move. Random effects start at zero, with
#' their SDs (or the diagonal of their Cholesky factor) at the scale of the
#' random-effect prior.
#'
#' The first path splits each unit's observations into `K` consecutive runs
#' of equal length, so every phase is entered; the others draw the number
#' and the times of each unit's moves at random from R's random number
#' generator.
#'
#' @param target Model description from [posterior_target()]
#' @param n_starts Number of starting points
#' @return List of `n_starts` parameter vectors
#' @keywords internal
initial_values <- function(target, n_starts) {
  n_starts <- as.integer(n_starts)
  if (length(n_starts) != 1L || is.na(n_starts) || n_starts < 1L) {
    stop("`n_starts` must be a positive integer.", call. = FALSE)
  }
  paths <- c(list(phase_paths(target, random = FALSE)),
             lapply(seq_len(n_starts - 1L), function(i) {
               phase_paths(target, random = TRUE)
             }))
  lapply(paths, start_from_paths, target = target)
}


#' Provisional phase of every observation
#'
#' @param target Model description from [posterior_target()]
#' @param random If `FALSE`, split each unit into `K` runs of equal length;
#'   if `TRUE`, draw the number of moves uniformly from those possible and
#'   their times uniformly without replacement
#' @return Integer vector of 0-based phases, one per observation
#' @keywords internal
phase_paths <- function(target, random) {
  K <- target$n_phases
  phase <- integer(length(target$y))
  for (i in seq_along(target$unit_start)) {
    rows <- (target$unit_start[i]:target$unit_end[i]) + 1L
    n <- length(rows)
    if (random) {
      n_moves <- sample.int(min(K - 1L, n - 1L) + 1L, 1L) - 1L
      move_at <- (seq_len(n - 1L))[sample.int(n - 1L, n_moves)]
      path <- vapply(seq_len(n) - 1L, function(s) sum(move_at <= s),
                     integer(1))
    } else {
      path <- floor(K * (seq_len(n) - 1L) / n)
      for (s in seq_len(n)[-1L]) path[s] <- min(path[s], path[s - 1L] + 1L)
    }
    phase[rows] <- as.integer(path)
  }
  phase
}


#' Parameter vector fitted to provisional phase paths
#'
#' @param phase 0-based phase of every observation, from [phase_paths()]
#' @param target Model description from [posterior_target()]
#' @return Parameter vector
#' @keywords internal
start_from_paths <- function(phase, target) {
  K <- target$n_phases
  y <- target$y
  theta <- numeric(target$n_params)
  family_names <- c("gaussian", "poisson", "binomial")

  for (k in seq_len(K)) {
    pred <- target$dynamics[[k]]
    w <- ifelse(phase == k - 1L, 1, 1e-3)
    family <- family_names[target$family[k] + 1L]
    beta <- ridge_glm(pred$X, y, w, family)
    theta[pred$offset + seq_along(beta)] <- beta
    log_sigma <- if (family == "gaussian") {
      resid <- y - drop(pred$X %*% beta)
      0.5 * log(sum(w * resid^2) / sum(w))
    } else {
      log(target$prior$sigma_params[1])
    }
    theta[target$sigma_offset + k] <- log_sigma
  }

  unit_first <- rep(FALSE, length(y))
  unit_first[target$unit_start + 1L] <- TRUE
  previous <- c(NA_integer_, phase[-length(phase)])
  for (j in seq_len(K - 1L)) {
    pred <- target$transitions[[j]]
    at_risk <- !unit_first & previous == j - 1L
    moved <- as.numeric(phase[at_risk] == j)
    beta <- ridge_glm(pred$X[at_risk, , drop = FALSE], moved,
                      rep(1, sum(at_risk)), "binomial")
    theta[pred$offset + seq_along(beta)] <- beta
  }

  re_scale <- target$prior$re_params[1]
  for (pred in c(target$transitions, target$dynamics)) {
    for (term in pred$re) {
      hyper <- term$offset + term$n_groups * term$n_coefs
      if (term$kind == 2L) {
        n_c <- term$n_coefs
        diagonal <- cumsum(c(1L, n_c - seq_len(n_c - 1L) + 1L))
        theta[hyper + diagonal] <- re_scale
      } else {
        theta[hyper + seq_len(if (term$kind == 0L) 1L else term$n_coefs)] <-
          log(re_scale)
      }
    }
  }
  theta
}


#' Weighted ridge regression with a canonical link
#'
#' Maximises `sum(w * loglik(y | X beta)) - lambda / 2 * sum(beta^2)` by
#' Newton steps with step halving, for the Gaussian (identity link), Poisson
#' (log link) and binomial (logit link, `y` a proportion of one trial)
#' families.
#'
#' @param X Design matrix
#' @param y Outcome
#' @param w Non-negative observation weights
#' @param family `"gaussian"`, `"poisson"` or `"binomial"`
#' @param lambda Ridge penalty
#' @return Coefficient vector
#' @keywords internal
ridge_glm <- function(X, y, w, family, lambda = 1e-3) {
  p <- ncol(X)
  objective <- function(beta) {
    eta <- drop(X %*% beta)
    ll <- switch(family,
      gaussian = -0.5 * (y - eta)^2,
      poisson = y * eta - exp(eta),
      binomial = y * eta - (pmax(eta, 0) + log1p(exp(-abs(eta))))
    )
    sum(w * ll) - 0.5 * lambda * sum(beta^2)
  }
  beta <- numeric(p)
  value <- objective(beta)
  for (iter in seq_len(100L)) {
    eta <- drop(X %*% beta)
    mu <- switch(family, gaussian = eta, poisson = exp(eta),
                 binomial = stats::plogis(eta))
    v <- switch(family, gaussian = rep(1, length(eta)), poisson = mu,
                binomial = mu * (1 - mu))
    grad <- drop(crossprod(X, w * (y - mu))) - lambda * beta
    H <- crossprod(X, X * (w * v)) + diag(lambda, p)
    step <- drop(solve(H, grad))
    t <- 1
    repeat {
      candidate <- beta + t * step
      new_value <- objective(candidate)
      if (is.finite(new_value) && new_value >= value - 1e-12) break
      t <- t / 2
      if (t < 1e-10) return(beta)
    }
    beta <- candidate
    converged <- abs(new_value - value) < 1e-10 * (1 + abs(value))
    value <- new_value
    if (converged) break
  }
  beta
}


#' Evaluate an expression with a seeded random number generator
#'
#' Seeds R's generator with `seed` for `expr` and restores the caller's
#' generator state afterwards.
#'
#' @param seed Integer seed
#' @param expr Expression to evaluate
#' @return The value of `expr`
#' @keywords internal
with_seed <- function(seed, expr) {
  env <- globalenv()
  had_seed <- exists(".Random.seed", envir = env, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = env, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = env)
    } else if (exists(".Random.seed", envir = env, inherits = FALSE)) {
      rm(".Random.seed", envir = env)
    }
  })
  set.seed(seed)
  expr
}
