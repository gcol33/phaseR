#' Fit a phase model
#'
#' @param model A `phase_model` object
#' @param data A data frame with columns: id, time, outcome(s), covariates
#' @param chains Number of MCMC chains (default 4)
#' @param iter Total iterations per chain (default 2000)
#' @param warmup Warmup iterations per chain (default iter/2)
#' @param backend Estimation backend: "hmc" (default), the No-U-Turn sampler
#'   with step-size and diagonal-metric adaptation during warmup, or
#'   "laplace", a Gaussian approximation at the posterior mode
#' @param cores Number of cores for parallel chains (default 1)
#' @param prior Prior specification from `prior()` (default: weakly
#'   informative); both backends evaluate it for every model family
#' @param seed Integer seed. Two fits with the same `seed`, data and settings
#'   return identical draws. With `seed = NULL` the seed is drawn from R's
#'   random number generator, so `set.seed()` before the call makes the fit
#'   reproducible; the seed used is stored in the returned object as `seed`.
#' @param ... Backend settings. For `backend = "hmc"`: `adapt_delta`, the target
#'   acceptance statistic of step-size adaptation (default 0.8; raise it
#'   towards 1 when there are divergent transitions), and `max_treedepth`
#'   (default 10). For `backend = "laplace"`: `n_samples` (default 1000)
#'   and `n_starts` (default 5), the number of data-informed starting points
#'   from which the posterior mode is searched; the highest mode found is
#'   kept.
#'
#' @return A `phase_fit` object
#' @export
#'
#' @examples
#' \dontrun{
#' model <- phase_model(
#'   phases("inactive", "active"),
#'   transition("inactive", "active", ~ x),
#'   dynamics("inactive", y ~ 1),
#'   dynamics("active", y ~ x)
#' )
#' dat <- sim_phaseR(n_units = 50, n_times = 5)
#' fit <- fit_phaseR(model, data = dat)
#' }
#'
fit_phaseR <- function(model,
                       data,
                       chains = 4,
                       iter = 2000,
                       warmup = floor(iter / 2),
                       backend = c("hmc", "laplace"),
                       cores = 1,
                       prior = NULL,
                       seed = NULL,
                       ...) {

  backend <- match.arg(backend)

  # Set up priors (use defaults if not specified)
  if (is.null(prior)) {
    prior <- prior()
  }

  # Validate
  validate_model_data(model, data)

  # Prepare data for C++
  stan_data <- prepare_stan_data(model, data, prior)

  seed <- resolve_seed(seed)

  # Dispatch to backend
  fit_raw <- switch(backend,
    hmc = fit_hmc(stan_data, chains, iter, warmup, cores, seed, ...),
    laplace = fit_laplace(stan_data, seed = seed, ...),
    stop("Unknown backend")
  )

  fit <- convert_to_phase_fit(fit_raw, model, data, backend, chains, iter, warmup, stan_data)
  fit$seed <- seed
  fit
}


#' Resolve the seed of a fit
#'
#' @param seed `NULL` or a single integer-valued number.
#' @return An integer seed. `NULL` is replaced by a draw from R's random
#'   number generator.
#' @keywords internal
resolve_seed <- function(seed) {
  if (is.null(seed)) {
    return(sample.int(.Machine$integer.max, 1L))
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != round(seed) || abs(seed) > .Machine$integer.max) {
    stop("`seed` must be NULL or a single integer.", call. = FALSE)
  }
  as.integer(seed)
}


#' Prepare data for C++ likelihood
#'
#' Builds the design of a phase model with `K` phases: one transition
#' linear predictor per consecutive pair of phases (from the `transition()`
#' between them, intercept-only if there is none) and one dynamics linear
#' predictor per phase (from its `dynamics()`, an intercept-only Gaussian
#' model of the response if there is none). Each may contain random-effect
#' terms, and each phase has its own outcome family.
#'
#' The parameter vector holds the transition coefficients, the dynamics
#' coefficients and the log residual SD of each phase, followed by the
#' random-effect parameters of each predictor in the same order.
#'
#' @param model A `phase_model`
#' @param data Data frame
#' @param prior Prior specification from [prior()]
#' @return List with the outcome, the unit ranges, the design matrices
#'   (`X_trans_list`, `X_dyn_list`), the random-effect structures
#'   (`re_trans`, `re_dyn`, from [build_multi_re()]), the family code of each
#'   phase, the parameter count and names, and the prior in C++ form
#' @keywords internal
prepare_stan_data <- function(model, data, prior = NULL) {

  if (is.null(prior)) prior <- prior()
  prior_cpp <- prior_to_cpp(prior)

  # Sort data by id, then time
  data <- data[order(data$id, data$time), ]

  unit_ids <- unique(data$id)
  n_units <- length(unit_ids)
  unit_start <- integer(n_units)
  unit_end <- integer(n_units)
  for (i in seq_along(unit_ids)) {
    idx <- which(data$id == unit_ids[i])
    unit_start[i] <- min(idx) - 1L  # 0-indexed for C++
    unit_end[i] <- max(idx) - 1L
  }

  K <- model$phases$n_phases
  phase_names <- model$phases$names
  response_var <- model$dynamics[[1]]$response

  for (tr in model$transitions) {
    if (match(tr$to, phase_names) != match(tr$from, phase_names) + 1L) {
      stop(sprintf(paste0("Transition %s -> %s skips a phase; units move ",
                          "through the phases in order, so transitions ",
                          "link consecutive phases."), tr$from, tr$to),
           call. = FALSE)
    }
  }
  trans_formula <- lapply(seq_len(K - 1L), function(j) {
    for (tr in model$transitions) {
      if (tr$from == phase_names[j]) return(tr$formula)
    }
    ~ 1
  })
  dyn_spec <- lapply(phase_names, function(p) {
    for (d in model$dynamics) {
      if (d$phase == p) return(list(formula = d$formula, family = d$family))
    }
    list(formula = stats::as.formula(paste(response_var, "~ 1")),
         family = stats::gaussian())
  })

  trans_parsed <- lapply(trans_formula, parse_formula_re)
  dyn_parsed <- lapply(dyn_spec, function(d) parse_formula_re(d$formula))

  interactions <- unique(unlist(lapply(c(trans_parsed, dyn_parsed),
                                       function(p) p$nested_interactions)))
  if (length(interactions) > 0) {
    data <- create_nested_interactions(data, interactions)
  }

  X_trans_list <- lapply(trans_parsed, function(p) {
    model.matrix(p$fixed_formula, data)
  })
  X_dyn_list <- lapply(dyn_parsed, function(p) {
    model.matrix(update(p$fixed_formula, NULL ~ .), data)
  })
  re_trans <- lapply(trans_parsed, function(p) build_multi_re(data, p$re_terms))
  re_dyn <- lapply(dyn_parsed, function(p) build_multi_re(data, p$re_terms))
  family <- vapply(dyn_spec, function(d) encode_family(d$family), integer(1))

  trans_label <- if (K == 2L) {
    "trans"
  } else {
    paste0("trans_", phase_names[-K], "_to_", phase_names[-1])
  }
  re_names <- function(re, label, transition) {
    unlist(lapply(names(re$re_list), function(re_name) {
      build_re_param_names(label, re_name, re$re_list[[re_name]], transition)
    }))
  }
  param_names <- c(
    unlist(lapply(seq_len(K - 1L), function(j) {
      paste0("beta_", trans_label[j], "_", colnames(X_trans_list[[j]]))
    })),
    unlist(lapply(seq_len(K), function(k) {
      paste0("beta_", phase_names[k], "_", colnames(X_dyn_list[[k]]))
    })),
    paste0("log_sigma_", phase_names),
    unlist(lapply(seq_len(K - 1L), function(j) {
      re_names(re_trans[[j]], trans_label[j], TRUE)
    })),
    unlist(lapply(seq_len(K), function(k) {
      re_names(re_dyn[[k]], phase_names[k], FALSE)
    }))
  )
  n_params <- sum(vapply(X_trans_list, ncol, integer(1))) +
    sum(vapply(X_dyn_list, ncol, integer(1))) + K +
    sum(vapply(c(re_trans, re_dyn), function(r) r$total_re_params, numeric(1)))
  stopifnot(length(param_names) == n_params)

  has_trans_re <- any(vapply(re_trans, function(r) r$n_re_terms > 0, logical(1)))
  has_dyn_re <- any(vapply(re_dyn, function(r) r$n_re_terms > 0, logical(1)))
  first_levels <- function(r) {
    if (r$n_re_terms > 0) r$re_list[[1]]$levels else character(0)
  }

  list(
    id = as.integer(factor(data$id)),
    time = data$time,
    y = data[[response_var]],
    n_units = n_units,
    unit_start = unit_start,
    unit_end = unit_end,
    k_phases = K,
    phase_names = phase_names,
    X_trans_list = X_trans_list,
    X_dyn_list = X_dyn_list,
    n_trans_coef = vapply(X_trans_list, ncol, integer(1)),
    n_dyn_coef = vapply(X_dyn_list, ncol, integer(1)),
    family = family,
    is_glm = any(family > 0L),
    re_trans = re_trans,
    re_dyn = re_dyn,
    has_re = has_trans_re || has_dyn_re,
    has_trans_re = has_trans_re,
    has_dyn_re = has_dyn_re,
    has_slopes = any(vapply(c(re_trans, re_dyn), function(r) r$has_slopes,
                            logical(1))),
    trans_re_levels = first_levels(re_trans[[1]]),
    dyn_re_levels = stats::setNames(lapply(re_dyn, first_levels), phase_names),
    n_params = as.integer(n_params),
    param_names = param_names,
    prior_cpp = prior_cpp
  )
}


#' Create interaction variables for nested RE
#'
#' Creates columns like "site:plot" from "site" and "plot" columns.
#'
#' @param data Data frame
#' @param interactions Character vector of interaction terms (e.g., "site:plot")
#' @return Modified data frame with new columns
#' @keywords internal
#'
create_nested_interactions <- function(data, interactions) {

  for (int_term in interactions) {
    if (int_term %in% names(data)) {
      next  # Already exists
    }

    # Split the interaction term
    parts <- strsplit(int_term, ":")[[1]]

    # Validate all parts exist
    missing <- parts[!parts %in% names(data)]
    if (length(missing) > 0) {
      stop(sprintf("Grouping variable(s) '%s' not found in data",
                   paste(missing, collapse = "', '")), call. = FALSE)
    }

    # Create the interaction variable using base::interaction
    int_args <- lapply(parts, function(p) data[[p]])
    int_args$drop <- TRUE
    int_args$sep <- ":"
    interaction_values <- do.call(base::interaction, int_args)
    data[[int_term]] <- interaction_values
  }

  data
}


#' Validate nested structure
#'
#' Checks that nested groups are properly nested (each level of b is within one level of a).
#'
#' @param data Data frame
#' @param outer Outer grouping variable name
#' @param inner Inner grouping variable name
#' @return TRUE if valid, otherwise throws an error
#' @keywords internal
#'
validate_nesting <- function(data, outer, inner) {

  # For each level of inner, there should be exactly one level of outer
  mapping <- unique(data[, c(outer, inner), drop = FALSE])
  inner_counts <- table(mapping[[inner]])

  if (any(inner_counts > 1)) {
    bad_levels <- names(inner_counts)[inner_counts > 1]
    stop(sprintf("Invalid nesting: '%s' levels appear in multiple '%s' groups: %s",
                 inner, outer, paste(bad_levels[1:min(3, length(bad_levels))], collapse = ", ")),
         call. = FALSE)
  }

  TRUE
}


#' Build RE parameter names
#'
#' Generates parameter names for RE terms including slopes
#'
#' @param component Label of the linear predictor (e.g. "trans", a phase
#'   name)
#' @param re_name Name of grouping variable
#' @param re_info RE info structure from build_multi_re
#' @param transition Whether the term belongs to a transition predictor
#'   (effects named `u_...`) or a dynamics predictor (`v_...`)
#' @return Character vector of parameter names
#' @keywords internal
#'
build_re_param_names <- function(component, re_name, re_info,
                                 transition = component == "trans") {

  n_coefs <- re_info$n_coefs
  levels <- re_info$levels

  coef_names <- c(if (re_info$has_intercept) "intercept", re_info$slope_vars)
  prefix <- paste0(if (transition) "u_" else "v_", component, "_", re_name)
  hyper <- paste0(component, "_re_", re_name)

  # Effects: one per group (and coefficient), level-major
  names <- if (n_coefs == 1) {
    paste0(prefix, "_", levels)
  } else {
    paste0(prefix, "_", rep(levels, each = n_coefs), "_",
           rep(coef_names, times = length(levels)))
  }

  # Covariance parameters
  if (n_coefs == 1) {
    names <- c(names, paste0("log_sigma_", hyper))
  } else if (re_info$correlated) {
    # SDs, then the unconstrained values of the correlation matrix's
    # Cholesky factor: atanh of the canonical partial correlation of each
    # pair (i, j), i > j, row by row
    pairs <- do.call(rbind, lapply(2:n_coefs, function(i) cbind(i, seq_len(i - 1))))
    names <- c(names, paste0("log_sigma_", hyper, "_", coef_names),
               paste0("atanh_cpc_", hyper, "_", pairs[, 1], "_", pairs[, 2]))
  } else {
    names <- c(names, paste0("log_sigma_", hyper, "_", coef_names))
  }

  names
}


#' Encode R family object to integer
#' @keywords internal
encode_family <- function(family) {
  if (is.null(family)) return(0L)
  fam <- family$family
  switch(fam,
    gaussian = 0L,
    poisson = 1L,
    binomial = 2L,
    0L  # default to gaussian
  )
}


#' Convert raw fit to phase_fit object
#' @keywords internal
convert_to_phase_fit <- function(fit_raw, model, data, backend, chains, iter, warmup, stan_data) {

  k_phases <- stan_data$k_phases

  fit_obj <- list(
    draws = fit_raw$draws,
    chain_id = fit_raw$chain_id,
    param_names = fit_raw$param_names,
    model = model,
    data = data,
    backend = backend,
    n_chains = chains,
    n_iter = iter,
    n_warmup = warmup,
    n_units = length(unique(data$id)),
    n_times = max(data$time),
    diagnostics = fit_raw$diagnostics,
    k_phases = k_phases
  )

  fit_obj$has_re <- stan_data$has_re
  fit_obj$has_trans_re <- stan_data$has_trans_re
  fit_obj$has_dyn_re <- stan_data$has_dyn_re
  fit_obj$trans_re_levels <- stan_data$trans_re_levels
  fit_obj$dyn_re_levels <- stan_data$dyn_re_levels
  fit_obj$phase_names <- stan_data$phase_names

  structure(fit_obj, class = "phase_fit")
}
