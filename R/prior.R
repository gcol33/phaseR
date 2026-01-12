#' Specify priors for phase model
#'
#' Creates a prior specification object with either simple SD parameters
#' or full distribution specifications using formula syntax.
#'
#' @param ... Prior specifications using formula syntax (e.g., `beta ~ normal(0, 2.5)`)
#'   or named SD parameters for backward compatibility (`beta_sd = 2.5`).
#' @param beta_sd Prior SD for regression coefficients (default 2.5, ignored if
#'   formula syntax is used)
#' @param sigma_sd Prior SD for log(sigma) parameters (default 1, ignored if
#'   formula syntax is used)
#' @param re_sd Prior SD for random effect hyperparameter (default 0.5, ignored if
#'   formula syntax is used)
#'
#' @return A `phase_prior` object
#' @export
#'
#' @details
#' The default priors are weakly informative:
#' \itemize{
#'   \item Regression coefficients: Normal(0, 2.5)
#'   \item Residual SD (sigma): Half-Normal(0, 1)
#'   \item Random effect SD: Half-Normal(0, 0.5)
#'   \item Correlation (LKJ): LKJ(2)
#' }
#'
#' ## Supported distributions
#'
#' \describe{
#'   \item{`normal(mean, sd)`}{Normal distribution}
#'   \item{`student_t(df, mean, scale)`}{Student-t distribution}
#'   \item{`cauchy(location, scale)`}{Cauchy distribution}
#'   \item{`half_normal(sd)`}{Half-Normal (positive values only)}
#'   \item{`half_cauchy(scale)`}{Half-Cauchy (positive values only)}
#'   \item{`lkj(eta)`}{LKJ prior for correlation matrices}
#' }
#'
#' ## Parameter classes
#'
#' \describe{
#'   \item{`beta`}{All regression coefficients (fixed effects)}
#'   \item{`sigma`}{Residual standard deviations}
#'   \item{`re`}{Random effect standard deviations}
#'   \item{`cor`}{Correlation matrices (for correlated random slopes)}
#' }
#'
#' @examples
#' # Default priors
#' prior()
#'
#' # Simple SD specification (backward compatible)
#' prior(beta_sd = 1)
#'
#' # Full distribution specification
#' prior(
#'   beta ~ normal(0, 2.5),
#'   sigma ~ half_cauchy(1),
#'   re ~ half_normal(0.5),
#'   cor ~ lkj(2)
#' )
#'
#' # Student-t for robust regression coefficients
#' prior(beta ~ student_t(3, 0, 2.5))
#'
prior <- function(..., beta_sd = 2.5, sigma_sd = 1, re_sd = 0.5) {

  dots <- list(...)

  # Check if using formula syntax or old SD syntax
  if (length(dots) > 0) {
    # New formula syntax
    specs <- parse_prior_specs(dots)
    return(build_prior_from_specs(specs))
  }


  # Old SD syntax (backward compatible)
  if (beta_sd <= 0 || sigma_sd <= 0 || re_sd <= 0) {
    stop("Prior SDs must be positive", call. = FALSE)
  }

  structure(
    list(
      beta = prior_normal(0, beta_sd),
      sigma = prior_half_normal(sigma_sd),
      re = prior_half_normal(re_sd),
      cor = prior_lkj(2),
      # Keep old fields for backward compatibility
      beta_sd = beta_sd,
      sigma_sd = sigma_sd,
      re_sd = re_sd
    ),
    class = "phase_prior"
  )
}


#' @export
print.phase_prior <- function(x, ...) {
  cat("phaseR prior specification\n")

  # Print each component
  if (!is.null(x$beta)) {
    cat("  Coefficients (beta): ", format_prior(x$beta), "\n", sep = "")
  }
  if (!is.null(x$sigma)) {
    cat("  Residual SD (sigma): ", format_prior(x$sigma), "\n", sep = "")
  }
  if (!is.null(x$re)) {
    cat("  RE SD:               ", format_prior(x$re), "\n", sep = "")
  }
  if (!is.null(x$cor)) {
    cat("  Correlation:         ", format_prior(x$cor), "\n", sep = "")
  }

  invisible(x)
}


#' Check if priors are default
#' @keywords internal
is_default_prior <- function(prior) {
  if (is.null(prior)) return(TRUE)

  # Check old-style fields
  if (!is.null(prior$beta_sd)) {
    return(prior$beta_sd == 2.5 && prior$sigma_sd == 1 && prior$re_sd == 0.5)
  }

  # Check new-style specs
  is_default_normal <- function(p, mean, sd) {
    p$type == "normal" && p$mean == mean && p$sd == sd
  }
  is_default_half_normal <- function(p, sd) {
    p$type == "half_normal" && p$sd == sd
  }
  is_default_lkj <- function(p, eta) {
    p$type == "lkj" && p$eta == eta
  }

  is_default_normal(prior$beta, 0, 2.5) &&
    is_default_half_normal(prior$sigma, 1) &&
    is_default_half_normal(prior$re, 0.5) &&
    is_default_lkj(prior$cor, 2)
}


# =============================================================================
# Distribution Constructors
# =============================================================================

#' Normal prior distribution
#'
#' @param mean Mean of the distribution
#' @param sd Standard deviation (must be positive)
#' @return A prior specification object
#' @export
#'
#' @examples
#' prior(beta ~ normal(0, 2.5))
#'
normal <- function(mean = 0, sd = 1) {
  prior_normal(mean, sd)
}

#' @keywords internal
prior_normal <- function(mean = 0, sd = 1) {
  if (sd <= 0) stop("sd must be positive", call. = FALSE)
  structure(
    list(type = "normal", mean = mean, sd = sd),
    class = "prior_spec"
  )
}


#' Student-t prior distribution
#'
#' @param df Degrees of freedom (must be positive)
#' @param mean Location parameter
#' @param scale Scale parameter (must be positive)
#' @return A prior specification object
#' @export
#'
#' @examples
#' # Robust prior for coefficients
#' prior(beta ~ student_t(3, 0, 2.5))
#'
student_t <- function(df = 3, mean = 0, scale = 2.5) {
  prior_student_t(df, mean, scale)
}

#' @keywords internal
prior_student_t <- function(df = 3, mean = 0, scale = 2.5) {
  if (df <= 0) stop("df must be positive", call. = FALSE)
  if (scale <= 0) stop("scale must be positive", call. = FALSE)
  structure(
    list(type = "student_t", df = df, mean = mean, scale = scale),
    class = "prior_spec"
  )
}


#' Cauchy prior distribution
#'
#' @param location Location parameter
#' @param scale Scale parameter (must be positive)
#' @return A prior specification object
#' @export
#'
#' @examples
#' prior(beta ~ cauchy(0, 2.5))
#'
cauchy <- function(location = 0, scale = 1) {
  prior_cauchy(location, scale)
}

#' @keywords internal
prior_cauchy <- function(location = 0, scale = 1) {

  if (scale <= 0) stop("scale must be positive", call. = FALSE)
  structure(
    list(type = "cauchy", location = location, scale = scale),
    class = "prior_spec"
  )
}


#' Half-Normal prior distribution
#'
#' A Normal(0, sd) distribution truncated to positive values.
#' Suitable for scale parameters like sigma.
#'
#' @param sd Standard deviation of the underlying normal (must be positive)
#' @return A prior specification object
#' @export
#'
#' @examples
#' prior(sigma ~ half_normal(1))
#'
half_normal <- function(sd = 1) {
  prior_half_normal(sd)
}

#' @keywords internal
prior_half_normal <- function(sd = 1) {
  if (sd <= 0) stop("sd must be positive", call. = FALSE)
  structure(
    list(type = "half_normal", sd = sd),
    class = "prior_spec"
  )
}


#' Half-Cauchy prior distribution
#'
#' A Cauchy(0, scale) distribution truncated to positive values.
#' Has heavier tails than half-normal.
#'
#' @param scale Scale parameter (must be positive)
#' @return A prior specification object
#' @export
#'
#' @examples
#' prior(sigma ~ half_cauchy(1))
#'
half_cauchy <- function(scale = 1) {
  prior_half_cauchy(scale)
}

#' @keywords internal
prior_half_cauchy <- function(scale = 1) {
  if (scale <= 0) stop("scale must be positive", call. = FALSE)
  structure(
    list(type = "half_cauchy", scale = scale),
    class = "prior_spec"
  )
}


#' LKJ prior for correlation matrices
#'
#' The LKJ(eta) prior places more mass on the identity matrix
#' as eta increases. eta = 1 is uniform over correlation matrices.
#'
#' @param eta Shape parameter (must be positive). eta = 1 is uniform,
#'   larger values favor weaker correlations.
#' @return A prior specification object
#' @export
#'
#' @examples
#' prior(cor ~ lkj(2))
#'
lkj <- function(eta = 2) {
  prior_lkj(eta)
}

#' @keywords internal
prior_lkj <- function(eta = 2) {
  if (eta <= 0) stop("eta must be positive", call. = FALSE)
  structure(
    list(type = "lkj", eta = eta),
    class = "prior_spec"
  )
}


# =============================================================================
# Prior Specification Parsing
# =============================================================================

#' Parse prior specifications from formula syntax
#' @keywords internal
parse_prior_specs <- function(dots) {

  specs <- list()

  for (i in seq_along(dots)) {
    spec <- dots[[i]]

    # Must be a formula
    if (!inherits(spec, "formula")) {
      stop("Prior specifications must use formula syntax: param ~ distribution(...)",
           call. = FALSE)
    }

    # Extract LHS (parameter class) and RHS (distribution)
    lhs <- as.character(spec[[2]])
    rhs <- spec[[3]]

    # Validate parameter class
    valid_classes <- c("beta", "sigma", "re", "cor")
    if (!lhs %in% valid_classes) {
      stop("Unknown parameter class '", lhs, "'. Valid classes: ",
           paste(valid_classes, collapse = ", "), call. = FALSE)
    }

    # Evaluate RHS to get prior_spec object
    prior_obj <- tryCatch(
      eval(rhs, envir = parent.frame(3)),
      error = function(e) {
        stop("Error evaluating prior for '", lhs, "': ", e$message,
             call. = FALSE)
      }
    )

    if (!inherits(prior_obj, "prior_spec")) {
      stop("Right-hand side of '", lhs, "' must be a distribution function ",
           "(e.g., normal(), half_cauchy())", call. = FALSE)
    }

    # Validate distribution is appropriate for parameter class
    validate_prior_for_class(lhs, prior_obj)

    specs[[lhs]] <- prior_obj
  }

  specs
}


#' Validate that a prior distribution is appropriate for a parameter class
#' @keywords internal
validate_prior_for_class <- function(param_class, prior_obj) {

  type <- prior_obj$type

  # sigma and re must use positive-constrained distributions
  if (param_class %in% c("sigma", "re")) {
    positive_types <- c("half_normal", "half_cauchy")
    if (!type %in% positive_types) {
      stop("Parameter class '", param_class, "' requires a positive distribution ",
           "(half_normal or half_cauchy), not '", type, "'", call. = FALSE)
    }
  }

  # cor must use LKJ
  if (param_class == "cor") {
    if (type != "lkj") {
      stop("Parameter class 'cor' requires lkj() prior, not '", type, "'",
           call. = FALSE)
    }
  }

  # beta can use normal, student_t, or cauchy
  if (param_class == "beta") {
    unbounded_types <- c("normal", "student_t", "cauchy")
    if (!type %in% unbounded_types) {
      stop("Parameter class 'beta' requires an unbounded distribution ",
           "(normal, student_t, or cauchy), not '", type, "'", call. = FALSE)
    }
  }

  invisible(TRUE)
}


#' Build phase_prior object from parsed specifications
#' @keywords internal
build_prior_from_specs <- function(specs) {

  # Start with defaults
  result <- list(
    beta = prior_normal(0, 2.5),
    sigma = prior_half_normal(1),
    re = prior_half_normal(0.5),
    cor = prior_lkj(2)
  )

  # Override with user specifications
  for (name in names(specs)) {
    result[[name]] <- specs[[name]]
  }

  # Add backward-compatible fields
  if (result$beta$type == "normal") {
    result$beta_sd <- result$beta$sd
  } else {
    result$beta_sd <- NA
  }

  if (result$sigma$type == "half_normal") {
    result$sigma_sd <- result$sigma$sd
  } else {
    result$sigma_sd <- NA
  }

  if (result$re$type == "half_normal") {
    result$re_sd <- result$re$sd
  } else {
    result$re_sd <- NA
  }

  structure(result, class = "phase_prior")
}


#' Format prior for printing
#' @keywords internal
format_prior <- function(prior_obj) {
  switch(prior_obj$type,
    "normal" = sprintf("Normal(%.2g, %.2g)", prior_obj$mean, prior_obj$sd),
    "student_t" = sprintf("Student-t(%.1f, %.2g, %.2g)",
                          prior_obj$df, prior_obj$mean, prior_obj$scale),
    "cauchy" = sprintf("Cauchy(%.2g, %.2g)", prior_obj$location, prior_obj$scale),
    "half_normal" = sprintf("Half-Normal(%.2g)", prior_obj$sd),
    "half_cauchy" = sprintf("Half-Cauchy(%.2g)", prior_obj$scale),
    "lkj" = sprintf("LKJ(%.1f)", prior_obj$eta),
    "unknown"
  )
}


# =============================================================================
# Convert priors to C++ format
# =============================================================================

#' Convert prior specification to C++ format
#'
#' Creates a list that can be passed to C++ with type codes and parameters.
#'
#' @param prior A phase_prior object
#' @return List with type codes and parameter vectors for each parameter class
#' @keywords internal
#'
prior_to_cpp <- function(prior) {

  if (is.null(prior)) {
    prior <- prior()  # Default priors
  }

  # Type codes for C++
  # 0 = normal, 1 = student_t, 2 = cauchy, 3 = half_normal, 4 = half_cauchy, 5 = lkj
  type_code <- function(type) {
    switch(type,
      "normal" = 0L,
      "student_t" = 1L,
      "cauchy" = 2L,
      "half_normal" = 3L,
      "half_cauchy" = 4L,
      "lkj" = 5L,
      stop("Unknown prior type: ", type)
    )
  }

  # Extract parameters as vector (padded to 3 elements for consistent C++ interface)
  params_vec <- function(prior_obj) {
    params <- switch(prior_obj$type,
      "normal" = c(prior_obj$mean, prior_obj$sd),
      "student_t" = c(prior_obj$df, prior_obj$mean, prior_obj$scale),
      "cauchy" = c(prior_obj$location, prior_obj$scale),
      "half_normal" = c(prior_obj$sd),
      "half_cauchy" = c(prior_obj$scale),
      "lkj" = c(prior_obj$eta)
    )
    # Pad to length 3
    c(params, rep(0, 3 - length(params)))
  }

  list(
    beta_type = type_code(prior$beta$type),
    beta_params = params_vec(prior$beta),
    sigma_type = type_code(prior$sigma$type),
    sigma_params = params_vec(prior$sigma),
    re_type = type_code(prior$re$type),
    re_params = params_vec(prior$re),
    cor_type = type_code(prior$cor$type),
    cor_params = params_vec(prior$cor)
  )
}
