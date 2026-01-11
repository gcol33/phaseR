#' Specify priors for phase model
#'
#' Creates a prior specification object. Currently supports Normal priors
#' on all coefficients.
#'
#' @param beta_sd Prior SD for regression coefficients (default 2.5)
#' @param sigma_sd Prior SD for log(sigma) parameters (default 1)
#' @param re_sd Prior SD for random effect hyperparameter (default 0.5)
#'
#' @return A `phase_prior` object
#' @export
#'
#' @details
#' The default priors are weakly informative:
#' \itemize{
#'   \item Regression coefficients: Normal(0, 2.5)
#'   \item Log residual SD: Normal(0, 1)
#'   \item Random effect SD: Half-Normal(0, 0.5)
#' }
#'
#' @examples
#' # Default priors
#' prior()
#'
#' # Tighter priors on coefficients
#' prior(beta_sd = 1)
#'
prior <- function(beta_sd = 2.5, sigma_sd = 1, re_sd = 0.5) {

  if (beta_sd <= 0 || sigma_sd <= 0 || re_sd <= 0) {
    stop("Prior SDs must be positive", call. = FALSE)
  }

  structure(
    list(
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
  cat(sprintf("  Coefficients: Normal(0, %.1f)\n", x$beta_sd))
  cat(sprintf("  Log(sigma):   Normal(0, %.1f)\n", x$sigma_sd))
  cat(sprintf("  RE sigma:     Half-Normal(0, %.1f)\n", x$re_sd))
  invisible(x)
}


#' Check if priors are default
#' @keywords internal
is_default_prior <- function(prior) {
  if (is.null(prior)) return(TRUE)
  prior$beta_sd == 2.5 && prior$sigma_sd == 1 && prior$re_sd == 0.5
}
