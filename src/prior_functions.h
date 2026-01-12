#ifndef PRIOR_FUNCTIONS_H
#define PRIOR_FUNCTIONS_H

#include <Rcpp.h>
#include <cmath>

// Prior type codes (must match R/prior.R)
// 0 = normal, 1 = student_t, 2 = cauchy, 3 = half_normal, 4 = half_cauchy, 5 = lkj

namespace prior {

// =============================================================================
// Constants
// =============================================================================

const double LOG_2PI = std::log(2.0 * M_PI);
const double LOG_PI = std::log(M_PI);


// =============================================================================
// Log-density functions
// =============================================================================

// Normal(mean, sd): log p(x) = -0.5*log(2*pi) - log(sd) - 0.5*((x-mean)/sd)^2
inline double log_normal(double x, double mean, double sd) {
  double z = (x - mean) / sd;
  return -0.5 * LOG_2PI - std::log(sd) - 0.5 * z * z;
}

// Student-t(df, mean, scale): uses standard t-distribution formula
// log p(x) = log(Gamma((df+1)/2)) - log(Gamma(df/2)) - 0.5*log(df*pi)
//            - log(scale) - ((df+1)/2)*log(1 + ((x-mean)/scale)^2/df)
inline double log_student_t(double x, double df, double mean, double scale) {
  double z = (x - mean) / scale;
  // Use lgamma for log-gamma function
  double log_const = std::lgamma((df + 1.0) / 2.0) - std::lgamma(df / 2.0)
                     - 0.5 * std::log(df * M_PI) - std::log(scale);
  return log_const - ((df + 1.0) / 2.0) * std::log(1.0 + z * z / df);
}

// Cauchy(location, scale): Student-t with df=1
// log p(x) = -log(pi) - log(scale) - log(1 + ((x-location)/scale)^2)
inline double log_cauchy(double x, double location, double scale) {
  double z = (x - location) / scale;
  return -LOG_PI - std::log(scale) - std::log(1.0 + z * z);
}

// Half-Normal(0, sd): Normal truncated to x > 0
// log p(x) = log(2) - 0.5*log(2*pi) - log(sd) - 0.5*(x/sd)^2, for x > 0
// For log(sigma) parameterization: x = log(sigma), so sigma = exp(x)
// We evaluate the half-normal on sigma, then add Jacobian log(sigma)
inline double log_half_normal(double log_sigma, double sd) {
  double sigma = std::exp(log_sigma);
  if (sigma <= 0) return -INFINITY;
  double z = sigma / sd;
  // Half-normal on sigma + Jacobian
  return std::log(2.0) - 0.5 * LOG_2PI - std::log(sd) - 0.5 * z * z + log_sigma;
}

// Half-Cauchy(0, scale): Cauchy truncated to x > 0
// log p(x) = log(2) - log(pi) - log(scale) - log(1 + (x/scale)^2), for x > 0
// For log(sigma) parameterization: same Jacobian adjustment
inline double log_half_cauchy(double log_sigma, double scale) {
  double sigma = std::exp(log_sigma);
  if (sigma <= 0) return -INFINITY;
  double z = sigma / scale;
  // Half-Cauchy on sigma + Jacobian
  return std::log(2.0) - LOG_PI - std::log(scale) - std::log(1.0 + z * z) + log_sigma;
}

// LKJ(eta) prior on correlation matrix
// For a single correlation rho (2x2 case):
// log p(rho) = (eta - 1) * log(1 - rho^2) + const
// For Cholesky factor L where rho = L[1,0], with constraint L[0,0]=1, L[1,1]=sqrt(1-L[1,0]^2)
// We store L[1,0] directly, so the Jacobian is already handled by Cholesky parameterization
inline double log_lkj_chol_2x2(double L_10, double eta) {
  // L_10 is the off-diagonal of lower Cholesky factor
  // rho = L_10 (since L[0,0] = 1 for correlation matrix Cholesky)
  double rho_sq = L_10 * L_10;
  if (rho_sq >= 1.0) return -INFINITY;
  // LKJ prior: (1 - rho^2)^(eta-1)
  return (eta - 1.0) * std::log(1.0 - rho_sq);
}

// General LKJ for larger correlation matrices
// For d-dimensional correlation matrix, the log-prior on Cholesky factor is:
// sum over k=2..d of (d - k + 2*eta - 2) * log(L[k,k])
// where L[k,k] = sqrt(1 - sum(L[k,j]^2 for j<k))
inline double log_lkj_chol(const Rcpp::NumericVector& L_lower, int d, double eta) {
  // L_lower contains the strictly lower triangular elements column by column
  // For d=2: L_lower = {L[1,0]} (length 1)
  // For d=3: L_lower = {L[1,0], L[2,0], L[2,1]} (length 3)

  if (d == 1) return 0.0;  // No correlation for 1D

  double log_prior = 0.0;
  int idx = 0;

  for (int k = 1; k < d; k++) {
    // Sum of squares of elements in row k (below diagonal)
    double sum_sq = 0.0;
    for (int j = 0; j < k; j++) {
      double L_kj = L_lower[idx++];
      sum_sq += L_kj * L_kj;
    }

    if (sum_sq >= 1.0) return -INFINITY;

    // L[k,k] = sqrt(1 - sum_sq)
    double L_kk = std::sqrt(1.0 - sum_sq);

    // LKJ contribution for this row
    log_prior += (d - k + 2.0 * eta - 2.0) * std::log(L_kk);
  }

  return log_prior;
}


// =============================================================================
// Dispatcher function
// =============================================================================

// Evaluate log-prior for an unbounded parameter (beta, etc.)
inline double log_prior_unbounded(double x, int type, const Rcpp::NumericVector& params) {
  switch (type) {
    case 0:  // normal(mean, sd)
      return log_normal(x, params[0], params[1]);
    case 1:  // student_t(df, mean, scale)
      return log_student_t(x, params[0], params[1], params[2]);
    case 2:  // cauchy(location, scale)
      return log_cauchy(x, params[0], params[1]);
    default:
      Rcpp::stop("Invalid prior type for unbounded parameter: %d", type);
      return -INFINITY;
  }
}

// Evaluate log-prior for a positive parameter (sigma, stored as log_sigma)
inline double log_prior_positive(double log_sigma, int type, const Rcpp::NumericVector& params) {
  switch (type) {
    case 3:  // half_normal(sd)
      return log_half_normal(log_sigma, params[0]);
    case 4:  // half_cauchy(scale)
      return log_half_cauchy(log_sigma, params[0]);
    default:
      Rcpp::stop("Invalid prior type for positive parameter: %d", type);
      return -INFINITY;
  }
}

// Evaluate LKJ log-prior for correlation Cholesky factor
inline double log_prior_cor(const Rcpp::NumericVector& L_lower, int d, int type,
                            const Rcpp::NumericVector& params) {
  if (type != 5) {
    Rcpp::stop("Invalid prior type for correlation: %d", type);
    return -INFINITY;
  }
  return log_lkj_chol(L_lower, d, params[0]);  // params[0] = eta
}

}  // namespace prior

#endif  // PRIOR_FUNCTIONS_H
