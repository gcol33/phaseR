#ifndef PRIOR_FUNCTIONS_H
#define PRIOR_FUNCTIONS_H

#include <Rcpp.h>
#include <cmath>
#include "scalar_math.h"

// Prior type codes (must match prior_to_cpp() in R/prior.R):
// 0 = normal, 1 = student_t, 2 = cauchy, 3 = half_normal, 4 = half_cauchy,
// 5 = lkj.
//
// Every log density returns its value and writes its derivative with
// respect to the argument into d.

namespace prior {

using phaser::kLog2;
using phaser::kLog2Pi;
using phaser::kLogPi;

// Normal(mean, sd).
template <class T>
T log_normal(const T& x, double mean, double sd, T& d) {
  T z = (x - mean) / sd;
  d = -z / sd;
  return -0.5 * kLog2Pi - std::log(sd) - 0.5 * z * z;
}

// Normal(0, exp(log_sd)) density of x; d_x and d_log_sd are its derivatives.
template <class T>
T log_normal_log_sd(const T& x, const T& log_sd, T& d_x, T& d_log_sd) {
  using std::exp;
  T inv_sd = exp(-log_sd);
  T z = x * inv_sd;
  d_x = -z * inv_sd;
  d_log_sd = z * z - 1.0;
  return -0.5 * kLog2Pi - log_sd - 0.5 * z * z;
}

// Student-t(df, mean, scale).
template <class T>
T log_student_t(const T& x, double df, double mean, double scale, T& d) {
  using std::log1p;
  T z = (x - mean) / scale;
  double log_const = std::lgamma((df + 1.0) / 2.0) - std::lgamma(df / 2.0) -
                     0.5 * std::log(df * M_PI) - std::log(scale);
  d = -(df + 1.0) * z / (scale * (df + z * z));
  return log_const - ((df + 1.0) / 2.0) * log1p(z * z / df);
}

// Cauchy(location, scale).
template <class T>
T log_cauchy(const T& x, double location, double scale, T& d) {
  using std::log1p;
  T z = (x - location) / scale;
  d = -2.0 * z / (scale * (1.0 + z * z));
  return -kLogPi - std::log(scale) - log1p(z * z);
}

// Half-Normal(0, sd) prior on sigma = exp(log_sigma), as a density of
// log_sigma: the Half-Normal density of sigma times the Jacobian sigma.
template <class T>
T log_half_normal_log_scale(const T& log_sigma, double sd, T& d) {
  using std::exp;
  T z = exp(log_sigma) / sd;
  d = 1.0 - z * z;
  return kLog2 - 0.5 * kLog2Pi - std::log(sd) - 0.5 * z * z + log_sigma;
}

// Half-Cauchy(0, scale) prior on sigma = exp(log_sigma), as a density of
// log_sigma.
template <class T>
T log_half_cauchy_log_scale(const T& log_sigma, double scale, T& d) {
  using std::exp;
  using std::log1p;
  T z = exp(log_sigma) / scale;
  T z2 = z * z;
  d = 1.0 - 2.0 * z2 / (1.0 + z2);
  return kLog2 - kLogPi - std::log(scale) - log1p(z2) + log_sigma;
}

// Prior of an unbounded parameter (a regression coefficient).
template <class T>
T log_prior_unbounded(const T& x, int type, const Rcpp::NumericVector& p,
                      T& d) {
  switch (type) {
    case 0: return log_normal(x, p[0], p[1], d);
    case 1: return log_student_t(x, p[0], p[1], p[2], d);
    case 2: return log_cauchy(x, p[0], p[1], d);
    default:
      Rcpp::stop("Invalid prior type for unbounded parameter: %d", type);
  }
}

// Prior of a positive parameter sigma, as a density of log(sigma).
template <class T>
T log_prior_positive(const T& log_sigma, int type,
                     const Rcpp::NumericVector& p, T& d) {
  switch (type) {
    case 3: return log_half_normal_log_scale(log_sigma, p[0], d);
    case 4: return log_half_cauchy_log_scale(log_sigma, p[0], d);
    default:
      Rcpp::stop("Invalid prior type for positive parameter: %d", type);
  }
}

// Positive-parameter prior on a real-valued scale whose sign is not
// identified (the diagonal of a Cholesky factor): the half distribution
// reflected to the whole real line, Normal(0, sd) or Cauchy(0, scale).
template <class T>
T log_prior_symmetric_scale(const T& x, int type,
                            const Rcpp::NumericVector& p, T& d) {
  switch (type) {
    case 3: return log_normal(x, 0.0, p[0], d);
    case 4: return log_cauchy(x, 0.0, p[0], d);
    default:
      Rcpp::stop("Invalid prior type for positive parameter: %d", type);
  }
}

// The prior of each parameter class, as built by prior_to_cpp() in R:
// regression coefficients (beta), residual SDs (sigma), random-effect SDs
// (re) and random-effect correlations (cor). SDs are parameterised by their
// logarithm, so sigma() and re_sd() include the Jacobian of exp().
class PriorSet {
public:
  explicit PriorSet(const Rcpp::List& p)
    : beta_type_(Rcpp::as<int>(p["beta_type"])),
      beta_params_(Rcpp::as<Rcpp::NumericVector>(p["beta_params"])),
      sigma_type_(Rcpp::as<int>(p["sigma_type"])),
      sigma_params_(Rcpp::as<Rcpp::NumericVector>(p["sigma_params"])),
      re_type_(Rcpp::as<int>(p["re_type"])),
      re_params_(Rcpp::as<Rcpp::NumericVector>(p["re_params"])),
      cor_type_(Rcpp::as<int>(p["cor_type"])),
      cor_params_(Rcpp::as<Rcpp::NumericVector>(p["cor_params"])) {
    if (cor_type_ != 5) {
      Rcpp::stop("Invalid prior type for correlation: %d", cor_type_);
    }
  }

  template <class T>
  T beta(const T& x, T& d) const {
    return log_prior_unbounded(x, beta_type_, beta_params_, d);
  }

  template <class T>
  T sigma(const T& log_sigma, T& d) const {
    return log_prior_positive(log_sigma, sigma_type_, sigma_params_, d);
  }

  template <class T>
  T re_sd(const T& log_sigma, T& d) const {
    return log_prior_positive(log_sigma, re_type_, re_params_, d);
  }

  template <class T>
  T re_scale(const T& x, T& d) const {
    return log_prior_symmetric_scale(x, re_type_, re_params_, d);
  }

  double cor_eta() const { return cor_params_[0]; }

private:
  int beta_type_;
  Rcpp::NumericVector beta_params_;
  int sigma_type_;
  Rcpp::NumericVector sigma_params_;
  int re_type_;
  Rcpp::NumericVector re_params_;
  int cor_type_;
  Rcpp::NumericVector cor_params_;
};

}  // namespace prior

#endif  // PRIOR_FUNCTIONS_H
