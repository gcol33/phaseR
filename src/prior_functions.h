#ifndef PRIOR_FUNCTIONS_H
#define PRIOR_FUNCTIONS_H

#include <Rcpp.h>
#include <cmath>
#include "scalar_math.h"

// Prior type codes (must match prior_to_cpp() in R/prior.R):
// 0 = normal, 1 = student_t, 2 = cauchy, 3 = half_normal, 4 = half_cauchy,
// 5 = lkj.

namespace prior {

using phaser::kLog2;
using phaser::kLog2Pi;
using phaser::kLogPi;

// Normal(mean, sd) log density.
template <class T>
T log_normal(const T& x, double mean, double sd) {
  using std::log;
  T z = (x - mean) / sd;
  return -0.5 * kLog2Pi - log(sd) - 0.5 * z * z;
}

// Normal(0, exp(log_sd)) log density of x.
template <class T>
T log_normal_log_sd(const T& x, const T& log_sd) {
  using std::exp;
  T z = x * exp(-log_sd);
  return -0.5 * kLog2Pi - log_sd - 0.5 * z * z;
}

// Student-t(df, mean, scale) log density.
template <class T>
T log_student_t(const T& x, double df, double mean, double scale) {
  using std::log;
  using std::log1p;
  T z = (x - mean) / scale;
  double log_const = std::lgamma((df + 1.0) / 2.0) - std::lgamma(df / 2.0) -
                     0.5 * std::log(df * M_PI) - std::log(scale);
  return log_const - ((df + 1.0) / 2.0) * log1p(z * z / df);
}

// Cauchy(location, scale) log density.
template <class T>
T log_cauchy(const T& x, double location, double scale) {
  using std::log1p;
  T z = (x - location) / scale;
  return -kLogPi - std::log(scale) - log1p(z * z);
}

// Half-Normal(0, sd) prior on sigma = exp(log_sigma), as a density of
// log_sigma: the Half-Normal density of sigma times the Jacobian sigma.
template <class T>
T log_half_normal_log_scale(const T& log_sigma, double sd) {
  using std::exp;
  T z = exp(log_sigma) / sd;
  return kLog2 - 0.5 * kLog2Pi - std::log(sd) - 0.5 * z * z + log_sigma;
}

// Half-Cauchy(0, scale) prior on sigma = exp(log_sigma), as a density of
// log_sigma.
template <class T>
T log_half_cauchy_log_scale(const T& log_sigma, double scale) {
  using std::exp;
  using std::log1p;
  T z = exp(log_sigma) / scale;
  return kLog2 - kLogPi - std::log(scale) - log1p(z * z) + log_sigma;
}

// Prior of an unbounded parameter (a regression coefficient).
template <class T>
T log_prior_unbounded(const T& x, int type, const Rcpp::NumericVector& p) {
  switch (type) {
    case 0: return log_normal(x, p[0], p[1]);
    case 1: return log_student_t(x, p[0], p[1], p[2]);
    case 2: return log_cauchy(x, p[0], p[1]);
    default:
      Rcpp::stop("Invalid prior type for unbounded parameter: %d", type);
  }
}

// Prior of a positive parameter sigma, as a density of log(sigma).
template <class T>
T log_prior_positive(const T& log_sigma, int type,
                     const Rcpp::NumericVector& p) {
  switch (type) {
    case 3: return log_half_normal_log_scale(log_sigma, p[0]);
    case 4: return log_half_cauchy_log_scale(log_sigma, p[0]);
    default:
      Rcpp::stop("Invalid prior type for positive parameter: %d", type);
  }
}

// Positive-parameter prior on a real-valued scale whose sign is not
// identified (the diagonal of a Cholesky factor): the half distribution
// reflected to the whole real line, Normal(0, sd) or Cauchy(0, scale).
template <class T>
T log_prior_symmetric_scale(const T& x, int type,
                            const Rcpp::NumericVector& p) {
  switch (type) {
    case 3: return log_normal(x, 0.0, p[0]);
    case 4: return log_cauchy(x, 0.0, p[0]);
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
  T beta(const T& x) const {
    return log_prior_unbounded(x, beta_type_, beta_params_);
  }

  template <class T>
  T sigma(const T& log_sigma) const {
    return log_prior_positive(log_sigma, sigma_type_, sigma_params_);
  }

  template <class T>
  T re_sd(const T& log_sigma) const {
    return log_prior_positive(log_sigma, re_type_, re_params_);
  }

  template <class T>
  T re_scale(const T& x) const {
    return log_prior_symmetric_scale(x, re_type_, re_params_);
  }

  double cor_eta() const { return cor_params_[0]; }

  // Scale of the random-effect SD prior.
  double re_scale_parameter() const { return re_params_[0]; }

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
