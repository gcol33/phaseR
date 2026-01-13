#include <Rcpp.h>
#include <cmath>
#include <vector>

using namespace Rcpp;

// Log-sum-exp
inline double log_sum_exp_glm(double a, double b) {
  if (a == -INFINITY) return b;
  if (b == -INFINITY) return a;
  if (a > b) {
    return a + std::log1p(std::exp(b - a));
  } else {
    return b + std::log1p(std::exp(a - b));
  }
}

// Inverse logit
inline double invlogit_glm(double x) {
  return 1.0 / (1.0 + std::exp(-x));
}

// Log factorial (using lgamma)
inline double log_factorial(int n) {
  return std::lgamma(n + 1);
}

// Log density of normal
inline double log_dnorm_glm(double y, double mu, double sigma) {
  double z = (y - mu) / sigma;
  return -0.5 * std::log(2.0 * M_PI) - std::log(sigma) - 0.5 * z * z;
}

// Log PMF of Poisson
inline double log_dpois(double y, double lambda) {
  if (lambda <= 0) return -INFINITY;
  int y_int = static_cast<int>(y);
  return y_int * std::log(lambda) - lambda - log_factorial(y_int);
}

// Log PMF of Binomial (for proportion data with known n)
inline double log_dbinom(double y, int n, double p) {
  if (p <= 0 || p >= 1) return -INFINITY;
  int k = static_cast<int>(y * n + 0.5);  // y is proportion, convert to count
  return log_factorial(n) - log_factorial(k) - log_factorial(n - k) +
         k * std::log(p) + (n - k) * std::log(1 - p);
}

// Emission log-likelihood based on family
// family: 0 = gaussian, 1 = poisson, 2 = binomial
inline double log_emission(double y, double mu, double sigma, int family, int n_trials = 1) {
  switch (family) {
    case 0:  // gaussian
      return log_dnorm_glm(y, mu, sigma);
    case 1:  // poisson (log link)
      return log_dpois(y, std::exp(mu));
    case 2:  // binomial (logit link)
      return log_dbinom(y, n_trials, invlogit_glm(mu));
    default:
      return log_dnorm_glm(y, mu, sigma);
  }
}


// [[Rcpp::export]]
double phase_log_likelihood_glm(
    NumericVector params,
    List data,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    int family_0,      // 0=gaussian, 1=poisson, 2=binomial
    int family_1,
    int n_trials_0,    // for binomial
    int n_trials_1
) {
  int idx = 0;

  NumericVector beta_trans(n_trans_coef);
  for (int j = 0; j < n_trans_coef; j++) {
    beta_trans[j] = params[idx++];
  }

  NumericVector beta_0(n_dyn_coef_0);
  for (int j = 0; j < n_dyn_coef_0; j++) {
    beta_0[j] = params[idx++];
  }

  NumericVector beta_1(n_dyn_coef_1);
  for (int j = 0; j < n_dyn_coef_1; j++) {
    beta_1[j] = params[idx++];
  }

  double log_sigma_0 = params[idx++];
  double log_sigma_1 = params[idx++];
  double sigma_0 = std::exp(log_sigma_0);
  double sigma_1 = std::exp(log_sigma_1);

  // Use eager deep copies to prevent R GC issues
  std::vector<double> y = Rcpp::as<std::vector<double>>(data["y"]);
  NumericMatrix X_trans = Rcpp::as<NumericMatrix>(data["X_trans"]);
  NumericMatrix X_dyn = Rcpp::as<NumericMatrix>(data["X_dyn"]);
  int n_units = Rcpp::as<int>(data["n_units"]);
  std::vector<int> unit_start = Rcpp::as<std::vector<int>>(data["unit_start"]);
  std::vector<int> unit_end = Rcpp::as<std::vector<int>>(data["unit_end"]);

  double total_ll = 0.0;

  for (int i = 0; i < n_units; i++) {
    int start = unit_start[i];
    int end = unit_end[i];

    // Forward algorithm
    double mu_0_init = 0.0;
    for (int j = 0; j < n_dyn_coef_0; j++) {
      mu_0_init += beta_0[j] * X_dyn(start, j);
    }

    double alpha_0 = log_emission(y[start], mu_0_init, sigma_0, family_0, n_trials_0);
    double alpha_1 = -INFINITY;

    for (int t = start + 1; t <= end; t++) {
      double eta_trans = 0.0;
      for (int j = 0; j < n_trans_coef; j++) {
        eta_trans += beta_trans[j] * X_trans(t, j);
      }
      double p_trans = invlogit_glm(eta_trans);

      double mu_0 = 0.0, mu_1 = 0.0;
      for (int j = 0; j < n_dyn_coef_0; j++) {
        mu_0 += beta_0[j] * X_dyn(t, j);
      }
      for (int j = 0; j < n_dyn_coef_1; j++) {
        mu_1 += beta_1[j] * X_dyn(t, j);
      }

      double log_emit_0 = log_emission(y[t], mu_0, sigma_0, family_0, n_trials_0);
      double log_emit_1 = log_emission(y[t], mu_1, sigma_1, family_1, n_trials_1);

      double new_alpha_0 = alpha_0 + std::log(1.0 - p_trans) + log_emit_0;
      double from_0 = alpha_0 + std::log(p_trans) + log_emit_1;
      double from_1 = alpha_1 + log_emit_1;
      double new_alpha_1 = log_sum_exp_glm(from_0, from_1);

      alpha_0 = new_alpha_0;
      alpha_1 = new_alpha_1;
    }

    total_ll += log_sum_exp_glm(alpha_0, alpha_1);
  }

  return total_ll;
}


// [[Rcpp::export]]
double phase_log_prior_glm(
    NumericVector params,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    int family_0,
    int family_1
) {
  double log_prior = 0.0;
  int idx = 0;

  // Transition coefficients
  for (int j = 0; j < n_trans_coef; j++) {
    double x = params[idx++];
    log_prior += -0.5 * std::log(2.0 * M_PI * 6.25) - 0.5 * x * x / 6.25;
  }

  // Dynamics coefficients
  for (int j = 0; j < n_dyn_coef_0; j++) {
    double x = params[idx++];
    log_prior += -0.5 * std::log(2.0 * M_PI * 6.25) - 0.5 * x * x / 6.25;
  }
  for (int j = 0; j < n_dyn_coef_1; j++) {
    double x = params[idx++];
    log_prior += -0.5 * std::log(2.0 * M_PI * 6.25) - 0.5 * x * x / 6.25;
  }

  // Sigma priors (only meaningful for Gaussian)
  double log_sigma_0 = params[idx++];
  double log_sigma_1 = params[idx++];
  if (family_0 == 0) {
    log_prior += -0.5 * std::log(2.0 * M_PI) - 0.5 * log_sigma_0 * log_sigma_0;
  }
  if (family_1 == 0) {
    log_prior += -0.5 * std::log(2.0 * M_PI) - 0.5 * log_sigma_1 * log_sigma_1;
  }

  return log_prior;
}


// [[Rcpp::export]]
double phase_log_posterior_glm(
    NumericVector params,
    List data,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    int family_0,
    int family_1,
    int n_trials_0,
    int n_trials_1
) {
  return phase_log_likelihood_glm(params, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1,
                                   family_0, family_1, n_trials_0, n_trials_1) +
         phase_log_prior_glm(params, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1,
                              family_0, family_1);
}
