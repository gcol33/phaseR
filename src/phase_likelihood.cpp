#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// Log-sum-exp for numerical stability
inline double log_sum_exp(double a, double b) {
  if (a > b) {
    return a + std::log1p(std::exp(b - a));
  } else {
    return b + std::log1p(std::exp(a - b));
  }
}

// Inverse logit
inline double invlogit(double x) {
  return 1.0 / (1.0 + std::exp(-x));
}

// Log density of normal
inline double log_dnorm(double y, double mu, double sigma) {
  double z = (y - mu) / sigma;
  return -0.5 * std::log(2.0 * M_PI) - std::log(sigma) - 0.5 * z * z;
}

// [[Rcpp::export]]
double phase_log_likelihood(
    NumericVector params,
    List data,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1
) {

  // Unpack parameters
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

  // Unpack data
  IntegerVector id = data["id"];
  NumericVector y = data["y"];
  NumericMatrix X_trans = data["X_trans"];
  NumericMatrix X_dyn = data["X_dyn"];

  int N = y.size();
  int n_units = data["n_units"];

  // Forward algorithm for each unit
  double total_ll = 0.0;

  IntegerVector unit_start = data["unit_start"];
  IntegerVector unit_end = data["unit_end"];

  for (int i = 0; i < n_units; i++) {

    int start = unit_start[i];
    int end = unit_end[i];

    // Forward probabilities (log scale)
    // alpha[0] = log P(S_t = 0, Y_{1:t})
    // alpha[1] = log P(S_t = 1, Y_{1:t})
    double alpha_0, alpha_1;

    // Initial time point: always in phase 0
    double mu_0_init = 0.0;
    for (int j = 0; j < n_dyn_coef_0; j++) {
      mu_0_init += beta_0[j] * X_dyn(start, j);
    }
    alpha_0 = log_dnorm(y[start], mu_0_init, sigma_0);
    alpha_1 = -INFINITY;  // Can't be in phase 1 at t=1

    // Forward pass
    for (int t = start + 1; t <= end; t++) {

      // Transition probability
      double eta_trans = 0.0;
      for (int j = 0; j < n_trans_coef; j++) {
        eta_trans += beta_trans[j] * X_trans(t, j);
      }
      double p_trans = invlogit(eta_trans);

      // Emission probabilities
      double mu_0 = 0.0, mu_1 = 0.0;
      for (int j = 0; j < n_dyn_coef_0; j++) {
        mu_0 += beta_0[j] * X_dyn(t, j);
      }
      for (int j = 0; j < n_dyn_coef_1; j++) {
        mu_1 += beta_1[j] * X_dyn(t, j);
      }

      double log_emit_0 = log_dnorm(y[t], mu_0, sigma_0);
      double log_emit_1 = log_dnorm(y[t], mu_1, sigma_1);

      // Update forward probabilities
      // New alpha_0: came from phase 0, stayed in phase 0
      double new_alpha_0 = alpha_0 + std::log(1.0 - p_trans) + log_emit_0;

      // New alpha_1: came from phase 0 and transitioned, or came from phase 1
      double from_0 = alpha_0 + std::log(p_trans) + log_emit_1;
      double from_1 = alpha_1 + log_emit_1;  // P(stay in 1) = 1
      double new_alpha_1 = log_sum_exp(from_0, from_1);

      alpha_0 = new_alpha_0;
      alpha_1 = new_alpha_1;
    }

    // Total likelihood for this unit
    total_ll += log_sum_exp(alpha_0, alpha_1);
  }

  return total_ll;
}

// [[Rcpp::export]]
double phase_log_prior(
    NumericVector params,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1
) {

  double log_prior = 0.0;
  int idx = 0;

  // beta_trans ~ Normal(0, 2.5)
  for (int j = 0; j < n_trans_coef; j++) {
    double x = params[idx++];
    log_prior += -0.5 * std::log(2.0 * M_PI * 6.25) - 0.5 * x * x / 6.25;
  }

  // beta_0 ~ Normal(0, 2.5)
  for (int j = 0; j < n_dyn_coef_0; j++) {
    double x = params[idx++];
    log_prior += -0.5 * std::log(2.0 * M_PI * 6.25) - 0.5 * x * x / 6.25;
  }

  // beta_1 ~ Normal(0, 2.5)
  for (int j = 0; j < n_dyn_coef_1; j++) {
    double x = params[idx++];
    log_prior += -0.5 * std::log(2.0 * M_PI * 6.25) - 0.5 * x * x / 6.25;
  }

  // log_sigma ~ Normal(0, 1)
  double log_sigma_0 = params[idx++];
  double log_sigma_1 = params[idx++];
  log_prior += -0.5 * std::log(2.0 * M_PI) - 0.5 * log_sigma_0 * log_sigma_0;
  log_prior += -0.5 * std::log(2.0 * M_PI) - 0.5 * log_sigma_1 * log_sigma_1;

  return log_prior;
}

// [[Rcpp::export]]
double phase_log_posterior(
    NumericVector params,
    List data,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1
) {
  return phase_log_likelihood(params, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1) +
         phase_log_prior(params, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1);
}
