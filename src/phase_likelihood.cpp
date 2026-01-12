#include <Rcpp.h>
#include <cmath>
#include "prior_functions.h"

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


// ============================================================================
// Custom Prior Versions
// ============================================================================

// [[Rcpp::export]]
double phase_log_prior_custom(
    NumericVector params,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    int beta_type,
    NumericVector beta_params,
    int sigma_type,
    NumericVector sigma_params
) {

  double log_prior = 0.0;
  int idx = 0;

  // beta_trans priors
  for (int j = 0; j < n_trans_coef; j++) {
    log_prior += prior::log_prior_unbounded(params[idx++], beta_type, beta_params);
  }

  // beta_0 priors
  for (int j = 0; j < n_dyn_coef_0; j++) {
    log_prior += prior::log_prior_unbounded(params[idx++], beta_type, beta_params);
  }

  // beta_1 priors
  for (int j = 0; j < n_dyn_coef_1; j++) {
    log_prior += prior::log_prior_unbounded(params[idx++], beta_type, beta_params);
  }

  // sigma priors (stored as log_sigma)
  double log_sigma_0 = params[idx++];
  double log_sigma_1 = params[idx++];
  log_prior += prior::log_prior_positive(log_sigma_0, sigma_type, sigma_params);
  log_prior += prior::log_prior_positive(log_sigma_1, sigma_type, sigma_params);

  return log_prior;
}

// [[Rcpp::export]]
double phase_log_posterior_custom(
    NumericVector params,
    List data,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    int beta_type,
    NumericVector beta_params,
    int sigma_type,
    NumericVector sigma_params
) {
  return phase_log_likelihood(params, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1) +
         phase_log_prior_custom(params, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1,
                                beta_type, beta_params, sigma_type, sigma_params);
}


// ============================================================================
// Random Effects Versions
// ============================================================================

// [[Rcpp::export]]
double phase_log_likelihood_re(
    NumericVector params,
    List data,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    int n_trans_re,
    int n_dyn_re_0,
    int n_dyn_re_1,
    bool has_trans_re,
    bool has_dyn_re_0,
    bool has_dyn_re_1
) {

  int idx = 0;

  // Fixed effects
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

  // Random effects
  NumericVector u_trans(n_trans_re);
  if (has_trans_re) {
    for (int g = 0; g < n_trans_re; g++) {
      u_trans[g] = params[idx++];
    }
    idx++;  // Skip log_sigma_trans_re (used only in prior)
  }

  NumericVector v_0(n_dyn_re_0);
  if (has_dyn_re_0) {
    for (int g = 0; g < n_dyn_re_0; g++) {
      v_0[g] = params[idx++];
    }
    idx++;  // Skip log_sigma_dyn_re_0
  }

  NumericVector v_1(n_dyn_re_1);
  if (has_dyn_re_1) {
    for (int g = 0; g < n_dyn_re_1; g++) {
      v_1[g] = params[idx++];
    }
    idx++;  // Skip log_sigma_dyn_re_1
  }

  double sigma_0 = std::exp(log_sigma_0);
  double sigma_1 = std::exp(log_sigma_1);

  // Data
  NumericVector y = data["y"];
  NumericMatrix X_trans = data["X_trans"];
  NumericMatrix X_dyn = data["X_dyn"];

  int n_units = data["n_units"];
  IntegerVector unit_start = data["unit_start"];
  IntegerVector unit_end = data["unit_end"];

  // RE indices (1-indexed from R)
  IntegerVector trans_re_idx;
  IntegerVector dyn_re_idx_0;
  IntegerVector dyn_re_idx_1;

  if (has_trans_re) trans_re_idx = data["trans_re_idx"];
  if (has_dyn_re_0) dyn_re_idx_0 = data["dyn_re_idx_0"];
  if (has_dyn_re_1) dyn_re_idx_1 = data["dyn_re_idx_1"];

  double total_ll = 0.0;

  for (int i = 0; i < n_units; i++) {

    int start = unit_start[i];
    int end = unit_end[i];

    // Get group indices for this unit (convert from 1-indexed to 0-indexed)
    int g_trans = has_trans_re ? trans_re_idx[start] - 1 : 0;
    int g_dyn_0 = has_dyn_re_0 ? dyn_re_idx_0[start] - 1 : 0;
    int g_dyn_1 = has_dyn_re_1 ? dyn_re_idx_1[start] - 1 : 0;

    double alpha_0, alpha_1;

    // Initial time point
    double mu_0_init = 0.0;
    for (int j = 0; j < n_dyn_coef_0; j++) {
      mu_0_init += beta_0[j] * X_dyn(start, j);
    }
    if (has_dyn_re_0) {
      mu_0_init += v_0[g_dyn_0];
    }

    alpha_0 = log_dnorm(y[start], mu_0_init, sigma_0);
    alpha_1 = -INFINITY;

    // Forward pass
    for (int t = start + 1; t <= end; t++) {

      double eta_trans = 0.0;
      for (int j = 0; j < n_trans_coef; j++) {
        eta_trans += beta_trans[j] * X_trans(t, j);
      }
      if (has_trans_re) {
        eta_trans += u_trans[g_trans];
      }
      double p_trans = invlogit(eta_trans);

      double mu_0 = 0.0, mu_1 = 0.0;
      for (int j = 0; j < n_dyn_coef_0; j++) {
        mu_0 += beta_0[j] * X_dyn(t, j);
      }
      for (int j = 0; j < n_dyn_coef_1; j++) {
        mu_1 += beta_1[j] * X_dyn(t, j);
      }
      if (has_dyn_re_0) {
        mu_0 += v_0[g_dyn_0];
      }
      if (has_dyn_re_1) {
        mu_1 += v_1[g_dyn_1];
      }

      double log_emit_0 = log_dnorm(y[t], mu_0, sigma_0);
      double log_emit_1 = log_dnorm(y[t], mu_1, sigma_1);

      double new_alpha_0 = alpha_0 + std::log(1.0 - p_trans) + log_emit_0;
      double from_0 = alpha_0 + std::log(p_trans) + log_emit_1;
      double from_1 = alpha_1 + log_emit_1;
      double new_alpha_1 = log_sum_exp(from_0, from_1);

      alpha_0 = new_alpha_0;
      alpha_1 = new_alpha_1;
    }

    total_ll += log_sum_exp(alpha_0, alpha_1);
  }

  return total_ll;
}


// [[Rcpp::export]]
double phase_log_prior_re(
    NumericVector params,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    int n_trans_re,
    int n_dyn_re_0,
    int n_dyn_re_1,
    bool has_trans_re,
    bool has_dyn_re_0,
    bool has_dyn_re_1
) {

  double log_prior = 0.0;
  int idx = 0;

  // Fixed effects priors
  for (int j = 0; j < n_trans_coef; j++) {
    double x = params[idx++];
    log_prior += -0.5 * std::log(2.0 * M_PI * 6.25) - 0.5 * x * x / 6.25;
  }
  for (int j = 0; j < n_dyn_coef_0; j++) {
    double x = params[idx++];
    log_prior += -0.5 * std::log(2.0 * M_PI * 6.25) - 0.5 * x * x / 6.25;
  }
  for (int j = 0; j < n_dyn_coef_1; j++) {
    double x = params[idx++];
    log_prior += -0.5 * std::log(2.0 * M_PI * 6.25) - 0.5 * x * x / 6.25;
  }

  double log_sigma_0 = params[idx++];
  double log_sigma_1 = params[idx++];
  log_prior += -0.5 * std::log(2.0 * M_PI) - 0.5 * log_sigma_0 * log_sigma_0;
  log_prior += -0.5 * std::log(2.0 * M_PI) - 0.5 * log_sigma_1 * log_sigma_1;

  // Random effects priors
  if (has_trans_re && n_trans_re > 0) {
    double log_sigma_re = params[idx + n_trans_re];
    double sigma_re = std::exp(log_sigma_re);
    double sigma_re_sq = sigma_re * sigma_re;
    for (int g = 0; g < n_trans_re; g++) {
      double u = params[idx++];
      log_prior += -0.5 * std::log(2.0 * M_PI) - std::log(sigma_re) - 0.5 * u * u / sigma_re_sq;
    }
    // Prior on log_sigma_re ~ Normal(0, 0.5)
    log_prior += -0.5 * std::log(2.0 * M_PI * 0.25) - 0.5 * log_sigma_re * log_sigma_re / 0.25;
    idx++;
  }

  if (has_dyn_re_0 && n_dyn_re_0 > 0) {
    double log_sigma_re = params[idx + n_dyn_re_0];
    double sigma_re = std::exp(log_sigma_re);
    double sigma_re_sq = sigma_re * sigma_re;
    for (int g = 0; g < n_dyn_re_0; g++) {
      double v = params[idx++];
      log_prior += -0.5 * std::log(2.0 * M_PI) - std::log(sigma_re) - 0.5 * v * v / sigma_re_sq;
    }
    log_prior += -0.5 * std::log(2.0 * M_PI * 0.25) - 0.5 * log_sigma_re * log_sigma_re / 0.25;
    idx++;
  }

  if (has_dyn_re_1 && n_dyn_re_1 > 0) {
    double log_sigma_re = params[idx + n_dyn_re_1];
    double sigma_re = std::exp(log_sigma_re);
    double sigma_re_sq = sigma_re * sigma_re;
    for (int g = 0; g < n_dyn_re_1; g++) {
      double v = params[idx++];
      log_prior += -0.5 * std::log(2.0 * M_PI) - std::log(sigma_re) - 0.5 * v * v / sigma_re_sq;
    }
    log_prior += -0.5 * std::log(2.0 * M_PI * 0.25) - 0.5 * log_sigma_re * log_sigma_re / 0.25;
    idx++;
  }

  return log_prior;
}


// [[Rcpp::export]]
double phase_log_posterior_re(
    NumericVector params,
    List data,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    int n_trans_re,
    int n_dyn_re_0,
    int n_dyn_re_1,
    bool has_trans_re,
    bool has_dyn_re_0,
    bool has_dyn_re_1
) {
  return phase_log_likelihood_re(params, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1,
                                  n_trans_re, n_dyn_re_0, n_dyn_re_1,
                                  has_trans_re, has_dyn_re_0, has_dyn_re_1) +
         phase_log_prior_re(params, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1,
                            n_trans_re, n_dyn_re_0, n_dyn_re_1,
                            has_trans_re, has_dyn_re_0, has_dyn_re_1);
}


// ============================================================================
// Custom Prior Versions for Random Effects
// ============================================================================

// [[Rcpp::export]]
double phase_log_prior_re_custom(
    NumericVector params,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    int n_trans_re,
    int n_dyn_re_0,
    int n_dyn_re_1,
    bool has_trans_re,
    bool has_dyn_re_0,
    bool has_dyn_re_1,
    int beta_type,
    NumericVector beta_params,
    int sigma_type,
    NumericVector sigma_params,
    int re_type,
    NumericVector re_params
) {

  double log_prior = 0.0;
  int idx = 0;

  // Fixed effects priors
  for (int j = 0; j < n_trans_coef; j++) {
    log_prior += prior::log_prior_unbounded(params[idx++], beta_type, beta_params);
  }
  for (int j = 0; j < n_dyn_coef_0; j++) {
    log_prior += prior::log_prior_unbounded(params[idx++], beta_type, beta_params);
  }
  for (int j = 0; j < n_dyn_coef_1; j++) {
    log_prior += prior::log_prior_unbounded(params[idx++], beta_type, beta_params);
  }

  // Residual sigma priors
  double log_sigma_0 = params[idx++];
  double log_sigma_1 = params[idx++];
  log_prior += prior::log_prior_positive(log_sigma_0, sigma_type, sigma_params);
  log_prior += prior::log_prior_positive(log_sigma_1, sigma_type, sigma_params);

  // Random effects priors
  if (has_trans_re && n_trans_re > 0) {
    double log_sigma_re = params[idx + n_trans_re];
    double sigma_re = std::exp(log_sigma_re);
    double sigma_re_sq = sigma_re * sigma_re;
    // RE values ~ Normal(0, sigma_re)
    for (int g = 0; g < n_trans_re; g++) {
      double u = params[idx++];
      log_prior += -0.5 * std::log(2.0 * M_PI) - std::log(sigma_re) - 0.5 * u * u / sigma_re_sq;
    }
    // Prior on RE SD
    log_prior += prior::log_prior_positive(log_sigma_re, re_type, re_params);
    idx++;
  }

  if (has_dyn_re_0 && n_dyn_re_0 > 0) {
    double log_sigma_re = params[idx + n_dyn_re_0];
    double sigma_re = std::exp(log_sigma_re);
    double sigma_re_sq = sigma_re * sigma_re;
    for (int g = 0; g < n_dyn_re_0; g++) {
      double v = params[idx++];
      log_prior += -0.5 * std::log(2.0 * M_PI) - std::log(sigma_re) - 0.5 * v * v / sigma_re_sq;
    }
    log_prior += prior::log_prior_positive(log_sigma_re, re_type, re_params);
    idx++;
  }

  if (has_dyn_re_1 && n_dyn_re_1 > 0) {
    double log_sigma_re = params[idx + n_dyn_re_1];
    double sigma_re = std::exp(log_sigma_re);
    double sigma_re_sq = sigma_re * sigma_re;
    for (int g = 0; g < n_dyn_re_1; g++) {
      double v = params[idx++];
      log_prior += -0.5 * std::log(2.0 * M_PI) - std::log(sigma_re) - 0.5 * v * v / sigma_re_sq;
    }
    log_prior += prior::log_prior_positive(log_sigma_re, re_type, re_params);
    idx++;
  }

  return log_prior;
}


// [[Rcpp::export]]
double phase_log_posterior_re_custom(
    NumericVector params,
    List data,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    int n_trans_re,
    int n_dyn_re_0,
    int n_dyn_re_1,
    bool has_trans_re,
    bool has_dyn_re_0,
    bool has_dyn_re_1,
    int beta_type,
    NumericVector beta_params,
    int sigma_type,
    NumericVector sigma_params,
    int re_type,
    NumericVector re_params
) {
  return phase_log_likelihood_re(params, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1,
                                  n_trans_re, n_dyn_re_0, n_dyn_re_1,
                                  has_trans_re, has_dyn_re_0, has_dyn_re_1) +
         phase_log_prior_re_custom(params, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1,
                                   n_trans_re, n_dyn_re_0, n_dyn_re_1,
                                   has_trans_re, has_dyn_re_0, has_dyn_re_1,
                                   beta_type, beta_params, sigma_type, sigma_params,
                                   re_type, re_params);
}
