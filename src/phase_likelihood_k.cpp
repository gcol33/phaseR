#include <Rcpp.h>
#include <cmath>
#include <vector>

using namespace Rcpp;

// Log-sum-exp for numerical stability
inline double log_sum_exp_k(double a, double b) {
  if (a == -INFINITY) return b;
  if (b == -INFINITY) return a;
  if (a > b) {
    return a + std::log1p(std::exp(b - a));
  } else {
    return b + std::log1p(std::exp(a - b));
  }
}

// Log-sum-exp for a vector
inline double log_sum_exp_vec(const std::vector<double>& v) {
  double max_val = -INFINITY;
  for (size_t i = 0; i < v.size(); i++) {
    if (v[i] > max_val) max_val = v[i];
  }
  if (max_val == -INFINITY) return -INFINITY;

  double sum = 0.0;
  for (size_t i = 0; i < v.size(); i++) {
    sum += std::exp(v[i] - max_val);
  }
  return max_val + std::log(sum);
}

// Inverse logit
inline double invlogit_k(double x) {
  return 1.0 / (1.0 + std::exp(-x));
}

// Log density of normal
inline double log_dnorm_k(double y, double mu, double sigma) {
  double z = (y - mu) / sigma;
  return -0.5 * std::log(2.0 * M_PI) - std::log(sigma) - 0.5 * z * z;
}


// [[Rcpp::export]]
double phase_log_likelihood_k(
    NumericVector params,
    List data,
    IntegerVector n_dyn_coef,    // length k: coefficients per phase
    IntegerVector n_trans_coef,  // length k-1: coefficients per transition
    int k_phases
) {
  // Unpack parameters
  int idx = 0;

  // Transition coefficients: one set per transition (k-1 transitions)
  std::vector<NumericVector> beta_trans(k_phases - 1);
  for (int tr = 0; tr < k_phases - 1; tr++) {
    beta_trans[tr] = NumericVector(n_trans_coef[tr]);
    for (int j = 0; j < n_trans_coef[tr]; j++) {
      beta_trans[tr][j] = params[idx++];
    }
  }

  // Dynamics coefficients: one set per phase
  std::vector<NumericVector> beta_dyn(k_phases);
  for (int p = 0; p < k_phases; p++) {
    beta_dyn[p] = NumericVector(n_dyn_coef[p]);
    for (int j = 0; j < n_dyn_coef[p]; j++) {
      beta_dyn[p][j] = params[idx++];
    }
  }

  // Sigmas: one per phase
  std::vector<double> sigma(k_phases);
  for (int p = 0; p < k_phases; p++) {
    double log_sigma = params[idx++];
    sigma[p] = std::exp(log_sigma);
  }

  // Data
  NumericVector y = data["y"];
  List X_trans_list = data["X_trans_list"];  // List of matrices
  List X_dyn_list = data["X_dyn_list"];      // List of matrices

  int n_units = data["n_units"];
  IntegerVector unit_start = data["unit_start"];
  IntegerVector unit_end = data["unit_end"];

  // Forward algorithm
  double total_ll = 0.0;

  for (int i = 0; i < n_units; i++) {
    int start = unit_start[i];
    int end = unit_end[i];

    // Forward probabilities (log scale): alpha[p] = log P(S_t = p, Y_{1:t})
    std::vector<double> alpha(k_phases);

    // Initial: always in phase 0
    NumericMatrix X_dyn_0 = X_dyn_list[0];
    double mu_init = 0.0;
    for (int j = 0; j < n_dyn_coef[0]; j++) {
      mu_init += beta_dyn[0][j] * X_dyn_0(start, j);
    }
    alpha[0] = log_dnorm_k(y[start], mu_init, sigma[0]);
    for (int p = 1; p < k_phases; p++) {
      alpha[p] = -INFINITY;
    }

    // Forward pass
    for (int t = start + 1; t <= end; t++) {
      std::vector<double> new_alpha(k_phases);

      for (int p = 0; p < k_phases; p++) {
        // Emission probability for phase p
        NumericMatrix X_dyn_p = X_dyn_list[p];
        double mu_p = 0.0;
        for (int j = 0; j < n_dyn_coef[p]; j++) {
          mu_p += beta_dyn[p][j] * X_dyn_p(t, j);
        }
        double log_emit_p = log_dnorm_k(y[t], mu_p, sigma[p]);

        // Contributions from previous states
        std::vector<double> contributions;

        // From same phase (stayed)
        if (p < k_phases - 1) {
          // Non-absorbing: P(stay) = 1 - P(trans to next)
          NumericMatrix X_trans_p = X_trans_list[p];
          double eta = 0.0;
          for (int j = 0; j < n_trans_coef[p]; j++) {
            eta += beta_trans[p][j] * X_trans_p(t, j);
          }
          double p_trans = invlogit_k(eta);
          contributions.push_back(alpha[p] + std::log(1.0 - p_trans) + log_emit_p);
        } else {
          // Absorbing: P(stay) = 1
          contributions.push_back(alpha[p] + log_emit_p);
        }

        // From previous phase (transitioned)
        if (p > 0) {
          NumericMatrix X_trans_prev = X_trans_list[p - 1];
          double eta = 0.0;
          for (int j = 0; j < n_trans_coef[p - 1]; j++) {
            eta += beta_trans[p - 1][j] * X_trans_prev(t, j);
          }
          double p_trans = invlogit_k(eta);
          contributions.push_back(alpha[p - 1] + std::log(p_trans) + log_emit_p);
        }

        new_alpha[p] = log_sum_exp_vec(contributions);
      }

      alpha = new_alpha;
    }

    // Total likelihood for this unit
    total_ll += log_sum_exp_vec(alpha);
  }

  return total_ll;
}


// [[Rcpp::export]]
double phase_log_prior_k(
    NumericVector params,
    IntegerVector n_dyn_coef,
    IntegerVector n_trans_coef,
    int k_phases
) {
  double log_prior = 0.0;
  int idx = 0;

  // Transition coefficients priors
  for (int tr = 0; tr < k_phases - 1; tr++) {
    for (int j = 0; j < n_trans_coef[tr]; j++) {
      double x = params[idx++];
      log_prior += -0.5 * std::log(2.0 * M_PI * 6.25) - 0.5 * x * x / 6.25;
    }
  }

  // Dynamics coefficients priors
  for (int p = 0; p < k_phases; p++) {
    for (int j = 0; j < n_dyn_coef[p]; j++) {
      double x = params[idx++];
      log_prior += -0.5 * std::log(2.0 * M_PI * 6.25) - 0.5 * x * x / 6.25;
    }
  }

  // Sigma priors
  for (int p = 0; p < k_phases; p++) {
    double log_sigma = params[idx++];
    log_prior += -0.5 * std::log(2.0 * M_PI) - 0.5 * log_sigma * log_sigma;
  }

  return log_prior;
}


// [[Rcpp::export]]
double phase_log_posterior_k(
    NumericVector params,
    List data,
    IntegerVector n_dyn_coef,
    IntegerVector n_trans_coef,
    int k_phases
) {
  return phase_log_likelihood_k(params, data, n_dyn_coef, n_trans_coef, k_phases) +
         phase_log_prior_k(params, n_dyn_coef, n_trans_coef, k_phases);
}
