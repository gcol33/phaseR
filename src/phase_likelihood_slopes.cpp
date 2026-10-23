#include <Rcpp.h>
#include <cmath>
#include <vector>
#include "phase_data.h"
#include "prior_functions.h"

using namespace Rcpp;

// Log-sum-exp for numerical stability
inline double log_sum_exp_s(double a, double b) {
  if (a > b) {
    return a + std::log1p(std::exp(b - a));
  } else {
    return b + std::log1p(std::exp(a - b));
  }
}

inline double invlogit_s(double x) {
  return 1.0 / (1.0 + std::exp(-x));
}

inline double log_dnorm_s(double y, double mu, double sigma) {
  double z = (y - mu) / sigma;
  return -0.5 * std::log(2.0 * M_PI) - std::log(sigma) - 0.5 * z * z;
}

// LKJ prior on correlation matrix (via Cholesky)
// This is the log-prior on the Cholesky factor L of a correlation matrix
// log p(L | eta) = (eta - 1) * sum of log diag(L) + adjustment
inline double log_lkj_cholesky(const std::vector<double>& L_vals, int dim, double eta = 1.0) {
  // L_vals contains lower triangle in column-major order
  // For eta = 1, this is uniform over correlation matrices
  // For eta > 1, favors identity matrix
  double log_p = 0.0;
  int idx = 0;
  for (int j = 0; j < dim; j++) {
    for (int i = j; i < dim; i++) {
      if (i == j) {
        // Diagonal element
        double L_jj = L_vals[idx];
        // Jacobian adjustment: (dim - j) * log(L_jj)
        log_p += (dim - j - 1 + 2.0 * (eta - 1.0)) * std::log(std::abs(L_jj));
      }
      idx++;
    }
  }
  return log_p;
}


// Compute RE contribution for a group using Z matrix
// re_coefs: vector of length n_coefs for this group
// Z_row: vector of length n_coefs (the row of Z matrix for this observation)
inline double compute_re_contribution(
    const std::vector<double>& re_coefs,
    const NumericVector& Z_row
) {
  double contrib = 0.0;
  int n = re_coefs.size();
  for (int k = 0; k < n; k++) {
    contrib += re_coefs[k] * Z_row[k];
  }
  return contrib;
}


// Transform raw eta values to actual RE coefficients using Cholesky
// L is lower triangular, stored as vector in column-major order
// eta ~ N(0, I), then L * eta ~ N(0, LL')
inline std::vector<double> apply_cholesky(
    const std::vector<double>& eta,
    const std::vector<double>& L_vals,
    int dim
) {
  std::vector<double> result(dim, 0.0);
  int L_idx = 0;

  // L is stored column-major: L11, L21, L31, ..., L22, L32, ...
  for (int j = 0; j < dim; j++) {
    for (int i = j; i < dim; i++) {
      result[i] += L_vals[L_idx] * eta[j];
      L_idx++;
    }
  }

  return result;
}


double phase_log_likelihood_slopes(
    NumericVector params,
    List data,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    List trans_re_info,
    List dyn_re_info_0,
    List dyn_re_info_1
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
  double sigma_0 = std::exp(log_sigma_0);
  double sigma_1 = std::exp(log_sigma_1);

  // Parse RE info for transition - use eager copies
  int n_trans_re_terms = Rcpp::as<int>(trans_re_info["n_terms"]);
  std::vector<int> trans_n_groups = Rcpp::as<std::vector<int>>(trans_re_info["n_groups_vec"]);
  std::vector<int> trans_n_coefs = Rcpp::as<std::vector<int>>(trans_re_info["n_coefs_vec"]);
  List trans_idx_list = Rcpp::as<List>(trans_re_info["idx_list"]);
  List trans_Z_list = Rcpp::as<List>(trans_re_info["Z_list"]);
  std::vector<int> trans_correlated = Rcpp::as<std::vector<int>>(trans_re_info["correlated_vec"]);

  // Extract transition RE values and covariance params
  std::vector<std::vector<std::vector<double>>> u_trans(n_trans_re_terms);

  for (int r = 0; r < n_trans_re_terms; r++) {
    int n_g = trans_n_groups[r];
    int n_c = trans_n_coefs[r];
    bool corr = trans_correlated[r];

    u_trans[r].resize(n_g);

    if (n_c == 1) {
      // Intercept only - simple case
      for (int g = 0; g < n_g; g++) {
        u_trans[r][g].resize(1);
        u_trans[r][g][0] = params[idx++];
      }
      idx++;  // skip log_sigma
    } else if (corr) {
      // Correlated: read raw eta values, then Cholesky params
      std::vector<std::vector<double>> eta_raw(n_g, std::vector<double>(n_c));
      for (int g = 0; g < n_g; g++) {
        for (int k = 0; k < n_c; k++) {
          eta_raw[g][k] = params[idx++];
        }
      }

      // Read Cholesky lower triangle
      int n_L = n_c * (n_c + 1) / 2;
      std::vector<double> L_vals(n_L);
      for (int l = 0; l < n_L; l++) {
        L_vals[l] = params[idx++];
      }

      // Transform each group's eta to actual RE
      for (int g = 0; g < n_g; g++) {
        u_trans[r][g] = apply_cholesky(eta_raw[g], L_vals, n_c);
      }
    } else {
      // Uncorrelated: read RE values and diagonal sigmas
      for (int g = 0; g < n_g; g++) {
        u_trans[r][g].resize(n_c);
        for (int k = 0; k < n_c; k++) {
          u_trans[r][g][k] = params[idx++];
        }
      }
      idx += n_c;  // skip log_sigmas
    }
  }

  // Parse RE info for dynamics 0 - use eager copies
  int n_dyn_re_terms_0 = Rcpp::as<int>(dyn_re_info_0["n_terms"]);
  std::vector<int> dyn_n_groups_0 = Rcpp::as<std::vector<int>>(dyn_re_info_0["n_groups_vec"]);
  std::vector<int> dyn_n_coefs_0 = Rcpp::as<std::vector<int>>(dyn_re_info_0["n_coefs_vec"]);
  List dyn_idx_list_0 = Rcpp::as<List>(dyn_re_info_0["idx_list"]);
  List dyn_Z_list_0 = Rcpp::as<List>(dyn_re_info_0["Z_list"]);
  std::vector<int> dyn_correlated_0 = Rcpp::as<std::vector<int>>(dyn_re_info_0["correlated_vec"]);

  std::vector<std::vector<std::vector<double>>> v_0(n_dyn_re_terms_0);

  for (int r = 0; r < n_dyn_re_terms_0; r++) {
    int n_g = dyn_n_groups_0[r];
    int n_c = dyn_n_coefs_0[r];
    bool corr = dyn_correlated_0[r];

    v_0[r].resize(n_g);

    if (n_c == 1) {
      for (int g = 0; g < n_g; g++) {
        v_0[r][g].resize(1);
        v_0[r][g][0] = params[idx++];
      }
      idx++;
    } else if (corr) {
      std::vector<std::vector<double>> eta_raw(n_g, std::vector<double>(n_c));
      for (int g = 0; g < n_g; g++) {
        for (int k = 0; k < n_c; k++) {
          eta_raw[g][k] = params[idx++];
        }
      }

      int n_L = n_c * (n_c + 1) / 2;
      std::vector<double> L_vals(n_L);
      for (int l = 0; l < n_L; l++) {
        L_vals[l] = params[idx++];
      }

      for (int g = 0; g < n_g; g++) {
        v_0[r][g] = apply_cholesky(eta_raw[g], L_vals, n_c);
      }
    } else {
      for (int g = 0; g < n_g; g++) {
        v_0[r][g].resize(n_c);
        for (int k = 0; k < n_c; k++) {
          v_0[r][g][k] = params[idx++];
        }
      }
      idx += n_c;
    }
  }

  // Parse RE info for dynamics 1 - use eager copies
  int n_dyn_re_terms_1 = Rcpp::as<int>(dyn_re_info_1["n_terms"]);
  std::vector<int> dyn_n_groups_1 = Rcpp::as<std::vector<int>>(dyn_re_info_1["n_groups_vec"]);
  std::vector<int> dyn_n_coefs_1 = Rcpp::as<std::vector<int>>(dyn_re_info_1["n_coefs_vec"]);
  List dyn_idx_list_1 = Rcpp::as<List>(dyn_re_info_1["idx_list"]);
  List dyn_Z_list_1 = Rcpp::as<List>(dyn_re_info_1["Z_list"]);
  std::vector<int> dyn_correlated_1 = Rcpp::as<std::vector<int>>(dyn_re_info_1["correlated_vec"]);

  std::vector<std::vector<std::vector<double>>> v_1(n_dyn_re_terms_1);

  for (int r = 0; r < n_dyn_re_terms_1; r++) {
    int n_g = dyn_n_groups_1[r];
    int n_c = dyn_n_coefs_1[r];
    bool corr = dyn_correlated_1[r];

    v_1[r].resize(n_g);

    if (n_c == 1) {
      for (int g = 0; g < n_g; g++) {
        v_1[r][g].resize(1);
        v_1[r][g][0] = params[idx++];
      }
      idx++;
    } else if (corr) {
      std::vector<std::vector<double>> eta_raw(n_g, std::vector<double>(n_c));
      for (int g = 0; g < n_g; g++) {
        for (int k = 0; k < n_c; k++) {
          eta_raw[g][k] = params[idx++];
        }
      }

      int n_L = n_c * (n_c + 1) / 2;
      std::vector<double> L_vals(n_L);
      for (int l = 0; l < n_L; l++) {
        L_vals[l] = params[idx++];
      }

      for (int g = 0; g < n_g; g++) {
        v_1[r][g] = apply_cholesky(eta_raw[g], L_vals, n_c);
      }
    } else {
      for (int g = 0; g < n_g; g++) {
        v_1[r][g].resize(n_c);
        for (int k = 0; k < n_c; k++) {
          v_1[r][g][k] = params[idx++];
        }
      }
      idx += n_c;
    }
  }

  // Observations and per-phase design matrices
  const phaser::TwoPhaseData d(data);
  const std::vector<double>& y = d.y;
  const int n_units = d.n_units;
  const std::vector<int>& unit_start = d.unit_start;
  const std::vector<int>& unit_end = d.unit_end;

  double total_ll = 0.0;

  for (int i = 0; i < n_units; i++) {

    int start = unit_start[i];
    int end = unit_end[i];

    double alpha_0, alpha_1;

    // Initial time point - compute mu_0
    double mu_0_init = 0.0;
    mu_0_init += phaser::linear_predictor(d.X_dyn_0, start, beta_0);

    // Add RE contributions for dynamics 0
    for (int r = 0; r < n_dyn_re_terms_0; r++) {
      IntegerVector re_idx = dyn_idx_list_0[r];
      int g = re_idx[start] - 1;
      int n_c = dyn_n_coefs_0[r];

      if (n_c == 1) {
        mu_0_init += v_0[r][g][0];
      } else {
        // Get Z row for this observation
        NumericMatrix Z_mat = dyn_Z_list_0[r];
        NumericVector Z_row = Z_mat(start, _);
        mu_0_init += compute_re_contribution(v_0[r][g], Z_row);
      }
    }

    alpha_0 = log_dnorm_s(y[start], mu_0_init, sigma_0);
    alpha_1 = -INFINITY;

    // Forward pass
    for (int t = start + 1; t <= end; t++) {

      // Transition linear predictor
      double eta_trans = 0.0;
      eta_trans += phaser::linear_predictor(d.X_trans, t, beta_trans);

      // Add transition RE
      for (int r = 0; r < n_trans_re_terms; r++) {
        IntegerVector re_idx = trans_idx_list[r];
        int g = re_idx[start] - 1;
        int n_c = trans_n_coefs[r];

        if (n_c == 1) {
          eta_trans += u_trans[r][g][0];
        } else {
          NumericMatrix Z_mat = trans_Z_list[r];
          NumericVector Z_row = Z_mat(t, _);
          eta_trans += compute_re_contribution(u_trans[r][g], Z_row);
        }
      }

      double p_trans = invlogit_s(eta_trans);

      // Dynamics linear predictors
      double mu_0 = 0.0, mu_1 = 0.0;
      mu_0 += phaser::linear_predictor(d.X_dyn_0, t, beta_0);
      mu_1 += phaser::linear_predictor(d.X_dyn_1, t, beta_1);

      // Add dynamics RE contributions
      for (int r = 0; r < n_dyn_re_terms_0; r++) {
        IntegerVector re_idx = dyn_idx_list_0[r];
        int g = re_idx[start] - 1;
        int n_c = dyn_n_coefs_0[r];

        if (n_c == 1) {
          mu_0 += v_0[r][g][0];
        } else {
          NumericMatrix Z_mat = dyn_Z_list_0[r];
          NumericVector Z_row = Z_mat(t, _);
          mu_0 += compute_re_contribution(v_0[r][g], Z_row);
        }
      }

      for (int r = 0; r < n_dyn_re_terms_1; r++) {
        IntegerVector re_idx = dyn_idx_list_1[r];
        int g = re_idx[start] - 1;
        int n_c = dyn_n_coefs_1[r];

        if (n_c == 1) {
          mu_1 += v_1[r][g][0];
        } else {
          NumericMatrix Z_mat = dyn_Z_list_1[r];
          NumericVector Z_row = Z_mat(t, _);
          mu_1 += compute_re_contribution(v_1[r][g], Z_row);
        }
      }

      double log_emit_0 = log_dnorm_s(y[t], mu_0, sigma_0);
      double log_emit_1 = log_dnorm_s(y[t], mu_1, sigma_1);

      double new_alpha_0 = alpha_0 + std::log(1.0 - p_trans) + log_emit_0;
      double from_0 = alpha_0 + std::log(p_trans) + log_emit_1;
      double from_1 = alpha_1 + log_emit_1;
      double new_alpha_1 = log_sum_exp_s(from_0, from_1);

      alpha_0 = new_alpha_0;
      alpha_1 = new_alpha_1;
    }

    total_ll += log_sum_exp_s(alpha_0, alpha_1);
  }

  return total_ll;
}


// Log prior of the random-effect terms of one model component, in the
// order the likelihood reads them. Intercept-only terms: effects
// ~ Normal(0, sigma_re), log(sigma_re) under the random-effect SD prior.
// Correlated terms: standardised effects ~ Normal(0, 1), an LKJ(eta) term
// on the Cholesky factor with eta from the correlation prior, and the
// random-effect SD prior reflected onto each diagonal element. Uncorrelated
// terms: effects ~ Normal(0, sigma_k), log(sigma_k) under the random-effect
// SD prior.
double log_prior_re_terms(const NumericVector& params, int& idx,
                          const List& re_info, const prior::PriorSet& pr) {
  int n_terms = Rcpp::as<int>(re_info["n_terms"]);
  std::vector<int> n_groups = Rcpp::as<std::vector<int>>(re_info["n_groups_vec"]);
  std::vector<int> n_coefs = Rcpp::as<std::vector<int>>(re_info["n_coefs_vec"]);
  std::vector<int> correlated = Rcpp::as<std::vector<int>>(re_info["correlated_vec"]);

  double log_prior = 0.0;
  for (int r = 0; r < n_terms; r++) {
    int n_g = n_groups[r];
    int n_c = n_coefs[r];

    if (n_c == 1) {
      log_prior += prior::log_prior_re_intercepts(params, idx, n_g, pr);
    } else if (correlated[r]) {
      for (int g = 0; g < n_g; g++) {
        for (int k = 0; k < n_c; k++) {
          log_prior += prior::log_normal(params[idx++], 0.0, 1.0);
        }
      }
      int n_L = n_c * (n_c + 1) / 2;
      std::vector<double> L_vals(n_L);
      for (int l = 0; l < n_L; l++) {
        L_vals[l] = params[idx++];
      }
      log_prior += log_lkj_cholesky(L_vals, n_c, pr.cor_eta());
      int l_idx = 0;
      for (int j = 0; j < n_c; j++) {
        log_prior += pr.re_scale(L_vals[l_idx]);
        l_idx += n_c - j;
      }
    } else {
      for (int g = 0; g < n_g; g++) {
        for (int k = 0; k < n_c; k++) {
          double log_sigma_k = params[idx + n_g * n_c + k];
          log_prior += prior::log_normal(params[idx++], 0.0, std::exp(log_sigma_k));
        }
      }
      for (int k = 0; k < n_c; k++) {
        log_prior += pr.re_sd(params[idx++]);
      }
    }
  }
  return log_prior;
}


double phase_log_prior_slopes(
    NumericVector params,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    List trans_re_info,
    List dyn_re_info_0,
    List dyn_re_info_1,
    const prior::PriorSet& pr
) {
  int idx = 0;
  double log_prior = prior::log_prior_fixed(
    params, idx, n_trans_coef + n_dyn_coef_0 + n_dyn_coef_1, 2, pr);
  log_prior += log_prior_re_terms(params, idx, trans_re_info, pr);
  log_prior += log_prior_re_terms(params, idx, dyn_re_info_0, pr);
  log_prior += log_prior_re_terms(params, idx, dyn_re_info_1, pr);
  return log_prior;
}


double phase_log_posterior_slopes(
    NumericVector params,
    List data,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    List trans_re_info,
    List dyn_re_info_0,
    List dyn_re_info_1,
    const prior::PriorSet& pr
) {
  return phase_log_likelihood_slopes(params, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1,
                                      trans_re_info, dyn_re_info_0, dyn_re_info_1) +
         phase_log_prior_slopes(params, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1,
                                 trans_re_info, dyn_re_info_0, dyn_re_info_1, pr);
}
