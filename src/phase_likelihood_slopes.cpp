#include <Rcpp.h>
#include <cmath>
#include <vector>

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


// [[Rcpp::export]]
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

  // Parse RE info for transition
  int n_trans_re_terms = trans_re_info["n_terms"];
  IntegerVector trans_n_groups = trans_re_info["n_groups_vec"];
  IntegerVector trans_n_coefs = trans_re_info["n_coefs_vec"];
  List trans_idx_list = trans_re_info["idx_list"];
  List trans_Z_list = trans_re_info["Z_list"];
  LogicalVector trans_correlated = trans_re_info["correlated_vec"];

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

  // Parse RE info for dynamics 0
  int n_dyn_re_terms_0 = dyn_re_info_0["n_terms"];
  IntegerVector dyn_n_groups_0 = dyn_re_info_0["n_groups_vec"];
  IntegerVector dyn_n_coefs_0 = dyn_re_info_0["n_coefs_vec"];
  List dyn_idx_list_0 = dyn_re_info_0["idx_list"];
  List dyn_Z_list_0 = dyn_re_info_0["Z_list"];
  LogicalVector dyn_correlated_0 = dyn_re_info_0["correlated_vec"];

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

  // Parse RE info for dynamics 1
  int n_dyn_re_terms_1 = dyn_re_info_1["n_terms"];
  IntegerVector dyn_n_groups_1 = dyn_re_info_1["n_groups_vec"];
  IntegerVector dyn_n_coefs_1 = dyn_re_info_1["n_coefs_vec"];
  List dyn_idx_list_1 = dyn_re_info_1["idx_list"];
  List dyn_Z_list_1 = dyn_re_info_1["Z_list"];
  LogicalVector dyn_correlated_1 = dyn_re_info_1["correlated_vec"];

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

  // Data
  NumericVector y = data["y"];
  NumericMatrix X_trans = data["X_trans"];
  NumericMatrix X_dyn = data["X_dyn"];

  int n_units = data["n_units"];
  IntegerVector unit_start = data["unit_start"];
  IntegerVector unit_end = data["unit_end"];

  double total_ll = 0.0;

  for (int i = 0; i < n_units; i++) {

    int start = unit_start[i];
    int end = unit_end[i];

    double alpha_0, alpha_1;

    // Initial time point - compute mu_0
    double mu_0_init = 0.0;
    for (int j = 0; j < n_dyn_coef_0; j++) {
      mu_0_init += beta_0[j] * X_dyn(start, j);
    }

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
      for (int j = 0; j < n_trans_coef; j++) {
        eta_trans += beta_trans[j] * X_trans(t, j);
      }

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
      for (int j = 0; j < n_dyn_coef_0; j++) {
        mu_0 += beta_0[j] * X_dyn(t, j);
      }
      for (int j = 0; j < n_dyn_coef_1; j++) {
        mu_1 += beta_1[j] * X_dyn(t, j);
      }

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


// [[Rcpp::export]]
double phase_log_prior_slopes(
    NumericVector params,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    List trans_re_info,
    List dyn_re_info_0,
    List dyn_re_info_1
) {

  double log_prior = 0.0;
  int idx = 0;

  // Fixed effects priors: Normal(0, 2.5)
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

  // Residual SD priors: Normal(0, 1) on log scale
  double log_sigma_0 = params[idx++];
  double log_sigma_1 = params[idx++];
  log_prior += -0.5 * std::log(2.0 * M_PI) - 0.5 * log_sigma_0 * log_sigma_0;
  log_prior += -0.5 * std::log(2.0 * M_PI) - 0.5 * log_sigma_1 * log_sigma_1;

  // Transition RE priors
  int n_trans_re_terms = trans_re_info["n_terms"];
  IntegerVector trans_n_groups = trans_re_info["n_groups_vec"];
  IntegerVector trans_n_coefs = trans_re_info["n_coefs_vec"];
  LogicalVector trans_correlated = trans_re_info["correlated_vec"];

  for (int r = 0; r < n_trans_re_terms; r++) {
    int n_g = trans_n_groups[r];
    int n_c = trans_n_coefs[r];
    bool corr = trans_correlated[r];

    if (n_c == 1) {
      // Simple intercept: Normal(0, sigma_re)
      double log_sigma_re = params[idx + n_g];
      double sigma_re = std::exp(log_sigma_re);
      double sigma_re_sq = sigma_re * sigma_re;

      for (int g = 0; g < n_g; g++) {
        double u = params[idx++];
        log_prior += -0.5 * std::log(2.0 * M_PI) - std::log(sigma_re) - 0.5 * u * u / sigma_re_sq;
      }
      // Half-Normal(0, 0.5) on log scale
      log_prior += -0.5 * std::log(2.0 * M_PI * 0.25) - 0.5 * log_sigma_re * log_sigma_re / 0.25;
      idx++;
    } else if (corr) {
      // Correlated: eta ~ N(0, I), L ~ LKJ
      for (int g = 0; g < n_g; g++) {
        for (int k = 0; k < n_c; k++) {
          double eta = params[idx++];
          log_prior += -0.5 * std::log(2.0 * M_PI) - 0.5 * eta * eta;
        }
      }

      // Cholesky prior (LKJ with eta=1 and half-normal on diagonal)
      int n_L = n_c * (n_c + 1) / 2;
      std::vector<double> L_vals(n_L);
      for (int l = 0; l < n_L; l++) {
        L_vals[l] = params[idx++];
      }

      // LKJ prior on correlation structure
      log_prior += log_lkj_cholesky(L_vals, n_c, 1.0);

      // Half-normal prior on diagonal elements (standard deviations)
      int l_idx = 0;
      for (int j = 0; j < n_c; j++) {
        for (int i = j; i < n_c; i++) {
          if (i == j) {
            // Diagonal: half-normal(0, 1) on positive scale
            double L_jj = L_vals[l_idx];
            log_prior += -0.5 * std::log(2.0 * M_PI) - 0.5 * L_jj * L_jj;
          }
          l_idx++;
        }
      }
    } else {
      // Uncorrelated: each coefficient has its own sigma
      for (int g = 0; g < n_g; g++) {
        for (int k = 0; k < n_c; k++) {
          double log_sigma_k = params[idx + n_g * n_c + k];
          double sigma_k = std::exp(log_sigma_k);
          double u = params[idx++];
          log_prior += -0.5 * std::log(2.0 * M_PI) - std::log(sigma_k) - 0.5 * u * u / (sigma_k * sigma_k);
        }
      }
      // Priors on log_sigmas
      for (int k = 0; k < n_c; k++) {
        double log_sigma_k = params[idx++];
        log_prior += -0.5 * std::log(2.0 * M_PI * 0.25) - 0.5 * log_sigma_k * log_sigma_k / 0.25;
      }
    }
  }

  // Dynamics 0 RE priors (same logic)
  int n_dyn_re_terms_0 = dyn_re_info_0["n_terms"];
  IntegerVector dyn_n_groups_0 = dyn_re_info_0["n_groups_vec"];
  IntegerVector dyn_n_coefs_0 = dyn_re_info_0["n_coefs_vec"];
  LogicalVector dyn_correlated_0 = dyn_re_info_0["correlated_vec"];

  for (int r = 0; r < n_dyn_re_terms_0; r++) {
    int n_g = dyn_n_groups_0[r];
    int n_c = dyn_n_coefs_0[r];
    bool corr = dyn_correlated_0[r];

    if (n_c == 1) {
      double log_sigma_re = params[idx + n_g];
      double sigma_re = std::exp(log_sigma_re);
      double sigma_re_sq = sigma_re * sigma_re;

      for (int g = 0; g < n_g; g++) {
        double v = params[idx++];
        log_prior += -0.5 * std::log(2.0 * M_PI) - std::log(sigma_re) - 0.5 * v * v / sigma_re_sq;
      }
      log_prior += -0.5 * std::log(2.0 * M_PI * 0.25) - 0.5 * log_sigma_re * log_sigma_re / 0.25;
      idx++;
    } else if (corr) {
      for (int g = 0; g < n_g; g++) {
        for (int k = 0; k < n_c; k++) {
          double eta = params[idx++];
          log_prior += -0.5 * std::log(2.0 * M_PI) - 0.5 * eta * eta;
        }
      }

      int n_L = n_c * (n_c + 1) / 2;
      std::vector<double> L_vals(n_L);
      for (int l = 0; l < n_L; l++) {
        L_vals[l] = params[idx++];
      }

      log_prior += log_lkj_cholesky(L_vals, n_c, 1.0);

      int l_idx = 0;
      for (int j = 0; j < n_c; j++) {
        for (int i = j; i < n_c; i++) {
          if (i == j) {
            double L_jj = L_vals[l_idx];
            log_prior += -0.5 * std::log(2.0 * M_PI) - 0.5 * L_jj * L_jj;
          }
          l_idx++;
        }
      }
    } else {
      for (int g = 0; g < n_g; g++) {
        for (int k = 0; k < n_c; k++) {
          double log_sigma_k = params[idx + n_g * n_c + k];
          double sigma_k = std::exp(log_sigma_k);
          double v = params[idx++];
          log_prior += -0.5 * std::log(2.0 * M_PI) - std::log(sigma_k) - 0.5 * v * v / (sigma_k * sigma_k);
        }
      }
      for (int k = 0; k < n_c; k++) {
        double log_sigma_k = params[idx++];
        log_prior += -0.5 * std::log(2.0 * M_PI * 0.25) - 0.5 * log_sigma_k * log_sigma_k / 0.25;
      }
    }
  }

  // Dynamics 1 RE priors
  int n_dyn_re_terms_1 = dyn_re_info_1["n_terms"];
  IntegerVector dyn_n_groups_1 = dyn_re_info_1["n_groups_vec"];
  IntegerVector dyn_n_coefs_1 = dyn_re_info_1["n_coefs_vec"];
  LogicalVector dyn_correlated_1 = dyn_re_info_1["correlated_vec"];

  for (int r = 0; r < n_dyn_re_terms_1; r++) {
    int n_g = dyn_n_groups_1[r];
    int n_c = dyn_n_coefs_1[r];
    bool corr = dyn_correlated_1[r];

    if (n_c == 1) {
      double log_sigma_re = params[idx + n_g];
      double sigma_re = std::exp(log_sigma_re);
      double sigma_re_sq = sigma_re * sigma_re;

      for (int g = 0; g < n_g; g++) {
        double v = params[idx++];
        log_prior += -0.5 * std::log(2.0 * M_PI) - std::log(sigma_re) - 0.5 * v * v / sigma_re_sq;
      }
      log_prior += -0.5 * std::log(2.0 * M_PI * 0.25) - 0.5 * log_sigma_re * log_sigma_re / 0.25;
      idx++;
    } else if (corr) {
      for (int g = 0; g < n_g; g++) {
        for (int k = 0; k < n_c; k++) {
          double eta = params[idx++];
          log_prior += -0.5 * std::log(2.0 * M_PI) - 0.5 * eta * eta;
        }
      }

      int n_L = n_c * (n_c + 1) / 2;
      std::vector<double> L_vals(n_L);
      for (int l = 0; l < n_L; l++) {
        L_vals[l] = params[idx++];
      }

      log_prior += log_lkj_cholesky(L_vals, n_c, 1.0);

      int l_idx = 0;
      for (int j = 0; j < n_c; j++) {
        for (int i = j; i < n_c; i++) {
          if (i == j) {
            double L_jj = L_vals[l_idx];
            log_prior += -0.5 * std::log(2.0 * M_PI) - 0.5 * L_jj * L_jj;
          }
          l_idx++;
        }
      }
    } else {
      for (int g = 0; g < n_g; g++) {
        for (int k = 0; k < n_c; k++) {
          double log_sigma_k = params[idx + n_g * n_c + k];
          double sigma_k = std::exp(log_sigma_k);
          double v = params[idx++];
          log_prior += -0.5 * std::log(2.0 * M_PI) - std::log(sigma_k) - 0.5 * v * v / (sigma_k * sigma_k);
        }
      }
      for (int k = 0; k < n_c; k++) {
        double log_sigma_k = params[idx++];
        log_prior += -0.5 * std::log(2.0 * M_PI * 0.25) - 0.5 * log_sigma_k * log_sigma_k / 0.25;
      }
    }
  }

  return log_prior;
}


// [[Rcpp::export]]
double phase_log_posterior_slopes(
    NumericVector params,
    List data,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    List trans_re_info,
    List dyn_re_info_0,
    List dyn_re_info_1
) {
  return phase_log_likelihood_slopes(params, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1,
                                      trans_re_info, dyn_re_info_0, dyn_re_info_1) +
         phase_log_prior_slopes(params, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1,
                                 trans_re_info, dyn_re_info_0, dyn_re_info_1);
}
