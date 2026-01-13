#include <Rcpp.h>
#include <cmath>
#include <vector>

using namespace Rcpp;

// Log-sum-exp for numerical stability
inline double log_sum_exp_m(double a, double b) {
  if (a > b) {
    return a + std::log1p(std::exp(b - a));
  } else {
    return b + std::log1p(std::exp(a - b));
  }
}

inline double invlogit_m(double x) {
  return 1.0 / (1.0 + std::exp(-x));
}

inline double log_dnorm_m(double y, double mu, double sigma) {
  double z = (y - mu) / sigma;
  return -0.5 * std::log(2.0 * M_PI) - std::log(sigma) - 0.5 * z * z;
}

// [[Rcpp::export]]
double phase_log_likelihood_multi_re(
    NumericVector params,
    List data,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    List trans_re_info,   // List with n_terms, n_groups_vec, idx_list
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

  // Parse multi-RE info for transition - use eager copies
  int n_trans_re_terms = Rcpp::as<int>(trans_re_info["n_terms"]);
  std::vector<int> trans_n_groups = Rcpp::as<std::vector<int>>(trans_re_info["n_groups_vec"]);
  List trans_idx_list = Rcpp::as<List>(trans_re_info["idx_list"]);

  // Extract transition RE values
  std::vector<NumericVector> u_trans(n_trans_re_terms);
  for (int r = 0; r < n_trans_re_terms; r++) {
    int n_g = trans_n_groups[r];
    NumericVector u_r(n_g);
    for (int g = 0; g < n_g; g++) {
      u_r[g] = params[idx++];
    }
    idx++;  // skip log_sigma for this RE term
    u_trans[r] = u_r;
  }

  // Parse multi-RE info for dynamics 0 - use eager copies
  int n_dyn_re_terms_0 = Rcpp::as<int>(dyn_re_info_0["n_terms"]);
  std::vector<int> dyn_n_groups_0 = Rcpp::as<std::vector<int>>(dyn_re_info_0["n_groups_vec"]);
  List dyn_idx_list_0 = Rcpp::as<List>(dyn_re_info_0["idx_list"]);

  std::vector<NumericVector> v_0(n_dyn_re_terms_0);
  for (int r = 0; r < n_dyn_re_terms_0; r++) {
    int n_g = dyn_n_groups_0[r];
    NumericVector v_r(n_g);
    for (int g = 0; g < n_g; g++) {
      v_r[g] = params[idx++];
    }
    idx++;  // skip log_sigma
    v_0[r] = v_r;
  }

  // Parse multi-RE info for dynamics 1 - use eager copies
  int n_dyn_re_terms_1 = Rcpp::as<int>(dyn_re_info_1["n_terms"]);
  std::vector<int> dyn_n_groups_1 = Rcpp::as<std::vector<int>>(dyn_re_info_1["n_groups_vec"]);
  List dyn_idx_list_1 = Rcpp::as<List>(dyn_re_info_1["idx_list"]);

  std::vector<NumericVector> v_1(n_dyn_re_terms_1);
  for (int r = 0; r < n_dyn_re_terms_1; r++) {
    int n_g = dyn_n_groups_1[r];
    NumericVector v_r(n_g);
    for (int g = 0; g < n_g; g++) {
      v_r[g] = params[idx++];
    }
    idx++;  // skip log_sigma
    v_1[r] = v_r;
  }

  // Data - use eager deep copies to prevent R GC issues
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

    double alpha_0, alpha_1;

    // Initial time point - compute mu_0
    double mu_0_init = 0.0;
    for (int j = 0; j < n_dyn_coef_0; j++) {
      mu_0_init += beta_0[j] * X_dyn(start, j);
    }
    // Add all dynamics RE terms for phase 0
    for (int r = 0; r < n_dyn_re_terms_0; r++) {
      IntegerVector re_idx = dyn_idx_list_0[r];
      int g = re_idx[start] - 1;  // Convert to 0-indexed
      mu_0_init += v_0[r][g];
    }

    alpha_0 = log_dnorm_m(y[start], mu_0_init, sigma_0);
    alpha_1 = -INFINITY;

    // Forward pass
    for (int t = start + 1; t <= end; t++) {

      // Transition linear predictor
      double eta_trans = 0.0;
      for (int j = 0; j < n_trans_coef; j++) {
        eta_trans += beta_trans[j] * X_trans(t, j);
      }
      // Add all transition RE terms
      for (int r = 0; r < n_trans_re_terms; r++) {
        IntegerVector re_idx = trans_idx_list[r];
        int g = re_idx[start] - 1;  // Use unit's first obs for group membership
        eta_trans += u_trans[r][g];
      }
      double p_trans = invlogit_m(eta_trans);

      // Dynamics linear predictors
      double mu_0 = 0.0, mu_1 = 0.0;
      for (int j = 0; j < n_dyn_coef_0; j++) {
        mu_0 += beta_0[j] * X_dyn(t, j);
      }
      for (int j = 0; j < n_dyn_coef_1; j++) {
        mu_1 += beta_1[j] * X_dyn(t, j);
      }

      // Add all dynamics RE terms
      for (int r = 0; r < n_dyn_re_terms_0; r++) {
        IntegerVector re_idx = dyn_idx_list_0[r];
        int g = re_idx[start] - 1;
        mu_0 += v_0[r][g];
      }
      for (int r = 0; r < n_dyn_re_terms_1; r++) {
        IntegerVector re_idx = dyn_idx_list_1[r];
        int g = re_idx[start] - 1;
        mu_1 += v_1[r][g];
      }

      double log_emit_0 = log_dnorm_m(y[t], mu_0, sigma_0);
      double log_emit_1 = log_dnorm_m(y[t], mu_1, sigma_1);

      double new_alpha_0 = alpha_0 + std::log(1.0 - p_trans) + log_emit_0;
      double from_0 = alpha_0 + std::log(p_trans) + log_emit_1;
      double from_1 = alpha_1 + log_emit_1;
      double new_alpha_1 = log_sum_exp_m(from_0, from_1);

      alpha_0 = new_alpha_0;
      alpha_1 = new_alpha_1;
    }

    total_ll += log_sum_exp_m(alpha_0, alpha_1);
  }

  return total_ll;
}


// [[Rcpp::export]]
double phase_log_prior_multi_re(
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

  for (int r = 0; r < n_trans_re_terms; r++) {
    int n_g = trans_n_groups[r];
    double log_sigma_re = params[idx + n_g];  // sigma is after the RE values
    double sigma_re = std::exp(log_sigma_re);
    double sigma_re_sq = sigma_re * sigma_re;

    // RE values: Normal(0, sigma_re)
    for (int g = 0; g < n_g; g++) {
      double u = params[idx++];
      log_prior += -0.5 * std::log(2.0 * M_PI) - std::log(sigma_re) - 0.5 * u * u / sigma_re_sq;
    }
    // Hyperprior on sigma_re: Half-Normal(0, 0.5) on log scale
    log_prior += -0.5 * std::log(2.0 * M_PI * 0.25) - 0.5 * log_sigma_re * log_sigma_re / 0.25;
    idx++;  // move past log_sigma_re
  }

  // Dynamics 0 RE priors
  int n_dyn_re_terms_0 = dyn_re_info_0["n_terms"];
  IntegerVector dyn_n_groups_0 = dyn_re_info_0["n_groups_vec"];

  for (int r = 0; r < n_dyn_re_terms_0; r++) {
    int n_g = dyn_n_groups_0[r];
    double log_sigma_re = params[idx + n_g];
    double sigma_re = std::exp(log_sigma_re);
    double sigma_re_sq = sigma_re * sigma_re;

    for (int g = 0; g < n_g; g++) {
      double v = params[idx++];
      log_prior += -0.5 * std::log(2.0 * M_PI) - std::log(sigma_re) - 0.5 * v * v / sigma_re_sq;
    }
    log_prior += -0.5 * std::log(2.0 * M_PI * 0.25) - 0.5 * log_sigma_re * log_sigma_re / 0.25;
    idx++;
  }

  // Dynamics 1 RE priors
  int n_dyn_re_terms_1 = dyn_re_info_1["n_terms"];
  IntegerVector dyn_n_groups_1 = dyn_re_info_1["n_groups_vec"];

  for (int r = 0; r < n_dyn_re_terms_1; r++) {
    int n_g = dyn_n_groups_1[r];
    double log_sigma_re = params[idx + n_g];
    double sigma_re = std::exp(log_sigma_re);
    double sigma_re_sq = sigma_re * sigma_re;

    for (int g = 0; g < n_g; g++) {
      double v = params[idx++];
      log_prior += -0.5 * std::log(2.0 * M_PI) - std::log(sigma_re) - 0.5 * v * v / sigma_re_sq;
    }
    log_prior += -0.5 * std::log(2.0 * M_PI * 0.25) - 0.5 * log_sigma_re * log_sigma_re / 0.25;
    idx++;
  }

  return log_prior;
}


// [[Rcpp::export]]
double phase_log_posterior_multi_re(
    NumericVector params,
    List data,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    List trans_re_info,
    List dyn_re_info_0,
    List dyn_re_info_1
) {
  return phase_log_likelihood_multi_re(params, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1,
                                        trans_re_info, dyn_re_info_0, dyn_re_info_1) +
         phase_log_prior_multi_re(params, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1,
                                   trans_re_info, dyn_re_info_0, dyn_re_info_1);
}
