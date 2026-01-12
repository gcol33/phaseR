#include <Rcpp.h>
#include <cmath>
#include <random>

using namespace Rcpp;

// Forward declaration
double phase_log_posterior_multi_re(
    NumericVector params, List data,
    int n_trans_coef, int n_dyn_coef_0, int n_dyn_coef_1,
    List trans_re_info, List dyn_re_info_0, List dyn_re_info_1);

// Numerical gradient for multi-RE model
NumericVector numerical_gradient_multi_re(
    NumericVector params,
    List data,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    List trans_re_info,
    List dyn_re_info_0,
    List dyn_re_info_1,
    double eps = 1e-6
) {
  int n = params.size();
  NumericVector grad(n);

  for (int i = 0; i < n; i++) {
    NumericVector params_plus = clone(params);
    NumericVector params_minus = clone(params);
    params_plus[i] += eps;
    params_minus[i] -= eps;

    double f_plus = phase_log_posterior_multi_re(params_plus, data, n_trans_coef,
                                                  n_dyn_coef_0, n_dyn_coef_1,
                                                  trans_re_info, dyn_re_info_0, dyn_re_info_1);
    double f_minus = phase_log_posterior_multi_re(params_minus, data, n_trans_coef,
                                                   n_dyn_coef_0, n_dyn_coef_1,
                                                   trans_re_info, dyn_re_info_0, dyn_re_info_1);
    grad[i] = (f_plus - f_minus) / (2.0 * eps);
  }

  return grad;
}

// Leapfrog step
void leapfrog_step_multi_re(
    NumericVector& q,
    NumericVector& p,
    List data,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    List trans_re_info,
    List dyn_re_info_0,
    List dyn_re_info_1,
    double step_size
) {
  int n = q.size();

  NumericVector grad = numerical_gradient_multi_re(q, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1,
                                                    trans_re_info, dyn_re_info_0, dyn_re_info_1);
  for (int i = 0; i < n; i++) {
    p[i] += 0.5 * step_size * grad[i];
  }

  for (int i = 0; i < n; i++) {
    q[i] += step_size * p[i];
  }

  grad = numerical_gradient_multi_re(q, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1,
                                      trans_re_info, dyn_re_info_0, dyn_re_info_1);
  for (int i = 0; i < n; i++) {
    p[i] += 0.5 * step_size * grad[i];
  }
}

// [[Rcpp::export]]
List phaseR_nuts_sampler_multi_re(
    List data,
    NumericVector init,
    int n_iter,
    int n_warmup,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    List trans_re_info,
    List dyn_re_info_0,
    List dyn_re_info_1,
    double target_accept = 0.8,
    int max_treedepth = 10
) {
  int n_params = init.size();
  int n_samples = n_iter - n_warmup;

  NumericMatrix draws(n_samples, n_params);
  NumericVector current = clone(init);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> normal(0, 1);
  std::uniform_real_distribution<> uniform(0, 1);

  double step_size = 0.1;
  double log_step_size = std::log(step_size);
  double H_bar = 0.0;
  double gamma = 0.05;
  double t0 = 10.0;
  double kappa = 0.75;
  double mu = std::log(10.0 * step_size);

  int n_divergent = 0;
  double total_accept = 0.0;
  int sample_idx = 0;

  for (int iter = 0; iter < n_iter; iter++) {

    NumericVector p(n_params);
    for (int i = 0; i < n_params; i++) {
      p[i] = normal(gen);
    }

    double current_H = -phase_log_posterior_multi_re(current, data, n_trans_coef,
                                                      n_dyn_coef_0, n_dyn_coef_1,
                                                      trans_re_info, dyn_re_info_0, dyn_re_info_1);
    for (int i = 0; i < n_params; i++) {
      current_H += 0.5 * p[i] * p[i];
    }

    NumericVector q = clone(current);
    NumericVector p_current = clone(p);

    int n_leapfrog = std::max(1, std::min(1024, (int)std::ceil(1.0 / step_size)));
    if (n_leapfrog > 100) n_leapfrog = 100;

    bool divergent = false;
    for (int l = 0; l < n_leapfrog; l++) {
      leapfrog_step_multi_re(q, p, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1,
                              trans_re_info, dyn_re_info_0, dyn_re_info_1, step_size);

      double proposed_H = -phase_log_posterior_multi_re(q, data, n_trans_coef,
                                                         n_dyn_coef_0, n_dyn_coef_1,
                                                         trans_re_info, dyn_re_info_0, dyn_re_info_1);
      for (int i = 0; i < n_params; i++) {
        proposed_H += 0.5 * p[i] * p[i];
      }

      if (std::abs(proposed_H - current_H) > 1000.0) {
        divergent = true;
        break;
      }
    }

    if (divergent) {
      n_divergent++;
      if (iter >= n_warmup) {
        for (int i = 0; i < n_params; i++) {
          draws(sample_idx, i) = current[i];
        }
        sample_idx++;
      }
      continue;
    }

    double proposed_H = -phase_log_posterior_multi_re(q, data, n_trans_coef,
                                                       n_dyn_coef_0, n_dyn_coef_1,
                                                       trans_re_info, dyn_re_info_0, dyn_re_info_1);
    for (int i = 0; i < n_params; i++) {
      proposed_H += 0.5 * p[i] * p[i];
    }

    double log_accept = current_H - proposed_H;
    double accept_prob = std::min(1.0, std::exp(log_accept));

    if (uniform(gen) < accept_prob) {
      current = q;
    }

    total_accept += accept_prob;

    // Dual averaging during warmup
    if (iter < n_warmup) {
      double w = 1.0 / (iter + t0);
      H_bar = (1.0 - w) * H_bar + w * (target_accept - accept_prob);
      log_step_size = mu - std::sqrt((double)iter) / gamma * H_bar;
      step_size = std::exp(log_step_size);
      step_size = std::max(1e-4, std::min(1.0, step_size));
    }

    if (iter >= n_warmup) {
      for (int i = 0; i < n_params; i++) {
        draws(sample_idx, i) = current[i];
      }
      sample_idx++;
    }
  }

  return List::create(
    Named("draws") = draws,
    Named("n_divergent") = n_divergent,
    Named("accept_prob") = total_accept / n_iter,
    Named("step_size") = step_size
  );
}
