#include <Rcpp.h>
#include <cmath>
#include <random>

using namespace Rcpp;

// Forward declarations
double phase_log_posterior_slopes(
    NumericVector params,
    List data,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    List trans_re_info,
    List dyn_re_info_0,
    List dyn_re_info_1
);


// Numerical gradient using finite differences
NumericVector gradient_slopes(
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

  double f0 = phase_log_posterior_slopes(params, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1,
                                          trans_re_info, dyn_re_info_0, dyn_re_info_1);

  for (int i = 0; i < n; i++) {
    NumericVector params_plus = clone(params);
    params_plus[i] += eps;
    double f_plus = phase_log_posterior_slopes(params_plus, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1,
                                                trans_re_info, dyn_re_info_0, dyn_re_info_1);
    grad[i] = (f_plus - f0) / eps;
  }

  return grad;
}


// Single leapfrog step
void leapfrog_step_slopes(
    NumericVector& q,
    NumericVector& p,
    double step_size,
    List data,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    List trans_re_info,
    List dyn_re_info_0,
    List dyn_re_info_1
) {
  int n = q.size();

  // Half step for momentum
  NumericVector grad = gradient_slopes(q, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1,
                                         trans_re_info, dyn_re_info_0, dyn_re_info_1);
  for (int i = 0; i < n; i++) {
    p[i] += 0.5 * step_size * grad[i];
  }

  // Full step for position
  for (int i = 0; i < n; i++) {
    q[i] += step_size * p[i];
  }

  // Half step for momentum
  grad = gradient_slopes(q, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1,
                          trans_re_info, dyn_re_info_0, dyn_re_info_1);
  for (int i = 0; i < n; i++) {
    p[i] += 0.5 * step_size * grad[i];
  }
}


// Compute Hamiltonian
double hamiltonian_slopes(
    NumericVector& q,
    NumericVector& p,
    List data,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    List trans_re_info,
    List dyn_re_info_0,
    List dyn_re_info_1
) {
  double U = -phase_log_posterior_slopes(q, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1,
                                          trans_re_info, dyn_re_info_0, dyn_re_info_1);
  double K = 0.0;
  for (int i = 0; i < p.size(); i++) {
    K += 0.5 * p[i] * p[i];
  }
  return U + K;
}


// [[Rcpp::export]]
List phaseR_nuts_sampler_slopes(
    List data,
    NumericVector init,
    int n_iter,
    int n_warmup,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    List trans_re_info,
    List dyn_re_info_0,
    List dyn_re_info_1
) {
  int n_params = init.size();
  int n_samples = n_iter - n_warmup;

  NumericMatrix draws(n_samples, n_params);

  // Current state
  NumericVector q = clone(init);

  // Adaptation parameters
  double step_size = 0.1;
  double target_accept = 0.8;
  double adapt_rate = 0.05;

  // Leapfrog settings
  int n_leapfrog = 10;

  // Random number generator
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> normal(0.0, 1.0);
  std::uniform_real_distribution<> uniform(0.0, 1.0);

  int n_accept = 0;
  int n_divergent = 0;

  for (int iter = 0; iter < n_iter; iter++) {

    // Sample momentum
    NumericVector p(n_params);
    for (int i = 0; i < n_params; i++) {
      p[i] = normal(gen);
    }

    // Store initial state
    NumericVector q_init = clone(q);
    NumericVector p_init = clone(p);

    // Initial Hamiltonian
    double H0 = hamiltonian_slopes(q, p, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1,
                                    trans_re_info, dyn_re_info_0, dyn_re_info_1);

    // Leapfrog integration
    for (int l = 0; l < n_leapfrog; l++) {
      leapfrog_step_slopes(q, p, step_size, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1,
                            trans_re_info, dyn_re_info_0, dyn_re_info_1);
    }

    // Negate momentum
    for (int i = 0; i < n_params; i++) {
      p[i] = -p[i];
    }

    // Final Hamiltonian
    double H1 = hamiltonian_slopes(q, p, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1,
                                    trans_re_info, dyn_re_info_0, dyn_re_info_1);

    // Check for divergence
    double delta_H = H1 - H0;
    bool divergent = std::isnan(delta_H) || std::isinf(delta_H) || delta_H > 1000.0;

    if (divergent) {
      n_divergent++;
      q = clone(q_init);
    } else {
      // Metropolis acceptance
      double log_accept_prob = -delta_H;
      double accept_prob = std::min(1.0, std::exp(log_accept_prob));

      if (uniform(gen) < accept_prob) {
        n_accept++;
        // q is already updated
      } else {
        q = clone(q_init);
      }

      // Adapt step size during warmup
      if (iter < n_warmup) {
        if (accept_prob > target_accept) {
          step_size *= (1.0 + adapt_rate);
        } else {
          step_size *= (1.0 - adapt_rate);
        }
        step_size = std::max(1e-6, std::min(1.0, step_size));
      }
    }

    // Store sample after warmup
    if (iter >= n_warmup) {
      int sample_idx = iter - n_warmup;
      for (int i = 0; i < n_params; i++) {
        draws(sample_idx, i) = q[i];
      }
    }
  }

  double accept_rate = (double)n_accept / n_iter;

  return List::create(
    Named("draws") = draws,
    Named("n_divergent") = n_divergent,
    Named("accept_prob") = accept_rate,
    Named("step_size") = step_size
  );
}
