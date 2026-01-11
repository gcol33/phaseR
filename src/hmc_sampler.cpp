#include <Rcpp.h>
#include <cmath>
#include <random>

using namespace Rcpp;

// Forward declaration
double phase_log_posterior(NumericVector params, List data,
                           int n_trans_coef, int n_dyn_coef_0, int n_dyn_coef_1);

// Numerical gradient
NumericVector numerical_gradient(
    NumericVector params,
    List data,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    double eps = 1e-6
) {

  int n = params.size();
  NumericVector grad(n);

  for (int i = 0; i < n; i++) {
    NumericVector params_plus = clone(params);
    NumericVector params_minus = clone(params);
    params_plus[i] += eps;
    params_minus[i] -= eps;

    double f_plus = phase_log_posterior(params_plus, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1);
    double f_minus = phase_log_posterior(params_minus, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1);

    grad[i] = (f_plus - f_minus) / (2.0 * eps);
  }

  return grad;
}

// Single leapfrog step
void leapfrog_step(
    NumericVector& q,
    NumericVector& p,
    List data,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    double step_size
) {

  int n = q.size();

  // Half step for momentum
  NumericVector grad = numerical_gradient(q, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1);
  for (int i = 0; i < n; i++) {
    p[i] += 0.5 * step_size * grad[i];
  }

  // Full step for position
  for (int i = 0; i < n; i++) {
    q[i] += step_size * p[i];
  }

  // Half step for momentum
  grad = numerical_gradient(q, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1);
  for (int i = 0; i < n; i++) {
    p[i] += 0.5 * step_size * grad[i];
  }
}

// [[Rcpp::export]]
List phaseR_nuts_sampler(
    List data,
    NumericVector init,
    int n_iter,
    int n_warmup,
    int n_trans_coef,
    int n_dyn_coef_0,
    int n_dyn_coef_1,
    double target_accept = 0.8,
    int max_treedepth = 10
) {

  int n_params = init.size();
  int n_samples = n_iter - n_warmup;

  NumericMatrix draws(n_samples, n_params);
  NumericVector accept_probs(n_iter);
  int n_divergent = 0;

  NumericVector q = clone(init);
  double step_size = 0.1;  // Will be adapted

  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> normal(0, 1);
  std::uniform_real_distribution<> uniform(0, 1);

  // Dual averaging parameters for step size adaptation
  double mu = std::log(10.0 * step_size);
  double log_step_size_bar = 0.0;
  double H_bar = 0.0;
  double gamma = 0.05;
  double t0 = 10.0;
  double kappa = 0.75;

  for (int iter = 0; iter < n_iter; iter++) {

    // Check for user interrupt
    if (iter % 100 == 0) {
      Rcpp::checkUserInterrupt();
    }

    // Sample momentum
    NumericVector p(n_params);
    for (int i = 0; i < n_params; i++) {
      p[i] = normal(gen);
    }

    NumericVector q_proposal = clone(q);
    NumericVector p_proposal = clone(p);

    // Current Hamiltonian
    double current_H = -phase_log_posterior(q, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1);
    for (int i = 0; i < n_params; i++) {
      current_H += 0.5 * p[i] * p[i];
    }

    // Leapfrog integration (simplified, not full NUTS)
    int n_steps = std::max(1, static_cast<int>(1.0 / step_size));
    n_steps = std::min(n_steps, 100);  // Cap at 100 steps

    for (int step = 0; step < n_steps; step++) {
      leapfrog_step(q_proposal, p_proposal, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1, step_size);
    }

    // Proposed Hamiltonian
    double proposed_H = -phase_log_posterior(q_proposal, data, n_trans_coef, n_dyn_coef_0, n_dyn_coef_1);
    for (int i = 0; i < n_params; i++) {
      proposed_H += 0.5 * p_proposal[i] * p_proposal[i];
    }

    // Accept/reject
    double log_accept = current_H - proposed_H;
    double accept_prob = std::min(1.0, std::exp(log_accept));
    accept_probs[iter] = accept_prob;

    // Check for divergence
    if (std::isnan(proposed_H) || proposed_H > 1e10) {
      n_divergent++;
      accept_prob = 0.0;
    }

    if (uniform(gen) < accept_prob) {
      q = q_proposal;
    }

    // Dual averaging step size adaptation during warmup
    if (iter < n_warmup) {
      double w = 1.0 / (iter + t0);
      H_bar = (1.0 - w) * H_bar + w * (target_accept - accept_prob);
      double log_step_size = mu - std::sqrt(static_cast<double>(iter + 1)) / gamma * H_bar;
      step_size = std::exp(log_step_size);
      // Bound step size
      step_size = std::max(1e-6, std::min(step_size, 1.0));
      double m = std::pow(static_cast<double>(iter + 1), -kappa);
      log_step_size_bar = m * log_step_size + (1.0 - m) * log_step_size_bar;
    } else {
      step_size = std::exp(log_step_size_bar);
    }

    // Store sample after warmup
    if (iter >= n_warmup) {
      for (int i = 0; i < n_params; i++) {
        draws(iter - n_warmup, i) = q[i];
      }
    }
  }

  double mean_accept = 0.0;
  for (int i = 0; i < n_iter; i++) {
    mean_accept += accept_probs[i];
  }
  mean_accept /= n_iter;

  return List::create(
    Named("draws") = draws,
    Named("n_divergent") = n_divergent,
    Named("accept_prob") = mean_accept,
    Named("step_size") = step_size
  );
}
