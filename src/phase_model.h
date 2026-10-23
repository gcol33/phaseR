#ifndef PHASER_PHASE_MODEL_H
#define PHASER_PHASE_MODEL_H

#include <Rcpp.h>
#include <vector>
#include "corr_cholesky.h"
#include "emission.h"
#include "phase_process.h"
#include "prior_functions.h"
#include "scalar_math.h"

namespace phaser {

// A random-effect term of one linear predictor: n_coefs coefficients per
// level of a grouping factor, multiplying the columns of Z (an intercept
// when Z is empty). Its parameters start at `offset`: the n_groups * n_coefs
// effects (level-major), followed by its covariance parameters:
//
//   kScalar    one coefficient, effects ~ Normal(0, sigma), then log(sigma);
//   kDiagonal  effects of coefficient j ~ Normal(0, sigma_j), then
//              log(sigma_1), ..., log(sigma_n);
//   kCholesky  correlated coefficients b = diag(sigma) L z with standardised
//              effects z ~ Normal(0, I), followed by log(sigma_1), ...,
//              log(sigma_n) and the n (n - 1) / 2 unconstrained values of
//              the Cholesky factor L of the correlation matrix
//              (corr_cholesky.h), which gets the LKJ prior.
struct RandomEffectTerm {
  enum Kind { kScalar = 0, kDiagonal = 1, kCholesky = 2 };

  int kind;
  int n_groups;
  int n_coefs;
  int offset;
  std::vector<int> group;  // level of each observation, 0-based; -1 for a
                           // level without an effect (new in predictions)
  std::vector<double> Z;   // n_obs x n_coefs, row-major; empty: intercept

  int n_effects() const { return n_groups * n_coefs; }
  int hyper_offset() const { return offset + n_effects(); }
  int n_hyper() const {
    if (kind == kScalar) return 1;
    if (kind == kDiagonal) return n_coefs;
    return n_coefs * (n_coefs + 1) / 2;
  }
};

// A linear predictor X beta + sum of random-effect terms, beta starting at
// `offset` in the parameter vector.
struct LinearPredictor {
  int n_coef;
  int offset;
  std::vector<double> X;  // n_obs x n_coef, row-major
  std::vector<RandomEffectTerm> re;
};

// Log posterior of a phase model with K ordered phases, built from the
// description returned by posterior_target() in R. Linear predictor c is the
// logit of moving from phase c to c + 1 for c < K - 1, and the mean (on the
// link scale) of phase c - (K - 1) otherwise. The log residual SD of phase k
// is parameter sigma_offset + k.
//
// The gradient is exact: the backward recursion of the phase process gives
// the derivatives of each unit's log likelihood with respect to its emission
// and transition terms, which the chain rule carries through the families,
// the linear predictors and the random-effect parameterisations. Every
// function is a template on the scalar type, so forward-mode scalars give
// higher derivatives through the same code.
class PhaseModel {
public:
  explicit PhaseModel(const Rcpp::List& target);

  int dim() const { return n_params_; }

  // Parameter indices of the random effects of every term (the effects, not
  // their covariance parameters), in increasing order.
  std::vector<int> random_effect_indices() const;

  // For each unit, the indices of the random effects its likelihood depends
  // on, in increasing order.
  std::vector<std::vector<int> > unit_random_effects() const;

  // Per observation t and phase k, in arrays indexed [t * K + k]:
  //   mean        mean of the outcome in phase k (response scale);
  //   marginal    p(phase_t = k) given the covariates alone;
  //   predictive  p(phase_t = k | y_0, ..., y_{t-1}) within the unit;
  //   filtered    p(phase_t = k | y_0, ..., y_t) within the unit;
  //   move        probability of moving from phase k to k + 1 between the
  //               previous observation and t (zero at a unit's first
  //               observation and for the last phase).
  struct PhaseSummaries {
    std::vector<double> mean, marginal, predictive, filtered, move;
  };
  void phase_summaries(const double* q, PhaseSummaries& out) const;

  // Per random-effect term, in the order of the predictors: the
  // coefficients of every level (level-major), the SD of each coefficient
  // and their correlation matrix (row-major).
  struct TermValues {
    std::vector<double> coefficients, sd, cor;
  };
  std::vector<TermValues> random_effect_values(const double* q) const;

  int n_phases() const { return K_; }
  int n_obs() const { return n_obs_; }

  // Log likelihood of each observation given the earlier observations of
  // its unit, log p(y_s | y_0, ..., y_{s-1}), marginalised over the latent
  // phases; out has one entry per observation. The entries of a unit sum to
  // its log likelihood.
  void pointwise_log_likelihood(const double* q, double* out) const;

  // Log likelihood of each unit.
  std::vector<double> unit_log_likelihoods(const double* q) const;

  // Log prior density of each random effect given its covariance
  // parameters, indexed like the parameter vector (zero elsewhere).
  std::vector<double> effect_log_priors(const double* q) const;

  // Prior SD of each random effect given its covariance parameters (1 for
  // the standardised effects of a correlated term), indexed like the
  // parameter vector (zero elsewhere).
  std::vector<double> effect_prior_sds(const double* q) const;

  // Log posterior at q (length dim()). Unless grad is null, writes its
  // gradient into grad (length dim()).
  template <class T>
  T log_posterior(const T* q, T* grad) const {
    std::vector<T> eta, d_eta;
    predictors(q, eta);
    if (grad != nullptr) {
      for (int i = 0; i < n_params_; i++) grad[i] = T(0.0);
      d_eta.assign(eta.size(), T(0.0));
    }
    T lp = log_likelihood(q, eta, grad != nullptr ? d_eta.data() : nullptr,
                          grad);
    if (grad != nullptr) predictors_adjoint(q, d_eta, grad);
    return lp + log_prior(q, grad);
  }

private:
  int n_params_;
  int K_;
  int n_obs_;
  int sigma_offset_;
  std::vector<int> unit_start_;
  std::vector<int> unit_end_;
  std::vector<int> family_;                     // per phase
  std::vector<std::vector<Observation> > obs_;  // per phase, per observation
  std::vector<LinearPredictor> lp_;             // K - 1 transitions, K phases
  prior::PriorSet prior_;
  int max_unit_size_;

  int n_predictors() const { return static_cast<int>(lp_.size()); }

  // Coefficients of the random-effect term: the effects themselves, or
  // diag(sigma) L z for a correlated term.
  template <class T>
  void term_coefficients(const RandomEffectTerm& term, const T* q,
                         std::vector<T>& b) const {
    using std::exp;
    const int n_c = term.n_coefs;
    b.assign(q + term.offset, q + term.offset + term.n_effects());
    if (term.kind != RandomEffectTerm::kCholesky) return;
    const T* log_sd = q + term.hyper_offset();
    CorrFactor<T> cf;
    corr_factor(log_sd + n_c, n_c, prior_.cor_eta(), false, cf);
    for (int g = 0; g < term.n_groups; g++) {
      const T* z = q + term.offset + g * n_c;
      T* bg = b.data() + g * n_c;
      for (int i = 0; i < n_c; i++) {
        T sum(0.0);
        for (int j = 0; j <= i; j++) sum = sum + cf.L[i * n_c + j] * z[j];
        bg[i] = exp(log_sd[i]) * sum;
      }
    }
  }

  // Adds to grad the derivatives with respect to the term's parameters,
  // given the derivatives d_b with respect to its coefficients.
  template <class T>
  void term_coefficients_adjoint(const RandomEffectTerm& term, const T* q,
                                 const std::vector<T>& d_b, T* grad) const {
    using std::exp;
    const int n_c = term.n_coefs;
    if (term.kind != RandomEffectTerm::kCholesky) {
      for (int e = 0; e < term.n_effects(); e++) {
        grad[term.offset + e] = grad[term.offset + e] + d_b[e];
      }
      return;
    }
    const T* log_sd = q + term.hyper_offset();
    T* d_log_sd = grad + term.hyper_offset();
    T* d_y = d_log_sd + n_c;
    const int n_y = n_c * (n_c - 1) / 2;
    CorrFactor<T> cf;
    corr_factor(log_sd + n_c, n_c, prior_.cor_eta(), true, cf);
    std::vector<T> d_L(n_c * n_c, T(0.0));
    for (int g = 0; g < term.n_groups; g++) {
      const T* z = q + term.offset + g * n_c;
      T* d_z = grad + term.offset + g * n_c;
      const T* d_bg = d_b.data() + g * n_c;
      for (int i = 0; i < n_c; i++) {
        const T sd = exp(log_sd[i]);
        T Lz(0.0);
        for (int j = 0; j <= i; j++) {
          Lz = Lz + cf.L[i * n_c + j] * z[j];
          d_z[j] = d_z[j] + sd * cf.L[i * n_c + j] * d_bg[i];
          d_L[i * n_c + j] = d_L[i * n_c + j] + sd * z[j] * d_bg[i];
        }
        d_log_sd[i] = d_log_sd[i] + sd * Lz * d_bg[i];
      }
    }
    for (int m = 0; m < n_y; m++) {
      for (int i = 0; i < n_c * n_c; i++) {
        d_y[m] = d_y[m] + d_L[i] * cf.dL[m * n_c * n_c + i];
      }
    }
  }

  // Linear predictors of every observation: eta[c * n_obs + t].
  template <class T>
  void predictors(const T* q, std::vector<T>& eta) const {
    eta.assign(static_cast<std::size_t>(n_predictors()) * n_obs_, T(0.0));
    std::vector<T> b;
    for (int c = 0; c < n_predictors(); c++) {
      const LinearPredictor& lp = lp_[c];
      T* eta_c = eta.data() + static_cast<std::size_t>(c) * n_obs_;
      const T* beta = q + lp.offset;
      for (int t = 0; t < n_obs_; t++) {
        const double* x = lp.X.data() + static_cast<std::size_t>(t) * lp.n_coef;
        T sum(0.0);
        for (int j = 0; j < lp.n_coef; j++) sum = sum + beta[j] * x[j];
        eta_c[t] = sum;
      }
      for (const RandomEffectTerm& term : lp.re) {
        term_coefficients(term, q, b);
        const int n_c = term.n_coefs;
        for (int t = 0; t < n_obs_; t++) {
          if (term.group[t] < 0) continue;
          const T* bg = b.data() + term.group[t] * n_c;
          if (term.Z.empty()) {
            eta_c[t] = eta_c[t] + bg[0];
          } else {
            const double* z = term.Z.data() + static_cast<std::size_t>(t) * n_c;
            for (int j = 0; j < n_c; j++) eta_c[t] = eta_c[t] + bg[j] * z[j];
          }
        }
      }
    }
  }

  // Adds to grad the derivatives with respect to the parameters, given the
  // derivatives d_eta with respect to the linear predictors.
  template <class T>
  void predictors_adjoint(const T* q, const std::vector<T>& d_eta,
                          T* grad) const {
    std::vector<T> d_b;
    for (int c = 0; c < n_predictors(); c++) {
      const LinearPredictor& lp = lp_[c];
      const T* d_eta_c = d_eta.data() + static_cast<std::size_t>(c) * n_obs_;
      T* d_beta = grad + lp.offset;
      for (int t = 0; t < n_obs_; t++) {
        const double* x = lp.X.data() + static_cast<std::size_t>(t) * lp.n_coef;
        for (int j = 0; j < lp.n_coef; j++) {
          d_beta[j] = d_beta[j] + d_eta_c[t] * x[j];
        }
      }
      for (const RandomEffectTerm& term : lp.re) {
        const int n_c = term.n_coefs;
        d_b.assign(term.n_effects(), T(0.0));
        for (int t = 0; t < n_obs_; t++) {
          if (term.group[t] < 0) continue;
          T* d_bg = d_b.data() + term.group[t] * n_c;
          if (term.Z.empty()) {
            d_bg[0] = d_bg[0] + d_eta_c[t];
          } else {
            const double* z = term.Z.data() + static_cast<std::size_t>(t) * n_c;
            for (int j = 0; j < n_c; j++) d_bg[j] = d_bg[j] + d_eta_c[t] * z[j];
          }
        }
        term_coefficients_adjoint(term, q, d_b, grad);
      }
    }
  }

  // Emission and transition terms of one unit, in the layout of
  // phase_process.h, and the derivatives of the log emission densities with
  // respect to the phase means (d_mu) and log residual SDs (d_log_sigma).
  template <class T>
  void unit_terms(int start, int n, const T* q, const std::vector<T>& eta,
                  T* log_emit, T* log_stay, T* log_move, T* d_mu,
                  T* d_log_sigma) const {
    const int K = K_;
    for (int s = 0; s < n; s++) {
      const int t = start + s;
      for (int k = 0; k < K; k++) {
        const T& mu = eta[static_cast<std::size_t>(K - 1 + k) * n_obs_ + t];
        log_emit[s * K + k] =
          log_emission(family_[k], obs_[k][t], mu, q[sigma_offset_ + k],
                       d_mu[s * K + k], d_log_sigma[s * K + k]);
      }
      if (s == 0) continue;
      for (int k = 0; k < K - 1; k++) {
        const T& x = eta[static_cast<std::size_t>(k) * n_obs_ + t];
        log_move[s * K + k] = log_inv_logit(x);
        log_stay[s * K + k] = log1m_inv_logit(x);
      }
    }
  }

  // Log likelihood. Unless d_eta is null, adds its derivatives with respect
  // to the linear predictors into d_eta and those with respect to the log
  // residual SDs into grad. A transition predictor x enters as
  // log inv_logit(x) (move) and log(1 - inv_logit(x)) (stay), with
  // derivatives 1 - inv_logit(x) and -inv_logit(x).
  template <class T>
  T log_likelihood(const T* q, const std::vector<T>& eta, T* d_eta,
                   T* grad) const {
    using std::exp;
    const int K = K_;
    const std::size_t size = static_cast<std::size_t>(max_unit_size_) * K;
    std::vector<T> log_emit(size), log_stay(size), log_move(size), alpha(size),
      d_mu(size), d_log_sigma(size);
    std::vector<T> beta, w_emit, w_stay, w_move;
    if (d_eta != nullptr) {
      beta.resize(size);
      w_emit.resize(size);
      w_stay.resize(size);
      w_move.resize(size);
    }
    T total(0.0);
    for (std::size_t i = 0; i < unit_start_.size(); i++) {
      const int start = unit_start_[i];
      const int n = unit_end_[i] - start + 1;
      unit_terms(start, n, q, eta, log_emit.data(), log_stay.data(),
                 log_move.data(), d_mu.data(), d_log_sigma.data());
      T ll = forward(n, K, log_emit.data(), log_stay.data(), log_move.data(),
                     alpha.data());
      total = total + ll;
      if (d_eta == nullptr) continue;

      backward(n, K, log_emit.data(), log_stay.data(), log_move.data(),
               alpha.data(), ll, beta.data(), w_emit.data(), w_stay.data(),
               w_move.data());
      for (int s = 0; s < n; s++) {
        const int t = start + s;
        for (int k = 0; k < K; k++) {
          const int j = s * K + k;
          T* d = d_eta + static_cast<std::size_t>(K - 1 + k) * n_obs_ + t;
          *d = *d + w_emit[j] * d_mu[j];
          grad[sigma_offset_ + k] =
            grad[sigma_offset_ + k] + w_emit[j] * d_log_sigma[j];
        }
        if (s == 0) continue;
        for (int k = 0; k < K - 1; k++) {
          const int j = s * K + k;
          T* d = d_eta + static_cast<std::size_t>(k) * n_obs_ + t;
          *d = *d + w_move[j] * exp(log_stay[j]) - w_stay[j] * exp(log_move[j]);
        }
      }
    }
    return total;
  }

  // Log prior; adds its gradient into grad unless grad is null.
  template <class T>
  T log_prior(const T* q, T* grad) const {
    T lp(0.0), d(0.0), d_sd(0.0);
    auto add = [grad](int i, const T& v) {
      if (grad != nullptr) grad[i] = grad[i] + v;
    };
    for (const LinearPredictor& pred : lp_) {
      for (int j = 0; j < pred.n_coef; j++) {
        lp = lp + prior_.beta(q[pred.offset + j], d);
        add(pred.offset + j, d);
      }
    }
    for (int k = 0; k < K_; k++) {
      lp = lp + prior_.sigma(q[sigma_offset_ + k], d);
      add(sigma_offset_ + k, d);
    }

    for (const LinearPredictor& pred : lp_) {
      for (const RandomEffectTerm& term : pred.re) {
        const int n_c = term.n_coefs;
        const int effects = term.offset;
        const int hyper = term.hyper_offset();
        if (term.kind == RandomEffectTerm::kCholesky) {
          for (int e = 0; e < term.n_effects(); e++) {
            lp = lp + prior::log_normal(q[effects + e], 0.0, 1.0, d);
            add(effects + e, d);
          }
          for (int j = 0; j < n_c; j++) {
            lp = lp + prior_.re_sd(q[hyper + j], d);
            add(hyper + j, d);
          }
          CorrFactor<T> cf;
          corr_factor(q + hyper + n_c, n_c, prior_.cor_eta(), grad != nullptr,
                      cf);
          lp = lp + cf.log_density;
          if (grad != nullptr) {
            for (int m = 0; m < n_c * (n_c - 1) / 2; m++) {
              add(hyper + n_c + m, cf.d_log_density[m]);
            }
          }
        } else {
          for (int e = 0; e < term.n_effects(); e++) {
            const int j = (term.kind == RandomEffectTerm::kScalar) ? 0 : e % n_c;
            lp = lp + prior::log_normal_log_sd(q[effects + e], q[hyper + j],
                                               d, d_sd);
            add(effects + e, d);
            add(hyper + j, d_sd);
          }
          for (int j = 0; j < term.n_hyper(); j++) {
            lp = lp + prior_.re_sd(q[hyper + j], d);
            add(hyper + j, d);
          }
        }
      }
    }
    return lp;
  }
};

}  // namespace phaser

#endif
