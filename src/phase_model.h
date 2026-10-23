#ifndef PHASER_PHASE_MODEL_H
#define PHASER_PHASE_MODEL_H

#include <Rcpp.h>
#include <vector>
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
//   kCholesky  standardised effects z ~ Normal(0, I) and the lower-triangular
//              factor L, packed by column, with coefficients L z.
struct RandomEffectTerm {
  enum Kind { kScalar = 0, kDiagonal = 1, kCholesky = 2 };

  int kind;
  int n_groups;
  int n_coefs;
  int offset;
  std::vector<int> group;  // level of each observation, 0-based
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
class PhaseModel {
public:
  explicit PhaseModel(const Rcpp::List& target);

  int dim() const { return n_params_; }

  template <class T>
  T log_posterior(const T* q) const {
    std::vector<T> eta;
    predictors(q, eta);
    return log_likelihood(q, eta) + log_prior(q);
  }

private:
  int n_params_;
  int K_;
  int n_obs_;
  int sigma_offset_;
  std::vector<int> unit_start_;
  std::vector<int> unit_end_;
  std::vector<int> family_;                 // per phase
  std::vector<std::vector<Observation> > obs_;  // per phase, per observation
  std::vector<LinearPredictor> lp_;         // K - 1 transitions, then K phases
  prior::PriorSet prior_;
  int max_unit_size_;

  int n_predictors() const { return static_cast<int>(lp_.size()); }

  // Coefficients of the random-effect term: the effects themselves, or
  // L z for a term with a Cholesky factor.
  template <class T>
  void term_coefficients(const RandomEffectTerm& term, const T* q,
                         std::vector<T>& b) const {
    const int n_c = term.n_coefs;
    b.assign(q + term.offset, q + term.offset + term.n_effects());
    if (term.kind != RandomEffectTerm::kCholesky) return;
    const T* L = q + term.hyper_offset();
    for (int g = 0; g < term.n_groups; g++) {
      const T* z = q + term.offset + g * n_c;
      T* bg = b.data() + g * n_c;
      for (int i = 0; i < n_c; i++) bg[i] = T(0.0);
      int l = 0;
      for (int j = 0; j < n_c; j++) {
        for (int i = j; i < n_c; i++) bg[i] = bg[i] + L[l++] * z[j];
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

  // Emission and transition terms of one unit, in the layout of
  // phase_process.h.
  template <class T>
  void unit_terms(int start, int n, const T* q, const std::vector<T>& eta,
                  T* log_emit, T* log_stay, T* log_move) const {
    const int K = K_;
    for (int s = 0; s < n; s++) {
      const int t = start + s;
      for (int k = 0; k < K; k++) {
        const T& mu = eta[static_cast<std::size_t>(K - 1 + k) * n_obs_ + t];
        log_emit[s * K + k] =
          log_emission(family_[k], obs_[k][t], mu, q[sigma_offset_ + k]);
      }
      if (s == 0) continue;
      for (int k = 0; k < K - 1; k++) {
        const T& x = eta[static_cast<std::size_t>(k) * n_obs_ + t];
        log_move[s * K + k] = log_inv_logit(x);
        log_stay[s * K + k] = log1m_inv_logit(x);
      }
    }
  }

  template <class T>
  T log_likelihood(const T* q, const std::vector<T>& eta) const {
    const std::size_t size = static_cast<std::size_t>(max_unit_size_) * K_;
    std::vector<T> log_emit(size), log_stay(size), log_move(size), alpha(size);
    T total(0.0);
    for (std::size_t i = 0; i < unit_start_.size(); i++) {
      const int start = unit_start_[i];
      const int n = unit_end_[i] - start + 1;
      unit_terms(start, n, q, eta, log_emit.data(), log_stay.data(),
                 log_move.data());
      total = total + forward(n, K_, log_emit.data(), log_stay.data(),
                              log_move.data(), alpha.data());
    }
    return total;
  }

  template <class T>
  T log_prior(const T* q) const {
    using std::abs;
    using std::log;
    T lp(0.0);
    for (const LinearPredictor& pred : lp_) {
      for (int j = 0; j < pred.n_coef; j++) {
        lp = lp + prior_.beta(q[pred.offset + j]);
      }
    }
    for (int k = 0; k < K_; k++) lp = lp + prior_.sigma(q[sigma_offset_ + k]);

    for (const LinearPredictor& pred : lp_) {
      for (const RandomEffectTerm& term : pred.re) {
        const int n_c = term.n_coefs;
        const T* effects = q + term.offset;
        const T* hyper = q + term.hyper_offset();
        if (term.kind == RandomEffectTerm::kCholesky) {
          for (int e = 0; e < term.n_effects(); e++) {
            lp = lp + prior::log_normal(effects[e], 0.0, 1.0);
          }
          const double eta = prior_.cor_eta();
          int l = 0;
          for (int j = 0; j < n_c; j++) {
            const T& L_jj = hyper[l];
            lp = lp + (n_c - j - 1 + 2.0 * (eta - 1.0)) * log(abs(L_jj)) +
                 prior_.re_scale(L_jj);
            l += n_c - j;
          }
        } else {
          for (int e = 0; e < term.n_effects(); e++) {
            const int j = (term.kind == RandomEffectTerm::kScalar) ? 0 : e % n_c;
            lp = lp + prior::log_normal_log_sd(effects[e], hyper[j]);
          }
          for (int j = 0; j < term.n_hyper(); j++) {
            lp = lp + prior_.re_sd(hyper[j]);
          }
        }
      }
    }
    return lp;
  }
};

}  // namespace phaser

#endif
