#include "phase_model.h"

#include <algorithm>
#include <string>

namespace phaser {

namespace {

// Row-major copy of an R matrix with n_rows rows.
std::vector<double> row_major(const Rcpp::NumericMatrix& M, int n_rows,
                              const std::string& what) {
  if (M.nrow() != n_rows) {
    Rcpp::stop("%s has %d rows for %d observations.", what, M.nrow(), n_rows);
  }
  const int n_cols = M.ncol();
  std::vector<double> out(static_cast<std::size_t>(n_rows) * n_cols);
  for (int t = 0; t < n_rows; t++) {
    for (int j = 0; j < n_cols; j++) {
      out[static_cast<std::size_t>(t) * n_cols + j] = M(t, j);
    }
  }
  return out;
}

// Marks parameters [offset, offset + n) as read by the model.
void claim(std::vector<int>& used, int offset, int n, const char* what) {
  if (offset < 0 || n < 0 || offset + n > static_cast<int>(used.size())) {
    Rcpp::stop("Parameters of %s lie outside the parameter vector.", what);
  }
  for (int i = offset; i < offset + n; i++) used[i]++;
}

LinearPredictor read_predictor(const Rcpp::List& spec, int n_obs,
                               std::vector<int>& used) {
  LinearPredictor lp;
  Rcpp::NumericMatrix X = spec["X"];
  lp.n_coef = X.ncol();
  lp.offset = Rcpp::as<int>(spec["offset"]);
  lp.X = row_major(X, n_obs, "A design matrix");
  claim(used, lp.offset, lp.n_coef, "a linear predictor");

  Rcpp::List terms = spec["re"];
  for (R_xlen_t r = 0; r < terms.size(); r++) {
    Rcpp::List t = terms[r];
    RandomEffectTerm term;
    term.kind = Rcpp::as<int>(t["kind"]);
    term.n_groups = Rcpp::as<int>(t["n_groups"]);
    term.n_coefs = Rcpp::as<int>(t["n_coefs"]);
    term.offset = Rcpp::as<int>(t["offset"]);
    term.group = Rcpp::as<std::vector<int> >(t["group"]);
    if (term.kind < 0 || term.kind > 2 ||
        (term.kind == RandomEffectTerm::kScalar && term.n_coefs != 1)) {
      Rcpp::stop("Invalid random-effect term.");
    }
    if (static_cast<int>(term.group.size()) != n_obs) {
      Rcpp::stop("A random-effect term has %d group indices for %d "
                 "observations.", static_cast<int>(term.group.size()), n_obs);
    }
    for (int g : term.group) {
      if (g < -1 || g >= term.n_groups) {
        Rcpp::stop("Random-effect group index out of range.");
      }
    }
    if (!Rf_isNull(t["Z"])) {
      Rcpp::NumericMatrix Z = t["Z"];
      if (Z.ncol() != term.n_coefs) {
        Rcpp::stop("A random-effect design matrix has %d columns for %d "
                   "coefficients.", Z.ncol(), term.n_coefs);
      }
      term.Z = row_major(Z, n_obs, "A random-effect design matrix");
    } else if (term.n_coefs != 1) {
      Rcpp::stop("A random-effect term with several coefficients needs "
                 "a design matrix.");
    }
    claim(used, term.offset, term.n_effects() + term.n_hyper(),
          "a random-effect term");
    lp.re.push_back(term);
  }
  return lp;
}

}  // namespace


PhaseModel::PhaseModel(const Rcpp::List& target)
  : n_params_(Rcpp::as<int>(target["n_params"])),
    K_(Rcpp::as<int>(target["n_phases"])),
    sigma_offset_(Rcpp::as<int>(target["sigma_offset"])),
    unit_start_(Rcpp::as<std::vector<int> >(target["unit_start"])),
    unit_end_(Rcpp::as<std::vector<int> >(target["unit_end"])),
    family_(Rcpp::as<std::vector<int> >(target["family"])),
    prior_(Rcpp::as<Rcpp::List>(target["prior"])),
    max_unit_size_(0) {
  const std::vector<double> y = Rcpp::as<std::vector<double> >(target["y"]);
  const std::vector<int> n_trials =
    Rcpp::as<std::vector<int> >(target["n_trials"]);
  n_obs_ = static_cast<int>(y.size());

  if (K_ < 2) Rcpp::stop("A phase model needs at least two phases.");
  if (static_cast<int>(family_.size()) != K_ ||
      static_cast<int>(n_trials.size()) != K_) {
    Rcpp::stop("`family` and `n_trials` need one entry per phase.");
  }
  if (unit_start_.size() != unit_end_.size()) {
    Rcpp::stop("`unit_start` and `unit_end` differ in length.");
  }
  for (std::size_t i = 0; i < unit_start_.size(); i++) {
    if (unit_start_[i] < 0 || unit_end_[i] >= n_obs_ ||
        unit_end_[i] < unit_start_[i]) {
      Rcpp::stop("Invalid observation range of unit %d.",
                 static_cast<int>(i) + 1);
    }
    max_unit_size_ = std::max(max_unit_size_, unit_end_[i] - unit_start_[i] + 1);
  }

  obs_.resize(K_);
  for (int k = 0; k < K_; k++) {
    obs_[k].reserve(n_obs_);
    for (int t = 0; t < n_obs_; t++) {
      obs_[k].push_back(Observation(y[t], family_[k], n_trials[k]));
    }
  }

  std::vector<int> used(n_params_, 0);
  Rcpp::List trans = target["transitions"];
  Rcpp::List dyn = target["dynamics"];
  if (trans.size() != K_ - 1 || dyn.size() != K_) {
    Rcpp::stop("A model with %d phases needs %d transition and %d dynamics "
               "predictors.", K_, K_ - 1, K_);
  }
  for (int c = 0; c < K_ - 1; c++) {
    lp_.push_back(read_predictor(trans[c], n_obs_, used));
  }
  for (int k = 0; k < K_; k++) {
    lp_.push_back(read_predictor(dyn[k], n_obs_, used));
  }
  claim(used, sigma_offset_, K_, "the residual SDs");
  for (int i = 0; i < n_params_; i++) {
    if (used[i] != 1) {
      Rcpp::stop("Parameter %d is read %d times by the model; the parameter "
                 "layout is inconsistent.", i + 1, used[i]);
    }
  }
}


void PhaseModel::pointwise_log_likelihood(const double* q, double* out) const {
  std::vector<double> eta;
  predictors(q, eta);
  const std::size_t size = static_cast<std::size_t>(max_unit_size_) * K_;
  std::vector<double> log_emit(size), log_stay(size), log_move(size),
    alpha(size), d_mu(size), d_log_sigma(size);
  for (std::size_t i = 0; i < unit_start_.size(); i++) {
    const int start = unit_start_[i];
    const int n = unit_end_[i] - start + 1;
    unit_terms(start, n, q, eta, log_emit.data(), log_stay.data(),
               log_move.data(), d_mu.data(), d_log_sigma.data());
    forward(n, K_, log_emit.data(), log_stay.data(), log_move.data(),
            alpha.data());
    double previous = 0.0;
    for (int s = 0; s < n; s++) {
      const double current = log_sum_exp(alpha.data() + s * K_,
                                         n_reachable(s, K_));
      out[start + s] = current - previous;
      previous = current;
    }
  }
}

void PhaseModel::phase_summaries(const double* q, PhaseSummaries& out) const {
  const int K = K_;
  const std::size_t n_cells = static_cast<std::size_t>(n_obs_) * K;
  out.mean.assign(n_cells, 0.0);
  out.marginal.assign(n_cells, 0.0);
  out.predictive.assign(n_cells, 0.0);
  out.filtered.assign(n_cells, 0.0);
  out.move.assign(n_cells, 0.0);

  std::vector<double> eta;
  predictors(q, eta);
  for (int t = 0; t < n_obs_; t++) {
    for (int k = 0; k < K; k++) {
      const double mu = eta[static_cast<std::size_t>(K - 1 + k) * n_obs_ + t];
      double m = mu;
      if (family_[k] == kPoisson) m = std::exp(mu);
      if (family_[k] == kBinomial) m = std::exp(log_inv_logit(mu));
      out.mean[t * K + k] = m;
    }
  }

  const std::size_t size = static_cast<std::size_t>(max_unit_size_) * K;
  std::vector<double> log_emit(size), log_stay(size), log_move(size),
    alpha(size), d_mu(size), d_log_sigma(size), marginal(K), next(K);
  for (std::size_t i = 0; i < unit_start_.size(); i++) {
    const int start = unit_start_[i];
    const int n = unit_end_[i] - start + 1;
    unit_terms(start, n, q, eta, log_emit.data(), log_stay.data(),
               log_move.data(), d_mu.data(), d_log_sigma.data());
    forward(n, K, log_emit.data(), log_stay.data(), log_move.data(),
            alpha.data());

    for (int s = 0; s < n; s++) {
      const int t = start + s;
      const int reach = n_reachable(s, K);
      const double log_norm = log_sum_exp(alpha.data() + s * K, reach);
      for (int k = 0; k < reach; k++) {
        out.filtered[t * K + k] = std::exp(alpha[s * K + k] - log_norm);
      }
      if (s == 0) {
        out.predictive[t * K] = 1.0;
        out.marginal[t * K] = 1.0;
        continue;
      }
      for (int k = 0; k < K - 1; k++) {
        out.move[t * K + k] = std::exp(log_move[s * K + k]);
      }
      // One step of the phase process from the previous observation's
      // filtered and marginal distributions.
      for (int k = 0; k < K; k++) {
        const double stay = (k < K - 1) ? std::exp(log_stay[s * K + k]) : 1.0;
        const double in_pred = out.filtered[(t - 1) * K + k] * stay +
          (k > 0 ? out.filtered[(t - 1) * K + k - 1] * out.move[t * K + k - 1] : 0.0);
        const double in_marg = out.marginal[(t - 1) * K + k] * stay +
          (k > 0 ? out.marginal[(t - 1) * K + k - 1] * out.move[t * K + k - 1] : 0.0);
        out.predictive[t * K + k] = in_pred;
        out.marginal[t * K + k] = in_marg;
      }
    }
  }
}

std::vector<int> PhaseModel::random_effect_indices() const {
  std::vector<int> out;
  for (const LinearPredictor& pred : lp_) {
    for (const RandomEffectTerm& term : pred.re) {
      for (int e = 0; e < term.n_effects(); e++) out.push_back(term.offset + e);
    }
  }
  std::sort(out.begin(), out.end());
  return out;
}


std::vector<std::vector<int> > PhaseModel::unit_random_effects() const {
  std::vector<std::vector<int> > out(unit_start_.size());
  for (std::size_t i = 0; i < unit_start_.size(); i++) {
    std::vector<int>& effects = out[i];
    for (const LinearPredictor& pred : lp_) {
      for (const RandomEffectTerm& term : pred.re) {
        for (int t = unit_start_[i]; t <= unit_end_[i]; t++) {
          if (term.group[t] < 0) continue;
          const int first = term.offset + term.group[t] * term.n_coefs;
          for (int j = 0; j < term.n_coefs; j++) effects.push_back(first + j);
        }
      }
    }
    std::sort(effects.begin(), effects.end());
    effects.erase(std::unique(effects.begin(), effects.end()), effects.end());
  }
  return out;
}

}  // namespace phaser
