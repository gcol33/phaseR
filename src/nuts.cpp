#include "nuts.h"

#include <cmath>
#include <limits>
#include <vector>

namespace phaser {

namespace {

const double kInf = std::numeric_limits<double>::infinity();

// Step of the central finite-difference gradient.
const double kGradientStep = 1e-6;

// Energy error above which a trajectory is declared divergent.
const double kMaxDeltaH = 1000.0;

typedef std::vector<double> Vec;

double log_sum_exp(double a, double b) {
  if (a == -kInf) return b;
  if (b == -kInf) return a;
  double m = (a > b) ? a : b;
  return m + std::log(std::exp(a - m) + std::exp(b - m));
}

double dot(const Vec& a, const Vec& b) {
  double s = 0.0;
  for (std::size_t i = 0; i < a.size(); i++) s += a[i] * b[i];
  return s;
}

Vec sum(const Vec& a, const Vec& b) {
  Vec out(a.size());
  for (std::size_t i = 0; i < a.size(); i++) out[i] = a[i] + b[i];
  return out;
}

void add_to(Vec& a, const Vec& b) {
  for (std::size_t i = 0; i < a.size(); i++) a[i] += b[i];
}

// No-U-turn criterion between the two ends of a trajectory with summed
// momentum rho: both ends still move along rho.
bool no_u_turn(const Vec& p_sharp_minus, const Vec& p_sharp_plus,
               const Vec& rho) {
  return dot(p_sharp_plus, rho) > 0.0 && dot(p_sharp_minus, rho) > 0.0;
}

// Position, momentum, log density and its gradient.
struct PhasePoint {
  Vec q;
  Vec p;
  Vec grad;
  double lp;
};

struct Transition {
  double accept_stat;
  int treedepth;
  int n_leapfrog;
  bool divergent;
  double energy;
};

// Hamiltonian dynamics with a diagonal Euclidean metric and the NUTS
// transition built on it.
class Nuts {
public:
  Nuts(LogDensity& target, Rng& rng, int max_treedepth)
    : target_(target), rng_(rng), max_treedepth_(max_treedepth),
      n_(target.dim()), inv_metric_(target.dim(), 1.0), divergent_(false) {}

  Vec& inv_metric() { return inv_metric_; }

  double hamiltonian(const PhasePoint& z) const {
    double k = 0.0;
    for (int i = 0; i < n_; i++) k += z.p[i] * z.p[i] * inv_metric_[i];
    return -z.lp + 0.5 * k;
  }

  void p_sharp(const Vec& p, Vec& out) const {
    out.resize(n_);
    for (int i = 0; i < n_; i++) out[i] = inv_metric_[i] * p[i];
  }

  void sample_momentum(PhasePoint& z) {
    for (int i = 0; i < n_; i++) {
      z.p[i] = rng_.normal() / std::sqrt(inv_metric_[i]);
    }
  }

  void leapfrog(PhasePoint& z, double epsilon) {
    for (int i = 0; i < n_; i++) z.p[i] += 0.5 * epsilon * z.grad[i];
    for (int i = 0; i < n_; i++) z.q[i] += epsilon * inv_metric_[i] * z.p[i];
    z.lp = target_.log_density_gradient(z.q, z.grad);
    for (int i = 0; i < n_; i++) z.p[i] += 0.5 * epsilon * z.grad[i];
  }

  double energy_after_step(PhasePoint& z, double epsilon) {
    leapfrog(z, epsilon);
    double h = hamiltonian(z);
    return std::isnan(h) ? kInf : h;
  }

  // Heuristic initial step size: double or halve epsilon until the
  // acceptance probability of a single leapfrog step crosses 0.8.
  double init_stepsize(const PhasePoint& z0, double epsilon) {
    const double log_target = std::log(0.8);
    PhasePoint z = z0;
    sample_momentum(z);
    double h0 = hamiltonian(z);
    double delta_h = h0 - energy_after_step(z, epsilon);
    int direction = (delta_h > log_target) ? 1 : -1;

    while (true) {
      z = z0;
      sample_momentum(z);
      h0 = hamiltonian(z);
      delta_h = h0 - energy_after_step(z, epsilon);
      if (direction == 1 && !(delta_h > log_target)) break;
      if (direction == -1 && !(delta_h < log_target)) break;
      epsilon = (direction == 1) ? 2.0 * epsilon : 0.5 * epsilon;
      if (epsilon > 1e7) {
        Rcpp::stop("NUTS step size diverged to infinity; the posterior "
                   "may be improper.");
      }
      if (epsilon == 0.0) {
        Rcpp::stop("NUTS step size collapsed to zero; the log density "
                   "is not finite near the current state.");
      }
    }
    return epsilon;
  }

  // One NUTS transition from z with step size epsilon. On return z holds
  // the selected state.
  Transition transition(PhasePoint& z_current, double epsilon) {
    divergent_ = false;

    PhasePoint z = z_current;
    sample_momentum(z);

    PhasePoint z_fwd = z;
    PhasePoint z_bck = z;
    PhasePoint z_sample = z;
    PhasePoint z_propose = z;

    // Momenta and sharp momenta at the ends of the forward and backward
    // subtrees: p_fwd_bck is the backward end of the forward subtree, etc.
    Vec p_fwd_fwd = z.p, p_fwd_bck = z.p, p_bck_fwd = z.p, p_bck_bck = z.p;
    Vec p_sharp0;
    p_sharp(z.p, p_sharp0);
    Vec p_sharp_fwd_fwd = p_sharp0, p_sharp_fwd_bck = p_sharp0;
    Vec p_sharp_bck_fwd = p_sharp0, p_sharp_bck_bck = p_sharp0;

    Vec rho = z.p;
    double log_sum_weight = 0.0;
    double h0 = hamiltonian(z);
    int n_leapfrog = 0;
    double sum_metro_prob = 0.0;
    int depth = 0;

    while (depth < max_treedepth_) {
      Vec rho_fwd(n_, 0.0);
      Vec rho_bck(n_, 0.0);
      double log_sum_weight_subtree = -kInf;
      bool valid_subtree;

      if (rng_.uniform() > 0.5) {
        // The current trajectory becomes the backward subtree.
        rho_bck = rho;
        p_bck_fwd = p_fwd_fwd;
        p_sharp_bck_fwd = p_sharp_fwd_fwd;
        z = z_fwd;
        valid_subtree = build_tree(depth, z, z_propose, p_sharp_fwd_bck,
                                   p_sharp_fwd_fwd, rho_fwd, p_fwd_bck,
                                   p_fwd_fwd, h0, epsilon, n_leapfrog,
                                   log_sum_weight_subtree, sum_metro_prob);
        z_fwd = z;
      } else {
        // The current trajectory becomes the forward subtree.
        rho_fwd = rho;
        p_fwd_bck = p_bck_bck;
        p_sharp_fwd_bck = p_sharp_bck_bck;
        z = z_bck;
        valid_subtree = build_tree(depth, z, z_propose, p_sharp_bck_fwd,
                                   p_sharp_bck_bck, rho_bck, p_bck_fwd,
                                   p_bck_bck, h0, -epsilon, n_leapfrog,
                                   log_sum_weight_subtree, sum_metro_prob);
        z_bck = z;
      }

      if (!valid_subtree) break;
      ++depth;

      // Biased progressive sampling: favour the new subtree.
      if (log_sum_weight_subtree > log_sum_weight) {
        z_sample = z_propose;
      } else if (rng_.uniform() <
                 std::exp(log_sum_weight_subtree - log_sum_weight)) {
        z_sample = z_propose;
      }
      log_sum_weight = log_sum_exp(log_sum_weight, log_sum_weight_subtree);

      rho = sum(rho_bck, rho_fwd);
      bool persist = no_u_turn(p_sharp_bck_bck, p_sharp_fwd_fwd, rho);
      Vec rho_extended = sum(rho_bck, p_fwd_bck);
      persist = persist &&
        no_u_turn(p_sharp_bck_bck, p_sharp_fwd_bck, rho_extended);
      rho_extended = sum(rho_fwd, p_bck_fwd);
      persist = persist &&
        no_u_turn(p_sharp_bck_fwd, p_sharp_fwd_fwd, rho_extended);
      if (!persist) break;
    }

    Transition t;
    t.accept_stat = (n_leapfrog > 0) ? sum_metro_prob / n_leapfrog : 0.0;
    t.treedepth = depth;
    t.n_leapfrog = n_leapfrog;
    t.divergent = divergent_;
    t.energy = hamiltonian(z_sample);

    z_current.q = z_sample.q;
    z_current.grad = z_sample.grad;
    z_current.lp = z_sample.lp;
    return t;
  }

private:
  // Builds a subtree of 2^depth leapfrog steps starting from z in the
  // direction of epsilon. z ends at the far edge of the subtree, z_propose
  // holds the state drawn from it, rho accumulates its summed momentum, and
  // p_beg / p_end (with their sharp versions) are the momenta at its two
  // ends. Returns false when the subtree diverged or contains a U-turn.
  bool build_tree(int depth, PhasePoint& z, PhasePoint& z_propose,
                  Vec& p_sharp_beg, Vec& p_sharp_end, Vec& rho,
                  Vec& p_beg, Vec& p_end, double h0, double epsilon,
                  int& n_leapfrog, double& log_sum_weight,
                  double& sum_metro_prob) {
    if (depth == 0) {
      double h = energy_after_step(z, epsilon);
      ++n_leapfrog;
      if (h - h0 > kMaxDeltaH) divergent_ = true;

      log_sum_weight = log_sum_exp(log_sum_weight, h0 - h);
      sum_metro_prob += (h0 - h > 0.0) ? 1.0 : std::exp(h0 - h);

      z_propose = z;
      p_sharp(z.p, p_sharp_beg);
      p_sharp_end = p_sharp_beg;
      add_to(rho, z.p);
      p_beg = z.p;
      p_end = p_beg;
      return !divergent_;
    }

    // First half.
    Vec p_init_end(n_), p_sharp_init_end(n_), rho_init(n_, 0.0);
    double log_sum_weight_init = -kInf;
    bool valid_init = build_tree(depth - 1, z, z_propose, p_sharp_beg,
                                 p_sharp_init_end, rho_init, p_beg,
                                 p_init_end, h0, epsilon, n_leapfrog,
                                 log_sum_weight_init, sum_metro_prob);
    if (!valid_init) return false;

    // Second half.
    PhasePoint z_propose_final = z;
    Vec p_final_beg(n_), p_sharp_final_beg(n_), rho_final(n_, 0.0);
    double log_sum_weight_final = -kInf;
    bool valid_final = build_tree(depth - 1, z, z_propose_final,
                                  p_sharp_final_beg, p_sharp_end, rho_final,
                                  p_final_beg, p_end, h0, epsilon, n_leapfrog,
                                  log_sum_weight_final, sum_metro_prob);
    if (!valid_final) return false;

    // Multinomial draw between the two halves.
    double log_sum_weight_subtree =
      log_sum_exp(log_sum_weight_init, log_sum_weight_final);
    log_sum_weight = log_sum_exp(log_sum_weight, log_sum_weight_subtree);
    if (log_sum_weight_final > log_sum_weight_subtree) {
      z_propose = z_propose_final;
    } else if (rng_.uniform() <
               std::exp(log_sum_weight_final - log_sum_weight_subtree)) {
      z_propose = z_propose_final;
    }

    Vec rho_subtree = sum(rho_init, rho_final);
    add_to(rho, rho_subtree);

    bool persist = no_u_turn(p_sharp_beg, p_sharp_end, rho_subtree);
    Vec rho_extended = sum(rho_init, p_final_beg);
    persist = persist &&
      no_u_turn(p_sharp_beg, p_sharp_final_beg, rho_extended);
    rho_extended = sum(rho_final, p_init_end);
    persist = persist &&
      no_u_turn(p_sharp_init_end, p_sharp_end, rho_extended);
    return persist;
  }

  LogDensity& target_;
  Rng& rng_;
  int max_treedepth_;
  int n_;
  Vec inv_metric_;
  bool divergent_;
};

// Dual averaging of log step size (Nesterov 2009; Hoffman and Gelman 2014,
// section 3.2) with Stan's constants.
class StepsizeAdaptation {
public:
  explicit StepsizeAdaptation(double delta)
    : delta_(delta), mu_(0.0), gamma_(0.05), kappa_(0.75), t0_(10.0) {
    restart();
  }

  void set_mu(double mu) { mu_ = mu; }

  void restart() {
    counter_ = 0.0;
    s_bar_ = 0.0;
    x_bar_ = 0.0;
  }

  double learn(double accept_stat) {
    counter_ += 1.0;
    if (accept_stat > 1.0) accept_stat = 1.0;
    double eta = 1.0 / (counter_ + t0_);
    s_bar_ = (1.0 - eta) * s_bar_ + eta * (delta_ - accept_stat);
    double x = mu_ - s_bar_ * std::sqrt(counter_) / gamma_;
    double x_eta = std::pow(counter_, -kappa_);
    x_bar_ = (1.0 - x_eta) * x_bar_ + x_eta * x;
    return std::exp(x);
  }

  double final_stepsize() const { return std::exp(x_bar_); }

private:
  double delta_, mu_, gamma_, kappa_, t0_;
  double counter_, s_bar_, x_bar_;
};

// Stan's warmup schedule for metric adaptation: an initial fast interval,
// a series of doubling slow windows in which the posterior variance is
// estimated, and a terminal fast interval.
class AdaptationWindows {
public:
  explicit AdaptationWindows(int num_warmup)
    : num_warmup_(num_warmup), init_buffer_(75), term_buffer_(50),
      base_window_(25), window_size_(0), next_window_(-1), counter_(0),
      enabled_(num_warmup >= 20) {
    if (!enabled_) return;
    if (init_buffer_ + base_window_ + term_buffer_ > num_warmup_) {
      init_buffer_ = static_cast<int>(0.15 * num_warmup_);
      term_buffer_ = static_cast<int>(0.1 * num_warmup_);
      base_window_ = num_warmup_ - (init_buffer_ + term_buffer_);
    }
    window_size_ = base_window_;
    next_window_ = init_buffer_ + window_size_ - 1;
  }

  bool enabled() const { return enabled_; }

  bool in_window() const {
    return counter_ >= init_buffer_ &&
      counter_ < num_warmup_ - term_buffer_ && counter_ != num_warmup_;
  }

  bool end_of_window() const {
    return counter_ == next_window_ && counter_ != num_warmup_;
  }

  void compute_next_window() {
    int last = num_warmup_ - term_buffer_ - 1;
    if (next_window_ == last) return;
    window_size_ *= 2;
    next_window_ = counter_ + window_size_;
    if (next_window_ != last && next_window_ + 2 * window_size_ >= last + 1) {
      next_window_ = last;
    }
  }

  void advance() { ++counter_; }

private:
  int num_warmup_, init_buffer_, term_buffer_, base_window_;
  int window_size_, next_window_, counter_;
  bool enabled_;
};

// Welford running variance.
class VarianceEstimator {
public:
  explicit VarianceEstimator(int n) : n_(n) { restart(); }

  void restart() {
    count_ = 0;
    mean_.assign(n_, 0.0);
    m2_.assign(n_, 0.0);
  }

  void add(const Vec& x) {
    ++count_;
    for (int i = 0; i < n_; i++) {
      double delta = x[i] - mean_[i];
      mean_[i] += delta / count_;
      m2_[i] += delta * (x[i] - mean_[i]);
    }
  }

  // Sample variance shrunk towards 1e-3, as in Stan.
  void regularized_variance(Vec& out) const {
    double n = static_cast<double>(count_);
    for (int i = 0; i < n_; i++) {
      double var = m2_[i] / (n - 1.0);
      out[i] = (n / (n + 5.0)) * var + 1e-3 * (5.0 / (n + 5.0));
    }
  }

private:
  int n_;
  int count_;
  Vec mean_, m2_;
};

}  // namespace


double LogDensity::log_density(const std::vector<double>& q) {
  for (int i = 0; i < dim_; i++) buffer_[i] = q[i];
  double lp = evaluate(buffer_);
  return std::isnan(lp) ? -kInf : lp;
}

double LogDensity::log_density_gradient(const std::vector<double>& q,
                                        std::vector<double>& grad) {
  grad.resize(dim_);
  for (int i = 0; i < dim_; i++) buffer_[i] = q[i];
  double lp = evaluate(buffer_);
  for (int i = 0; i < dim_; i++) {
    buffer_[i] = q[i] + kGradientStep;
    double f_plus = evaluate(buffer_);
    buffer_[i] = q[i] - kGradientStep;
    double f_minus = evaluate(buffer_);
    buffer_[i] = q[i];
    grad[i] = (f_plus - f_minus) / (2.0 * kGradientStep);
  }
  return std::isnan(lp) ? -kInf : lp;
}


Rcpp::List run_nuts(LogDensity& target, Rng& rng,
                    const NutsSettings& settings) {
  const int n = target.dim();
  const int n_iter = settings.n_iter;
  const int n_warmup = settings.n_warmup;
  const int n_samples = n_iter - n_warmup;

  // Initial state: retry until the log density and its gradient are finite.
  PhasePoint z;
  z.p.assign(n, 0.0);
  bool initialized = false;
  for (int attempt = 0; attempt < 100 && !initialized; attempt++) {
    Rcpp::NumericVector init = rng.initial_values(n);
    z.q.assign(init.begin(), init.end());
    z.lp = target.log_density_gradient(z.q, z.grad);
    initialized = std::isfinite(z.lp);
    for (int i = 0; i < n && initialized; i++) {
      initialized = std::isfinite(z.grad[i]);
    }
  }
  if (!initialized) {
    Rcpp::stop("NUTS could not find initial values with a finite log "
               "density and gradient after 100 attempts.");
  }

  Nuts sampler(target, rng, settings.max_treedepth);
  StepsizeAdaptation stepsize_adaptation(settings.adapt_delta);
  AdaptationWindows windows(n_warmup);
  VarianceEstimator variance(n);

  double epsilon = sampler.init_stepsize(z, 1.0);
  stepsize_adaptation.set_mu(std::log(10.0 * epsilon));

  Rcpp::NumericMatrix draws(n_samples, n);
  Rcpp::NumericVector accept_stat(n_samples), energy(n_samples), lp(n_samples);
  Rcpp::IntegerVector treedepth(n_samples), n_leapfrog(n_samples);
  Rcpp::LogicalVector divergent(n_samples);
  int n_divergent_warmup = 0;

  for (int iter = 0; iter < n_iter; iter++) {
    if (iter % 20 == 0) Rcpp::checkUserInterrupt();

    Transition t = sampler.transition(z, epsilon);

    if (iter < n_warmup) {
      if (t.divergent) ++n_divergent_warmup;
      epsilon = stepsize_adaptation.learn(t.accept_stat);

      if (windows.enabled()) {
        bool update = false;
        if (windows.in_window()) variance.add(z.q);
        if (windows.end_of_window()) {
          windows.compute_next_window();
          variance.regularized_variance(sampler.inv_metric());
          variance.restart();
          update = true;
        }
        windows.advance();
        if (update) {
          epsilon = sampler.init_stepsize(z, epsilon);
          stepsize_adaptation.set_mu(std::log(10.0 * epsilon));
          stepsize_adaptation.restart();
        }
      }

      if (iter == n_warmup - 1) {
        epsilon = stepsize_adaptation.final_stepsize();
      }
    } else {
      int s = iter - n_warmup;
      for (int i = 0; i < n; i++) draws(s, i) = z.q[i];
      accept_stat[s] = t.accept_stat;
      treedepth[s] = t.treedepth;
      n_leapfrog[s] = t.n_leapfrog;
      divergent[s] = t.divergent;
      energy[s] = t.energy;
      lp[s] = z.lp;
    }
  }

  int n_divergent = 0;
  double mean_accept = 0.0;
  for (int s = 0; s < n_samples; s++) {
    n_divergent += divergent[s] ? 1 : 0;
    mean_accept += accept_stat[s];
  }
  if (n_samples > 0) mean_accept /= n_samples;

  return Rcpp::List::create(
    Rcpp::Named("draws") = draws,
    Rcpp::Named("n_divergent") = n_divergent,
    Rcpp::Named("n_divergent_warmup") = n_divergent_warmup,
    Rcpp::Named("accept_prob") = mean_accept,
    Rcpp::Named("step_size") = epsilon,
    Rcpp::Named("inv_metric") = Rcpp::wrap(sampler.inv_metric()),
    Rcpp::Named("accept_stat") = accept_stat,
    Rcpp::Named("treedepth") = treedepth,
    Rcpp::Named("n_leapfrog") = n_leapfrog,
    Rcpp::Named("divergent") = divergent,
    Rcpp::Named("energy") = energy,
    Rcpp::Named("lp") = lp
  );
}

}  // namespace phaser
