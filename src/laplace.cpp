#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "block_precision.h"
#include "dual.h"
#include "phase_model.h"

namespace phaser {

namespace {

typedef Dual<double> D1;
typedef Dual<D1> D2;

// Greedy colouring of the random effects such that no two effects of one
// colour share a structurally nonzero row of the Hessian d2f/du2. The Hessian
// couples two effects when some unit's likelihood depends on both (and each
// effect with itself through its prior); rows[a] lists the effects coupled
// with effect a. Seeding all effects of one colour at once then recovers
// their Hessian columns from a single Hessian-vector product.
std::vector<std::vector<int> > colour_effects(
    const std::vector<std::vector<int> >& rows) {
  const int n = static_cast<int>(rows.size());
  std::vector<int> colour(n, -1);
  std::vector<int> forbidden_at(n + 1, -1);
  int n_colours = 0;
  for (int a = 0; a < n; a++) {
    for (int r : rows[a]) {
      for (int b : rows[r]) {
        if (colour[b] >= 0) forbidden_at[colour[b]] = a;
      }
    }
    int c = 0;
    while (forbidden_at[c] == a) c++;
    colour[a] = c;
    n_colours = std::max(n_colours, c + 1);
  }
  std::vector<std::vector<int> > out(n_colours);
  for (int a = 0; a < n; a++) out[colour[a]].push_back(a);
  return out;
}

}  // namespace


// Laplace approximation to the marginal posterior of the parameters theta of
// a phase model, with the random effects u integrated out:
//
//   log p(theta | y) ~= f(theta, u*) + (n_u / 2) log(2 pi) - log det(A) / 2,
//
// where f is the log joint posterior, u*(theta) maximises f over u, and
// A = -d2f/du2 at u*. Given theta the effects have a Gaussian prior, so this
// is the approximation of Tierney and Kadane (1986) used by INLA and TMB.
// Its gradient is exact. With phi(q) = -log det(-d2f/du2) / 2 and r the
// derivatives of f + phi with respect to u,
//
//   d/dtheta = d(f + phi)/dtheta + r^T du*/dtheta,
//   du*/dtheta = A^{-1} d2f/du dtheta,
//
// so the last term is d2f/dtheta du s with s = A^{-1} r: one
// Hessian-vector product. The derivative of phi with respect to every
// parameter is tr(A^{-1} d(d2f/du2)) / 2; on the Hessian's sparsity pattern
// it is a sum over the colours c of the effects of w_c^T H'' c, the third
// derivative of f along c and w_c (the entries of A^{-1} in c's columns),
// one pass with nested forward-mode dual numbers through the exact gradient
// of the log posterior per colour.
class LaplaceMarginal {
public:
  explicit LaplaceMarginal(const Rcpp::List& target)
    : model_(target), have_value_(false), have_gradient_(false),
      inner_converged_(true) {
    u_index_ = model_.random_effect_indices();
    const int n = model_.dim();
    std::vector<int> u_position(n, -1);
    for (std::size_t r = 0; r < u_index_.size(); r++) {
      u_position[u_index_[r]] = static_cast<int>(r);
    }
    for (int i = 0; i < n; i++) {
      if (u_position[i] < 0) theta_index_.push_back(i);
    }

    // Effects coupled in the Hessian: those sharing a unit, and each with
    // itself.
    const int n_u = n_u_count();
    rows_.assign(n_u, std::vector<int>());
    for (int a = 0; a < n_u; a++) rows_[a].push_back(a);
    for (const std::vector<int>& effects : model_.unit_random_effects()) {
      for (int a : effects) {
        for (int b : effects) rows_[u_position[a]].push_back(u_position[b]);
      }
    }
    for (std::vector<int>& r : rows_) {
      std::sort(r.begin(), r.end());
      r.erase(std::unique(r.begin(), r.end()), r.end());
    }
    colours_ = colour_effects(rows_);
    precision_ = BlockPrecision(rows_);
    for (const std::vector<int>& effects : model_.unit_random_effects()) {
      unit_block_.push_back(effects.empty()
                              ? -1 : precision_.block_of(u_position[effects[0]]));
    }
    u_.assign(n_u, 0.0);
  }

  int n_theta() const { return static_cast<int>(theta_index_.size()); }
  int n_u_count() const { return static_cast<int>(u_index_.size()); }
  const std::vector<int>& theta_index() const { return theta_index_; }
  const std::vector<int>& u_index() const { return u_index_; }
  bool inner_converged() const { return inner_converged_; }

  // Laplace log marginal posterior at theta; its gradient into grad unless
  // grad is null.
  double evaluate(const std::vector<double>& theta, std::vector<double>* grad) {
    if (!(have_value_ && theta == theta_)) {
      theta_ = theta;
      find_mode();
      have_value_ = true;
      have_gradient_ = false;
    }
    if (grad != nullptr) {
      if (!have_gradient_) {
        compute_gradient();
        have_gradient_ = true;
      }
      *grad = gradient_;
    }
    return value_;
  }

  // Number of blocks of effects whose Newton iterations from u = 0 and from
  // the effects shifted by +2 and -2 prior SDs end at different maxima at
  // theta: the effects given theta have several modes there, and a Laplace
  // approximation around one of them misses the others.
  int multimodal_blocks(const std::vector<double>& theta) {
    const int n_u = n_u_count();
    if (n_u == 0) return 0;
    const int n_blocks = precision_.n_blocks();
    std::vector<double> base = full(theta, std::vector<double>(n_u, 0.0));
    const std::vector<double> sd = model_.effect_prior_sds(base.data());
    std::vector<std::vector<double> > ends;
    std::vector<double> f_block;
    std::vector<char> converged;
    for (double shift : {0.0, 2.0, -2.0}) {
      std::vector<double> q = base;
      for (int r = 0; r < n_u; r++) q[u_index_[r]] = shift * sd[u_index_[r]];
      if (newton_blocks(q, f_block, converged)) ends.push_back(q);
    }
    have_value_ = false;
    int count = 0;
    for (int b = 0; b < n_blocks; b++) {
      bool differ = false;
      for (std::size_t k = 1; k < ends.size() && !differ; k++) {
        for (int r : precision_.members(b)) {
          const double x0 = ends[0][u_index_[r]];
          if (std::fabs(ends[k][u_index_[r]] - x0) > 1e-4 * (1.0 + std::fabs(x0))) {
            differ = true;
            break;
          }
        }
      }
      if (differ) count++;
    }
    return count;
  }

  // Hessian of the log marginal posterior (row-major), by central
  // differences of its exact gradient.
  std::vector<double> hessian(const std::vector<double>& theta) {
    const int p = n_theta();
    std::vector<double> H(static_cast<std::size_t>(p) * p);
    std::vector<double> g_plus, g_minus;
    for (int j = 0; j < p; j++) {
      const double h = 1e-5 * std::max(1.0, std::fabs(theta[j]));
      std::vector<double> step = theta;
      step[j] = theta[j] + h;
      evaluate(step, &g_plus);
      step[j] = theta[j] - h;
      evaluate(step, &g_minus);
      for (int i = 0; i < p; i++) H[i * p + j] = (g_plus[i] - g_minus[i]) / (2.0 * h);
    }
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < i; j++) {
        const double m = 0.5 * (H[i * p + j] + H[j * p + i]);
        H[i * p + j] = m;
        H[j * p + i] = m;
      }
    }
    evaluate(theta, nullptr);
    return H;
  }

  // At theta: the mode u*, the diagonal of A^{-1} (the variances of the
  // effects given theta) and du*/dtheta = A^{-1} d2f/du dtheta
  // (n_u x n_theta, row-major).
  void conditional(const std::vector<double>& theta, std::vector<double>& mode,
                   std::vector<double>& variance, std::vector<double>& slope) {
    evaluate(theta, nullptr);
    mode = u_;
    variance.resize(n_u_count());
    for (int r = 0; r < n_u_count(); r++) variance[r] = precision_.inverse_at(r, r);
    const int p = n_theta();
    const int n_u = n_u_count();
    slope.assign(static_cast<std::size_t>(n_u) * p, 0.0);
    std::vector<D1> qd(q_.size()), gd(q_.size());
    std::vector<double> column(n_u);
    for (int j = 0; j < p; j++) {
      for (std::size_t i = 0; i < q_.size(); i++) qd[i] = D1(q_[i], 0.0);
      qd[theta_index_[j]].d = 1.0;
      model_.log_posterior(qd.data(), gd.data());
      for (int r = 0; r < n_u; r++) column[r] = gd[u_index_[r]].d;
      precision_.solve(column.data());
      for (int r = 0; r < n_u; r++) slope[r * p + j] = column[r];
    }
  }

  // Draws of the effects' deviation from their mode given theta, A^{-1/2} z,
  // for the rows of z (each of length n_u).
  void effect_noise(const std::vector<double>& theta, std::vector<double>& z,
                    int n_rows) {
    evaluate(theta, nullptr);
    for (int i = 0; i < n_rows; i++) {
      precision_.solve_factor_transpose(z.data() + static_cast<std::size_t>(i) * n_u_count());
    }
  }

private:
  PhaseModel model_;
  std::vector<int> theta_index_;
  std::vector<int> u_index_;
  std::vector<std::vector<int> > rows_;
  std::vector<std::vector<int> > colours_;

  std::vector<double> theta_;
  std::vector<double> u_;          // mode of the effects at theta_
  std::vector<double> q_;          // full parameter vector (theta_, u_)
  std::vector<double> g_;          // gradient of f at q_
  BlockPrecision precision_;       // A = -d2f/du2 at q_, factorised
  std::vector<int> unit_block_;    // block of each unit's effects, -1: none
  std::vector<double> gradient_;
  double value_;
  bool have_value_;
  bool have_gradient_;
  bool inner_converged_;

  std::vector<double> full(const std::vector<double>& theta,
                           const std::vector<double>& u) const {
    std::vector<double> q(model_.dim());
    for (std::size_t j = 0; j < theta_index_.size(); j++) q[theta_index_[j]] = theta[j];
    for (std::size_t r = 0; r < u_index_.size(); r++) q[u_index_[r]] = u[r];
    return q;
  }

  // A = -d2f/du2 at q, from one Hessian-vector product per colour.
  void precision_at(const std::vector<double>& q) {
    precision_.set_zero();
    std::vector<D1> qd(q.size()), gd(q.size());
    for (const std::vector<int>& colour : colours_) {
      for (std::size_t i = 0; i < q.size(); i++) qd[i] = D1(q[i], 0.0);
      for (int a : colour) qd[u_index_[a]].d = 1.0;
      model_.log_posterior(qd.data(), gd.data());
      for (int a : colour) {
        for (int r : rows_[a]) precision_.at(r, a) = -gd[u_index_[r]].d;
      }
    }
  }

  // Log posterior of each block of effects: the log likelihood of the units
  // whose effects lie in the block plus the effects' prior. Given theta, f is
  // the sum of these and terms that do not depend on the effects, so the
  // blocks can be maximised separately.
  std::vector<double> block_values(const std::vector<double>& q) const {
    std::vector<double> out(precision_.n_blocks(), 0.0);
    const std::vector<double> unit_ll = model_.unit_log_likelihoods(q.data());
    for (std::size_t i = 0; i < unit_ll.size(); i++) {
      if (unit_block_[i] >= 0) out[unit_block_[i]] += unit_ll[i];
    }
    const std::vector<double> priors = model_.effect_log_priors(q.data());
    for (int r = 0; r < n_u_count(); r++) {
      out[precision_.block_of(r)] += priors[u_index_[r]];
    }
    return out;
  }

  // Newton's method with backtracking from the effects in q, run in every
  // block of effects with its own damping, step length and stopping
  // decision, so the mode reached in one block does not depend on the
  // others. On return q holds the effects reached, f_block the value of each
  // block and converged whether each block stopped at a strict maximum.
  // Returns false if the log posterior or its Hessian in the effects is not
  // finite along the way.
  bool newton_blocks(std::vector<double>& q, std::vector<double>& f_block,
                     std::vector<char>& converged) {
    const int n_u = n_u_count();
    const int n_blocks = precision_.n_blocks();
    std::vector<double> g(q.size());
    if (!std::isfinite(model_.log_posterior(q.data(), g.data()))) return false;
    f_block = block_values(q);
    converged.assign(n_blocks, 0);
    std::vector<double> delta(n_u, 0.0), decrement(n_blocks, 0.0),
      step(n_blocks, 1.0), last_decrement(n_blocks, HUGE_VAL);
    std::vector<char> done(n_blocks, 0), pending(n_blocks, 0),
      damped(n_blocks, 0);
    std::vector<double> q_new;
    for (int iter = 0; iter < 200; iter++) {
      precision_at(q);
      int n_active = 0;
      for (int b = 0; b < n_blocks; b++) {
        if (done[b]) continue;
        const double max_diag = precision_.max_abs_diagonal(b);
        double damping = 0.0;
        while (std::isfinite(max_diag) && damping < 1e12 * (1.0 + max_diag) &&
               !precision_.factor_block(b, damping)) {
          damping = std::max(1e-8 * (1.0 + max_diag), 10.0 * damping);
        }
        if (!std::isfinite(max_diag) || damping >= 1e12 * (1.0 + max_diag)) {
          return false;
        }
        damped[b] = damping > 0.0;
        decrement[b] = 0.0;
        for (int r : precision_.members(b)) delta[r] = g[u_index_[r]];
        precision_.solve_block(b, delta.data());
        for (int r : precision_.members(b)) decrement[b] += g[u_index_[r]] * delta[r];
        // The Newton decrement is twice the remaining ascent of the local
        // quadratic model and falls quadratically near the mode. Iterate
        // until it vanishes or, once small, stops falling: the block is
        // then at the resolution of double precision.
        const double scale = 1.0 + std::fabs(f_block[b]);
        if (!damped[b] &&
            (decrement[b] < 1e-24 * scale ||
             (decrement[b] < 1e-12 * scale && decrement[b] > 0.25 * last_decrement[b]))) {
          done[b] = 1;
          converged[b] = 1;
          continue;
        }
        last_decrement[b] = damped[b] ? HUGE_VAL : decrement[b];
        step[b] = 1.0;
        pending[b] = 1;
        n_active++;
      }
      if (n_active == 0) break;

      // Backtracking in every active block at once.
      q_new = q;
      for (int halving = 0; halving < 60 && n_active > 0; halving++) {
        for (int b = 0; b < n_blocks; b++) {
          if (!pending[b]) continue;
          for (int r : precision_.members(b)) {
            q_new[u_index_[r]] = q[u_index_[r]] + step[b] * delta[r];
          }
        }
        const std::vector<double> f_new = block_values(q_new);
        for (int b = 0; b < n_blocks; b++) {
          if (!pending[b]) continue;
          if (std::isfinite(f_new[b]) &&
              f_new[b] >= f_block[b] + 1e-4 * step[b] * decrement[b]) {
            pending[b] = 0;
            f_block[b] = f_new[b];
            n_active--;
          } else {
            step[b] *= 0.5;
          }
        }
      }
      for (int b = 0; b < n_blocks; b++) {
        if (!pending[b]) continue;
        // No ascent along the Newton direction: the block is at the limit
        // of floating-point resolution.
        for (int r : precision_.members(b)) q_new[u_index_[r]] = q[u_index_[r]];
        pending[b] = 0;
        done[b] = 1;
        converged[b] = !damped[b];
      }
      q.swap(q_new);
      if (!std::isfinite(model_.log_posterior(q.data(), g.data()))) return false;
    }
    return true;
  }

  // Maximises f over u at theta_ by Newton's method from u = 0, the prior
  // mean of the effects, so that the mode reached, and with it the
  // approximation, is a function of theta alone. Sets u_, q_, g_,
  // precision_ and value_.
  void find_mode() {
    const int n_u = n_u_count();
    inner_converged_ = true;
    q_ = full(theta_, std::vector<double>(n_u, 0.0));
    g_.assign(q_.size(), 0.0);
    if (n_u > 0) {
      std::vector<double> f_block;
      std::vector<char> converged;
      if (!newton_blocks(q_, f_block, converged)) {
        value_ = -HUGE_VAL;
        return;
      }
      for (char c : converged) {
        if (!c) inner_converged_ = false;
      }
      for (int r = 0; r < n_u; r++) u_[r] = q_[u_index_[r]];
    }
    const double f = model_.log_posterior(q_.data(), g_.data());
    if (!std::isfinite(f)) {
      value_ = -HUGE_VAL;
      return;
    }
    if (n_u > 0) {
      precision_at(q_);
      if (!precision_.factor(0.0)) {
        // No strict maximum in the effects: the approximation is not
        // defined at these parameters.
        value_ = -HUGE_VAL;
        return;
      }
    }
    value_ = f + 0.5 * n_u * kLog2Pi - (n_u > 0 ? 0.5 * precision_.log_det() : 0.0);
  }

  void compute_gradient() {
    const int p = n_theta();
    const int n_u = n_u_count();
    const std::size_t n = q_.size();
    if (!std::isfinite(value_)) {
      gradient_.assign(p, NAN);
      return;
    }
    gradient_.assign(p, 0.0);
    for (int j = 0; j < p; j++) gradient_[j] = g_[theta_index_[j]];
    if (n_u == 0) return;

    // Derivatives of phi with respect to every parameter.
    std::vector<double> d_phi(n, 0.0);
    std::vector<D2> qdd(n), gdd(n);
    for (const std::vector<int>& colour : colours_) {
      for (std::size_t i = 0; i < n; i++) qdd[i] = D2(D1(q_[i], 0.0), D1(0.0, 0.0));
      for (int a : colour) {
        qdd[u_index_[a]].v.d = 1.0;
        for (int r : rows_[a]) qdd[u_index_[r]].d.v = precision_.inverse_at(r, a);
      }
      model_.log_posterior(qdd.data(), gdd.data());
      for (std::size_t i = 0; i < n; i++) d_phi[i] += 0.5 * gdd[i].d.d;
    }

    // s = A^{-1} r, then d2f/dtheta du s by one Hessian-vector product.
    std::vector<double> s(n_u);
    for (int r = 0; r < n_u; r++) s[r] = g_[u_index_[r]] + d_phi[u_index_[r]];
    precision_.solve(s.data());
    std::vector<D1> qd(n), gd(n);
    for (std::size_t i = 0; i < n; i++) qd[i] = D1(q_[i], 0.0);
    for (int r = 0; r < n_u; r++) qd[u_index_[r]].d = s[r];
    model_.log_posterior(qd.data(), gd.data());
    for (int j = 0; j < p; j++) {
      gradient_[j] += d_phi[theta_index_[j]] + gd[theta_index_[j]].d;
    }
  }
};

}  // namespace phaser


using namespace Rcpp;

namespace {

std::vector<double> theta_vector(const phaser::LaplaceMarginal& lap,
                                 const NumericVector& theta) {
  if (theta.size() != lap.n_theta()) {
    stop("`theta` has length %d; the model has %d parameters outside the "
         "random effects.", static_cast<int>(theta.size()), lap.n_theta());
  }
  return as<std::vector<double> >(theta);
}

NumericMatrix row_major_matrix(const std::vector<double>& x, int n_row,
                               int n_col) {
  NumericMatrix M(n_row, n_col);
  for (int i = 0; i < n_row; i++) {
    for (int j = 0; j < n_col; j++) M(i, j) = x[static_cast<std::size_t>(i) * n_col + j];
  }
  return M;
}

IntegerVector one_based(const std::vector<int>& index) {
  IntegerVector out(index.size());
  for (std::size_t i = 0; i < index.size(); i++) out[i] = index[i] + 1;
  return out;
}

}  // namespace


// [[Rcpp::export]]
SEXP phaseR_laplace(List target) {
  return XPtr<phaser::LaplaceMarginal>(new phaser::LaplaceMarginal(target), true);
}


// [[Rcpp::export]]
List phaseR_laplace_layout(SEXP laplace) {
  XPtr<phaser::LaplaceMarginal> lap(laplace);
  return List::create(Named("theta") = one_based(lap->theta_index()),
                      Named("u") = one_based(lap->u_index()));
}


// [[Rcpp::export]]
List phaseR_laplace_objective(SEXP laplace, NumericVector theta) {
  XPtr<phaser::LaplaceMarginal> lap(laplace);
  std::vector<double> grad;
  double value = lap->evaluate(theta_vector(*lap, theta), &grad);
  return List::create(Named("value") = value, Named("gradient") = wrap(grad),
                      Named("converged") = lap->inner_converged());
}


// [[Rcpp::export]]
NumericMatrix phaseR_laplace_hessian(SEXP laplace, NumericVector theta) {
  XPtr<phaser::LaplaceMarginal> lap(laplace);
  const int p = lap->n_theta();
  return row_major_matrix(lap->hessian(theta_vector(*lap, theta)), p, p);
}


// [[Rcpp::export]]
List phaseR_laplace_conditional(SEXP laplace, NumericVector theta) {
  XPtr<phaser::LaplaceMarginal> lap(laplace);
  std::vector<double> mode, variance, slope;
  lap->conditional(theta_vector(*lap, theta), mode, variance, slope);
  const int n_u = lap->n_u_count();
  return List::create(
    Named("mode") = wrap(mode),
    Named("variance") = wrap(variance),
    Named("slope") = row_major_matrix(slope, n_u, lap->n_theta()));
}


// [[Rcpp::export]]
NumericMatrix phaseR_laplace_effect_noise(SEXP laplace, NumericVector theta,
                                          NumericMatrix z) {
  XPtr<phaser::LaplaceMarginal> lap(laplace);
  const int n_u = lap->n_u_count();
  if (z.ncol() != n_u) stop("`z` needs one column per random effect.");
  const int n = z.nrow();
  std::vector<double> x(static_cast<std::size_t>(n) * n_u);
  for (int i = 0; i < n; i++) {
    for (int r = 0; r < n_u; r++) x[static_cast<std::size_t>(i) * n_u + r] = z(i, r);
  }
  lap->effect_noise(theta_vector(*lap, theta), x, n);
  return row_major_matrix(x, n, n_u);
}


// [[Rcpp::export]]
int phaseR_laplace_multimodal(SEXP laplace, NumericVector theta) {
  XPtr<phaser::LaplaceMarginal> lap(laplace);
  return lap->multimodal_blocks(theta_vector(*lap, theta));
}
