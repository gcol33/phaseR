#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "dense_cholesky.h"
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
// Its gradient is exact:
//
//   d/dtheta_j = df/dtheta_j + tr(A^{-1} D_j H) / 2,
//
// with H = d2f/du2 and D_j the derivative along (e_j, du*/dtheta_j),
// du*/dtheta = A^{-1} d2f/du dtheta (the envelope theorem removes the
// df/du term). Second and third derivatives come from forward-mode dual
// numbers run through the exact gradient of the log posterior.
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

  // At theta: the mode u*, the precision A of the effects given theta, and
  // du*/dtheta (n_u x n_theta, row-major).
  void conditional(const std::vector<double>& theta, std::vector<double>& mode,
                   std::vector<double>& precision, std::vector<double>& slope) {
    std::vector<double> grad;
    evaluate(theta, &grad);
    mode = u_;
    precision = A_;
    slope = du_dtheta_;
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
  std::vector<double> A_;          // -d2f/du2 at q_
  std::vector<double> L_;          // Cholesky factor of A_
  std::vector<double> du_dtheta_;  // n_u x n_theta
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
  void precision_at(const std::vector<double>& q, std::vector<double>& A) const {
    const int n_u = n_u_count();
    A.assign(static_cast<std::size_t>(n_u) * n_u, 0.0);
    std::vector<D1> qd(q.size()), gd(q.size());
    for (const std::vector<int>& colour : colours_) {
      for (std::size_t i = 0; i < q.size(); i++) qd[i] = D1(q[i], 0.0);
      for (int a : colour) qd[u_index_[a]].d = 1.0;
      model_.log_posterior(qd.data(), gd.data());
      for (int a : colour) {
        for (int r : rows_[a]) A[r * n_u + a] = -gd[u_index_[r]].d;
      }
    }
  }

  // Maximises f over u at theta_ by Newton's method with backtracking,
  // starting from u = 0, the prior mean of the effects, so that the mode
  // reached, and with it the approximation, is a function of theta alone;
  // sets u_, q_, g_, A_, L_ and value_.
  void find_mode() {
    const int n_u = n_u_count();
    q_ = full(theta_, std::vector<double>(n_u, 0.0));
    g_.assign(q_.size(), 0.0);
    double f = model_.log_posterior(q_.data(), g_.data());
    inner_converged_ = true;

    if (n_u > 0) {
      inner_converged_ = false;
      std::vector<double> delta(n_u), q_new(q_.size()), g_new(q_.size());
      double last_decrement = HUGE_VAL;
      for (int iter = 0; iter < 200; iter++) {
        precision_at(q_, A_);
        L_ = A_;
        double damping = 0.0;
        while (!cholesky(L_, n_u)) {
          double max_diag = 0.0;
          for (int i = 0; i < n_u; i++) max_diag = std::max(max_diag, std::fabs(A_[i * n_u + i]));
          damping = std::max(1e-8 * (1.0 + max_diag), 10.0 * damping);
          L_ = A_;
          for (int i = 0; i < n_u; i++) L_[i * n_u + i] += damping;
        }
        for (int r = 0; r < n_u; r++) delta[r] = g_[u_index_[r]];
        cholesky_solve(L_, n_u, delta.data());
        double decrement = 0.0;
        for (int r = 0; r < n_u; r++) decrement += g_[u_index_[r]] * delta[r];
        // The Newton decrement is twice the remaining ascent of the local
        // quadratic model and falls quadratically near the mode. Iterate
        // until it vanishes or, once small, stops falling: the iterate is
        // then at the resolution of double precision.
        const double scale = 1.0 + std::fabs(f);
        if (damping == 0.0 &&
            (decrement < 1e-24 * scale ||
             (decrement < 1e-12 * scale && decrement > 0.25 * last_decrement))) {
          inner_converged_ = true;
          break;
        }
        last_decrement = damping == 0.0 ? decrement : HUGE_VAL;
        double t = 1.0;
        double f_new = f;
        bool accepted = false;
        for (int halving = 0; halving < 60; halving++) {
          q_new = q_;
          for (int r = 0; r < n_u; r++) q_new[u_index_[r]] += t * delta[r];
          f_new = model_.log_posterior(q_new.data(), g_new.data());
          if (std::isfinite(f_new) && f_new >= f + 1e-4 * t * decrement) {
            accepted = true;
            break;
          }
          t *= 0.5;
        }
        if (!accepted) {
          // No ascent along the Newton direction: the gradient is at the
          // limit of floating-point resolution.
          inner_converged_ = damping == 0.0;
          break;
        }
        q_.swap(q_new);
        g_.swap(g_new);
        f = f_new;
      }
      for (int r = 0; r < n_u; r++) u_[r] = q_[u_index_[r]];
      precision_at(q_, A_);
      L_ = A_;
      if (!cholesky(L_, n_u)) {
        Rcpp::stop("The log posterior has no strict maximum in the random "
                   "effects at these parameters; the Laplace approximation "
                   "is not defined. Use backend = \"hmc\".");
      }
    }
    value_ = f + 0.5 * n_u * kLog2Pi - 0.5 * cholesky_log_det(L_, n_u);
  }

  // Total derivative of the log marginal along (e_j, du*/dtheta_j). The
  // df/du term vanishes at the exact mode and is kept for the residual of
  // the Newton iteration.
  void compute_gradient() {
    const int p = n_theta();
    const int n_u = n_u_count();
    const std::size_t n = q_.size();
    gradient_.assign(p, 0.0);
    for (int j = 0; j < p; j++) gradient_[j] = g_[theta_index_[j]];
    du_dtheta_.assign(static_cast<std::size_t>(n_u) * p, 0.0);
    if (n_u == 0) return;

    // du*/dtheta_j = A^{-1} d2f/du dtheta_j.
    std::vector<D1> qd(n), gd(n);
    std::vector<double> column(n_u);
    for (int j = 0; j < p; j++) {
      for (std::size_t i = 0; i < n; i++) qd[i] = D1(q_[i], 0.0);
      qd[theta_index_[j]].d = 1.0;
      model_.log_posterior(qd.data(), gd.data());
      for (int r = 0; r < n_u; r++) column[r] = gd[u_index_[r]].d;
      cholesky_solve(L_, n_u, column.data());
      for (int r = 0; r < n_u; r++) {
        du_dtheta_[r * p + j] = column[r];
        gradient_[j] += g_[u_index_[r]] * column[r];
      }
    }

    // tr(A^{-1} D_j H) / 2: third derivatives along the total direction
    // (e_j, du*/dtheta_j) and the colour seeds, contracted with A^{-1} on
    // the Hessian's sparsity pattern.
    const std::vector<double> A_inv = cholesky_inverse(L_, n_u);
    std::vector<D2> qdd(n), gdd(n);
    for (int j = 0; j < p; j++) {
      double trace = 0.0;
      for (const std::vector<int>& colour : colours_) {
        for (std::size_t i = 0; i < n; i++) qdd[i] = D2(D1(q_[i], 0.0), D1(0.0, 0.0));
        qdd[theta_index_[j]].d.v = 1.0;
        for (int r = 0; r < n_u; r++) qdd[u_index_[r]].d.v = du_dtheta_[r * p + j];
        for (int a : colour) qdd[u_index_[a]].v.d = 1.0;
        model_.log_posterior(qdd.data(), gdd.data());
        for (int a : colour) {
          for (int r : rows_[a]) trace += A_inv[r * n_u + a] * gdd[u_index_[r]].d.d;
        }
      }
      gradient_[j] += 0.5 * trace;
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
  std::vector<double> mode, precision, slope;
  lap->conditional(theta_vector(*lap, theta), mode, precision, slope);
  const int n_u = lap->n_u_count();
  return List::create(
    Named("mode") = wrap(mode),
    Named("precision") = row_major_matrix(precision, n_u, n_u),
    Named("slope") = row_major_matrix(slope, n_u, lap->n_theta()));
}
