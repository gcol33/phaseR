#ifndef PHASER_NUTS_H
#define PHASER_NUTS_H

#include <Rcpp.h>
#include <vector>
#include "rng.h"

namespace phaser {

// Log density of a target distribution on an unconstrained parameter space,
// known up to an additive constant. Subclasses implement evaluate(); the
// gradient is obtained by central finite differences of evaluate(), so every
// model only has to provide its log posterior. A NaN log density is treated
// as -Inf, i.e. a point outside the support.
class LogDensity {
public:
  explicit LogDensity(int dim) : dim_(dim), buffer_(dim) {}
  virtual ~LogDensity() {}

  int dim() const { return dim_; }

  double log_density(const std::vector<double>& q);

  // Log density at q; writes its gradient into grad.
  double log_density_gradient(const std::vector<double>& q,
                              std::vector<double>& grad);

protected:
  virtual double evaluate(const Rcpp::NumericVector& q) = 0;

private:
  int dim_;
  Rcpp::NumericVector buffer_;
};

struct NutsSettings {
  int n_iter;
  int n_warmup;
  int max_treedepth;
  double adapt_delta;
};

// No-U-Turn sampler (Hoffman and Gelman 2014) in the multinomial form used by
// Stan (Betancourt 2017): trajectories are built by repeated doubling in a
// random direction, states are drawn from the trajectory with probability
// proportional to exp(-H), and doubling stops at a U-turn of any subtree
// (generalised criterion on the summed momentum), at a divergence
// (H - H0 > 1000) or at max_treedepth. Warmup adapts the step size by dual
// averaging towards an average acceptance statistic of adapt_delta and a
// diagonal inverse metric from the posterior variance in Stan's expanding
// windows (75 / 25-doubling / 50 iterations).
//
// Returns the post-warmup draws and per-iteration diagnostics.
Rcpp::List run_nuts(LogDensity& target, Rng& rng, const NutsSettings& settings);

}  // namespace phaser

#endif
