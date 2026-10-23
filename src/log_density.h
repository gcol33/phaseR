#ifndef PHASER_LOG_DENSITY_H
#define PHASER_LOG_DENSITY_H

#include <cmath>
#include <limits>

namespace phaser {

// Log density of a target distribution on an unconstrained parameter space,
// known up to an additive constant, with its exact gradient. Subclasses
// implement evaluate(); a NaN log density is reported as -Inf, i.e. a point
// outside the support.
class LogDensity {
public:
  explicit LogDensity(int dim) : dim_(dim) {}
  virtual ~LogDensity() {}

  int dim() const { return dim_; }

  double log_density(const double* q) const {
    return finite_or_minus_inf(evaluate(q, nullptr));
  }

  // Log density at q; writes its gradient into grad (length dim()).
  double log_density_gradient(const double* q, double* grad) const {
    return finite_or_minus_inf(evaluate(q, grad));
  }

protected:
  // Log density at q, and its gradient into grad unless grad is null.
  virtual double evaluate(const double* q, double* grad) const = 0;

private:
  int dim_;

  static double finite_or_minus_inf(double lp) {
    return std::isnan(lp) ? -std::numeric_limits<double>::infinity() : lp;
  }
};

}  // namespace phaser

#endif
