#ifndef PHASER_EMISSION_H
#define PHASER_EMISSION_H

#include <cmath>
#include "scalar_math.h"

namespace phaser {

// Outcome families of the phase dynamics (codes as in encode_family() in R).
enum Family { kGaussian = 0, kPoisson = 1, kBinomial = 2 };

// One observation of a phase's outcome: its value and, for the binomial
// family, the number of trials (the outcome is the proportion of successes).
// The parts of the log density that do not depend on the parameters are
// computed once.
struct Observation {
  double y;
  int count;           // Poisson count, or binomial successes
  int n_trials;
  double log_const;    // -log(y!) (Poisson) or log(choose(n, k)) (binomial)

  Observation(double y_, int family, int n_trials_)
    : y(y_), count(0), n_trials(n_trials_), log_const(0.0) {
    if (family == kPoisson) {
      count = static_cast<int>(y);
      log_const = -std::lgamma(count + 1.0);
    } else if (family == kBinomial) {
      count = static_cast<int>(y * n_trials + 0.5);
      log_const = std::lgamma(n_trials + 1.0) - std::lgamma(count + 1.0) -
                  std::lgamma(n_trials - count + 1.0);
    }
  }
};

// Log density of an observation given the linear predictor mu of its phase:
// Normal(mu, exp(log_sigma)) for the Gaussian family, Poisson(exp(mu)) and
// Binomial(n_trials, inv_logit(mu)). log_sigma enters only the Gaussian
// family. Writes the derivatives with respect to mu and log_sigma into
// d_mu and d_log_sigma.
template <class T>
T log_emission(int family, const Observation& obs, const T& mu,
               const T& log_sigma, T& d_mu, T& d_log_sigma) {
  using std::exp;
  switch (family) {
    case kPoisson: {
      T rate = exp(mu);
      d_mu = obs.count - rate;
      d_log_sigma = T(0.0);
      return obs.count * mu - rate + obs.log_const;
    }
    case kBinomial: {
      T log_p = log_inv_logit(mu);
      d_mu = obs.count - obs.n_trials * exp(log_p);
      d_log_sigma = T(0.0);
      return obs.log_const + obs.count * log_p +
             (obs.n_trials - obs.count) * log1m_inv_logit(mu);
    }
    default: {
      T inv_sigma = exp(-log_sigma);
      T z = (obs.y - mu) * inv_sigma;
      d_mu = z * inv_sigma;
      d_log_sigma = z * z - 1.0;
      return -0.5 * kLog2Pi - log_sigma - 0.5 * z * z;
    }
  }
}

}  // namespace phaser

#endif
