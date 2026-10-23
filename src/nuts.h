#ifndef PHASER_NUTS_H
#define PHASER_NUTS_H

#include <Rcpp.h>
#include <vector>
#include "log_density.h"
#include "rng.h"

namespace phaser {

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
Rcpp::List run_nuts(const LogDensity& target, Rng& rng, const NutsSettings& settings);

}  // namespace phaser

#endif
