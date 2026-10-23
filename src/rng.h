#ifndef PHASER_RNG_H
#define PHASER_RNG_H

#include <Rcpp.h>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

namespace phaser {

// Random number stream for one MCMC chain.
//
// The 64-bit Mersenne Twister is seeded from the pair (seed, chain) through
// std::seed_seq, so chains that share a user seed run on decorrelated streams
// and a chain's draws do not depend on how many chains run or in what order.
// Uniform and normal variates are computed here from the raw engine output
// rather than through std:: distributions, whose algorithms are left to the
// standard library implementation, so the stream for a given (seed, chain)
// does not change with the standard library the package is built against.
class Rng {
public:
  Rng(int seed, int chain) : has_spare_(false), spare_(0.0) {
    std::seed_seq seq{static_cast<std::uint32_t>(seed),
                      static_cast<std::uint32_t>(chain)};
    engine_.seed(seq);
  }

  // Uniform on [0, 1) with 53 random bits.
  double uniform() {
    return static_cast<double>(engine_() >> 11) * 0x1.0p-53;
  }

  // Standard normal by the Marsaglia polar method.
  double normal() {
    if (has_spare_) {
      has_spare_ = false;
      return spare_;
    }
    double u, v, s;
    do {
      u = 2.0 * uniform() - 1.0;
      v = 2.0 * uniform() - 1.0;
      s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);
    double scale = std::sqrt(-2.0 * std::log(s) / s);
    spare_ = v * scale;
    has_spare_ = true;
    return u * scale;
  }

  // Initial values: `centre` plus independent Normal(0, 0.1^2) noise per
  // coordinate.
  std::vector<double> initial_values(const std::vector<double>& centre) {
    std::vector<double> q(centre.size());
    for (std::size_t i = 0; i < q.size(); i++) {
      q[i] = centre[i] + 0.1 * normal();
    }
    return q;
  }

private:
  std::mt19937_64 engine_;
  bool has_spare_;
  double spare_;
};

}  // namespace phaser

#endif
