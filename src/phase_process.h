#ifndef PHASER_PHASE_PROCESS_H
#define PHASER_PHASE_PROCESS_H

#include "scalar_math.h"

namespace phaser {

// Latent phase process of one unit observed at n times, with K ordered
// phases: the unit is in phase 0 at its first observation; between
// consecutive observations it either stays in its phase or moves to the
// next one, and the last phase is absorbing. The arrays are indexed by
// observation s = 0, ..., n - 1 and phase k as [s * K + k]:
//
//   log_emit  log density of observation s given phase k;
//   log_stay  log probability of staying in phase k between observations
//             s - 1 and s (s >= 1, k < K - 1);
//   log_move  log probability of moving from phase k to phase k + 1 between
//             observations s - 1 and s (s >= 1, k < K - 1).
//
// The emission and transition terms are computed by the model from its
// families and linear predictors; the recursion only combines them.

// Phases reachable at observation s.
inline int n_reachable(int s, int K) { return (s + 1 < K) ? s + 1 : K; }

// Forward recursion: alpha[s * K + k] = log p(y_0, ..., y_s, phase_s = k)
// for the reachable phases k < n_reachable(s, K). Returns the log
// likelihood of the unit's observations, log sum_k exp(alpha[(n-1) * K + k]).
template <class T>
T forward(int n, int K, const T* log_emit, const T* log_stay,
          const T* log_move, T* alpha) {
  alpha[0] = log_emit[0];
  for (int s = 1; s < n; s++) {
    const T* prev = alpha + (s - 1) * K;
    T* cur = alpha + s * K;
    const int reach_prev = n_reachable(s - 1, K);
    const int reach = n_reachable(s, K);
    for (int k = 0; k < reach; k++) {
      const bool from_same = k < reach_prev;
      const bool from_previous = k > 0;
      T in_same, in_previous;
      if (from_same) {
        in_same = (k < K - 1) ? T(prev[k] + log_stay[s * K + k]) : prev[k];
      }
      if (from_previous) {
        in_previous = prev[k - 1] + log_move[s * K + k - 1];
      }
      T log_in;
      if (from_same && from_previous) {
        log_in = log_sum_exp(in_same, in_previous);
      } else if (from_same) {
        log_in = in_same;
      } else {
        log_in = in_previous;
      }
      cur[k] = log_in + log_emit[s * K + k];
    }
  }
  return log_sum_exp(alpha + (n - 1) * K, n_reachable(n - 1, K));
}

}  // namespace phaser

#endif
